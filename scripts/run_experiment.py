"""Build the dataset, split by CD-HIT cluster, train/evaluate models, write results/."""
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.metrics import precision_recall_curve, roc_curve
from sklearn.model_selection import StratifiedKFold

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from amp.data import balanced_by_length, load_negatives, load_positives, sample_fragments
from amp.evaluate import bootstrap_ci, score
from amp.features import featurize
from amp.models import zoo
from amp.split import cdhit_clusters, grouped_folds, purge_folds, similarity_edges

SEED = 0
IDENTITY = 0.4
RAW, RES, PROC = ROOT / "data" / "raw", ROOT / "results", ROOT / "data" / "processed"
RES.mkdir(exist_ok=True); PROC.mkdir(parents=True, exist_ok=True)


NEG_PEP, NEG_FRAG = "neg_secreted_peptides", "neg_secreted_fragments"


def build():
    rng = np.random.default_rng(SEED)
    pos = load_positives(RAW)
    sec = load_negatives(RAW / "uniprot_secreted.fasta", NEG_PEP)
    sec = sec[~sec.seq.isin(set(pos.seq))]
    # Whole short Swiss-Prot entries starting with the initiator Met are mostly unprocessed precursors,
    # whereas AMP records are mature peptides; drop them so a leading-M shortcut does not separate classes.
    n_sec_before_m = len(sec)
    sec = sec[~sec.seq.str.startswith("M")]
    prot = load_negatives(RAW / "uniprot_secreted_proteins.fasta", "protein", max_len=800)
    # (1) real secreted peptides, exact length balance: both classes subsampled per length bin
    pos_m, sec_m = balanced_by_length(pos, sec, rng)
    assert (np.sort(pos_m.seq.str.len().values) == np.sort(sec_m.seq.str.len().values)).all()
    # (2) random windows of secreted proteins, one per positive, exactly length matched
    frag = sample_fragments(prot, pos.seq.str.len().values, rng)
    frag = frag[~frag.seq.isin(set(pos.seq)) & ~frag.seq.isin(set(sec.seq))]
    df = pd.concat([pos.assign(kind="AMP"), sec_m.assign(kind=NEG_PEP), frag.assign(kind=NEG_FRAG)],
                   ignore_index=True)
    df["amp_matched_peptides"] = df.seq.isin(set(pos_m.seq)) & (df.kind == "AMP")
    df["parent"] = df.id.str.split(":").str[0].where(df.kind == NEG_FRAG, "")
    df["cdhit_cluster"] = cdhit_clusters(df.seq.tolist(), IDENTITY)
    df["cluster"] = df.cdhit_cluster
    edges, graph = similarity_edges(df.seq.tolist())
    np.save(PROC / "similarity_edges.npy", edges)
    sizes = df.cluster.value_counts()
    stats = {"positives_raw": pos.attrs["counts"],
             "pool_secreted_peptides_unique": len(sec), "pool_secreted_proteins": len(prot),
             "n_pos_all": int((df.kind == "AMP").sum()),
             "n_pos_matched_to_peptides": int(df.amp_matched_peptides.sum()),
             "n_neg_secreted_peptides": int((df.kind == NEG_PEP).sum()),
             "n_neg_secreted_fragments": int((df.kind == NEG_FRAG).sum()),
             "secreted_peptides_dropped_for_leading_M": n_sec_before_m - len(sec),
             "n_cdhit_clusters": int(df.cdhit_cluster.nunique()), "cdhit_identity": IDENTITY,
             "largest_cdhit_cluster": int(sizes.iloc[0]), "similarity_graph": graph}
    return df, stats


def amp_mask(df, neg_kind):
    """AMP rows used together with a given negative type (peptide negatives need a length-balanced subset)."""
    m = (df.kind == "AMP").values
    return m & df.amp_matched_peptides.values if neg_kind == NEG_PEP else m


def oof_predictions(df, folds, train_kind, test_kind, model_name, feats):
    """Train on AMP+train_kind rows of each fold's train part, predict AMP+test_kind rows of its test part."""
    kind_name, make = zoo(SEED)[model_name]
    X, k = feats[kind_name], df.kind.values
    mtr, mte = amp_mask(df, train_kind), amp_mask(df, test_kind)
    out_idx, out_p = [], []
    for tr, te in folds:
        tr = tr[mtr[tr] | (k[tr] == train_kind)]
        te = te[mte[te] | (k[te] == test_kind)]
        m = make().fit(X[tr], df.label.values[tr])
        out_idx.append(te); out_p.append(m.predict_proba(X[te])[:, 1])
    return np.concatenate(out_idx), np.concatenate(out_p)


def main():
    df, stats = build()
    df.to_csv(PROC / "dataset.csv", index=False)
    print(json.dumps(stats, indent=1))
    y = df.label.values
    kinds = {k: featurize(df.seq, k) for k in ["length", "phys", "aac+phys", "dipep+aac+phys"]}

    # Both splits stratify on (class, negative type) so train/test composition is comparable.
    strat = (df.kind.astype("category").cat.codes.values * 2 + df.amp_matched_peptides.values)
    cluster_folds, purged = purge_folds(grouped_folds(strat, df.cluster.values, 5, SEED),
                                        np.load(PROC / "similarity_edges.npy"))
    stats["purged_test_rows_per_fold"] = [int(x) for x in purged]
    stats["test_rows_per_fold_after_purge"] = [int(len(te)) for _, te in cluster_folds]
    random_folds = list(StratifiedKFold(5, shuffle=True, random_state=SEED).split(np.zeros(len(df)), strat))

    rows, curves = [], {}
    settings = [("secreted peptides", NEG_PEP, NEG_PEP),
                ("secreted protein fragments", NEG_FRAG, NEG_FRAG),
                ("train fragments -> test peptides", NEG_FRAG, NEG_PEP),
                ("train peptides -> test fragments", NEG_PEP, NEG_FRAG)]
    for split_name, folds in [("cluster (CD-HIT 40% + purge)", cluster_folds), ("random", random_folds)]:
        for setting, trk, tek in settings:
            if split_name == "random" and "->" in setting:
                continue
            for model in zoo(SEED):
                idx, p = oof_predictions(df, folds, trk, tek, model, kinds)
                yy = y[idx]
                s = score(yy, p)
                (rlo, rhi), (alo, ahi) = bootstrap_ci(yy, p, df.cluster.values[idx], 200, SEED)
                rows.append(dict(split=split_name, negatives=setting, model=model, n_test=len(yy), prevalence=float(yy.mean()),
                                 roc_auc=s["roc_auc"], roc_lo=rlo, roc_hi=rhi,
                                 pr_auc=s["pr_auc"], pr_lo=alo, pr_hi=ahi))
                curves[(split_name, setting, model)] = (yy, p)
                print(rows[-1], flush=True)
    tab = pd.DataFrame(rows)
    tab.to_csv(RES / "evaluation_table.csv", index=False)

    (RES / "dataset_stats.json").write_text(json.dumps(stats, indent=1))

    # --- figure
    best = "dipeptide+AAC+physchem (HGB)"
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.4))
    for lab, (sp, st) in {"cluster split": ("cluster (CD-HIT 40% + purge)", "secreted peptides"),
                          "random split": ("random", "secreted peptides")}.items():
        yy, p = curves[(sp, st, best)]
        f, t, _ = roc_curve(yy, p); ax[0].plot(f, t, label=f"{lab} (AUC {score(yy,p)['roc_auc']:.3f})")
        pr, rc, _ = precision_recall_curve(yy, p); ax[1].plot(rc, pr, label=f"{lab} (AP {score(yy,p)['pr_auc']:.3f})")
    ax[0].plot([0, 1], [0, 1], "k:"); ax[0].set(xlabel="FPR", ylabel="TPR", title=f"ROC, secreted-peptide negatives\n{best}")
    ax[1].set(xlabel="Recall", ylabel="Precision", title="Precision-recall")
    ax[0].legend(); ax[1].legend()
    sub = tab[tab.split == "cluster (CD-HIT 40% + purge)"]
    piv = sub.pivot(index="model", columns="negatives", values="roc_auc")
    piv = piv[[x[0] for x in settings]]
    piv.plot.barh(ax=ax[2]); ax[2].set(xlabel="ROC-AUC (cluster split)", xlim=(0.4, 1), title="Effect of negative set")
    ax[2].legend(fontsize=6, loc="upper center", bbox_to_anchor=(0.5, -0.18), ncol=2); fig.tight_layout(); fig.savefig(RES / "evaluation.png", dpi=150)

    md = tab.assign(roc=lambda d: d.apply(lambda r: f"{r.roc_auc:.3f} [{r.roc_lo:.3f}-{r.roc_hi:.3f}]", axis=1),
                    pr=lambda d: d.apply(lambda r: f"{r.pr_auc:.3f} [{r.pr_lo:.3f}-{r.pr_hi:.3f}]", axis=1)
                    )[["split", "negatives", "model", "n_test", "prevalence", "roc", "pr"]].rename(columns={"roc": "ROC-AUC [95% CI]", "pr": "PR-AUC (avg. precision) [95% CI]"})
    (RES / "evaluation_table.md").write_text(md.to_markdown(index=False) if hasattr(md, "to_markdown") else md.to_string())


if __name__ == "__main__":
    main()
