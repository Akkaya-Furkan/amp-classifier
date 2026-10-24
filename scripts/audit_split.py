"""Leakage audit of the CD-HIT split: nearest-neighbour identity of test sequences to the training set.

Compares random folds, CD-HIT-cluster-only folds and CD-HIT folds with similarity-edge purging (the split used), with a null
(test sequences shuffled: keeps composition, destroys order); 100 sampled test sequences per fold, all 5 folds.
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.model_selection import StratifiedKFold

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from amp.audit import nn_identity_many
from amp.split import grouped_folds, purge_folds

SEED, N = 0, 100


def main():
    df = pd.read_csv(ROOT / "data" / "processed" / "dataset.csv")
    strat = df.kind.astype("category").cat.codes.values * 2 + df.amp_matched_peptides.values
    rand = list(StratifiedKFold(5, shuffle=True, random_state=SEED).split(np.zeros(len(df)), strat))
    cd = grouped_folds(strat, df.cdhit_cluster.values, 5, SEED)
    grp, _ = purge_folds(cd, np.load(ROOT / "data" / "processed" / "similarity_edges.npy"))
    folds = {}
    for i in range(5):
        folds[f"random/fold{i}"], folds[f"cdhit_only/fold{i}"], folds[f"cdhit_purged/fold{i}"] = rand[i], cd[i], grp[i]
    rng = np.random.default_rng(SEED)
    out = {}
    for name, (tr, te) in folds.items():
        train = df.seq.values[tr].tolist()
        q = rng.choice(te, N, replace=False)
        qs = df.seq.values[q].tolist()
        shuf = ["".join(rng.permutation(list(s))) for s in qs]
        for tag, queries in [(name, qs), (name + "_shuffled_null", shuf)]:
            ids = nn_identity_many(queries, train)
            for j, den in enumerate(["over_shorter", "over_longer"]):
                v = ids[:, j]
                out[f"{tag}/{den}"] = dict(median=float(np.median(v)), p90=float(np.percentile(v, 90)),
                                           frac_ge_0_4=float((v >= 0.4).mean()), frac_ge_0_6=float((v >= 0.6).mean()),
                                           frac_ge_0_8=float((v >= 0.8).mean()))
            print(tag, out[f"{tag}/over_longer"], flush=True)
    summary = {}
    for kind in ["random", "cdhit_only", "cdhit_purged"]:
        for suffix in ["", "_shuffled_null"]:
            for den in ["over_shorter", "over_longer"]:
                rows = [out[f"{kind}/fold{i}{suffix}/{den}"] for i in range(5)]
                summary[f"{kind}{suffix}/{den}"] = {k: float(np.mean([r[k] for r in rows])) for k in rows[0]}
    out["summary_mean_over_folds"] = summary
    (ROOT / "results" / "split_audit.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
