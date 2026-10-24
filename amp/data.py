"""Build the labelled peptide table from raw downloads."""
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd

STD = re.compile(r"^[ACDEFGHIKLMNPQRSTVWY]+$")
MIN_LEN, MAX_LEN = 10, 100


def read_fasta(path):
    out, name, buf = [], None, []
    for line in Path(path).read_text().splitlines():
        if line.startswith(">"):
            if name is not None:
                out.append((name, "".join(buf)))
            name, buf = line[1:], []
        else:
            buf.append(line.strip())
    if name is not None:
        out.append((name, "".join(buf)))
    return out


def clean(seq, max_len=MAX_LEN):
    seq = seq.strip().upper() if seq else ""
    return seq if STD.match(seq) and MIN_LEN <= len(seq) <= max_len else None


def candidate_dbaasp(raw):
    """DBAASP monomers eligible as positives before activity filtering: standard upper-case residues
    (lower case = D-amino acids), length 10-100, no N-terminal modification."""
    out = []
    for r in json.loads((Path(raw) / "dbaasp.json").read_text()):
        if r.get("complexity") != "monomer" or r.get("nTerminus"):
            continue
        s = r.get("sequence") or ""
        if s == s.upper() and clean(s):
            out.append(r)
    return out


MIC_MAX_UG_ML = 32.0


def _conc(x):
    m = re.match(r"^\s*([0-9]*\.?[0-9]+)", str(x or ""))
    return None if (m is None or str(x).strip().startswith((">", "<"))) else float(m.group(1))


def dbaasp_active(rec, length):
    """True if any MIC record shows activity at <= 32 ug/ml (uM converted with 110 Da/residue)."""
    for measure, _species, conc, unit in rec["acts"]:
        c = _conc(conc)
        if not str(measure or "").upper().startswith("MIC") or c is None:
            continue
        if unit == "µM":
            c = c * 0.110 * length
        elif unit != "µg/ml":
            continue
        if c <= MIC_MAX_UG_ML:
            return True
    return False


def load_positives(raw):
    """APD3 + activity-confirmed DBAASP monomers, standard amino acids only, exact-sequence dedup.

    DBAASP records are kept only if DBAASP lists a MIC <= 32 ug/ml (see `dbaasp_active`); records with
    no usable MIC, or only weaker ones, are dropped (not used as negatives)."""
    act_path = Path(raw) / "dbaasp_activity.json"
    activity = json.loads(act_path.read_text()) if act_path.exists() else None
    recs = []
    for name, s in read_fasta(Path(raw) / "apd3.fasta"):
        # APD3 lower-case does not occur; D-amino/modified entries are not marked, so
        # only the standard-alphabet filter applies.
        recs.append(("APD3", name.split("|")[0], s))
    n_cand = n_active = n_nodata = 0
    for r in candidate_dbaasp(raw):
        n_cand += 1
        if activity is not None:
            a = activity.get(str(r["id"]))
            if a is None:
                n_nodata += 1
                continue
            if not dbaasp_active(a, len(r["sequence"])):
                continue
        n_active += 1
        recs.append(("DBAASP", r.get("dbaaspId"), r["sequence"]))
    df = pd.DataFrame(recs, columns=["source", "id", "seq"])
    df["seq"] = df["seq"].map(clean)
    n_raw = len(df)
    df = df.dropna(subset=["seq"])
    n_clean = len(df)
    df = df.groupby("seq", as_index=False).agg(source=("source", lambda x: "+".join(sorted(set(x)))),
                                               id=("id", "first"))
    df["label"] = 1
    df.attrs["counts"] = {"raw": n_raw, "after_filter": n_clean, "unique": len(df),
                          "dbaasp_candidates": n_cand, "dbaasp_activity_missing": n_nodata,
                          "dbaasp_mic_active": n_active}
    return df


def load_negatives(path, tag, max_len=MAX_LEN):
    recs = [(tag, name.split("|")[1] if "|" in name else name, clean(s, max_len)) for name, s in read_fasta(path)]
    df = pd.DataFrame(recs, columns=["source", "id", "seq"]).dropna(subset=["seq"])
    df = df.drop_duplicates("seq")
    df["label"] = 0
    return df


def length_matched_sample(neg, pos_lengths, ratio, rng, bins=None):
    """Sample negatives so their length histogram follows the positives' (ratio x).

    Matching is done in length bins; if a bin lacks negatives, all of them are taken
    (the shortfall is reported by the caller through the resulting size).
    """
    bins = bins if bins is not None else np.array([10, 15, 20, 25, 30, 40, 50, 60, 80, 101])
    pb = np.digitize(pos_lengths, bins)
    nb = np.digitize(neg["seq"].str.len().values, bins)
    picks = []
    for b in np.unique(pb):
        want = int(round((pb == b).sum() * ratio))
        pool = np.where(nb == b)[0]
        take = min(want, len(pool))
        picks.extend(rng.choice(pool, take, replace=False))
    return neg.iloc[sorted(picks)]


def balanced_by_length(pos, neg, rng, bins=None):
    """Subsample BOTH classes at each exact length to the smaller count (exact length balance)."""
    bins = bins if bins is not None else np.arange(10, 102)
    pb = np.digitize(pos["seq"].str.len().values, bins)
    nb = np.digitize(neg["seq"].str.len().values, bins)
    pi, ni = [], []
    for b in np.unique(np.concatenate([pb, nb])):
        P, N = np.where(pb == b)[0], np.where(nb == b)[0]
        k = min(len(P), len(N))
        pi.extend(rng.choice(P, k, replace=False)); ni.extend(rng.choice(N, k, replace=False))
    return pos.iloc[sorted(pi)], neg.iloc[sorted(ni)]


def sample_fragments(proteins, lengths, rng, skip=30):
    """One random window per requested length, each from a different protein (so no parent protein
    contributes two windows to different CV folds); first `skip` residues are excluded to avoid
    signal peptides. Returns a negatives frame with label 0."""
    seqs = proteins["seq"].tolist()
    ids = proteins["id"].tolist()
    recs, used = [], set()
    for L in lengths:
        for _ in range(200):
            j = int(rng.integers(len(seqs)))
            if j not in used and len(seqs[j]) >= skip + L + 10:
                used.add(j)
                st = int(rng.integers(skip, len(seqs[j]) - L + 1))
                frag = seqs[j][st:st + L]
                if STD.match(frag):
                    recs.append(("fragment", f"{ids[j]}:{st}-{st+L}", frag))
                    break
    df = pd.DataFrame(recs, columns=["source", "id", "seq"]).drop_duplicates("seq")
    df["label"] = 0
    return df
