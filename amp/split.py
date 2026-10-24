"""Similarity-aware splitting with CD-HIT clusters."""
import subprocess
import tempfile
from pathlib import Path

import numpy as np
from collections import defaultdict
from sklearn.model_selection import StratifiedGroupKFold


def cdhit_clusters(seqs, identity=0.4, word=2, threads=4):
    """Cluster sequences with cd-hit; return an array of cluster ids aligned to `seqs`."""
    with tempfile.TemporaryDirectory() as d:
        fa, out = Path(d) / "in.fa", Path(d) / "out"
        fa.write_text("".join(f">s{i}\n{s}\n" for i, s in enumerate(seqs)))
        subprocess.run(["cd-hit", "-i", str(fa), "-o", str(out), "-c", str(identity),
                        "-n", str(word), "-l", "5", "-T", str(threads), "-M", "4000", "-d", "0"],
                       check=True, capture_output=True)
        return parse_clstr(Path(str(out) + ".clstr").read_text(), len(seqs))


def parse_clstr(text, n):
    cid = np.full(n, -1)
    cur = -1
    for line in text.splitlines():
        if line.startswith(">Cluster"):
            cur = int(line.split()[1])
        elif line.strip():
            idx = int(line.split(">s")[1].split("...")[0])
            cid[idx] = cur
    if (cid < 0).any():
        raise ValueError("unassigned sequences in cd-hit output")
    return cid


def grouped_folds(y, groups, n_splits=5, seed=0):
    """Stratified folds in which no cluster spans train and test."""
    cv = StratifiedGroupKFold(n_splits=n_splits, shuffle=True, random_state=seed)
    return list(cv.split(np.zeros(len(y)), y, groups))


def similarity_edges(seqs, k=5, max_df=60, min_ident_shorter=0.8, min_ident_longer=0.5):
    """Pairs of sequences that are near-duplicates / strong local homologs that CD-HIT may place in
    different clusters.

    Candidates: pairs sharing an exact k-mer (k-mers present in more than `max_df` sequences are
    skipped as low-complexity, which bounds the pair count). An edge is kept if the best local alignment
    (BLOSUM62) has identical-position identity >= `min_ident_shorter` of the shorter sequence
    (containment / near-duplicate) or >= `min_ident_longer` of the longer one.
    Returns (edges array of shape (m, 2), info dict).
    """
    from .audit import make_aligner
    index = defaultdict(list)
    for i, s in enumerate(seqs):
        for km in {s[j:j + k] for j in range(len(s) - k + 1)}:
            index[km].append(i)
    pairs = set()
    for ids in index.values():
        if 1 < len(ids) <= max_df:
            for a in range(len(ids)):
                for b in range(a + 1, len(ids)):
                    pairs.add((ids[a], ids[b]))
    al = make_aligner()
    edges = []
    for a, b in sorted(pairs):
        sa, sb = seqs[a], seqs[b]
        if al.score(sa, sb) <= 0:
            continue
        ident = al.align(sa, sb)[0].counts().identities
        if ident / min(len(sa), len(sb)) >= min_ident_shorter or ident / max(len(sa), len(sb)) >= min_ident_longer:
            edges.append((a, b))
    return np.array(edges, dtype=int).reshape(-1, 2), {"candidate_pairs": len(pairs), "edges": len(edges)}


def purge_folds(folds, edges):
    """Remove from each fold's test part every sequence with a similarity edge to that fold's training
    part. Single-linkage grouping of this data collapses into one giant component, so the folds are
    built from CD-HIT clusters and the residual cross-fold similarity is purged instead.
    Returns (purged folds, list of numbers of purged test rows)."""
    out, dropped = [], []
    for tr, te in folds:
        in_tr = np.zeros(int(max(tr.max(), te.max())) + 1, dtype=bool)
        in_tr[tr] = True
        bad = set()
        for a, b in edges:
            if a < len(in_tr) and b < len(in_tr):
                if in_tr[a] and not in_tr[b]:
                    bad.add(b)
                elif in_tr[b] and not in_tr[a]:
                    bad.add(a)
        keep = np.array([i for i in te if i not in bad], dtype=int)
        out.append((tr, keep)); dropped.append(len(te) - len(keep))
    return out, dropped
