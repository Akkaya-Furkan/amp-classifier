"""Leakage audit: nearest-neighbour identity between test and training sequences."""
import numpy as np
from Bio.Align import PairwiseAligner, substitution_matrices


def make_aligner():
    a = PairwiseAligner()
    a.mode = "local"
    a.substitution_matrix = substitution_matrices.load("BLOSUM62")
    a.open_gap_score, a.extend_gap_score = -11, -1
    return a


def nn_identity(query, train, aligner=None, top=20):
    """Nearest-neighbour identity of `query` in `train`: the 20 best-scoring local alignments are
    re-examined and the maximum is returned for each denominator,
    as (identity/shorter length, identity/longer length). Zeros if nothing aligns."""
    al = aligner or make_aligner()
    scored = sorted(((al.score(query, t), t) for t in train), key=lambda x: -x[0])[:top]
    best = [0.0, 0.0]
    for sc, t in scored:
        if sc <= 0:
            break
        ident = al.align(query, t)[0].counts().identities
        best[0] = max(best[0], ident / min(len(query), len(t)))
        best[1] = max(best[1], ident / max(len(query), len(t)))
    return tuple(best)


def nn_identity_many(queries, train):
    al = make_aligner()
    return np.array([nn_identity(q, train, al) for q in queries])  # shape (n, 2)
