"""Sequence features: amino-acid composition and physicochemical descriptors."""
import math

import numpy as np

AA = "ACDEFGHIKLMNPQRSTVWY"
# Kyte-Doolittle hydropathy
KD = dict(zip(AA, [1.8, 2.5, -3.5, -3.5, 2.8, -0.4, -3.2, 4.5, -3.9, 3.8,
                   1.9, -3.5, -1.6, -3.5, -4.5, -0.8, -0.7, 4.2, -0.9, -1.3]))
# pKa values (EMBOSS-like)
PK_POS = {"K": 10.8, "R": 12.5, "H": 6.5}
PK_NEG = {"D": 3.9, "E": 4.1, "C": 8.5, "Y": 10.1}
PK_NTERM, PK_CTERM = 8.6, 3.6
PHYS_NAMES = ["length", "net_charge_7.4", "pI", "mean_hydropathy", "hydrophobic_moment",
              "frac_hydrophobic", "frac_cationic", "frac_anionic", "frac_aromatic", "frac_gly_pro"]


def net_charge(seq, ph=7.4):
    q = 1 / (1 + 10 ** (ph - PK_NTERM)) - 1 / (1 + 10 ** (PK_CTERM - ph))
    for a in seq:
        if a in PK_POS:
            q += 1 / (1 + 10 ** (ph - PK_POS[a]))
        elif a in PK_NEG:
            q -= 1 / (1 + 10 ** (PK_NEG[a] - ph))
    return q


def isoelectric_point(seq):
    lo, hi = 0.0, 14.0
    for _ in range(40):
        mid = (lo + hi) / 2
        if net_charge(seq, mid) > 0:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


def hydrophobic_moment(seq, angle_deg=100.0):
    ang = math.radians(angle_deg)
    s = sum(KD[a] * math.sin(i * ang) for i, a in enumerate(seq))
    c = sum(KD[a] * math.cos(i * ang) for i, a in enumerate(seq))
    return math.hypot(s, c) / len(seq)


def composition(seq):
    n = len(seq)
    return np.array([seq.count(a) / n for a in AA])


def dipeptides(seq):
    idx = {a: i for i, a in enumerate(AA)}
    v = np.zeros(400)
    for x, y in zip(seq, seq[1:]):
        v[idx[x] * 20 + idx[y]] += 1
    return v / max(len(seq) - 1, 1)


def physchem(seq):
    n = len(seq)
    f = lambda s: sum(seq.count(a) for a in s) / n
    return np.array([n, net_charge(seq), isoelectric_point(seq),
                     sum(KD[a] for a in seq) / n, hydrophobic_moment(seq),
                     f("AILMFVWC"), f("KR"), f("DE"), f("FWY"), f("GP")])


def featurize(seqs, kind="aac+phys"):
    """kind: 'length', 'phys', 'aac', 'aac+phys', 'dipep+aac+phys'."""
    rows = []
    for s in seqs:
        parts = []
        if kind == "length":
            parts = [np.array([len(s)])]
        else:
            if "phys" in kind:
                parts.append(physchem(s))
            if "aac" in kind:
                parts.append(composition(s))
            if "dipep" in kind:
                parts.append(dipeptides(s))
        rows.append(np.concatenate(parts))
    return np.vstack(rows)
