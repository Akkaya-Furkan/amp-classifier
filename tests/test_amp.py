import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from amp.data import clean, length_matched_sample
from amp.features import composition, dipeptides, featurize, hydrophobic_moment, isoelectric_point, net_charge
from amp.split import cdhit_clusters, grouped_folds, parse_clstr
import pandas as pd


def test_clean_filters_nonstandard_and_length():
    assert clean("KWKLFKKIEK") == "KWKLFKKIEK"
    assert clean("KWKLFKKIEX") is None
    assert clean("KWK") is None
    assert clean("kwklfkkiek") == "KWKLFKKIEK"


def test_charge_signs_and_pi():
    assert net_charge("KKKKKKKKKK") > 5
    assert net_charge("DDDDDDDDDD") < -5
    assert isoelectric_point("KKKKKKKKKK") > 10 > isoelectric_point("DDDDDDDDDD") + 5


def test_composition_sums_to_one():
    assert composition("ACDEFGHIKL").sum() == pytest.approx(1)
    assert dipeptides("ACDEFGHIKL").sum() == pytest.approx(1)


def test_hydrophobic_moment_amphipathic_higher():
    helix = "LKKLLKKLLKKLLKKL"  # period-3.6 amphipathic pattern
    flat = "LLLLLLLLLLLLLLLL"
    assert hydrophobic_moment(helix) > hydrophobic_moment(flat)


def test_featurize_shapes():
    s = ["KWKLFKKIEKVGQNIRDGIIKAGPAVAVVGQATQIAK", "GLFDIVKKVVGALGSL"]
    assert featurize(s, "phys").shape == (2, 10)
    assert featurize(s, "aac+phys").shape == (2, 30)
    assert featurize(s, "dipep+aac+phys").shape == (2, 430)


def test_parse_clstr():
    txt = ">Cluster 0\n0\t10aa, >s0... *\n1\t9aa, >s2... at 50.00%\n>Cluster 1\n0\t12aa, >s1... *\n"
    assert parse_clstr(txt, 3).tolist() == [0, 1, 0]


def test_length_matching_respects_ratio():
    rng = np.random.default_rng(0)
    neg = pd.DataFrame({"seq": ["A" * n for n in rng.integers(10, 50, 2000)]})
    pos_len = np.array([12] * 20 + [30] * 20)
    out = length_matched_sample(neg, pos_len, 1.0, rng)
    assert abs(len(out) - 40) <= 2


def test_cdhit_and_grouped_folds_have_no_cluster_leak():
    pytest.importorskip("Bio")
    import shutil
    if shutil.which("cd-hit") is None:
        pytest.skip("cd-hit not installed")
    rng = np.random.default_rng(1)
    base = ["".join(rng.choice(list("ACDEFGHIKLMNPQRSTVWY"), 30)) for _ in range(40)]
    seqs = base + [b[:28] + "AA" for b in base]  # near-duplicates must cluster together
    cid = cdhit_clusters(seqs, 0.9, word=5)
    assert all(cid[i] == cid[i + 40] for i in range(40))
    y = np.array([0, 1] * 40)
    for tr, te in grouped_folds(y, cid, 4):
        assert not set(cid[tr]) & set(cid[te])


def test_nn_identity_and_no_alignment():
    pytest.importorskip("Bio")
    from amp.audit import nn_identity
    assert nn_identity("KWKLFKKIEK", ["KWKLFKKIEK", "GGGGGGGGGG"]) [0] == pytest.approx(1.0)
    assert nn_identity("WWWWWWWWWW", ["DDDDDDDDDD"]) == (0.0, 0.0)


def test_cdhit_keeps_ten_residue_sequences():
    import shutil
    if shutil.which("cd-hit") is None:
        pytest.skip("cd-hit not installed")
    seqs = ["ACDEFGHIKL", "MNPQRSTVWY", "KWKLFKKIEKVGQ"]
    assert (cdhit_clusters(seqs, 0.4) >= 0).all()


def test_similarity_edges_and_purge():
    pytest.importorskip("Bio")
    from amp.split import purge_folds, similarity_edges
    rng = np.random.default_rng(3)
    mk = lambda n: "".join(rng.choice(list("ACDEFGHIKLMNPQRSTVWY"), n))
    a = mk(25)
    seqs = [a, a[:12] + "W" + a[13:], a[5:20], mk(25), mk(25), mk(18)]  # 0,1,2 related
    edges, _ = similarity_edges(seqs)
    assert {(0, 1), (0, 2)} <= {tuple(e) for e in edges}
    folds = [(np.array([0, 3, 4]), np.array([1, 2, 5]))]
    (tr, te), = purge_folds(folds, edges)[0]
    assert te.tolist() == [5]


def test_dbaasp_active_rule():
    from amp.data import dbaasp_active
    assert dbaasp_active({"acts": [["MIC", "x", "8", "µg/ml"]]}, 20)
    assert not dbaasp_active({"acts": [["MIC", "x", ">64", "µg/ml"]]}, 20)
    assert not dbaasp_active({"acts": [["IC50", "x", "1", "µg/ml"]]}, 20)
    assert dbaasp_active({"acts": [["MIC", "x", "5", "µM"]]}, 20)       # 5 uM * 2.2 kDa = 11 ug/ml
    assert not dbaasp_active({"acts": [["MIC", "x", "50", "µM"]]}, 20)  # 110 ug/ml


def test_balanced_by_length_exact():
    from amp.data import balanced_by_length
    rng = np.random.default_rng(0)
    pos = pd.DataFrame({"seq": ["A" * n for n in [10] * 5 + [20] * 3]})
    neg = pd.DataFrame({"seq": ["C" * n for n in [10] * 2 + [20] * 9 + [30] * 4]})
    p, n = balanced_by_length(pos, neg, rng)
    assert sorted(p.seq.str.len()) == sorted(n.seq.str.len()) == [10, 10, 20, 20, 20]


def test_group_bootstrap_ci_contains_estimate():
    from amp.evaluate import bootstrap_ci, score
    rng = np.random.default_rng(0)
    y = rng.integers(0, 2, 600); p = y * 0.5 + rng.random(600)
    groups = np.repeat(np.arange(100), 6)
    (lo, hi), _ = bootstrap_ci(y, p, groups, 100)
    assert lo <= score(y, p)["roc_auc"] <= hi


def test_fragments_use_distinct_parents():
    from amp.data import sample_fragments
    rng = np.random.default_rng(0)
    prot = pd.DataFrame({"id": [f"P{i}" for i in range(50)],
                         "seq": ["".join(rng.choice(list("ACDEFGHIKLMNPQRSTVWY"), 150)) for _ in range(50)]})
    fr = sample_fragments(prot, [20] * 40, rng)
    parents = fr.id.str.split(":").str[0]
    assert len(fr) == 40 and parents.is_unique
    assert (fr.seq.str.len() == 20).all()
