"""Cross-validated evaluation helpers."""
import numpy as np
from sklearn.metrics import average_precision_score, roc_auc_score


def score(y, p):
    """ROC-AUC and average precision (reported as PR-AUC)."""
    return {"roc_auc": roc_auc_score(y, p), "pr_auc": average_precision_score(y, p)}


def bootstrap_ci(y, p, groups, n=200, seed=0):
    """95% CI from a bootstrap that resamples whole similarity groups (not single rows),
    so related sequences are not treated as independent. Conditional on the fitted predictions."""
    rng = np.random.default_rng(seed)
    order = np.argsort(groups, kind="stable")
    g_sorted = groups[order]
    uniq, start = np.unique(g_sorted, return_index=True)
    bounds = list(zip(start, list(start[1:]) + [len(g_sorted)]))
    r, a = [], []
    for _ in range(n):
        pick = rng.integers(0, len(uniq), len(uniq))
        idx = order[np.concatenate([np.arange(*bounds[i]) for i in pick])]
        if len(np.unique(y[idx])) < 2:
            continue
        s = score(y[idx], p[idx])
        r.append(s["roc_auc"]); a.append(s["pr_auc"])
    return np.percentile(r, [2.5, 97.5]), np.percentile(a, [2.5, 97.5])
