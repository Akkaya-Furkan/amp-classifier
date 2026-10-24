"""Model zoo: (feature kind, estimator) pairs."""
from sklearn.ensemble import HistGradientBoostingClassifier, RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler


def zoo(seed=0):
    lr = lambda: make_pipeline(StandardScaler(), LogisticRegression(max_iter=2000, C=1.0))
    return {
        "length_only (LR)": ("length", lr),
        "length_only (HGB)": ("length", lambda: HistGradientBoostingClassifier(random_state=seed, max_iter=100)),
        "physchem (LR)": ("phys", lr),
        "AAC+physchem (LR)": ("aac+phys", lr),
        "AAC+physchem (RF)": ("aac+phys", lambda: RandomForestClassifier(
            n_estimators=300, n_jobs=-1, random_state=seed, min_samples_leaf=2)),
        "dipeptide+AAC+physchem (HGB)": ("dipep+aac+phys", lambda: HistGradientBoostingClassifier(
            random_state=seed, max_iter=300, learning_rate=0.08)),
    }
