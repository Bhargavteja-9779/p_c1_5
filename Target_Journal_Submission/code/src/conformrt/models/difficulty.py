"""Meta-learned difficulty (scale) model σ(x; system, calibrants) for normalised conformal scores.

Trained on residuals obtained by SIMULATING the calibrant protocol on training-fold chromatographic systems whose base
scores come from base models that never saw those systems (inner cross-fitting). Target: log(|residual| + 0.01 min).
Feature groups (used for the ablation study):
  compound  — RDKit physico-chemical descriptors, elemental counts
  ion       — counts of acidic/basic (ionisable) groups + formal charge
  novelty   — max Tanimoto similarity of the compound to the base-model training set; exact-seen flag
  geometry  — position of the compound's base score relative to the calibrant scores (relative position, extrapolation
              distance, gap to nearest calibrant), calibrant RT span, number of calibrants
  system    — mobile-phase pH class, organic modifier, column USP class, column length/particle/temperature/flow,
              gradient duration, dead time
"""
from __future__ import annotations
import numpy as np
import pandas as pd

GROUPS = {
    "compound": ["MolLogP", "TPSA", "MolWt", "NumHDonors", "NumHAcceptors", "NumRotatableBonds", "FractionCSP3",
                 "NumAromaticRings", "RingCount", "HeavyAtomCount", "MolMR", "LabuteASA", "n_N", "n_O", "n_S", "n_P",
                 "n_F", "n_Cl", "n_Br", "n_I", "has_Si"],
    "ion": ["ion_carboxylic_acid", "ion_sulfonic_acid", "ion_phosphate_acid", "ion_phenol", "ion_tetrazole_or_acidic_NH",
            "ion_aliphatic_amine_prim", "ion_aliphatic_amine_sec", "ion_aliphatic_amine_tert", "ion_quaternary_N",
            "ion_amidine_guanidine", "ion_pyridine_like_N", "ion_imidazole_like", "n_acidic_groups", "n_basic_groups",
            "formal_charge"],
    "novelty": ["maxsim", "seen"],
    "geometry": ["pos", "extrap", "gap", "cal_span", "n_cal"],
    "system": ["ph_class", "modifier", "usp", "col_len", "particle", "temp", "flow", "grad_dur", "t0"],
}
CATEGORICAL = ["ph_class", "modifier", "usp"]
CATS = {"ph_class": ["acidic", "buffered", "basic", "none"],
        "modifier": ["acn", "meoh", "mixed", "acn+other", "meoh+other", "mixed+other", "other"],
        "usp": ["L1", "L11", "L7", "L10", "L109", "L43", "other"]}


def encode_system(meta_row: dict) -> dict:
    d = {k: meta_row.get(k) for k in ["ph_class", "modifier", "col_len", "particle", "temp", "flow", "grad_dur", "t0"]}
    d["usp"] = meta_row.get("usp") if meta_row.get("usp") in CATS["usp"] else "other"
    for c in CATEGORICAL:
        v = d[c] if d[c] in CATS[c] else CATS[c][-1]
        d[c] = CATS[c].index(v)
    return d


def geometry(s_query: np.ndarray, s_cal: np.ndarray, y_cal: np.ndarray) -> dict:
    lo, hi = float(np.min(s_cal)), float(np.max(s_cal))
    span = max(hi - lo, 1e-9)
    pos = (s_query - lo) / span
    extrap = np.maximum(0, np.maximum(lo - s_query, s_query - hi)) / span
    gap = np.min(np.abs(s_query[:, None] - np.asarray(s_cal)[None, :]), axis=1) / span
    return {"pos": pos, "extrap": extrap, "gap": gap, "cal_span": float(np.max(y_cal) - np.min(y_cal)), "n_cal": len(s_cal)}


def geometry_loo(s_query: np.ndarray, s_cal: np.ndarray, y_cal: np.ndarray) -> dict:
    """Geometry features for each leave-one-out calibrant set: arrays of shape (n, m)."""
    n = len(s_cal); m = len(s_query)
    order = np.argsort(s_cal); sc = s_cal[order]; yc = y_cal[order]
    lo_all, hi_all = sc[0], sc[-1]
    out = {k: np.empty((n, m)) for k in ["pos", "extrap", "gap"]}
    span_y = np.empty(n)
    d = np.abs(s_query[None, :] - sc[:, None])  # (n, m) distances sorted by calibrant score
    # nearest and second-nearest distance per query (for LOO gap)
    part = np.partition(d, 1, axis=0) if n > 1 else np.vstack([d, d])
    d1, d2 = part[0], part[1]
    argmin = np.argmin(d, axis=0)
    ys = np.sort(y_cal)
    for jj in range(n):
        i = order[jj]
        lo = sc[1] if jj == 0 else lo_all
        hi = sc[-2] if jj == n - 1 else hi_all
        span = max(hi - lo, 1e-9)
        out["pos"][i] = (s_query - lo) / span
        out["extrap"][i] = np.maximum(0, np.maximum(lo - s_query, s_query - hi)) / span
        g = np.where(argmin == jj, d2, d1)
        out["gap"][i] = g / span
        yr = np.delete(yc, jj)
        span_y[i] = yr.max() - yr.min()
    out["cal_span"] = span_y
    out["n_cal"] = n - 1
    return out


class DifficultyModel:
    def __init__(self, groups=("compound", "ion", "novelty", "geometry", "system"), n_estimators=400, learning_rate=0.05,
                 num_leaves=31, min_child_samples=100, n_jobs=4, seed=0):
        self.groups = list(groups)
        self.cols = [c for g in self.groups for c in GROUPS[g]]
        self.kw = dict(n_estimators=n_estimators, learning_rate=learning_rate, num_leaves=num_leaves,
                       min_child_samples=min_child_samples, subsample=0.8, subsample_freq=1, colsample_bytree=0.8,
                       n_jobs=n_jobs, random_state=seed, verbose=-1)

    def fit(self, df: pd.DataFrame, target: np.ndarray):
        import lightgbm as lgb
        cat_idx = [i for i, c in enumerate(self.cols) if c in CATEGORICAL]
        X = df[self.cols].to_numpy(dtype=np.float32)
        self.model = lgb.LGBMRegressor(**self.kw).fit(X, target, categorical_feature=cat_idx or "auto")
        return self

    def sigma(self, X) -> np.ndarray:
        """X: DataFrame (columns by name) or ndarray with columns in `self.cols` order."""
        if isinstance(X, pd.DataFrame):
            X = X[self.cols].to_numpy(dtype=np.float32)
        return np.exp(self.model.predict(X)) + 1e-3
