"""Base retention models that output a system-independent retention score s(x).

* RetentionOrderModel: pairwise learning-to-rank (XGBoost `rank:pairwise`) with query groups = chromatographic
  systems. It learns a score whose ORDER agrees with the elution order inside each training system, which is the
  quantity that transfers between RP systems (cf. retention-order index of Kretschmer et al., Nat. Methods 2026);
  absolute RTs on the target system are recovered by the calibrant mapping.
* SMRTRegressor: gradient-boosted regression of RT on the METLIN SMRT system (Domingo-Almenara et al. 2019), the
  common "pre-train on SMRT, then project" baseline.
"""
from __future__ import annotations
import numpy as np
import pandas as pd


class FeatureStore:
    """Holds the compound feature matrix (aligned with compounds.parquet rows) and an ik14 -> row index map."""

    def __init__(self, processed_dir: str):
        comp = pd.read_parquet(f"{processed_dir}/compounds.parquet")
        z = np.load(f"{processed_dir}/features.npz", allow_pickle=True)
        self.comp = comp
        self.index = pd.Series(np.arange(len(comp)), index=comp["ik14"].values)
        self.desc = z["desc"].astype(np.float32)
        self.desc_names = list(z["desc_names"])
        self.morgan = z["morgan"]
        self.bits = z["bits"]
        self.ion = z["ion"].astype(np.float32)
        self.ion_names = list(z["ion_names"])
        self.elem = z["elem"].astype(np.float32)
        self.elem_names = list(z["elem_names"])

    def rows(self, ik14) -> np.ndarray:
        return self.index.loc[list(ik14)].values

    def base_matrix(self, rows: np.ndarray) -> np.ndarray:
        return np.hstack([self.desc[rows], self.morgan[rows].astype(np.float32)])

    def physchem(self, rows: np.ndarray) -> pd.DataFrame:
        keep = ["MolLogP", "TPSA", "MolWt", "NumHDonors", "NumHAcceptors", "NumRotatableBonds", "FractionCSP3",
                "NumAromaticRings", "RingCount", "HeavyAtomCount", "MolMR", "LabuteASA"]
        idx = [self.desc_names.index(k) for k in keep]
        df = pd.DataFrame(self.desc[rows][:, idx], columns=keep)
        for j, n in enumerate(self.ion_names):
            df["ion_" + n] = self.ion[rows, j]
        for j, n in enumerate(self.elem_names):
            df[n] = self.elem[rows, j]
        df["n_acidic_groups"] = df[["ion_carboxylic_acid", "ion_sulfonic_acid", "ion_phosphate_acid", "ion_tetrazole_or_acidic_NH"]].sum(1)
        df["n_basic_groups"] = df[["ion_aliphatic_amine_prim", "ion_aliphatic_amine_sec", "ion_aliphatic_amine_tert",
                                   "ion_quaternary_N", "ion_amidine_guanidine"]].sum(1)
        return df


class RetentionOrderModel:
    def __init__(self, n_estimators=600, eta=0.05, max_depth=8, subsample=0.8, colsample_bytree=0.5, min_child_weight=5,
                 n_jobs=4, seed=0):
        self.params = dict(objective="rank:pairwise", eta=eta, max_depth=max_depth, subsample=subsample,
                           colsample_bytree=colsample_bytree, min_child_weight=min_child_weight, tree_method="hist",
                           nthread=n_jobs, seed=seed, lambdarank_pair_method="mean", lambdarank_num_pair_per_sample=8,
                           eval_metric="ndcg")
        self.n_estimators = n_estimators
        self.booster = None

    def fit(self, fs: FeatureStore, train: pd.DataFrame):
        """train: long table with columns system, ik14, rt (one row per system x compound)."""
        import xgboost as xgb
        train = train.sort_values(["system", "rt"]).reset_index(drop=True)
        X = fs.base_matrix(fs.rows(train.ik14))
        # labels: within-system RT rank scaled to [0, 30] integer relevance-free ordering (pairwise only uses order)
        y = train.groupby("system")["rt"].rank(method="average").values.astype(np.float32)
        qid = pd.factorize(train.system)[0]
        d = xgb.DMatrix(X, label=y, qid=qid)
        self.booster = xgb.train(self.params, d, num_boost_round=self.n_estimators)
        return self

    def predict(self, fs: FeatureStore, rows: np.ndarray) -> np.ndarray:
        import xgboost as xgb
        return self.booster.predict(xgb.DMatrix(fs.base_matrix(rows)))


class SMRTRegressor:
    def __init__(self, n_estimators=1500, learning_rate=0.03, num_leaves=127, feature_fraction=0.5, n_jobs=4, seed=0):
        self.kw = dict(n_estimators=n_estimators, learning_rate=learning_rate, num_leaves=num_leaves,
                       colsample_bytree=feature_fraction, subsample=0.8, subsample_freq=1, n_jobs=n_jobs,
                       random_state=seed, verbose=-1)
        self.model = None

    def fit(self, fs: FeatureStore, smrt: pd.DataFrame):
        import lightgbm as lgb
        X = fs.base_matrix(fs.rows(smrt.ik14))
        self.model = lgb.LGBMRegressor(**self.kw).fit(X, smrt.rt.values)
        return self

    def predict(self, fs: FeatureStore, rows: np.ndarray) -> np.ndarray:
        return self.model.predict(fs.base_matrix(rows))
