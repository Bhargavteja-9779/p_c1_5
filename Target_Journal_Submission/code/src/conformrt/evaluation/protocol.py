"""Leave-source-out protocol, fold-wise model training (with inner cross-fitting) and simulation of the calibrant
protocol on training systems (difficulty-model training data)."""
from __future__ import annotations
import os, math
import numpy as np
import pandas as pd
from joblib import Parallel, delayed

from ..models.base_models import FeatureStore, RetentionOrderModel, SMRTRegressor
from ..models.mapping import MonotoneMap
from ..models.difficulty import DifficultyModel, GROUPS, encode_system, geometry

FOLD_NAMES = {0: "BGC", 1: "Folberth", 2: "PredRet", 3: "Souihi+Stoffel", 4: "Other sources"}


def fold_of_source(src: str) -> int:
    s = str(src)
    if "BGC" in s:
        return 0
    if "Folberth" in s:
        return 1
    if "PredRet" in s:
        return 2
    if "Souihi" in s or "Stoffel" in s:
        return 3
    return 4


class Data:
    def __init__(self, processed_dir: str, method_type: str = "RP"):
        self.fs = FeatureStore(processed_dir)
        rt = pd.read_parquet(f"{processed_dir}/rt_long.parquet")
        meta = pd.read_csv(f"{processed_dir}/systems.csv", dtype={"system": str})
        rt["system"] = rt["system"].astype(str)
        self.meta_all = meta.set_index("system")
        meta = meta[meta.method_type == method_type].copy()
        meta["fold"] = meta.source.map(fold_of_source)
        self.meta = meta.set_index("system")
        self.rt = rt[rt.system.isin(self.meta.index)].copy()
        smrt = self.meta[self.meta.source.str.contains("Domingo", na=False)]
        self.smrt_id = smrt.index[0] if len(smrt) else None
        # compound-level helpers
        comp = self.fs.comp
        self.formula = comp.set_index("ik14")["formula"]
        self.n_systems = comp.set_index("ik14")["n_systems"]
        self.superclass = comp.set_index("ik14")["superclass"]
        self._bitvects = None
        self._freq = {}

    def fold_frequency(self, exclude) -> pd.Series:
        """Number of systems of this separation mode outside the excluded source groups in which each compound was measured:
        the public frequency available to a laboratory whose own source group is held out ('common' calibrant rule)."""
        key = frozenset(int(x) for x in exclude)
        if key not in self._freq:
            keep = self.meta.index[~self.meta.fold.isin(key)]
            self._freq[key] = self.rt[self.rt.system.isin(keep)].groupby("ik14").system.nunique()
        return self._freq[key]

    def targets(self, min_compounds: int = 60) -> pd.DataFrame:
        m = self.meta[(self.meta.n_compounds >= min_compounds)]
        if self.smrt_id is not None:
            m = m.drop(index=self.smrt_id, errors="ignore")
        return m

    def training_table(self, exclude_folds, smrt_n: int = 15000, seed: int = 0) -> pd.DataFrame:
        keep = self.meta[~self.meta.fold.isin(exclude_folds)].index
        tr = self.rt[self.rt.system.isin(keep)]
        if self.smrt_id is not None:
            sm = self.rt[self.rt.system == self.smrt_id]
            if len(sm) > smrt_n:
                sm = sm.sample(smrt_n, random_state=seed)
            tr = pd.concat([tr[tr.system != self.smrt_id], sm])
        tr = tr[tr.groupby("system").ik14.transform("size") >= 3]
        return tr

    # ---------------- novelty (applicability-domain) features ----------------
    def bitvects(self):
        if self._bitvects is None:
            from rdkit import DataStructs
            bv = []
            for row in self.fs.bits:
                v = DataStructs.ExplicitBitVect(1024)
                on = np.flatnonzero(np.unpackbits(row))
                v.SetBitsFromList(on.tolist())
                bv.append(v)
            self._bitvects = bv
        return self._bitvects


_BV = None
_REF = None


def _maxsim_chunk(qidx):
    from rdkit import DataStructs
    return [max(DataStructs.BulkTanimotoSimilarity(_BV[q], _REF)) if _REF else 0.0 for q in qidx]


def max_similarity(data: Data, query_rows: np.ndarray, ref_rows: np.ndarray, n_jobs: int = 4) -> np.ndarray:
    """Max Tanimoto (Morgan r=2, 1024 bits) of each query compound to the reference (training) compounds."""
    global _BV, _REF
    _BV = data.bitvects(); _REF = [_BV[r] for r in ref_rows]
    chunks = np.array_split(np.asarray(query_rows), max(1, n_jobs * 8))
    if n_jobs > 1:
        import multiprocessing as mp
        ctx = mp.get_context("fork")
        with ctx.Pool(n_jobs) as pool:
            res = pool.map(_maxsim_chunk, chunks)
    else:
        res = [_maxsim_chunk(c) for c in chunks]
    return np.array([v for r in res for v in r], dtype=np.float32)


# ---------------- fold-wise base models ----------------
class FoldModels:
    """Trains and caches base models. key = frozenset of excluded folds."""

    def __init__(self, data: Data, cache_dir: str, base_kw: dict, smrt_n: int, n_jobs: int, log=None):
        self.data, self.cache_dir, self.base_kw, self.smrt_n, self.n_jobs, self.log = data, cache_dir, base_kw, smrt_n, n_jobs, log
        os.makedirs(cache_dir, exist_ok=True)
        self.models = {}

    def get(self, exclude: frozenset):
        if exclude in self.models:
            return self.models[exclude]
        import joblib
        tag = "excl_" + "_".join(str(f) for f in sorted(exclude)) + f"_t{self.base_kw.get('n_estimators')}"
        path = os.path.join(self.cache_dir, tag + ".joblib")
        tr = self.data.training_table(exclude, self.smrt_n)
        train_ik = np.unique(tr.ik14.values)
        if os.path.exists(path):
            model = joblib.load(path)
        else:
            if self.log: self.log.info(f"training base model excluding folds {sorted(exclude)}: {len(tr)} rows, {tr.system.nunique()} systems")
            model = RetentionOrderModel(n_jobs=self.n_jobs, **self.base_kw).fit(self.data.fs, tr)
            joblib.dump(model, path)
        self.models[exclude] = (model, train_ik)
        return self.models[exclude]


def simulate_calibration_rows(data: Data, systems, scores: dict, novelty: dict, n_list=(10, 20, 40), draws=3,
                              min_compounds=30, seed=0) -> pd.DataFrame:
    """Simulate the calibrant protocol on (training) systems whose base scores are out-of-system.
    Returns one row per held-out compound per draw with σ-model features and target log(|r|+0.01)."""
    rows = []
    fs = data.fs
    for sid in systems:
        tab = data.rt[data.rt.system == sid]
        N = len(tab)
        if N < min_compounds:
            continue
        ik = tab.ik14.values; y = tab.rt.values; s = scores[sid]
        frows = fs.rows(ik)
        pc = fs.physchem(frows)
        sysf = encode_system(data.meta.loc[sid].to_dict())
        ms, seen = novelty[sid]
        for n in n_list:
            if N < n + 10:
                continue
            for r in range(draws):
                rng = np.random.default_rng([int(sid), n, r, seed])
                cal = rng.choice(N, n, replace=False); te = np.setdiff1d(np.arange(N), cal)
                mp = MonotoneMap().fit(s[cal], y[cal])
                res = y[te] - mp.predict(s[te])
                g = geometry(s[te], s[cal], y[cal])
                d = pc.iloc[te].reset_index(drop=True).copy()
                d["maxsim"] = ms[te]; d["seen"] = seen[te]
                for k in ["pos", "extrap", "gap"]:
                    d[k] = g[k]
                d["cal_span"] = g["cal_span"]; d["n_cal"] = n
                for k, v in sysf.items():
                    d[k] = v
                d["target"] = np.log(np.abs(res) + 0.01)
                d["system"] = sid; d["n"] = n; d["draw"] = r
                rows.append(d)
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()
