"""User-facing interface: retention-time windows for candidate structures on a laboratory's own LC system.

    from conformrt.api import load_bundle, predict_windows
    bundle = load_bundle("models/conformrt_rp_bundle.joblib")       # built by scripts/train_deployable_models.py
    windows = predict_windows(calibrants, candidates, bundle, alpha=0.10, method="jackknife+")

calibrants: DataFrame with columns `smiles`, `rt` (min) measured on the target system (at least 1/alpha - 1 rows; the
            guarantee requires calibrants that are exchangeable with the analytes, see the article).
candidates: DataFrame with columns `id`, `smiles` (and optionally `rt_obs`, the observed RT of the feature).
method:     "jackknife+" (recommended default) or "normalised" (jackknife+ scaled by the meta-learned difficulty model;
            requires `system_meta`, a dict with the routinely reported chromatographic metadata, see encode_system()).
Returns the candidates with columns `rt_pred`, `lo`, `hi` and, if `rt_obs` is given, `keep` (rt_obs inside the window).
Distribution-free guarantee: P(t in window) >= 1 - 2*alpha for exchangeable calibrants (empirically close to 1 - alpha)."""
from __future__ import annotations
import numpy as np
import pandas as pd

from .preprocessing.features import featurize, tanimoto_max
from .models.mapping import MonotoneMap
from .models.difficulty import encode_system, geometry_loo
from .evaluation.intervals import jackknife_plus

PHYSCHEM = ["MolLogP", "TPSA", "MolWt", "NumHDonors", "NumHAcceptors", "NumRotatableBonds", "FractionCSP3", "NumAromaticRings",
            "RingCount", "HeavyAtomCount", "MolMR", "LabuteASA"]


def load_bundle(path: str) -> dict:
    import joblib
    return joblib.load(path)


def _features(smiles: pd.Series, bundle: dict, n_jobs: int = 1):
    F = featurize(smiles.reset_index(drop=True), n_jobs=n_jobs)
    if not F["ok"].all():
        bad = smiles.reset_index(drop=True)[~F["ok"]].tolist()
        raise ValueError(f"could not parse SMILES: {bad[:5]}")
    import xgboost as xgb
    X = np.hstack([F["desc"], F["morgan"].astype(np.float32)])
    score = bundle["base_model"].predict(xgb.DMatrix(X))
    names = list(F["desc_names"])
    C = pd.DataFrame(F["desc"][:, [names.index(k) for k in PHYSCHEM]], columns=PHYSCHEM)
    for j, n in enumerate(F["ion_names"]):
        C["ion_" + n] = F["ion"][:, j].astype(np.float32)
    for j, n in enumerate(F["elem_names"]):
        C[n] = F["elem"][:, j].astype(np.float32)
    C["n_acidic_groups"] = C[["ion_carboxylic_acid", "ion_sulfonic_acid", "ion_phosphate_acid", "ion_tetrazole_or_acidic_NH"]].sum(1)
    C["n_basic_groups"] = C[["ion_aliphatic_amine_prim", "ion_aliphatic_amine_sec", "ion_aliphatic_amine_tert", "ion_quaternary_N",
                             "ion_amidine_guanidine"]].sum(1)
    from rdkit import Chem
    from .preprocessing.features import mol_from_smiles
    ik = np.array([Chem.MolToInchiKey(mol_from_smiles(x))[:14] for x in smiles.tolist()])
    seen = np.isin(ik, bundle["train_ik"])  # exact match (InChIKey connectivity block) with a base-model training compound
    maxsim = tanimoto_max(F["bits"], bundle["train_bits"])
    C["maxsim"] = np.where(seen, 1.0, maxsim); C["seen"] = seen.astype(np.float32)
    return score, C


def predict_windows(calibrants: pd.DataFrame, candidates: pd.DataFrame, bundle: dict, alpha: float = 0.10,
                    method: str = "jackknife+", system_meta: dict | None = None, n_jobs: int = 1, chunk: int = 20000) -> pd.DataFrame:
    n = len(calibrants)
    if n < 1 / alpha - 1:
        raise ValueError(f"{n} calibrants give an infinite window at alpha = {alpha}: at least {int(np.ceil(1 / alpha - 1))} are needed")
    s_cal, C_cal = _features(calibrants["smiles"], bundle, n_jobs)
    y_cal = calibrants["rt"].to_numpy(float)
    loo = [MonotoneMap().fit(np.delete(s_cal, i), np.delete(y_cal, i)) for i in range(n)]
    mu_cal = np.array([loo[i].predict(s_cal[i:i + 1])[0] for i in range(n)])
    res = np.abs(y_cal - mu_cal)
    full = MonotoneMap().fit(s_cal, y_cal)
    if method == "normalised":
        if system_meta is None:
            raise ValueError("method='normalised' needs system_meta (chromatographic metadata of the target system)")
        from .evaluation.system_eval import _sigma_matrix, _sigma_cal_loo
        sysf = encode_system(system_meta); sig_model = bundle["difficulty_model"]
        sig_cal = _sigma_cal_loo(sig_model, C_cal.reset_index(drop=True), s_cal, y_cal, sysf)
        R = res / sig_cal
    elif method == "jackknife+":
        R = res
    else:
        raise ValueError("method must be 'jackknife+' or 'normalised'")
    out = []
    for start in range(0, len(candidates), chunk):  # chunks keep memory bounded for large candidate lists
        cand = candidates.iloc[start:start + chunk].reset_index(drop=True)
        s_q, C_q = _features(cand["smiles"], bundle, n_jobs)
        mu = np.vstack([m.predict(s_q) for m in loo])
        sig = _sigma_matrix(sig_model, C_q, s_q, s_cal, y_cal, sysf, loo=True) if method == "normalised" else None
        lo, hi = jackknife_plus(mu, R, sig, alpha)
        o = cand.copy(); o["rt_pred"] = full.predict(s_q); o["lo"] = lo; o["hi"] = hi
        if "rt_obs" in o:
            o["keep"] = (o.rt_obs >= o.lo) & (o.rt_obs <= o.hi)
        out.append(o)
    return pd.concat(out, ignore_index=True)
