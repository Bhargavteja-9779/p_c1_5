"""Selection-aware system-level calibration (SASC; EXPERIMENT_PLAN Addendum D).

A known calibrant-selection rule is replayed on other chromatographic systems (exchangeable units). For each system the
coverage of (normalised) jackknife+ windows whose score terms σ_{-i}(x)R_i are multiplied by c is computed on a grid of c;
c_ρ is the smallest c for which the AVERAGE coverage over training systems reaches 1 − α.
"""
from __future__ import annotations
import numpy as np
from ..models.mapping import MonotoneMap
from .intervals import k_upper, k_lower

GRID = np.round(np.arange(0.5, 6.0001, 0.02), 4)


def loo_quantities(s_cal, y_cal, s_te, sig_model_fn=None):
    """LOO maps; returns mu_loo_te (n,m), R (n,), sig_loo_te (n,m) or None. sig_model_fn(s_query, s_cal, y_cal) -> (sl_te, sl_cal)."""
    n = len(s_cal)
    loo = [MonotoneMap().fit(np.delete(s_cal, i), np.delete(y_cal, i)) for i in range(n)]
    mu_c = np.array([loo[i].predict(s_cal[i:i + 1])[0] for i in range(n)])
    mu_t = np.vstack([m.predict(s_te) for m in loo])
    res = np.abs(y_cal - mu_c)
    if sig_model_fn is None:
        return mu_t, res, None
    sl_te, sl_cal = sig_model_fn(s_te, s_cal, y_cal)
    return mu_t, res / sl_cal, sl_te


def windows(mu_t, R, sl_te, alpha, c=1.0):
    n = mu_t.shape[0]; ku, kl = k_upper(n, alpha), k_lower(n, alpha)
    S = (1.0 if sl_te is None else sl_te) * R[:, None] * c
    hi = np.full(mu_t.shape[1], np.inf) if ku > n else np.partition(mu_t + S, ku - 1, axis=0)[ku - 1]
    lo = np.full(mu_t.shape[1], -np.inf) if kl < 1 else np.partition(mu_t - S, kl - 1, axis=0)[kl - 1]
    return lo, hi


def coverage_curve(mu_t, R, sl_te, y_te, alpha, grid=GRID):
    out = np.empty(len(grid))
    for g, c in enumerate(grid):
        lo, hi = windows(mu_t, R, sl_te, alpha, c)
        out[g] = np.mean((lo <= y_te) & (y_te <= hi))
    return out


def choose_c(curves: np.ndarray, alpha: float, grid=GRID) -> float:
    """curves: (n_systems, len(grid)) system-level coverage curves (draw-averaged). Smallest c with mean coverage >= 1-α."""
    mean = np.nanmean(curves, axis=0)
    ok = np.flatnonzero(mean >= 1 - alpha)
    return float(grid[ok[0]]) if len(ok) else float(grid[-1])
