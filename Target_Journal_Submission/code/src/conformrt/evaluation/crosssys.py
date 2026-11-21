"""Cross-system conformal calibration (XSC) of retention-time windows.

A laboratory's calibrants are used only to fit the monotone map ĝ and a laboratory-specific error scale τ̂; the width
multiplier c is calibrated on hundreds of *other* chromatographic systems on which the laboratory's calibrant protocol
(number n, selection rule ρ) is replayed.  Systems — not compounds — are the exchangeable units.

Window:   C(x) = [ ĝ(s(x)) − c·τ̂·σ(x),  ĝ(s(x)) + c·τ̂·σ(x) ],
          τ̂ = sqrt( mean_i ( r_i / σ_{-i}(x_i) )² ),  r_i = t_i − ĝ_{-i}(s_i)   (leave-one-out residuals),
          σ ≡ 1 ('flat') or the meta-learned difficulty model ('σ').
Pivot:    u(x) = |t(x) − ĝ(s(x))| / (τ̂·σ(x)); per system s and calibrant draw d the coverage at c is
          cov_{s,d}(c) = mean_x 1{u(x) ≤ c}  and the smallest multiplier giving coverage ≥ 1−ε is c_{s,d}(ε).

Two calibration targets:
  * marginal ('avg'):  smallest c with  (1/S) Σ_s cov_s(c) ≥ 1 − α             (average retention over laboratories);
  * laboratory-level ('lab', (ε, δ)):  smallest c with  (1/S) Σ_s F_s(c) ≥ (1 − δ)(1 + 1/S),
    F_s(c) = fraction of draws on system s with c_{s,d}(ε) ≤ c.  For a new system exchangeable with the S training
    systems, P{ coverage of the laboratory's own windows ≥ 1 − ε } ≥ 1 − δ (exact when one draw per system is used:
    c = ⌈(1−δ)(S+1)⌉-th smallest c_s; the draw-averaged CDF is the smoothed version, see Supplementary S7).
Both work for any n ≥ 3 calibrants (no n ≥ 1/α − 1 requirement), because the quantile is taken over systems.
"""
from __future__ import annotations
import math
import numpy as np
from ..models.mapping import MonotoneMap

C_GRID = np.round(np.exp(np.linspace(np.log(0.05), np.log(40.0), 500)), 5)
EPS_LIST = (0.05, 0.10, 0.20)


def loo_residuals(s_cal, y_cal):
    n = len(s_cal)
    loo = [MonotoneMap().fit(np.delete(s_cal, i), np.delete(y_cal, i)) for i in range(n)]
    mu_c = np.array([loo[i].predict(s_cal[i:i + 1])[0] for i in range(n)])
    return y_cal - mu_c, loo


def scale_pivot(s_cal, y_cal, s_te, y_te=None, sig_te=None, sig_cal_loo=None, loo_res=None, shrink=0.0):
    """Returns (mu_te, tau, u_te or None). sig_te: σ(x) with the full calibrant geometry; sig_cal_loo: σ_{-i}(x_i)."""
    if loo_res is None:
        loo_res, _ = loo_residuals(s_cal, y_cal)
    z = loo_res / (1.0 if sig_cal_loo is None else sig_cal_loo)
    # shrink=k>0: empirical-Bayes shrinkage of the scale towards the meta-learned prediction (τ = 1 when σ is used),
    # i.e. k pseudo-calibrants with normalised residual 1 (stabilises τ̂ for very few calibrants)
    tau = float(np.sqrt((np.sum(z ** 2) + shrink) / (len(z) + shrink))) + 1e-6
    mu = MonotoneMap().fit(s_cal, y_cal).predict(s_te)
    scale = tau * (1.0 if sig_te is None else sig_te)
    u = None if y_te is None else np.abs(y_te - mu) / scale
    return mu, tau, u


def pivot_summary(u: np.ndarray, grid=C_GRID, eps_list=EPS_LIST):
    """Coverage curve on the grid and the per-draw required multipliers c(ε)."""
    us = np.sort(u)
    curve = np.searchsorted(us, grid, side="right") / len(us)
    m = len(us)
    req = {e: float(us[min(m - 1, math.ceil((1 - e) * m) - 1)]) for e in eps_list}
    return curve.astype(np.float32), req


def choose_avg(curves: np.ndarray, alpha: float, grid=C_GRID) -> float:
    """curves: (S, G) draw-averaged coverage curves per system."""
    S = curves.shape[0]
    target = 1 - alpha
    ok = np.flatnonzero(np.nanmean(curves, 0) >= target)
    return float(grid[ok[0]]) if len(ok) else float("inf")


def choose_lab(req_by_system: list, eps: float, delta: float, grid=C_GRID) -> float:
    """req_by_system: list over systems of arrays of c_{s,d}(ε) over draws. Draw-averaged system-level CDF."""
    S = len(req_by_system)
    F = np.mean([np.searchsorted(np.sort(r), grid, side="right") / len(r) for r in req_by_system], 0)
    target = (1 - delta) * (1 + 1 / S)
    if target > 1:
        return float("inf")
    ok = np.flatnonzero(F >= target)
    return float(grid[ok[0]]) if len(ok) else float("inf")


def choose_lab_exact(req_first_draw: np.ndarray, delta: float) -> float:
    """One draw per system: ⌈(1−δ)(S+1)⌉-th smallest required multiplier (exact finite-sample level 1−δ)."""
    S = len(req_first_draw); k = math.ceil((1 - delta) * (S + 1))
    return float("inf") if k > S else float(np.sort(req_first_draw)[k - 1])
