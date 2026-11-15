"""Prediction-interval constructions used in the study.

Notation: calibrants (s_i, y_i), i=1..n, measured on the target system; s = base retention score; y = RT (min).
`fit_map` is a factory returning a fresh mapping object (fit/predict).

Jackknife+ (Barber, Candès, Ramdas & Tibshirani, Ann. Stat. 2021):
    C(x) = [ q^-_α{ μ_{-i}(x) − σ_{-i}(x) R_i },  q^+_α{ μ_{-i}(x) + σ_{-i}(x) R_i } ],
    R_i = |y_i − μ_{-i}(s_i)| / σ_{-i}(x_i),
q^+_α{v} = ⌈(1−α)(n+1)⌉-th smallest of v_1..v_n (+∞ if that index exceeds n); q^-_α{v} = ⌊α(n+1)⌋-th smallest
(−∞ if the index is 0). With σ ≡ 1 this is the original jackknife+. The normalised version keeps the guarantee
P(Y ∈ C(X)) ≥ 1 − 2α under exchangeability of calibrants and test compound because σ_{-i} is computed by the same
symmetric procedure from the leave-one-out calibrant set (the proof of Barber et al., Thm 1, only uses that the
leave-two-out score matrix is exchangeable; see Supplementary S2 of the manuscript).
"""
from __future__ import annotations
import math
import numpy as np
from scipy import stats


def k_upper(n: int, alpha: float) -> int:
    return math.ceil((1 - alpha) * (n + 1))


def k_lower(n: int, alpha: float) -> int:
    return math.floor(alpha * (n + 1))


def jackknife_plus(mu_loo_test: np.ndarray, R: np.ndarray, sig_loo_test: np.ndarray | None, alpha: float):
    """mu_loo_test, sig_loo_test: (n, m) arrays (LOO model i evaluated at test point j); R: (n,) scores."""
    n = mu_loo_test.shape[0]
    sig = 1.0 if sig_loo_test is None else sig_loo_test
    lo_v = mu_loo_test - sig * R[:, None]
    hi_v = mu_loo_test + sig * R[:, None]
    ku, kl = k_upper(n, alpha), k_lower(n, alpha)
    m = mu_loo_test.shape[1]
    hi = np.full(m, np.inf) if ku > n else np.partition(hi_v, ku - 1, axis=0)[ku - 1]
    lo = np.full(m, -np.inf) if kl < 1 else np.partition(lo_v, kl - 1, axis=0)[kl - 1]
    return lo, hi


def split_conformal(mu_test: np.ndarray, resid_cal: np.ndarray, alpha: float, sig_test=None, sig_cal=None):
    n = len(resid_cal)
    sc = np.abs(resid_cal) / (1.0 if sig_cal is None else sig_cal)
    k = k_upper(n, alpha)
    q = np.inf if k > n else np.sort(sc)[k - 1]
    w = q * (1.0 if sig_test is None else sig_test)
    return mu_test - w, mu_test + w


def gaussian_rmse(mu_test: np.ndarray, loo_resid: np.ndarray, alpha: float):
    z = stats.norm.ppf(1 - alpha / 2)
    hw = z * np.sqrt(np.mean(loo_resid ** 2))
    return mu_test - hw, mu_test + hw


def winkler(lo, hi, y, alpha):
    lo = np.asarray(lo); hi = np.asarray(hi); y = np.asarray(y)
    w = hi - lo
    w = w + (2 / alpha) * np.clip(lo - y, 0, None) + (2 / alpha) * np.clip(y - hi, 0, None)
    return w


def loo_fits(fit_map, s_cal, y_cal):
    """Return list of n mapping objects fitted without calibrant i, and the full fit."""
    n = len(s_cal)
    idx = np.arange(n)
    loo = [fit_map().fit(s_cal[idx != i], y_cal[idx != i]) for i in range(n)]
    full = fit_map().fit(s_cal, y_cal)
    return loo, full
