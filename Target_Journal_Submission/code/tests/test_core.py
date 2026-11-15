"""Unit tests: run with `python -m pytest -q tests` from the code/ directory."""
import os, sys, math
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from conformrt.evaluation.intervals import jackknife_plus, split_conformal, k_upper, k_lower, winkler
from conformrt.models.mapping import MonotoneMap, PSplineMap, LinearMap
from conformrt.models.difficulty import geometry, geometry_loo


def _sim(n, m, rng, hetero=True):
    s = rng.uniform(0, 10, n + m)
    sd = 0.2 + (0.3 * s if hetero else 0.0)
    y = 1 + 2 * s + 0.08 * s ** 2 + rng.normal(0, sd)
    return s, y, sd


def test_quantile_indices():
    assert k_upper(19, 0.1) == 18 and k_upper(8, 0.1) == 9  # 9 > 8 -> infinite interval
    assert k_lower(19, 0.1) == 2 and k_lower(8, 0.1) == 0


def test_jackknife_plus_infinite_when_too_few():
    mu = np.zeros((8, 5)); lo, hi = jackknife_plus(mu, np.ones(8), None, 0.1)
    assert np.isinf(hi).all() and np.isinf(lo).all()


def test_monotone_map_is_monotone_and_interpolates():
    rng = np.random.default_rng(1); s = np.sort(rng.uniform(0, 5, 15)); y = np.exp(0.4 * s) + rng.normal(0, 0.05, 15)
    m = MonotoneMap().fit(s, y); g = np.linspace(-2, 7, 400); p = m.predict(g)
    assert np.all(np.diff(p) >= -1e-9)


def test_geometry_loo_matches_bruteforce():
    rng = np.random.default_rng(2); s_cal = rng.uniform(0, 1, 12); y_cal = rng.uniform(0, 10, 12); q = rng.uniform(-0.2, 1.2, 30)
    G = geometry_loo(q, s_cal, y_cal)
    for i in range(12):
        keep = np.arange(12) != i
        g = geometry(q, s_cal[keep], y_cal[keep])
        for k in ["pos", "extrap", "gap"]:
            assert np.allclose(G[k][i], g[k])
        assert np.isclose(G["cal_span"][i], g["cal_span"])


def _coverage(method, reps=300, n=20, m=200, alpha=0.1, seed=0):
    rng = np.random.default_rng(seed); cov = []
    for r in range(reps):
        s, y, sd = _sim(n, m, rng)
        sc, yc, st, yt = s[:n], y[:n], s[n:], y[n:]
        loo = [MonotoneMap().fit(np.delete(sc, i), np.delete(yc, i)) for i in range(n)]
        mu_c = np.array([loo[i].predict(sc[i:i + 1])[0] for i in range(n)]); mu_t = np.vstack([l.predict(st) for l in loo])
        if method == "jk+":
            lo, hi = jackknife_plus(mu_t, np.abs(yc - mu_c), None, alpha)
        else:  # normalised by a FIXED known scale function (true sd shape) -> guarantee must hold
            sig_c = 0.2 + 0.3 * sc; sig_t = np.broadcast_to(0.2 + 0.3 * st, mu_t.shape)
            lo, hi = jackknife_plus(mu_t, np.abs(yc - mu_c) / sig_c, sig_t, alpha)
        cov.append(np.mean((lo <= yt) & (yt <= hi)))
    return float(np.mean(cov))


def test_jackknife_plus_coverage_exchangeable():
    c = _coverage("jk+"); assert c >= 0.87, c          # theory: >= 1-2a = 0.8 worst case; ~1-a in practice


def test_normalised_jackknife_plus_coverage_exchangeable():
    c = _coverage("njk+"); assert c >= 0.87, c


def test_split_conformal_exact():
    rng = np.random.default_rng(3); cov = []
    for r in range(2000):
        res = rng.normal(size=19); yt = rng.normal(size=1)
        lo, hi = split_conformal(np.zeros(1), res, 0.1); cov.append((lo[0] <= yt[0] <= hi[0]))
    assert abs(np.mean(cov) - 0.9) < 0.02


def test_winkler_penalises_misses():
    assert np.isclose(winkler([0], [1], [2], 0.1)[0], 1 + 20 * 1)


def test_pspline_interval_contains_most():
    rng = np.random.default_rng(4); s = rng.uniform(0, 10, 200); y = 3 + s + rng.normal(0, 0.5, 200)
    m = PSplineMap().fit(s, y); lo, hi = m.predict_interval(s, 0.1)
    assert 0.85 < np.mean((lo <= y) & (y <= hi)) < 0.95
