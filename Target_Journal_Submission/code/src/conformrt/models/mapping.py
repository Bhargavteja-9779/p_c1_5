"""Maps from a system-independent base score s to retention time on the target system, fitted on calibrants.

All mappings implement fit(s, y) -> self and predict(s) -> rt. Some also provide a parametric prediction interval
(`predict_interval`) used by the classical baselines.
"""
from __future__ import annotations
import numpy as np
from scipy.interpolate import PchipInterpolator, BSpline
from scipy import stats
from sklearn.isotonic import IsotonicRegression


class MonotoneMap:
    """Isotonic regression (increasing) followed by PCHIP interpolation through the isotonic block means; linear
    extrapolation outside the calibrant score range using the slope of a robust (Theil–Sen) line through calibrants."""

    def fit(self, s, y):
        s = np.asarray(s, float); y = np.asarray(y, float)
        o = np.argsort(s); s, y = s[o], y[o]
        iso = IsotonicRegression(increasing=True, out_of_bounds="clip").fit(s, y)
        yf = iso.predict(s)
        # collapse ties in s and constant blocks to strictly increasing knots
        ks, ky = [], []
        for si, yi in zip(s, yf):
            if ks and si - ks[-1] < 1e-12:
                ky[-1] = 0.5 * (ky[-1] + yi)
            else:
                ks.append(si); ky.append(yi)
        ks, ky = np.array(ks), np.array(ky)
        self.lo, self.hi = ks[0], ks[-1]
        self.f = PchipInterpolator(ks, ky, extrapolate=False) if len(ks) >= 2 else None
        self.c0 = ky[0] if len(ky) else float(np.mean(y))
        try:
            slope = stats.theilslopes(y, s)[0]
        except Exception:
            slope = 0.0
        self.slope = max(slope, 0.0)
        self.ylo, self.yhi = ky[0], ky[-1]
        return self

    def predict(self, s):
        s = np.asarray(s, float)
        if self.f is None:
            return np.full_like(s, self.c0)
        out = self.f(np.clip(s, self.lo, self.hi))
        out = np.where(s < self.lo, self.ylo + self.slope * (s - self.lo), out)
        out = np.where(s > self.hi, self.yhi + self.slope * (s - self.hi), out)
        return out


class LinearMap:
    def fit(self, s, y):
        s = np.asarray(s, float); y = np.asarray(y, float)
        X = np.c_[np.ones_like(s), s]
        self.beta, *_ = np.linalg.lstsq(X, y, rcond=None)
        self.XtX_inv = np.linalg.pinv(X.T @ X)
        r = y - X @ self.beta; self.dof = max(len(y) - 2, 1); self.sigma = np.sqrt((r @ r) / self.dof)
        return self

    def predict(self, s):
        s = np.asarray(s, float); return self.beta[0] + self.beta[1] * s

    def predict_interval(self, s, alpha):
        s = np.asarray(s, float); X = np.c_[np.ones_like(s), s]
        h = np.einsum("ij,jk,ik->i", X, self.XtX_inv, X)
        t = stats.t.ppf(1 - alpha / 2, self.dof); m = X @ self.beta; hw = t * self.sigma * np.sqrt(1 + h)
        return m - hw, m + hw


class PSplineMap:
    """Penalised cubic B-spline regression (Eilers & Marx P-spline; GAM-like) with smoothing parameter chosen by GCV,
    and the classical Gaussian prediction interval  yhat ± t_{1-α/2, n-edf} σ̂ sqrt(1 + h(x))."""

    def __init__(self, n_knots=8, lambdas=np.logspace(-4, 3, 15)):
        self.n_knots, self.lambdas = n_knots, lambdas

    def _basis(self, s):
        k = 3
        inner = np.linspace(self.lo, self.hi, self.n_knots)
        t = np.r_[[self.lo] * k, inner, [self.hi] * k]
        nb = len(t) - k - 1
        sc = np.clip(s, self.lo, self.hi)
        B = BSpline.design_matrix(sc, t, k).toarray() if hasattr(BSpline, "design_matrix") else None
        return B, nb

    def fit(self, s, y):
        s = np.asarray(s, float); y = np.asarray(y, float)
        self.lo, self.hi = s.min(), s.max() + 1e-9
        B, nb = self._basis(s)
        D = np.diff(np.eye(nb), 2, axis=0); P = D.T @ D
        best = None
        for lam in self.lambdas:
            A = B.T @ B + lam * P
            Ai = np.linalg.pinv(A)
            H = B @ Ai @ B.T; edf = np.trace(H)
            r = y - H @ y; n = len(y)
            gcv = n * (r @ r) / max(n - edf, 1e-3) ** 2
            if best is None or gcv < best[0]:
                best = (gcv, lam, Ai, edf, r)
        _, self.lam, self.Ai, self.edf, r = best
        self.coef = self.Ai @ B.T @ y
        self.dof = max(len(y) - self.edf, 1.0)
        self.sigma = np.sqrt((r @ r) / self.dof)
        # linear extrapolation slope from end segments (Theil-Sen on all calibrants)
        try:
            self.slope = max(stats.theilslopes(y, s)[0], 0.0)
        except Exception:
            self.slope = 0.0
        self.ylo = float(self._inside(np.array([self.lo]))[0]); self.yhi = float(self._inside(np.array([self.hi]))[0])
        return self

    def _inside(self, s):
        B, _ = self._basis(s); return B @ self.coef

    def predict(self, s):
        s = np.asarray(s, float); out = self._inside(s)
        out = np.where(s < self.lo, self.ylo + self.slope * (s - self.lo), out)
        out = np.where(s > self.hi, self.yhi + self.slope * (s - self.hi), out)
        return out

    def predict_interval(self, s, alpha):
        s = np.asarray(s, float); B, _ = self._basis(s)
        h = np.einsum("ij,jk,ik->i", B, self.Ai, B)
        m = self.predict(s); t = stats.t.ppf(1 - alpha / 2, self.dof)
        hw = t * self.sigma * np.sqrt(1 + h)
        return m - hw, m + hw


class GPMap:
    """Gaussian-process regression of RT on the base score (RBF + linear (DotProduct) + white noise kernel,
    hyper-parameters by marginal likelihood). Predictive Gaussian interval mean ± z σ_pred."""

    def fit(self, s, y):
        from sklearn.gaussian_process import GaussianProcessRegressor
        from sklearn.gaussian_process.kernels import RBF, DotProduct, WhiteKernel, ConstantKernel as C
        s = np.asarray(s, float).reshape(-1, 1); y = np.asarray(y, float)
        self.mu_s, self.sd_s = s.mean(), s.std() + 1e-9
        self.mu_y = y.mean(); self.sd_y = y.std() + 1e-9
        k = C(1.0, (1e-3, 1e3)) * RBF(1.0, (1e-2, 1e2)) + C(0.1, (1e-4, 1e2)) * DotProduct(1.0, (1e-3, 1e2)) + WhiteKernel(0.05, (1e-6, 1.0))
        self.gp = GaussianProcessRegressor(kernel=k, normalize_y=False, n_restarts_optimizer=1, random_state=0)
        self.gp.fit((s - self.mu_s) / self.sd_s, (y - self.mu_y) / self.sd_y)
        return self

    def predict(self, s):
        s = np.asarray(s, float).reshape(-1, 1)
        return self.gp.predict((s - self.mu_s) / self.sd_s) * self.sd_y + self.mu_y

    def predict_interval(self, s, alpha):
        s = np.asarray(s, float).reshape(-1, 1)
        m, sd = self.gp.predict((s - self.mu_s) / self.sd_s, return_std=True)
        z = stats.norm.ppf(1 - alpha / 2)
        m = m * self.sd_y + self.mu_y; sd = sd * self.sd_y
        return m - z * sd, m + z * sd
