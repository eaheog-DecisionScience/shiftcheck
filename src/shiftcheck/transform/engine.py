"""Affine vs nonlinear-monotone quantile-map compatibility. Layer 2.

For feature j the empirical quantile map is q_T(p) vs q_S(p). We fit a robust
line (Theil-Sen) and measure the RMS deviation of the Q-Q curve from it, in
units of the target's robust scale. Bootstrap intervals are compared with a
user tolerance `tol` (an equivalence-style decision):

  upper CI < floor + tol  -> affine_compatible
  lower CI > floor + tol  -> nonlinear_monotone
  otherwise               -> inconclusive

`floor` is the gap expected from sampling noise alone under an exactly affine
map at these sample sizes (estimated by resampling the source and pushing an
independent resample through the fitted affine map); `tol` therefore bounds the
*excess* nonlinearity, in target standard deviations.

This characterizes the marginal relationship only; it does not identify the
mechanism that generated the target.
"""
import numpy as np
from scipy.stats import theilslopes
from shiftcheck.uncertainty.resampling import bootstrap_estimates

P_GRID = np.linspace(0.10, 0.90, 17)


def _qq(A, B):
    return np.quantile(A[:, 0], P_GRID), np.quantile(B[:, 0], P_GRID)


def _scale(qt):
    s = np.std(qt)
    return s if s > 0 else 1.0


def _fit(A, B):
    qs, qt = _qq(A, B)
    slope, intercept, _, _ = theilslopes(qt, qs)
    resid = qt - (intercept + slope * qs)
    gap = float(np.sqrt(np.mean(resid ** 2)) / _scale(qt))
    return np.array([gap, slope, intercept])


def _noise_floor(A, B, slope, intercept, n_rep, rng):
    """Median gap under an exactly affine target of the same sizes."""
    vals = np.empty(n_rep)
    for b in range(n_rep):
        As = A[rng.integers(0, len(A), len(A))]
        Ap = A[rng.integers(0, len(A), len(B))]
        vals[b] = _fit(As, intercept + slope * Ap)[0]
    return float(np.median(vals))


def analyze_transform(S, T, names, alpha, n_boot, rng, tol=0.1,
                      paired=False, block_length=None):
    out = []
    if S.shape[1] == 0:
        return {"features": [], "tolerance": tol}
    for j in range(S.shape[1]):
        A, B = S[:, [j]], T[:, [j]]
        qs, qt = _qq(A, B)
        if np.std(qs) == 0 or np.std(qt) == 0:
            out.append({"feature": names[j], "verdict": "not_applicable",
                        "reason": "constant feature"})
            continue
        gap, slope, intercept = bootstrap_estimates(_fit, A, B, n_boot, alpha, rng,
                                                    paired=paired, block_length=block_length)
        floor = _noise_floor(A, B, slope.value, intercept.value, max(50, n_boot // 3), rng)
        if slope.lower <= 0:
            verdict = "no_positive_affine_relation"
        elif gap.upper < floor + tol:
            verdict = "affine_compatible"
        elif gap.lower > floor + tol:
            verdict = "nonlinear_monotone"
        else:
            verdict = "inconclusive"
        s_, i_ = slope.value, intercept.value
        out.append({
            "feature": names[j], "verdict": verdict, "tolerance": tol, "noise_floor": floor,
            "nonlinearity_gap": gap.to_dict(), "affine_slope": slope.to_dict(),
            "affine_intercept": intercept.to_dict(),
            "quantile_curve": {"p": P_GRID.tolist(), "source": qs.tolist(), "target": qt.tolist(),
                               "affine_fit": (i_ + s_ * qs).tolist(), "scale": float(_scale(qt))},
        })
    return {"features": out, "tolerance": tol}
