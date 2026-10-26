"""Synthetic shift generator with known ground truth. Layer 1 (depends on core only).

Source: Gaussian copula (equicorrelation `rho`) with Gamma marginals.
Target: a composition of the requested shifts, applied to an independent draw.
"""
import numpy as np
from scipy.stats import norm, gamma, t as student_t, lognorm

KINDS = {"affine", "monotone", "copula", "location"}
MARGINALS = {
    "gamma": lambda u, shape: gamma.ppf(u, a=shape),
    "normal": lambda u, shape: norm.ppf(u),
    "t3": lambda u, shape: student_t.ppf(u, 3),
    "lognormal": lambda u, shape: lognorm.ppf(u, s=0.75),
}
POSITIVE_ONLY = {"gamma", "lognormal"}  # required by the x**(1+strength) monotone shift


def _gaussian_copula(n, d, rho, rng):
    C = np.full((d, d), rho)
    np.fill_diagonal(C, 1.0)
    Z = rng.multivariate_normal(np.zeros(d), C, size=n)
    return norm.cdf(Z)


def make_shift_pair(n_source=1000, n_target=1000, d=3, kinds=(), strength=0.5,
                    rho=0.3, rho_target=0.8, shape=3.0, marginal="gamma", seed=None):
    """Return (source, target, truth).

    kinds: subset of {"affine", "monotone", "copula"}; empty means no shift.
      affine   : x -> a + b x with a = strength, b = 1 + strength
      monotone : x -> x ** (1 + strength)   (strictly increasing, nonlinear)
      copula   : dependence changes from `rho` to `rho_target`, marginals untouched
      location : x -> x + strength * s0, where s0 is the source's *population* robust
                 scale (IQR/1.349), so the true normalized Wasserstein distance equals
                 `strength` exactly
    `marginal` is one of gamma | normal | t3 | lognormal.
    """
    kinds = set(kinds)
    if not kinds <= KINDS:
        raise ValueError(f"unknown kinds: {kinds - KINDS}")
    if marginal not in MARGINALS:
        raise ValueError(f"unknown marginal: {marginal}")
    if "monotone" in kinds and marginal not in POSITIVE_ONLY:
        raise ValueError("'monotone' (x**(1+strength)) needs a positive-support marginal")
    ppf = MARGINALS[marginal]
    s0 = float((ppf(0.75, shape) - ppf(0.25, shape)) / 1.349)
    rng = np.random.default_rng(seed)
    U_s = _gaussian_copula(n_source, d, rho, rng)
    U_t = _gaussian_copula(n_target, d, rho_target if "copula" in kinds else rho, rng)
    S = ppf(U_s, shape)
    T = ppf(U_t, shape)
    if "monotone" in kinds:
        T = T ** (1.0 + strength)
    if "affine" in kinds:
        T = strength + (1.0 + strength) * T
    if "location" in kinds:
        T = T + strength * s0
    truth = {"kinds": sorted(kinds), "strength": strength, "rho": rho, "marginal": marginal,
             "population_robust_scale": s0,
             "rho_target": rho_target if "copula" in kinds else rho,
             "marginal_shift": bool(kinds & {"affine", "monotone", "location"}),
             "nonlinear_marginal": "monotone" in kinds,
             "copula_shift": "copula" in kinds}
    return S, T, truth
