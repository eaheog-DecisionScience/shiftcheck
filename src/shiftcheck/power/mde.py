"""Simulation-based detectability: what shift sizes could these sample sizes reveal? Layer 2.

Marginal:   location shifts of delta source-SDs, detected by feature-wise KS at
            alpha/d (Bonferroni; conservative relative to Holm). Both samples
            are bootstrap resamples of the observed source.
Dependence: both samples are generated from the source's fitted Gaussian copula
            with empirical source marginals; the target's correlation matrix is
            blended toward all-ones, R_t = (1-w) R + w J. Detection uses the
            pairwise Spearman max-|diff| permutation screen (cheap) rather than
            the copula CvM test, so this is the power of the screen.
Effect size on the dependence side = mean |dSpearman| over pairs implied by R_t vs R.

This is power against *these* alternatives, a guide to resolution, not a
guarantee for other shift types.
"""
import numpy as np
from scipy.stats import ks_2samp, norm
from shiftcheck.core.ranks import pseudo_obs
from shiftcheck.core.scale import robust_scale
from shiftcheck.uncertainty.resampling import permutation_null

DELTAS = (0.02, 0.04, 0.06, 0.08, 0.1, 0.15, 0.2, 0.3, 0.5)
WEIGHTS = (0.05, 0.1, 0.2, 0.3, 0.5)


def _first_crossing(x, power, target):
    """Smallest effect with power >= target, by linear interpolation; None if never reached."""
    for k, pw in enumerate(power):
        if pw >= target:
            if k == 0:
                return float(x[0])
            x0, x1, p0, p1 = x[k - 1], x[k], power[k - 1], pw
            return float(x0 + (target - p0) * (x1 - x0) / (p1 - p0)) if p1 > p0 else float(x1)
    return None


def marginal_mde(S, n_source, n_target, alpha, rng, n_reps=60, deltas=DELTAS, target_power=0.8):
    d = S.shape[1]
    a = alpha / d
    feats = []
    for j in range(d):
        sc = robust_scale(S[:, j])
        pw = []
        for delta in deltas:
            rej = 0
            for _ in range(n_reps):
                x = S[rng.integers(0, len(S), n_source), j]
                y = S[rng.integers(0, len(S), n_target), j] + delta * sc
                rej += ks_2samp(x, y).pvalue < a
            pw.append(rej / n_reps)
        feats.append({"feature_index": j, "power": pw,
                      "mde_location_shift_sd": _first_crossing(list(deltas), pw, target_power)})
    return {"shift_family": "location shift, in source robust SD", "deltas": list(deltas),
            "target_power": target_power, "n_reps": n_reps, "features": feats}


def dependence_mde(S, n_source, n_target, alpha, rng, n_reps=30, n_perm=99,
                   weights=WEIGHTS, target_power=0.8):
    d = S.shape[1]
    if d < 2:
        return {"applicable": False}
    R = np.corrcoef(norm.ppf(pseudo_obs(S)).T)
    J = np.ones((d, d))
    iu = np.triu_indices(d, 1)
    rho = lambda r: 6.0 / np.pi * np.arcsin(r / 2.0)

    def gen(Rm, n):
        U = norm.cdf(rng.multivariate_normal(np.zeros(d), Rm, size=n, check_valid="ignore"))
        return np.column_stack([np.quantile(S[:, j], U[:, j]) for j in range(d)])

    def stat(A, B):
        return np.corrcoef(pseudo_obs(B).T)[iu] - np.corrcoef(pseudo_obs(A).T)[iu]

    effects, power = [], []
    for w in weights:
        Rt = (1 - w) * R + w * J
        effects.append(float(np.mean(np.abs(rho(Rt[iu]) - rho(R[iu])))))
        rej = 0
        for _ in range(n_reps):
            obs, null = permutation_null(stat, gen(R, n_source), gen(Rt, n_target), n_perm, rng)
            p = (1 + np.sum(np.max(np.abs(null), axis=1) >= np.max(np.abs(obs)))) / (n_perm + 1)
            rej += p < alpha
        power.append(rej / n_reps)
    return {"applicable": True, "shift_family": "correlation blended toward all-ones",
            "effect_unit": "mean |dSpearman| over pairs", "effects": effects, "power": power,
            "target_power": target_power, "n_reps": n_reps, "n_perm": n_perm,
            "mde_mean_abs_dspearman": _first_crossing(effects, power, target_power)}
