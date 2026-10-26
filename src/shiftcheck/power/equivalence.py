"""Equivalence checks: can "not detected" be upgraded to "negligibly different"? Layer 2.

A non-significant test is not evidence of no shift. A property is called
`equivalent_within_tolerance` only if a one-sided upper confidence bound for its
effect size lies below a user tolerance delta (the TOST logic: the effect is
shown to lie inside (-delta, delta)).

  marginal   : max over features of the upper bound of the normalized Wasserstein
               distance (source robust SDs)           < delta_marginal
  dependence : upper bound of max |dSpearman| over feature pairs < delta_dependence

Three constructions of the upper bound are computed from the *same* bootstrap
draws (source and target resampled separately), so they can be compared on
equal footing in the benchmark suite:

  conservative : (1 - alpha/2) percentile       (a (1-alpha) two-sided interval)
  tost         : (1 - alpha) percentile         (a (1-2*alpha) interval, the usual TOST level)
  basic        : 2*estimate - alpha percentile  (pivotal; corrects the first-order upward
                                                 bias of these non-negative distance statistics)

`interval` selects which one drives the decision; `upper_bounds_all` reports the
three. For these, `smallest_certifiable_delta` is the chosen upper bound.

`interval="boundary"` tests H0: effect >= delta directly (boundary-calibrated test,
intersection-union over features / pairs, certify iff max p < alpha):

  marginal   : per feature, simulate source/target pairs *at* the boundary (both drawn
               from the observed source's marginal, the target shifted by exactly
               delta robust SDs, so the true Wasserstein distance equals delta) and
               ask how often the statistic is as small as observed.
  dependence : per pair, simulate Gaussian-copula pairs whose Spearman correlations
               differ by exactly delta (direction taken from the observed difference)
               and do the same with |dSpearman|.

Validity at the boundary is by construction for location-type marginal shifts and
Gaussian-type dependence changes; robustness to other shift types is checked in the
benchmark (shiftcheck.benchmarks.boundary). `smallest_certifiable_delta` is the smallest
delta at which the test would certify (found by bisection with common random numbers).
"""
import numpy as np
from scipy.stats import wasserstein_distance, rankdata
from shiftcheck.core.ranks import pseudo_obs
from shiftcheck.core.scale import robust_scale

METHODS = ("conservative", "tost", "basic")


def upper_bounds(obs, draws, alpha):
    """Upper bounds per output column from bootstrap `draws` (B x k) around estimate `obs` (k,)."""
    obs, draws = np.atleast_1d(obs), np.atleast_2d(draws)
    return {"conservative": np.nanquantile(draws, 1 - alpha / 2, axis=0),
            "tost": np.nanquantile(draws, 1 - alpha, axis=0),
            "basic": 2 * obs - np.nanquantile(draws, alpha, axis=0)}


def equivalence_bounds(S, T, n_boot, alpha, rng, parts=("marginal", "dependence")):
    """Estimates and all three upper bounds for the requested effect sizes."""
    d = S.shape[1]
    scale = np.array([robust_scale(S[:, j]) for j in range(d)])
    iu = np.triu_indices(d, 1) if d > 1 else None
    do_m, do_d = "marginal" in parts, "dependence" in parts and d >= 2

    def w_stat(A, B):
        return np.array([wasserstein_distance(A[:, j], B[:, j]) for j in range(d)]) / scale

    def dep_stat(A, B):
        return np.array([np.max(np.abs(np.corrcoef(pseudo_obs(B).T)[iu] - np.corrcoef(pseudo_obs(A).T)[iu]))])

    w_obs = w_stat(S, T) if do_m else None
    dep_obs = dep_stat(S, T) if do_d else None
    w_draws = np.empty((n_boot, d)) if do_m else None
    dep_draws = np.empty((n_boot, 1)) if do_d else None
    for b in range(n_boot):
        Sb = S[rng.integers(0, len(S), len(S))]
        Tb = T[rng.integers(0, len(T), len(T))]
        if do_m:
            w_draws[b] = w_stat(Sb, Tb)
        if do_d:
            dep_draws[b] = dep_stat(Sb, Tb)
    out = {}
    if do_m:
        wb = upper_bounds(w_obs, w_draws, alpha)
        out["marginal"] = {"estimates": w_obs.tolist(), "upper": {m: float(np.max(wb[m])) for m in METHODS}}
    if do_d:
        db = upper_bounds(dep_obs, dep_draws, alpha)
        out["dependence"] = {"estimate": float(dep_obs[0]), "upper": {m: float(db[m][0]) for m in METHODS}}
    return out


def _grid_quantiles(x_sorted, m):
    """Sorted sample -> m equally spaced quantile-grid values (identity if len == m)."""
    k = len(x_sorted)
    if k == m:
        return x_sorted
    return np.interp((np.arange(m) + 0.5) / m * k - 0.5, np.arange(k), x_sorted)


def _bisect_smallest(p_fn, alpha, hi, iters=14):
    """Smallest delta in [0, hi] with p_fn(delta) < alpha (p_fn non-increasing); nan if none."""
    if p_fn(0.0) < alpha:
        return 0.0
    if p_fn(hi) >= alpha:
        return float("nan")
    lo = 0.0
    for _ in range(iters):
        mid = 0.5 * (lo + hi)
        lo, hi = (lo, mid) if p_fn(mid) < alpha else (mid, hi)
    return float(hi)


def _marginal_boundary(S, T, delta, alpha, rng, n_sim, find_smallest):
    d, m = S.shape[1], min(len(S), len(T))
    per_feature, deltas = [], []
    for j in range(d):
        sc = robust_scale(S[:, j])
        obs = np.mean(np.abs(_grid_quantiles(np.sort(S[:, j]), m) - _grid_quantiles(np.sort(T[:, j]), m))) / sc
        col = S[:, j]
        D = np.empty((n_sim, m))
        for b in range(n_sim):
            a = np.sort(col[rng.integers(0, len(col), len(S))])
            c = np.sort(col[rng.integers(0, len(col), len(T))])
            D[b] = (_grid_quantiles(a, m) - _grid_quantiles(c, m)) / sc

        def p_at(dl, D=D, obs=obs):
            return (1 + np.sum(np.mean(np.abs(D - dl), axis=1) <= obs)) / (n_sim + 1)

        per_feature.append(float(p_at(delta)))
        if find_smallest:
            deltas.append(_bisect_smallest(p_at, alpha, hi=max(3.0, 3 * obs)))
    out = {"p_per_feature": per_feature, "p": float(max(per_feature))}
    if find_smallest:
        out["smallest_certifiable_delta"] = float(np.max(deltas)) if not np.any(np.isnan(deltas)) else float("nan")
    return out


def _sp(x, y):
    """Row-wise Spearman correlation of (B, n) arrays."""
    cx = rankdata(x, axis=1)
    cy = rankdata(y, axis=1)
    cx = cx - cx.mean(axis=1, keepdims=True)
    cy = cy - cy.mean(axis=1, keepdims=True)
    return (cx * cy).sum(1) / np.sqrt((cx ** 2).sum(1) * (cy ** 2).sum(1))


def _dependence_boundary(S, T, delta, alpha, rng, n_sim, find_smallest, n_sim_search=149):
    d = S.shape[1]
    US, UT = pseudo_obs(S), pseudo_obs(T)
    CS, CT = np.corrcoef(US.T), np.corrcoef(UT.T)
    ps, pairs = [], []
    deltas = []
    for i, j in zip(*np.triu_indices(d, 1)):
        rho_s, diff = CS[i, j], CT[i, j] - CS[i, j]
        obs, sgn = abs(diff), (1.0 if diff >= 0 else -1.0)
        if abs(rho_s + sgn * delta) > 0.98:
            sgn = -sgn

        def build(B):
            Z1, Z2 = rng.standard_normal((2, B, len(S))), rng.standard_normal((2, B, len(T)))
            r_s = 2 * np.sin(np.pi * rho_s / 6)
            sp_s = _sp(Z1[0], r_s * Z1[0] + np.sqrt(1 - r_s ** 2) * Z1[1])
            return Z2, sp_s

        def p_from(Z2, sp_s, dl):
            rho_t = np.clip(rho_s + sgn * dl, -0.999, 0.999)
            r_t = 2 * np.sin(np.pi * rho_t / 6)
            sp_t = _sp(Z2[0], r_t * Z2[0] + np.sqrt(1 - r_t ** 2) * Z2[1])
            return (1 + np.sum(np.abs(sp_t - sp_s) <= obs)) / (len(sp_s) + 1)

        Z2, sp_s = build(n_sim)
        ps.append(float(p_from(Z2, sp_s, delta)))
        pairs.append((int(i), int(j)))
        if find_smallest:
            Z2b, sp_sb = build(n_sim_search)
            deltas.append(_bisect_smallest(lambda dl: p_from(Z2b, sp_sb, dl), alpha, hi=0.9, iters=10))
    out = {"p_per_pair": ps, "pairs": pairs, "p": float(max(ps))}
    if find_smallest:
        out["smallest_certifiable_delta"] = float(np.max(deltas)) if not np.any(np.isnan(deltas)) else float("nan")
    return out


def boundary_pvalues(S, T, delta_marginal, delta_dependence, alpha, rng, n_sim=499, find_smallest=True):
    """Boundary-calibrated equivalence p-values (intersection-union over features / pairs)."""
    out = {"marginal": _marginal_boundary(S, T, delta_marginal, alpha, rng, n_sim, find_smallest)}
    if S.shape[1] >= 2:
        out["dependence"] = _dependence_boundary(S, T, delta_dependence, alpha, rng, n_sim, find_smallest)
    return out


def analyze_equivalence(S, T, names, alpha, n_boot, rng, delta_marginal=0.1,
                        delta_dependence=0.1, interval="conservative", find_smallest=True):
    if S.shape[1] == 0:
        return {
            "interval": interval,
            "marginal": {"applicable": False, "equivalent": False, "delta": delta_marginal},
            "dependence": {"applicable": False, "equivalent": False, "delta": delta_dependence}
        }
    if interval == "boundary":
        bp = boundary_pvalues(S, T, delta_marginal, delta_dependence, alpha, rng, n_sim=n_boot,
                              find_smallest=find_smallest)
        out = {"interval": "boundary", "n_sim": n_boot,
               "marginal": {"delta": delta_marginal, "unit": "source robust SD", "p_value": bp["marginal"]["p"],
                            "p_per_feature": bp["marginal"]["p_per_feature"],
                            "smallest_certifiable_delta": bp["marginal"].get("smallest_certifiable_delta"),
                            "equivalent": bool(bp["marginal"]["p"] < alpha)}}
        if "dependence" not in bp:
            out["dependence"] = {"applicable": False}
        else:
            out["dependence"] = {"applicable": True, "delta": delta_dependence, "unit": "Spearman",
                                 "p_value": bp["dependence"]["p"], "p_per_pair": bp["dependence"]["p_per_pair"],
                                 "smallest_certifiable_delta": bp["dependence"].get("smallest_certifiable_delta"),
                                 "equivalent": bool(bp["dependence"]["p"] < alpha)}
        return out
    if interval not in METHODS:
        raise ValueError(f"interval must be one of {METHODS} or 'boundary'")
    b = equivalence_bounds(S, T, n_boot, alpha, rng)
    up = b["marginal"]["upper"]
    out = {"interval": interval, "n_boot": n_boot,
           "marginal": {"delta": delta_marginal, "unit": "source robust SD",
                        "max_wasserstein_estimate": float(max(b["marginal"]["estimates"])),
                        "upper_bounds_all": up, "smallest_certifiable_delta": up[interval],
                        "equivalent": bool(up[interval] < delta_marginal)}}
    if "dependence" not in b:
        out["dependence"] = {"applicable": False}
        return out
    dp = b["dependence"]["upper"]
    out["dependence"] = {"applicable": True, "delta": delta_dependence, "unit": "Spearman",
                         "max_abs_spearman_diff_estimate": b["dependence"]["estimate"],
                         "upper_bounds_all": dp, "smallest_certifiable_delta": dp[interval],
                         "equivalent": bool(dp[interval] < delta_dependence)}
    return out
