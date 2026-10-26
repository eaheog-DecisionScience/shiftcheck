"""Marginal-vs-dependence discrepancy decomposition. Layer 2 (core, uncertainty).

Four distributions are built from the two samples (A = source, B = target),
using within-sample Hazen plotting positions U = (rank - 0.5)/n and per-feature
Hazen quantile functions Q:

    A = (source marginals, source copula)      observed
    B = (target marginals, target copula)      observed
    M = Q_B(U_A) = (target marginals, source copula)   counterfactual
    C = Q_A(U_B) = (source marginals, target copula)   counterfactual

Cross-fitting. A counterfactual that reuses A's coordinate values (C) is not
independent of A, so pairwise distances between them run short and the
U-statistic energy distance is biased by about -d/n (the first version had
exactly this bias, found by the benchmark suite). Each sample is therefore
split into two halves; counterfactuals are built from one half of each sample
and compared only against the *other* halves, then the roles are swapped and
the two folds averaged. The bootstrap resamples within each half so that
duplicates never carry atoms across folds.

With D the squared energy distance (U-statistic) on features standardized by
the source's robust scale, the two orders of attribution give

    marginal   = [ D(A,M) + D(C,B) ] / 2     (marginals change first / second)
    dependence = [ D(A,C) + D(M,B) ] / 2     (copula changes first / second)
    interaction = D(A,B) - marginal - dependence

Averaging both orders is a two-player Shapley attribution. Energy distance is
not additive, so the interaction term is reported, not hidden. The share
marginal / (marginal + dependence) uses the clamped-at-zero contributions.
Total discrepancy is tested by label permutation; contributions and the share
by a percentile bootstrap that resamples source and target separately.
"""
import numpy as np
from scipy.spatial.distance import cdist
from scipy.stats import rankdata
from shiftcheck.core.scale import robust_scale
from shiftcheck.core.types import Estimate

UNIT = "squared energy distance, features standardized by source robust scale"


def _self(X):
    n = len(X)
    return cdist(X, X).sum() / (n * (n - 1))


def _cross(X, Y):
    return cdist(X, Y).mean()


def _plotting_positions(X):
    """Hazen positions (rank - 0.5) / n, so Q_X(U_X) reproduces X exactly."""
    n = len(X)
    return np.column_stack([(rankdata(X[:, j]) - 0.5) / n for j in range(X.shape[1])])


def _requantile(U, ref):
    return np.column_stack([np.quantile(ref[:, j], U[:, j], method="hazen") for j in range(ref.shape[1])])


def _fold(A_cf, B_cf, A_ev, B_ev):
    """Counterfactuals from (A_cf, B_cf); distances evaluated against (A_ev, B_ev)."""
    M = _requantile(_plotting_positions(A_cf), B_cf)
    C = _requantile(_plotting_positions(B_cf), A_cf)
    s = {"A": _self(A_ev), "B": _self(B_ev), "M": _self(M), "C": _self(C)}
    X = {"A": A_ev, "B": B_ev, "M": M, "C": C}

    def e(u, v):
        return 2 * _cross(X[u], X[v]) - s[u] - s[v]

    return np.array([e("A", "B"),
                     0.5 * (e("A", "M") + e("C", "B")),
                     0.5 * (e("A", "C") + e("M", "B"))])


def _components(A1, A2, B1, B2):
    total, marg, dep = 0.5 * (_fold(A1, B1, A2, B2) + _fold(A2, B2, A1, B1))
    pm, pd = max(marg, 0.0), max(dep, 0.0)
    share = pm / (pm + pd) if pm + pd > 0 else np.nan
    return np.array([total, marg, dep, total - marg - dep, share])


def _bootstrap(halves, n_boot, alpha, rng):
    obs = _components(*halves)
    draws = np.empty((n_boot, obs.size))
    for b in range(n_boot):
        draws[b] = _components(*[h[rng.integers(0, len(h), len(h))] for h in halves])
    lo, hi = np.nanquantile(draws, [alpha / 2, 1 - alpha / 2], axis=0)
    return [Estimate(float(o), float(l), float(h), "percentile bootstrap within cross-fit halves", n_boot)
            for o, l, h in zip(obs, lo, hi)]


def _permutation_total(S, T, n_perm, rng):
    X = np.vstack([S, T])
    D = cdist(X, X)
    n, m = len(S), len(T)

    def energy(a, b):
        saa = D[np.ix_(a, a)].sum() / (n * (n - 1))
        sbb = D[np.ix_(b, b)].sum() / (m * (m - 1))
        return 2 * D[np.ix_(a, b)].mean() - saa - sbb

    idx = np.arange(n + m)
    obs = energy(idx[:n], idx[n:])
    null = np.empty(n_perm)
    for k in range(n_perm):
        p = rng.permutation(n + m)
        null[k] = energy(p[:n], p[n:])
    return obs, float((1 + np.sum(null >= obs)) / (n_perm + 1))


def analyze_decomposition(S, T, names, alpha, n_boot, n_perm, rng, max_n=800):
    d = S.shape[1]
    if d < 2:
        return {"applicable": False, "reason": "need at least 2 features to separate dependence"}
    subsampled = None
    if len(S) > max_n or len(T) > max_n:
        subsampled = {"source": min(len(S), max_n), "target": min(len(T), max_n)}
        S = S[rng.choice(len(S), min(len(S), max_n), replace=False)]
        T = T[rng.choice(len(T), min(len(T), max_n), replace=False)]
    center = np.median(S, axis=0)
    scale = np.array([robust_scale(S[:, j]) for j in range(d)])
    A, B = (S - center) / scale, (T - center) / scale

    _, p_total = _permutation_total(A, B, n_perm, rng)
    ia, ib = rng.permutation(len(A)), rng.permutation(len(B))
    halves = (A[ia[: len(A) // 2]], A[ia[len(A) // 2:]], B[ib[: len(B) // 2]], B[ib[len(B) // 2:]])
    total, marg, dep, inter, share = _bootstrap(halves, n_boot, alpha, rng)
    if p_total >= alpha:
        verdict = "not_detected"
    elif share.lower > 0.5:
        verdict = "marginal_dominated"
    elif share.upper < 0.5:
        verdict = "dependence_dominated"
    else:
        verdict = "mixed_or_unclear"
    return {"applicable": True, "unit": UNIT, "verdict": verdict, "total_p": p_total,
            "total": total.to_dict(), "marginal": marg.to_dict(), "dependence": dep.to_dict(),
            "interaction": inter.to_dict(), "marginal_share": share.to_dict(),
            "subsampled_to": subsampled, "n_perm": n_perm}
