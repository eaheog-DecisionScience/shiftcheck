"""Dependence-structure diagnostics for continuous and categorical data. Layer 2.

Continuous features are assessed on pseudo-observations (ranks / (n + 1)), ensuring
marginal shifts do not leak into dependence comparisons.
Categorical features are assessed via Cramer's V contingency associations.
Supports IID, paired data (sign-flip permutation), and time series (circular shift).
"""
import numpy as np
from scipy.stats import chi2_contingency
from shiftcheck.core.ranks import pseudo_obs
from shiftcheck.uncertainty.resampling import permutation_null, p_value


def _ecdf_at(U, grid, chunk=64):
    out = np.empty(len(grid))
    for k in range(0, len(grid), chunk):
        g = grid[k:k + chunk]
        out[k:k + chunk] = np.all(U[:, None, :] <= g[None, :, :], axis=2).mean(axis=0)
    return out


def _cramers_v(x, y):
    """Compute Cramer's V association between two categorical series."""
    cats_x, cats_y = sorted(set(x)), sorted(set(y))
    r, k = len(cats_x), len(cats_y)
    if r < 2 or k < 2:
        return 0.0
    table = np.zeros((r, k), dtype=int)
    map_x = {c: idx for idx, c in enumerate(cats_x)}
    map_y = {c: idx for idx, c in enumerate(cats_y)}
    for xi, yi in zip(x, y):
        table[map_x[xi], map_y[yi]] += 1
    try:
        chi2 = chi2_contingency(table)[0]
        n = len(x)
        min_dim = min(r - 1, k - 1)
        return float(np.sqrt(chi2 / (n * min_dim))) if min_dim > 0 and n > 0 else 0.0
    except Exception:
        return 0.0


def analyze_dependence(S, T, names, alpha, n_perm, rng, n_grid=256,
                       paired=False, block_length=None,
                       S_cat=None, T_cat=None, cat_names=None):
    """Continuous copula and categorical association dependence diagnostics."""
    d = S.shape[1] if S is not None else 0
    cat_pairs = []
    if S_cat is not None and S_cat.shape[1] >= 2 and cat_names is not None:
        iu_cat = np.triu_indices(S_cat.shape[1], 1)
        for i, j in zip(*iu_cat):
            v_s = _cramers_v(S_cat[:, i].astype(str), S_cat[:, j].astype(str))
            v_t = _cramers_v(T_cat[:, i].astype(str), T_cat[:, j].astype(str))
            cat_pairs.append({
                "pair": [cat_names[i], cat_names[j]],
                "cramers_v_source": float(v_s),
                "cramers_v_target": float(v_t),
                "difference": float(v_t - v_s),
            })

    if d < 2:
        return {
            "applicable": False,
            "reason": "need at least 2 continuous features for copula test",
            "categorical_pairs": cat_pairs,
        }

    grid = rng.random((n_grid, d))
    iu = np.triu_indices(d, 1)

    def stat(A, B):
        UA, UB = pseudo_obs(A), pseudo_obs(B)
        cvm = np.mean((_ecdf_at(UA, grid) - _ecdf_at(UB, grid)) ** 2)
        diff = np.corrcoef(UB.T)[iu] - np.corrcoef(UA.T)[iu]
        return np.concatenate([[cvm], diff])

    obs, null = permutation_null(stat, S, T, n_perm, rng, paired=paired, block_length=block_length)
    p_global = p_value(obs[0], null[:, 0])
    max_null = np.max(np.abs(null[:, 1:]), axis=1)
    US, UT = pseudo_obs(S), pseudo_obs(T)
    rs, rt = np.corrcoef(US.T), np.corrcoef(UT.T)
    pairs = []
    for k, (i, j) in enumerate(zip(*iu)):
        adj = float((1 + np.sum(max_null >= abs(obs[1 + k]))) / (n_perm + 1))
        pairs.append({
            "pair": [names[i], names[j]],
            "spearman_source": float(rs[i, j]),
            "spearman_target": float(rt[i, j]),
            "difference": float(obs[1 + k]),
            "p_adjusted": adj,
            "detected": bool(adj < alpha),
        })

    return {
        "applicable": True,
        "copula_cvm_statistic": float(obs[0]),
        "copula_p": p_global,
        "detected": bool(p_global < alpha),
        "n_perm": n_perm,
        "n_grid": n_grid,
        "pairs": pairs,
        "spearman_difference_matrix": (rt - rs).tolist(),
        "features": list(names),
        "categorical_pairs": cat_pairs,
    }
