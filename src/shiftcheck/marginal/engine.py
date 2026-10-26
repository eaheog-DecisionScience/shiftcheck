"""Per-feature marginal diagnostics for continuous and categorical data. Layer 2 (core, uncertainty)."""
import numpy as np
from scipy.stats import ks_2samp, wasserstein_distance, chi2_contingency
from shiftcheck.uncertainty.resampling import bootstrap_estimates
from shiftcheck.uncertainty.multiplicity import adjust

QUANTILE_LEVELS = (0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95)


def _robust_scale(x):
    q75, q25 = np.quantile(x, [0.75, 0.25])
    s = (q75 - q25) / 1.349
    return s if s > 0 else (np.std(x) or 1.0)


def _lag1_autocorr(x):
    """Compute lag-1 serial autocorrelation."""
    if len(x) < 3:
        return 0.0
    x_centered = x - np.mean(x)
    var = np.var(x)
    if var <= 0:
        return 0.0
    return float(np.mean(x_centered[:-1] * x_centered[1:]) / var)


def analyze_marginals(S, T, names, alpha, n_boot, rng, multiplicity="holm",
                      paired=False, block_length=None):
    """Continuous marginal shift analysis supporting IID, paired, and time-series data."""
    d = S.shape[1]
    if d == 0:
        return {"features": [], "multiplicity": multiplicity, "applicable": False, "any_detected": False}

    ks_p, rows = [], []
    for j in range(d):
        s, t = S[:, j], T[:, j]
        scale = _robust_scale(s)
        ks = ks_2samp(s, t)
        ks_p.append(ks.pvalue)

        def w_norm(A, B, scale=scale):
            return wasserstein_distance(A[:, 0], B[:, 0]) / scale

        w = bootstrap_estimates(w_norm, S[:, [j]], T[:, [j]], n_boot, alpha, rng,
                                paired=paired, block_length=block_length)[0]
        q_shift = (np.quantile(t, QUANTILE_LEVELS) - np.quantile(s, QUANTILE_LEVELS)) / scale

        row = {
            "feature": names[j],
            "ks_statistic": float(ks.statistic),
            "ks_p": float(ks.pvalue),
            "wasserstein_normalized": w.to_dict(),
            "quantile_shift_in_source_sd": dict(zip(map(str, QUANTILE_LEVELS), map(float, q_shift))),
            "outside_source_range_frac": float(np.mean((t < s.min()) | (t > s.max()))),
        }

        # Paired data metrics
        if paired and len(s) == len(t):
            diff = t - s
            wilc_stat, wilc_p = None, None
            try:
                from scipy.stats import wilcoxon
                if np.all(diff == 0):
                    wilc_stat, wilc_p = 0.0, 1.0
                else:
                    res = wilcoxon(t, s)
                    wilc_stat, wilc_p = float(res.statistic), float(res.pvalue)
            except Exception:
                wilc_stat, wilc_p = None, None

            row["paired"] = {
                "mean_difference": float(np.mean(diff)),
                "median_difference": float(np.median(diff)),
                "positive_shift_frac": float(np.mean(diff > 0)),
                "wilcoxon_statistic": wilc_stat,
                "wilcoxon_p": wilc_p,
            }

        # Time series autocorrelation metrics
        if block_length is not None and block_length > 1:
            ac_s = _lag1_autocorr(s)
            ac_t = _lag1_autocorr(t)
            row["autocorrelation_lag1"] = {
                "source": ac_s,
                "target": ac_t,
                "difference": float(ac_t - ac_s),
            }

        rows.append(row)

    adj = adjust(ks_p, multiplicity)
    for r, a in zip(rows, adj):
        r["ks_p_adjusted"] = float(a)
        r["detected"] = bool(a < alpha)

    return {
        "features": rows,
        "multiplicity": multiplicity,
        "applicable": True,
        "any_detected": bool(any(r["detected"] for r in rows)),
    }


def analyze_categorical_marginals(S_cat, T_cat, names, alpha, n_boot, rng,
                                  paired=False, block_length=None, multiplicity="holm"):
    """Categorical marginal shift analysis: Chi-square of homogeneity and Total Variation Distance."""
    if S_cat is None or S_cat.shape[1] == 0:
        return {"features": [], "applicable": False, "any_detected": False}

    d_cat = S_cat.shape[1]
    chi2_p, rows = [], []
    for j in range(d_cat):
        s = S_cat[:, j].astype(str)
        t = T_cat[:, j].astype(str)
        cats = sorted(set(s) | set(t))
        s_counts = {c: int(np.sum(s == c)) for c in cats}
        t_counts = {c: int(np.sum(t == c)) for c in cats}
        table = np.array([[s_counts[c] for c in cats], [t_counts[c] for c in cats]])

        # Chi-squared test of homogeneity
        try:
            res = chi2_contingency(table)
            chi2_stat, pval = float(res.statistic), float(res.pvalue)
        except Exception:
            chi2_stat, pval = 0.0, 1.0
        chi2_p.append(pval)

        # Total Variation Distance (TVD) with bootstrap CI
        def tvd_stat(A, B, cats=cats):
            a_vals = A[:, 0].astype(str)
            b_vals = B[:, 0].astype(str)
            p_a = np.array([np.mean(a_vals == c) for c in cats])
            p_b = np.array([np.mean(b_vals == c) for c in cats])
            return 0.5 * np.sum(np.abs(p_a - p_b))

        tvd = bootstrap_estimates(tvd_stat, S_cat[:, [j]], T_cat[:, [j]], n_boot, alpha, rng,
                                  paired=paired, block_length=block_length)[0]

        shifts = {
            c: {
                "source_prop": float(s_counts[c] / len(s)),
                "target_prop": float(t_counts[c] / len(t)),
                "diff": float(t_counts[c] / len(t) - s_counts[c] / len(s)),
            }
            for c in cats
        }

        row = {
            "feature": names[j],
            "chi2_statistic": chi2_stat,
            "chi2_p": pval,
            "tvd": tvd.to_dict(),
            "categories": cats,
            "category_shifts": shifts,
        }

        # Paired categorical data
        if paired and len(s) == len(t):
            agreement = float(np.mean(s == t))
            row["paired"] = {
                "category_agreement_frac": agreement,
                "changed_category_frac": float(1.0 - agreement),
            }

        rows.append(row)

    adj = adjust(chi2_p, multiplicity) if chi2_p else []
    for r, a in zip(rows, adj):
        r["chi2_p_adjusted"] = float(a)
        r["detected"] = bool(a < alpha)

    return {
        "features": rows,
        "multiplicity": multiplicity,
        "applicable": True,
        "any_detected": bool(any(r["detected"] for r in rows)),
    }
