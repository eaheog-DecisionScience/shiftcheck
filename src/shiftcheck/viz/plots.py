"""Plots that consume report dictionaries only (no recomputation). Layer 1.

Takes `report.to_dict()` (or a DiagnosticReport via .to_dict()) so that
visuals can never disagree with the stored results.
"""
import numpy as np


def _plt():
    try:
        import matplotlib.pyplot as plt
    except ImportError as e:  # pragma: no cover
        raise ImportError("install shiftcheck[plots] for plotting") from e
    return plt


def _d(report):
    return report.to_dict() if hasattr(report, "to_dict") else report


def plot_quantile_map(report, feature=0, ax=None):
    """Q-Q curve with the robust affine fit: is the marginal change affine?"""
    plt = _plt()
    r = _d(report)["transform_compatibility"]["features"][feature]
    q = r["quantile_curve"]
    ax = ax or plt.subplots(figsize=(4.5, 4.2))[1]
    ax.plot(q["source"], q["target"], "o-", ms=3, label="empirical quantile map")
    ax.plot(q["source"], q["affine_fit"], "--", label="robust affine fit")
    g = r["nonlinearity_gap"]
    ax.set(xlabel="source quantile", ylabel="target quantile",
           title=f"{r['feature']}: gap {g['value']:.2f} [{g['lower']:.2f}, {g['upper']:.2f}]\n{r['verdict']}")
    ax.legend(frameon=False, fontsize=8)
    return ax


def plot_quantile_shift(report, ax=None):
    """Quantile shift (in source SDs) per feature: where in the distribution did it move?"""
    plt = _plt()
    feats = _d(report)["marginals"]["features"]
    ax = ax or plt.subplots(figsize=(5.5, 3.8))[1]
    for f in feats:
        qs = f["quantile_shift_in_source_sd"]
        ax.plot([float(k) for k in qs], list(qs.values()), "o-", ms=3, label=f["feature"])
    ax.axhline(0, color="grey", lw=0.8)
    ax.set(xlabel="quantile level p", ylabel="(Q_T(p) - Q_S(p)) / source sd", title="Quantile shift profile")
    ax.legend(frameon=False, fontsize=8)
    return ax


def plot_dependence_diff(report, ax=None):
    """Heatmap of Spearman(target) - Spearman(source) with detected pairs outlined."""
    plt = _plt()
    dep = _d(report)["dependence"]
    if not dep.get("applicable"):
        raise ValueError("dependence diagnostics not applicable (need >= 2 features)")
    M = np.array(dep["spearman_difference_matrix"])
    names = dep["features"]
    ax = ax or plt.subplots(figsize=(4.5, 4))[1]
    v = max(0.1, np.abs(M).max())
    im = ax.imshow(M, cmap="RdBu_r", vmin=-v, vmax=v)
    ax.figure.colorbar(im, ax=ax, label="Δ Spearman")
    ax.set_xticks(range(len(names)), names)
    ax.set_yticks(range(len(names)), names)
    idx = {n: i for i, n in enumerate(names)}
    for p in dep["pairs"]:
        i, j = idx[p["pair"][0]], idx[p["pair"][1]]
        for a, b in ((i, j), (j, i)):
            ax.text(b, a, f"{M[a, b]:+.2f}" + ("*" if p["detected"] else ""), ha="center", va="center", fontsize=8,
                    color="white" if abs(M[a, b]) > 0.6 * v else "black")
    ax.set_title(f"Dependence change (copula p={dep['copula_p']:.3g}); * = adj-p < alpha")
    return ax


def plot_decomposition(report, ax=None):
    """Marginal vs dependence vs interaction contributions to the total discrepancy, with CIs."""
    plt = _plt()
    dec = _d(report)["decomposition"]
    if not dec.get("applicable"):
        raise ValueError("decomposition not applicable (need >= 2 features)")
    keys = ["marginal", "dependence", "interaction", "total"]
    vals = np.array([dec[k]["value"] for k in keys])
    lo = vals - np.array([dec[k]["lower"] for k in keys])
    hi = np.array([dec[k]["upper"] for k in keys]) - vals
    ax = ax or plt.subplots(figsize=(4.8, 3.8))[1]
    ax.bar(keys, vals, yerr=[np.maximum(lo, 0), np.maximum(hi, 0)], capsize=4,
           color=["#1f77b4", "#d62728", "#7f7f7f", "#2ca02c"])
    ax.axhline(0, color="grey", lw=0.8)
    sh = dec["marginal_share"]
    sub = (f"\nmarginal share {sh['value']:.2f} [{sh['lower']:.2f}, {sh['upper']:.2f}]"
           if dec["verdict"] != "not_detected" else "")
    ax.set(ylabel="squared energy distance", title=f"{dec['verdict']}{sub}")
    return ax


def plot_power_curves(report, ax=None):
    """Detection power vs location-shift size per feature: how small a shift could we have seen?"""
    plt = _plt()
    pw = _d(report)["power"]
    if not pw:
        raise ValueError("run diagnose(..., power=True) first")
    m = pw["marginal"]
    ax = ax or plt.subplots(figsize=(5.2, 3.8))[1]
    for f in m["features"]:
        ax.plot(m["deltas"], f["power"], "o-", ms=3, label=f"feature {f['feature_index']}")
    ax.axhline(m["target_power"], color="grey", ls="--", lw=0.8)
    ax.set(xlabel="location shift (source robust SD)", ylabel="detection power", ylim=(0, 1.02),
           title="Marginal detectability at these sample sizes")
    ax.legend(frameon=False, fontsize=8)
    return ax
