"""Unified report and the `diagnose` entry point. Layer 3 (all engines)."""
import json
from dataclasses import dataclass, field
import numpy as np
from shiftcheck.core.validate import prepare
from shiftcheck.marginal.engine import analyze_marginals, analyze_categorical_marginals
from shiftcheck.transform.engine import analyze_transform
from shiftcheck.dependence.engine import analyze_dependence
from shiftcheck.decomposition.engine import analyze_decomposition
from shiftcheck.power.equivalence import analyze_equivalence
from shiftcheck.power.mde import marginal_mde, dependence_mde

LIMITATIONS = [
    "Diagnostics concern the feature distribution P(X); they cannot establish whether P(Y|X) is stable.",
    "Compatibility with an affine/monotone quantile map does not identify the process that generated the target.",
    "'inconclusive' means the data cannot tell 'no shift' from 'small shift'; only 'equivalent_within_tolerance' supports 'negligible'.",
    "Decomposition shares depend on the chosen distance (energy) and on counterfactual constructions from empirical quantile maps.",
    "Assumes exchangeable observations unless time_series or paired mode is explicitly specified.",
]


def _jsonable(o):
    if isinstance(o, (np.floating, np.integer)):
        return o.item()
    if isinstance(o, np.ndarray):
        return o.tolist()
    raise TypeError(f"not serializable: {type(o)}")


def _fmt(x, nd=3):
    return "n/a" if x is None or (isinstance(x, float) and np.isnan(x)) else f"{x:.{nd}f}"


def _status(detected, equivalent):
    return "detected" if detected else ("equivalent_within_tolerance" if equivalent else "inconclusive")


@dataclass
class DiagnosticReport:
    alpha: float
    n_source: int
    n_target: int
    features: list
    marginals: dict
    transform_compatibility: dict
    dependence: dict
    decomposition: dict
    equivalence: dict
    categorical_marginals: dict = None
    paired: bool = False
    time_series: bool = False
    block_length: int = None
    power: dict = None
    notes: list = field(default_factory=list)
    limitations: list = field(default_factory=lambda: list(LIMITATIONS))

    def profile(self):
        """Independent per-property statuses (never a single exclusive label)."""
        verdicts = [f["verdict"] for f in self.transform_compatibility["features"]]
        if "nonlinear_monotone" in verdicts:
            tf = "nonlinear_monotone_detected"
        elif all(v == "affine_compatible" for v in verdicts) and verdicts:
            tf = "affine_compatible"
        else:
            tf = "inconclusive"
        dep, eq = self.dependence, self.equivalence
        dec = self.decomposition
        p = {
            "marginal_difference": _status(self.marginals["any_detected"], eq["marginal"].get("equivalent", False)),
            "marginal_form": tf,
            "dependence_difference": ("not_applicable" if not dep.get("applicable") else
                                      _status(dep["detected"], eq["dependence"].get("equivalent", False))),
            "shift_attribution": dec["verdict"] if dec.get("applicable") else "not_applicable",
            "outcome_mechanism": "not_assessable_from_features",
        }
        if self.categorical_marginals and self.categorical_marginals.get("applicable"):
            p["categorical_difference"] = "detected" if self.categorical_marginals["any_detected"] else "inconclusive"
        if self.paired:
            p["sampling_design"] = "paired"
        elif self.time_series:
            p["sampling_design"] = f"time_series(block_length={self.block_length})"
        return p

    def summary(self):
        p = self.profile()
        if self.paired:
            hdr = f"shiftcheck report  (paired n={self.n_source}, alpha={self.alpha})"
        elif self.time_series:
            hdr = f"shiftcheck report  (time series n_source={self.n_source}, n_target={self.n_target}, block_length={self.block_length}, alpha={self.alpha})"
        else:
            hdr = f"shiftcheck report  (source n={self.n_source}, target n={self.n_target}, alpha={self.alpha})"

        lines = [hdr, ""]
        lines += [f"  {k:24s} {v}" for k, v in p.items()]
        lines.append("")

        if self.marginals.get("applicable", True) and self.marginals["features"]:
            lines.append("  Continuous Marginals:")
            for f in self.marginals["features"]:
                w = f["wasserstein_normalized"]
                lines.append(f"  [{f['feature']}] KS adj-p={f['ks_p_adjusted']:.3g}  "
                             f"W/sd={w['value']:.3f} [{w['lower']:.3f}, {w['upper']:.3f}]")
                if "paired" in f:
                    lines.append(f"    paired: mean diff={f['paired']['mean_difference']:+.3f}, "
                                 f"Wilcoxon p={_fmt(f['paired']['wilcoxon_p'])}, "
                                 f"frac(T>S)={f['paired']['positive_shift_frac']:.1%}")
                if "autocorrelation_lag1" in f:
                    ac = f["autocorrelation_lag1"]
                    lines.append(f"    autocorr lag-1: source={ac['source']:+.3f}, "
                                 f"target={ac['target']:+.3f}, diff={ac['difference']:+.3f}")

        if self.categorical_marginals and self.categorical_marginals.get("applicable"):
            lines.append("")
            lines.append("  Categorical Marginals:")
            for cf in self.categorical_marginals["features"]:
                tvd = cf["tvd"]
                lines.append(f"  [{cf['feature']}] Chi2 adj-p={cf['chi2_p_adjusted']:.3g}  "
                             f"TVD={tvd['value']:.3f} [{tvd['lower']:.3f}, {tvd['upper']:.3f}]")
                if "paired" in cf:
                    lines.append(f"    paired: category agreement={cf['paired']['category_agreement_frac']:.1%}, "
                                 f"changed={cf['paired']['changed_category_frac']:.1%}")

        for f in self.transform_compatibility["features"]:
            if "nonlinearity_gap" in f:
                g = f["nonlinearity_gap"]
                lines.append(f"  [{f['feature']}] qq-gap={g['value']:.3f} [{g['lower']:.3f}, {g['upper']:.3f}]"
                             f" (floor {f['noise_floor']:.3f} + tol {f['tolerance']}) -> {f['verdict']}")

        if self.dependence.get("applicable"):
            lines.append(f"  copula CvM p={self.dependence['copula_p']:.3g}")
            for pr in self.dependence["pairs"]:
                lines.append(f"    {pr['pair'][0]}~{pr['pair'][1]}: dSpearman={pr['difference']:+.3f} adj-p={pr['p_adjusted']:.3g}")

        if self.dependence.get("categorical_pairs"):
            lines.append("  Categorical Associations (Cramer's V):")
            for cp in self.dependence["categorical_pairs"]:
                lines.append(f"    {cp['pair'][0]}~{cp['pair'][1]}: Cramer's V source={cp['cramers_v_source']:.3f}, "
                             f"target={cp['cramers_v_target']:.3f}, diff={cp['difference']:+.3f}")

        dec = self.decomposition
        if dec.get("applicable"):
            c = lambda k: f"{dec[k]['value']:.4f} [{dec[k]['lower']:.4f}, {dec[k]['upper']:.4f}]"
            lines += ["", f"  discrepancy ({dec['unit']}), permutation p={dec['total_p']:.3g}",
                      f"    total       {c('total')}", f"    marginal    {c('marginal')}",
                      f"    dependence  {c('dependence')}", f"    interaction {c('interaction')}"]
            if dec["verdict"] != "not_detected":
                lines.append(f"    marginal share {dec['marginal_share']['value']:.2f} "
                             f"[{dec['marginal_share']['lower']:.2f}, {dec['marginal_share']['upper']:.2f}]")

        eq = self.equivalence
        if eq["marginal"].get("applicable", True):
            lines += ["", f"  equivalence ({eq['interval']}): marginal smallest certifiable tolerance="
                          f"{_fmt(eq['marginal'].get('smallest_certifiable_delta'))} W/sd (chosen {eq['marginal'].get('delta')})"
                          + (f", p={eq['marginal']['p_value']:.3g}" if "p_value" in eq["marginal"] else "")]
        if eq["dependence"].get("applicable"):
            lines.append(f"               dependence smallest certifiable tolerance="
                         f"{_fmt(eq['dependence'].get('smallest_certifiable_delta'))} max|dSpearman| (chosen {eq['dependence'].get('delta')})"
                         + (f", p={eq['dependence']['p_value']:.3g}" if "p_value" in eq["dependence"] else ""))

        if self.power:
            mm = [f["mde_location_shift_sd"] for f in self.power["marginal"]["features"]]
            fmt = lambda v: "n/a (> grid)" if v is None else f"{v:.3f}"
            lines += ["", "  detectability at these sample sizes (80% power):",
                      "    location shift in SD per feature: " + ", ".join(fmt(v) for v in mm)]
            if self.power["dependence"].get("applicable"):
                lines.append("    mean |dSpearman| screen: " + fmt(self.power["dependence"]["mde_mean_abs_dspearman"]))

        lines += ["", "Limitations:"] + [f"  - {s}" for s in self.limitations]
        lines += [f"Notes: {n}" for n in self.notes]
        return "\n".join(lines)

    def to_dict(self):
        keys = ("alpha", "n_source", "n_target", "features", "marginals", "transform_compatibility",
                "dependence", "decomposition", "equivalence", "categorical_marginals",
                "paired", "time_series", "block_length", "power", "notes", "limitations")
        d = {k: getattr(self, k) for k in keys}
        d["profile"] = self.profile()
        return d

    def to_json(self, path=None, indent=2):
        s = json.dumps(self.to_dict(), indent=indent, default=_jsonable)
        if path:
            with open(path, "w") as fh:
                fh.write(s)
        return s


def diagnose(source, target, alpha=0.05, n_boot=300, n_perm=300, tol=0.1,
             multiplicity="holm", feature_names=None, seed=None,
             delta_marginal=0.1, delta_dependence=0.1, max_n_decomp=800,
             power=False, power_reps=30, equivalence_interval="conservative", n_boot_equivalence=500,
             paired=False, time_series=False, block_length=None, categorical_features=None):
    rng = np.random.default_rng(seed)
    prep = prepare(
        source, target, feature_names=feature_names,
        paired=paired, time_series=time_series, block_length=block_length,
        categorical_features=categorical_features
    )
    S, T, names = prep.source, prep.target, prep.names
    S_cat, T_cat, cat_names = prep.categorical_source, prep.categorical_target, prep.categorical_names
    paired_flag = prep.paired
    block_len = prep.block_length

    marginals = analyze_marginals(
        S, T, names, alpha, n_boot, rng, multiplicity=multiplicity,
        paired=paired_flag, block_length=block_len
    )

    cat_marginals = analyze_categorical_marginals(
        S_cat, T_cat, cat_names, alpha, n_boot, rng,
        paired=paired_flag, block_length=block_len, multiplicity=multiplicity
    )

    rep = DiagnosticReport(
        alpha=alpha,
        n_source=len(S) if len(S) > 0 else (len(S_cat) if S_cat is not None else 0),
        n_target=len(T) if len(T) > 0 else (len(T_cat) if T_cat is not None else 0),
        features=names + cat_names,
        marginals=marginals,
        transform_compatibility=analyze_transform(
            S, T, names, alpha, n_boot, rng, tol, paired=paired_flag, block_length=block_len
        ),
        dependence=analyze_dependence(
            S, T, names, alpha, n_perm, rng, paired=paired_flag, block_length=block_len,
            S_cat=S_cat, T_cat=T_cat, cat_names=cat_names
        ),
        decomposition=analyze_decomposition(S, T, names, alpha, n_boot, n_perm, rng, max_n_decomp),
        equivalence=analyze_equivalence(
            S, T, names, alpha, n_boot_equivalence, rng,
            delta_marginal, delta_dependence, equivalence_interval
        ),
        categorical_marginals=cat_marginals,
        paired=paired_flag,
        time_series=prep.time_series,
        block_length=block_len,
        notes=prep.notes,
    )

    if power and len(S) > 0:
        rep.power = {
            "marginal": marginal_mde(S, len(S), len(T), alpha, rng, n_reps=2 * power_reps),
            "dependence": dependence_mde(S, len(S), len(T), alpha, rng, n_reps=power_reps),
        }
    return rep
