"""Markdown + figure report for benchmark results. Layer 4."""
import os
import re
import numpy as np

TRANSFORM_TOL = 0.1  # default excess-nonlinearity tolerance of diagnose()

EXPECT_NOTES = {
    "typeI": "No shift is present, so every detection rate should be <= alpha (flag: Wilson lower bound above alpha).",
    "power": "Informational: detection rates should rise with effect size and sample size.",
    "coverage": "Nominal 95% intervals (flag: Wilson upper bound below 0.95). The normalized-Wasserstein interval is expected to behave worst near a true effect of 0 because the statistic is non-negative and noise-inflated.",
    "transform": "Each scenario states the true (noise-free) gap of the simulated map. If the true gap is below the 0.1 tolerance, `nonlinear_monotone` is an error; if it is at or above, `affine_compatible` is an error (flag: Wilson lower bound above alpha). Beyond tolerance, the `nonlinear_monotone` rate is power.",
    "equiv_marginal": "Once the true effect >= the tolerance, certifying equivalence is an error and should occur <= alpha (flag: lower bound above alpha). Below the tolerance the rate is power to certify.",
    "equiv_dependence": "As for the marginal case, with the effect measured as max |dSpearman|.",
    "equiv_compare": "Three upper-bound constructions on identical bootstrap draws (conservative = 97.5th percentile, tost = 95th percentile, basic = pivotal 2*estimate - 5th percentile). A construction is *valid* if its certification rate is <= alpha whenever the true effect >= the tolerance (the boundary row, true effect = 0.1, is the stress test). Among valid constructions, higher certification rates below the tolerance mean a more practical test.",
    "decomp": "No-shift row: total-discrepancy test is a type-I check. Other rows: marginal share should be high for marginal-only, low for copula-only, and in between for both; `interaction` shows non-additivity.",
}


def _fmt(m):
    if "rate" in m:
        return f"{m['rate']:.2f} [{m['lo']:.2f},{m['hi']:.2f}]"
    return f"{m['mean']:.3g} (sd {m['sd']:.2g})"


def _row_metrics(metrics):
    parts = []
    for k, m in metrics.items():
        if "rate" in m or "mean" in m:
            parts.append(f"{k}={_fmt(m)}")
        else:
            cats = ", ".join(f"{c}={_fmt(v)}" for c, v in m.items() if c != "n")
            parts.append(f"{k}: {cats}")
    return "; ".join(parts)


def _true_effect(exp, params):
    if exp == "equiv_marginal":
        return params["strength"]
    if exp == "equiv_dependence":
        return 6 / np.pi * (np.arcsin(params["rho_target"] / 2) - np.arcsin(params["rho"] / 2))
    return None


def _flags(out):
    a, delta = out["alpha"], out["delta"]
    flags = []
    for exp, cells in out["experiments"].items():
        for tag, c in cells.items():
            M = c["metrics"]
            if exp == "typeI":
                for k, m in M.items():
                    if m["lo"] > a:
                        flags.append(f"typeI [{tag}] {k} = {m['rate']:.2f} [{m['lo']:.2f},{m['hi']:.2f}] exceeds alpha")
            elif exp == "coverage":
                for k, m in M.items():
                    if k.endswith("covers_truth") and m["hi"] < 0.95:
                        flags.append(f"coverage [{tag}] {k} = {m['rate']:.2f} [{m['lo']:.2f},{m['hi']:.2f}] below 0.95")
            elif exp == "transform":
                g = c["params"].get("pop_gap")
                if g is not None and g < TRANSFORM_TOL:      # truth within tolerance -> nonlinear is an error
                    m = M["verdict"].get("nonlinear_monotone")
                    if m and m["lo"] > a:
                        flags.append(f"transform [{tag}] `nonlinear_monotone` = {m['rate']:.2f} [{m['lo']:.2f},{m['hi']:.2f}] "
                                     f"although true gap {g:.3f} < tolerance")
                elif g is not None:                           # truth beyond tolerance -> affine_compatible is an error
                    m = M["verdict"].get("affine_compatible")
                    if m and m["lo"] > a:
                        flags.append(f"transform [{tag}] false `affine_compatible` = {m['rate']:.2f} [{m['lo']:.2f},{m['hi']:.2f}] "
                                     f"although true gap {g:.3f} >= tolerance")
            elif exp in ("equiv_marginal", "equiv_dependence"):
                eff = _true_effect(exp, c["params"])
                m = M["certified_equivalent"]
                if eff >= delta - 1e-9 and m["lo"] > a:
                    flags.append(f"{exp} [{tag}] certified equivalent at true effect {eff:.3f} >= tolerance: "
                                 f"{m['rate']:.2f} [{m['lo']:.2f},{m['hi']:.2f}]")
            elif exp == "equiv_compare":
                w = c["params"]["strength"]
                g = 6 / np.pi * (np.arcsin(c["params"]["rho_target"] / 2) - np.arcsin(c["params"]["rho"] / 2))
                for prop, eff in (("marginal", w), ("dependence", g)):
                    if eff < delta - 1e-9:
                        continue
                    for meth in ("conservative", "tost", "basic"):
                        m = M[f"{prop}_certified_{meth}"]
                        if m["lo"] > a:
                            flags.append(f"equiv_compare [{tag}] {prop}/{meth} certifies {m['rate']:.3f} "
                                         f"[{m['lo']:.3f},{m['hi']:.3f}] at true effect {eff:.3f} >= tolerance (invalid)")
            elif exp == "decomp" and tag.startswith("no shift"):
                m = M["total_detected"]
                if m["lo"] > a:
                    flags.append(f"decomp [{tag}] total test false-detection {m['rate']:.2f} [{m['lo']:.2f},{m['hi']:.2f}]")
    return flags


def write_report(out, out_dir):
    L = [f"# shiftcheck benchmark report", "",
         f"level `{out['level']}` · {out['n_tasks']} replicates · {out['seconds']:.0f}s · alpha={out['alpha']} · equivalence tolerance={out['delta']}",
         "", "Rates are shown as `rate [Wilson 95% interval]`. With few replicates intervals are wide; "
         "an unflagged cell means *no significant departure was found*, not that the property is proven.", ""]
    fl = _flags(out)
    L += ["## Flags (rates significantly off expectation)", ""]
    L += [f"- ⚠ {f}" for f in fl] if fl else ["- none"]
    L.append("")
    for exp, cells in out["experiments"].items():
        if exp == "equiv_compare":
            L += [f"## {exp}", "", EXPECT_NOTES[exp], "",
                  "Certification rate [Wilson 95%] by construction; mean upper bound in brackets.", "",
                  "| scenario | reps | marginal: conservative | tost | basic | dependence: conservative | tost | basic |",
                  "|---|---|---|---|---|---|---|---|"]
            for tag, c in cells.items():
                M = c["metrics"]
                cell = lambda p, m: f"{_fmt(M[f'{p}_certified_{m}'])} ({M[f'{p}_upper_{m}']['mean']:.3f})"
                L.append(f"| {tag} | {c['reps']} | " + " | ".join(cell(p, m) for p in ("marginal", "dependence")
                                                                      for m in ("conservative", "tost", "basic")) + " |")
            L.append("")
            continue
        L += [f"## {exp}", "", EXPECT_NOTES.get(exp, ""), "", "| scenario | reps | results |", "|---|---|---|"]
        for tag, c in cells.items():
            L.append(f"| {tag} | {c['reps']} | {_row_metrics(c['metrics'])} |")
        L.append("")
    with open(os.path.join(out_dir, "report.md"), "w") as fh:
        fh.write("\n".join(L))
    try:
        _figure(out, os.path.join(out_dir, "benchmark_figures.png"))
    except ImportError:  # matplotlib optional
        pass


def _figure(out, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    E = out["experiments"]
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.2))

    def curve(a, exp, pat, xs_fn, metric, label_fn):
        series = {}
        for tag, c in E.get(exp, {}).items():
            m = re.match(pat, tag)
            if not m:
                continue
            x, n = xs_fn(m, c), c["params"]["n"]
            mm = c["metrics"][metric]
            series.setdefault(n, []).append((x, mm["rate"], mm["lo"], mm["hi"]))
        for n, pts in sorted(series.items()):
            pts.sort()
            x, r, lo, hi = zip(*pts)
            a.errorbar(x, r, yerr=[np.maximum(np.array(r) - lo, 0), np.maximum(np.array(hi) - r, 0)], marker="o", ms=3, capsize=2, label=f"n={n}")

    curve(ax[0], "power", r"location ([\d.]+) SD", lambda m, c: float(m.group(1)), "marginal_detected", None)
    ax[0].set(xlabel="location shift (source robust SD)", ylabel="detection rate", title="Marginal screen power", ylim=(-0.02, 1.02))
    curve(ax[1], "power", r"copula rho 0.3->([\d.]+)", lambda m, c: 6 / np.pi * (np.arcsin(float(m.group(1)) / 2) - np.arcsin(0.15)),
          "copula_detected", None)
    ax[1].set(xlabel="true dSpearman", ylabel="detection rate", title="Copula test power", ylim=(-0.02, 1.02))
    curve(ax[2], "equiv_marginal", r"true W=([\d.]+)", lambda m, c: float(m.group(1)), "certified_equivalent", None)
    ax[2].axvline(out["delta"], color="grey", ls="--", lw=0.8)
    ax[2].axhline(out["alpha"], color="grey", ls=":", lw=0.8)
    ax[2].set(xlabel="true normalized Wasserstein", ylabel="rate certified 'equivalent'",
              title=f"Equivalence certification (tolerance {out['delta']})", ylim=(-0.02, 1.02))
    for a in ax:
        a.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=120)
