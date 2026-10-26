"""Benchmark suite with known ground truth. Layer 4: may use every other layer; nothing imports it.

Experiments (each cell is a scenario repeated over independent seeds):

  typeI      no shift: false-detection rates of the marginal screen (Holm) and the
             dependence tests, over sample sizes, dimensions and marginal families.
  power      detection rates vs location-shift size and vs copula change.
  coverage   bootstrap-interval coverage where the truth is known exactly
             (normalized Wasserstein for a pure location shift; affine slope/intercept).
  transform  verdict rates (affine_compatible / nonlinear_monotone / inconclusive).
  equiv_*    rate of certifying "equivalent within tolerance" as the true effect
             crosses the tolerance (must be <= alpha once the true effect >= tolerance).
  decomp     decomposition share and verdicts for marginal-only / copula-only / both,
             plus type-I error of the total-discrepancy permutation test.

All rates carry Wilson 95% intervals; flags mark rates *significantly* off their
expectation (so Monte Carlo noise is not reported as a defect).
"""
import json
import os
import time
import zlib
from concurrent.futures import ProcessPoolExecutor
import numpy as np

from shiftcheck.core.validate import prepare
from shiftcheck.simulate.generator import make_shift_pair
from shiftcheck.marginal.engine import analyze_marginals
from shiftcheck.transform.engine import analyze_transform
from shiftcheck.dependence.engine import analyze_dependence
from shiftcheck.decomposition.engine import analyze_decomposition
from shiftcheck.power.equivalence import analyze_equivalence, equivalence_bounds, METHODS

ALPHA = 0.05
DELTA = 0.1  # equivalence tolerance used throughout


def wilson(k, n, z=1.96):
    if n == 0:
        return (float("nan"),) * 3
    p = k / n
    den = 1 + z * z / n
    c = (p + z * z / (2 * n)) / den
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return p, max(0.0, c - h), min(1.0, c + h)


def _pop_gap(strength, kind="monotone", shape=3.0):
    """Noise-free nonlinearity gap of the simulated map, with the transform engine's own statistic."""
    from scipy.stats import gamma, theilslopes
    from shiftcheck.transform.engine import P_GRID
    qs = gamma.ppf(P_GRID, a=shape)
    qt = qs ** (1 + strength) if kind == "monotone" else strength + (1 + strength) * qs
    slope, icpt, _, _ = theilslopes(qt, qs)
    return float(np.sqrt(np.mean((qt - (icpt + slope * qs)) ** 2)) / np.std(qt))


def _rho_for_spearman_gap(rho, gap):
    """Gaussian-copula rho_target whose Spearman correlation exceeds that of `rho` by `gap`."""
    sp = 6 / np.pi * np.arcsin(rho / 2) + gap
    return float(2 * np.sin(np.pi * sp / 6))


# ----------------------------------------------------------------------------- workers

def _data(p, seed):
    return make_shift_pair(n_source=p["n"], n_target=p["n"], d=p.get("d", 3),
                           kinds=p.get("kinds", ()), strength=p.get("strength", 0.0),
                           rho=p.get("rho", 0.3), rho_target=p.get("rho_target", 0.8),
                           marginal=p.get("marginal", "gamma"), seed=seed)


def _work(task):
    exp, params, seed = task["exp"], task["params"], task["seed"]
    rng = np.random.default_rng(seed + 1)
    S, T, truth = _data(params, seed)
    prep = prepare(S, T)
    S, T, names = prep.source, prep.target, prep.names

    if exp in ("typeI", "power"):
        m = analyze_marginals(S, T, names, ALPHA, 10, rng)
        if exp == "power" and params["kinds"] == ["location"]:   # dependence test irrelevant here
            return {"marginal_detected": m["any_detected"]}
        dep = analyze_dependence(S, T, names, ALPHA, 99, rng)
        return {"marginal_detected": m["any_detected"], "copula_detected": dep["detected"],
                "any_pair_detected": any(q["detected"] for q in dep["pairs"])}
    if exp == "coverage":
        out = {}
        if "location" in params["kinds"]:
            w = analyze_marginals(S, T, names, ALPHA, 200, rng)["features"][0]["wasserstein_normalized"]
            out["W_covers_truth"] = w["lower"] <= params["strength"] <= w["upper"]
            out["W_estimate_bias"] = w["value"] - params["strength"]
        if "affine" in params["kinds"]:
            f = analyze_transform(S, T, names, ALPHA, 100, rng)["features"][0]
            s = params["strength"]
            out["slope_covers_truth"] = f["affine_slope"]["lower"] <= 1 + s <= f["affine_slope"]["upper"]
            out["intercept_covers_truth"] = f["affine_intercept"]["lower"] <= s <= f["affine_intercept"]["upper"]
        return out
    if exp == "transform":
        f = analyze_transform(S, T, names, ALPHA, 100, rng)["features"][0]
        return {"verdict": f["verdict"]}
    if exp in ("equiv_marginal", "equiv_dependence"):
        eq = analyze_equivalence(S, T, names, ALPHA, 100, rng, DELTA, DELTA)
        key = "marginal" if exp == "equiv_marginal" else "dependence"
        return {"certified_equivalent": eq[key]["equivalent"],
                "smallest_certifiable_delta": (eq["marginal"]["smallest_certifiable_delta"] if key == "marginal"
                                               else eq["dependence"]["smallest_certifiable_delta"])}
    if exp == "equiv_compare":
        b = equivalence_bounds(S, T, 100, ALPHA, rng)
        out = {}
        for m in METHODS:
            out[f"marginal_certified_{m}"] = b["marginal"]["upper"][m] < DELTA
            out[f"dependence_certified_{m}"] = b["dependence"]["upper"][m] < DELTA
            out[f"marginal_upper_{m}"] = b["marginal"]["upper"][m]
            out[f"dependence_upper_{m}"] = b["dependence"]["upper"][m]
        return out
    if exp == "decomp":
        d = analyze_decomposition(S, T, names, ALPHA, 40, 99, rng, max_n=10 ** 9)
        return {"total_detected": d["total_p"] < ALPHA, "verdict": d["verdict"],
                "marginal_share": d["marginal_share"]["value"],
                "marginal": d["marginal"]["value"], "dependence": d["dependence"]["value"],
                "interaction": d["interaction"]["value"]}
    raise ValueError(exp)


# ----------------------------------------------------------------------------- configs

def _grid(level):
    L = {"smoke": dict(reps=6, ns=(100, 300), ns_big=(300,), decomp_ns=(150,), cmp=6, cmp_big=4),
         "quick": dict(reps=30, ns=(100, 300, 1000), ns_big=(300, 1000, 3000), decomp_ns=(200, 400), cmp=300, cmp_big=120),
         "full": dict(reps=300, ns=(100, 300, 1000, 3000), ns_big=(300, 1000, 3000, 10000), decomp_ns=(200, 400, 800),
                      cmp=1000, cmp_big=400)}[level]
    reps, ns, nb, dn = L["reps"], L["ns"], L["ns_big"], L["decomp_ns"]
    cfg = []

    def add(exp, tag, params, r=reps):
        cfg.append({"exp": exp, "tag": tag, "params": params, "reps": r})

    for n in ns:
        for d in (2, 5):
            add("typeI", f"gamma d={d} n={n}", {"n": n, "d": d})
    for m in ("normal", "t3", "lognormal"):
        add("typeI", f"{m} d=3 n=300", {"n": 300, "d": 3, "marginal": m})
    for n in ns:
        for dl in (0.05, 0.1, 0.2, 0.3):
            add("power", f"location {dl} SD n={n}", {"n": n, "kinds": ["location"], "strength": dl})
        for rt in (0.5, 0.7, 0.9):
            add("power", f"copula rho 0.3->{rt} n={n}", {"n": n, "kinds": ["copula"], "rho_target": rt})
    n = 1000
    for dl in (0.1, 0.3):
        add("coverage", f"location {dl} SD n={n}", {"n": n, "kinds": ["location"], "strength": dl})
    for s in (0.2, 0.4):
        add("coverage", f"affine s={s} n={n}", {"n": n, "kinds": ["affine"], "strength": s})
    for n in nb:
        add("transform", f"affine s=0.4 (true gap 0) n={n}", {"n": n, "kinds": ["affine"], "strength": 0.4, "pop_gap": 0.0})
        for s in (0.1, 0.25, 0.5, 1.0):
            g = _pop_gap(s)
            add("transform", f"monotone s={s} (true gap {g:.3f}) n={n}",
                {"n": n, "kinds": ["monotone"], "strength": s, "pop_gap": g})
    for n in nb:
        for dl in (0.0, 0.05, 0.1, 0.15):
            add("equiv_marginal", f"true W={dl} n={n}", {"n": n, "kinds": ["location"], "strength": dl})
        for g in (0.0, 0.05, 0.1, 0.15):
            add("equiv_dependence", f"true dSpearman={g} n={n}",
                {"n": n, "kinds": ["copula"], "rho": 0.3, "rho_target": _rho_for_spearman_gap(0.3, g)})
    # Interval-construction comparison: same bootstrap draws, three upper bounds, marginal + dependence
    # effects set simultaneously (they are independent: location shift vs copula change).
    def cmp_cell(n, e, r):
        add("equiv_compare", f"true W={e}, dSpearman={e} n={n}",
            {"n": n, "kinds": ["location", "copula"], "strength": e, "rho": 0.3,
             "rho_target": _rho_for_spearman_gap(0.3, e)}, r=r)

    for n in (1000, 3000):
        for e in (0.0, 0.05, 0.1, 0.15):
            cmp_cell(n, e, L["cmp"])
    cmp_cell(300, 0.0, L["cmp"])
    for e in (0.0, 0.1):
        cmp_cell(10000, e, L["cmp_big"])
    for n in dn:
        add("decomp", f"no shift n={n}", {"n": n}, r=max(6, reps))
        add("decomp", f"marginal only (monotone 0.4) n={n}", {"n": n, "kinds": ["monotone"], "strength": 0.4}, r=max(6, reps))
        add("decomp", f"copula only (0.3->0.8) n={n}", {"n": n, "kinds": ["copula"], "rho_target": 0.8}, r=max(6, reps))
        add("decomp", f"both n={n}", {"n": n, "kinds": ["monotone", "copula"], "strength": 0.4, "rho_target": 0.8}, r=max(6, reps))
    return cfg


# ----------------------------------------------------------------------------- aggregation

def _aggregate(values):
    out = {}
    for k in values[0]:
        col = [v[k] for v in values]
        if isinstance(col[0], (bool, np.bool_)):
            p, lo, hi = wilson(int(np.sum(col)), len(col))
            out[k] = {"rate": p, "lo": lo, "hi": hi, "n": len(col)}
        elif isinstance(col[0], str):
            out[k] = {c: dict(zip(("rate", "lo", "hi"), wilson(col.count(c), len(col))))
                      for c in sorted(set(col))}
            out[k]["n"] = len(col)
        else:
            a = np.asarray(col, dtype=float)
            out[k] = {"mean": float(np.nanmean(a)), "sd": float(np.nanstd(a)), "n": int(np.sum(~np.isnan(a)))}
    return out


def run(level="quick", jobs=1, out_dir="benchmark_results", base_seed=20261023, progress=True, only=None):
    os.makedirs(out_dir, exist_ok=True)
    cfg = [c for c in _grid(level) if only is None or c["exp"] in only]
    tasks = []
    for ci, c in enumerate(cfg):
        for r in range(c["reps"]):
            tasks.append({"exp": c["exp"], "params": c["params"], "ci": ci, "rep": r})
    # Seed depends only on (base_seed, experiment, scenario, replicate), never on which other
    # scenarios are run, so partial re-runs reproduce exactly the same data.
    for t in tasks:
        key = zlib.crc32(f"{t['exp']}|{cfg[t['ci']]['tag']}".encode())
        t["seed"] = int(np.random.SeedSequence([base_seed, key, t["rep"]]).generate_state(1, dtype=np.uint32)[0])
    t0 = time.time()
    results = [None] * len(tasks)
    if jobs > 1:
        with ProcessPoolExecutor(jobs) as ex:
            for i, res in enumerate(ex.map(_work, tasks, chunksize=4)):
                results[i] = res
                if progress and i % 100 == 0:
                    print(f"  {i}/{len(tasks)}  {time.time() - t0:.0f}s", flush=True)
    else:
        for i, t in enumerate(tasks):
            results[i] = _work(t)
            if progress and i % 100 == 0:
                print(f"  {i}/{len(tasks)}  {time.time() - t0:.0f}s", flush=True)
    by_cfg = {}
    for t, r in zip(tasks, results):
        by_cfg.setdefault(t["ci"], []).append(r)
    out = {"level": level, "alpha": ALPHA, "delta": DELTA, "seconds": time.time() - t0,
           "n_tasks": len(tasks), "experiments": {}}
    for ci, c in enumerate(cfg):
        out["experiments"].setdefault(c["exp"], {})[c["tag"]] = {
            "params": c["params"], "reps": c["reps"], "metrics": _aggregate(by_cfg[ci])}
    with open(os.path.join(out_dir, "results.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    from shiftcheck.benchmarks.report import write_report
    write_report(out, out_dir)
    return out
