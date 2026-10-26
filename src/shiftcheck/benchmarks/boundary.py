"""Boundary validity and certification power of the equivalence constructions. Layer 4.

The decisive validity check for an equivalence test is a SINGLE statistic sitting
exactly at the tolerance (true effect == delta): a valid test certifies "equivalent"
with probability <= alpha there. (Putting every feature/pair at the tolerance at once
hides invalid constructions, because all must certify together.)

Cells
  validity (true effect = delta):
    marginal location d=1        marginal scale d=1 (different shape than the null simulation)
    dependence Gaussian d=2      dependence t-copula(nu=4) d=2 (different shape than the null simulation)
  power (true effect below delta): same families, effects 0 and delta/2.

Methods: conservative and tost (bootstrap percentile bounds, n_boot draws) and
boundary (boundary-calibrated test, n_sim simulations per statistic).
Rates carry Wilson 95% intervals; a method is flagged *invalid* only if the lower
interval end at the boundary exceeds alpha.
"""
import json
import os
import time
import zlib
import numpy as np
from scipy.stats import gamma, norm, t as student_t, rankdata

from shiftcheck.simulate.generator import make_shift_pair
from shiftcheck.power.equivalence import equivalence_bounds, _marginal_boundary, _dependence_boundary
from shiftcheck.benchmarks.runner import wilson, _rho_for_spearman_gap

ALPHA, DELTA, RHO0, NU = 0.05, 0.1, 0.3, 4


# --------------------------------------------------------------------------- data generators
def _gen_location(e):
    return lambda n, seed: make_shift_pair(n, n, d=1, kinds=["location"], strength=e, seed=seed)[:2]


def _scale_c(e, shape=3.0):
    """c such that T = m + (1+c)(X-m) has true normalized Wasserstein distance e from X."""
    u = (np.arange(2_000_000) + 0.5) / 2_000_000
    q = gamma.ppf(u, a=shape)
    m = float(np.median(q))
    s0 = float((gamma.ppf(0.75, a=shape) - gamma.ppf(0.25, a=shape)) / 1.349)
    return m, e * s0 / float(np.mean(np.abs(q - m)))


def _gen_scale(e):
    m, c = _scale_c(e)

    def gen(n, seed):
        S, T, _ = make_shift_pair(n, n, d=1, kinds=[], seed=seed)
        return S, m + (1 + c) * (T - m)
    return gen


def _gen_gauss_dep(e):
    rt = _rho_for_spearman_gap(RHO0, e)
    return lambda n, seed: make_shift_pair(n, n, d=2, kinds=["copula"], rho=RHO0, rho_target=rt, seed=seed)[:2]


def _tcopula(n, r, rng):
    Z = rng.multivariate_normal([0, 0], [[1, r], [r, 1]], size=n)
    return student_t.cdf(Z / np.sqrt(rng.chisquare(NU, size=(n, 1)) / NU), NU)


def _t_spearman(r, n=400_000, seed=5):
    U = _tcopula(n, r, np.random.default_rng(seed))
    return float(np.corrcoef(rankdata(U[:, 0]), rankdata(U[:, 1]))[0, 1])


def _t_rho_for_gap(gap):
    target = _t_spearman(RHO0) + gap
    lo, hi = RHO0, 0.99
    for _ in range(22):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if _t_spearman(mid) < target else (lo, mid)
    return 0.5 * (lo + hi)


def _gen_t_dep(e):
    rt = _t_rho_for_gap(e) if e > 0 else RHO0

    def gen(n, seed):
        rng = np.random.default_rng(seed)
        return _tcopula(n, RHO0, rng), _tcopula(n, rt, rng)
    return gen


# --------------------------------------------------------------------------- cells
def _cells(scale):
    r = lambda k: max(4, int(k * scale))
    cells = []
    for n in (1000, 3000):
        cells.append(("validity", "marginal location (single feature, true W = 0.1)", "marginal", _gen_location(DELTA), n, r(300)))
        cells.append(("validity", "dependence Gaussian (single pair, true dSpearman = 0.1)", "dependence", _gen_gauss_dep(DELTA), n, r(300)))
    for n in (1000, 3000):
        cells.append(("validity", "marginal SCALE shift (single feature, true W = 0.1)", "marginal", _gen_scale(DELTA), n, r(150)))
        cells.append(("validity", "dependence t-copula nu=4 (single pair, true dSpearman = 0.1)", "dependence", _gen_t_dep(DELTA), n, r(150)))
    for e in (0.0, 0.05):
        for n in (1000, 3000):
            cells.append(("power", f"marginal location true W = {e}", "marginal", _gen_location(e), n, r(100)))
            cells.append(("power", f"dependence Gaussian true dSpearman = {e}", "dependence", _gen_gauss_dep(e), n, r(100)))
    return cells


def _one(prop, S, T, n_boot, n_sim, rng):
    if prop == "marginal":
        b = equivalence_bounds(S, T, n_boot, ALPHA, rng, parts=("marginal",))["marginal"]["upper"]
        p = _marginal_boundary(S, T, DELTA, ALPHA, rng, n_sim, False)["p"]
    else:
        b = equivalence_bounds(S, T, n_boot, ALPHA, rng, parts=("dependence",))["dependence"]["upper"]
        p = _dependence_boundary(S, T, DELTA, ALPHA, rng, n_sim, False)["p"]
    return {"conservative": b["conservative"] < DELTA, "tost": b["tost"] < DELTA, "boundary": p < ALPHA}


def _md(results, meta):
    L = ["# Boundary validity and certification power of equivalence constructions", "",
         f"alpha={ALPHA}, tolerance={DELTA}, bootstrap draws={meta['n_boot']}, boundary simulations={meta['n_sim']}, "
         f"scale={meta['scale']}. Rates are certification rates `rate [Wilson 95%]`.", "",
         "**validity** rows: the true effect equals the tolerance, so a valid method must certify <= alpha.",
         "**power** rows: the true effect is below the tolerance, so higher is better (among valid methods).", "",
         "| type | scenario | n | reps | conservative | tost | boundary | verdict |", "|---|---|---|---|---|---|---|---|"]
    for c in results:
        cells, verdict = [], ""
        for m in ("conservative", "tost", "boundary"):
            p, lo, hi = c["rates"][m]
            cells.append(f"{p:.3f} [{lo:.3f},{hi:.3f}]")
        if c["type"] == "validity":
            bad = [m for m in ("conservative", "tost", "boundary") if c["rates"][m][1] > ALPHA]
            verdict = ("INVALID: " + ", ".join(bad)) if bad else "no method significantly > alpha"
        L.append(f"| {c['type']} | {c['scenario']} | {c['n']} | {c['reps']} | " + " | ".join(cells) + f" | {verdict} |")
    return "\n".join(L)


def run(out_dir="boundary_results", n_boot=500, n_sim=499, scale=1.0, seed=424242, progress=True):
    os.makedirs(out_dir, exist_ok=True)
    meta = {"n_boot": n_boot, "n_sim": n_sim, "scale": scale}
    results, t0 = [], time.time()
    for typ, scen, prop, gen, n, reps in _cells(scale):
        hits = {"conservative": 0, "tost": 0, "boundary": 0}
        key = zlib.crc32(f"{typ}|{scen}|{n}".encode())
        for rep in range(reps):
            sd = int(np.random.SeedSequence([seed, key, rep]).generate_state(1, dtype=np.uint32)[0])
            S, T = gen(n, sd)
            res = _one(prop, S, T, n_boot, n_sim, np.random.default_rng(sd + 1))
            for m in hits:
                hits[m] += bool(res[m])
        results.append({"type": typ, "scenario": scen, "n": n, "reps": reps + 1 if reps else 0,
                        "hits": hits, "rates": {m: wilson(hits[m], reps) for m in hits}})
        results[-1]["reps"] = reps
        with open(os.path.join(out_dir, "boundary_results.json"), "w") as fh:
            json.dump({"meta": meta, "results": results}, fh, indent=1)
        with open(os.path.join(out_dir, "boundary_report.md"), "w") as fh:
            fh.write(_md(results, meta))
        if progress:
            print(f"[{time.time() - t0:.0f}s] {typ} | {scen} | n={n} | " +
                  ", ".join(f"{m}={hits[m]}/{reps}" for m in hits), flush=True)
    return results


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="boundary_results")
    ap.add_argument("--n-boot", type=int, default=500)
    ap.add_argument("--n-sim", type=int, default=499)
    ap.add_argument("--scale", type=float, default=1.0, help="multiply replicate counts (1.0 ~ 35-45 min on one core)")
    a = ap.parse_args()
    run(a.out, a.n_boot, a.n_sim, a.scale)
