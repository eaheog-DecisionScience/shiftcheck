import numpy as np
from shiftcheck import diagnose, make_shift_pair
from shiftcheck.power.equivalence import analyze_equivalence, boundary_pvalues, equivalence_bounds


def _pair(e, n=1500, d=2, seed=1):
    S, T, _ = make_shift_pair(n, n, d=d, kinds=["location"], strength=e, seed=seed)
    return S, T


def test_bound_constructions_are_ordered():
    S, T = _pair(0.0)
    b = equivalence_bounds(S, T, 200, 0.05, np.random.default_rng(0))
    u = b["marginal"]["upper"]
    assert u["conservative"] >= u["tost"] > 0


def test_boundary_certifies_null_and_rejects_large_effect():
    rng = np.random.default_rng(0)
    S, T = _pair(0.0, n=3000)
    assert boundary_pvalues(S, T, 0.1, 0.1, 0.05, rng, n_sim=199, find_smallest=False)["marginal"]["p"] < 0.05
    S, T = _pair(0.3, n=1500)
    assert boundary_pvalues(S, T, 0.1, 0.1, 0.05, rng, n_sim=199, find_smallest=False)["marginal"]["p"] > 0.5


def test_boundary_smallest_certifiable_is_below_conservative_bound():
    rng = np.random.default_rng(1)
    S, T = _pair(0.0, n=2000)
    eq_b = analyze_equivalence(S, T, ["a", "b"], 0.05, 199, rng, interval="boundary")
    eq_c = analyze_equivalence(S, T, ["a", "b"], 0.05, 199, rng, interval="conservative")
    assert eq_b["marginal"]["smallest_certifiable_delta"] < eq_c["marginal"]["smallest_certifiable_delta"]


def test_diagnose_boundary_end_to_end_and_summary():
    S, T, _ = make_shift_pair(800, 800, d=3, kinds=["location"], strength=0.0, seed=4)
    r = diagnose(S, T, seed=4, n_boot=40, n_perm=40, max_n_decomp=200, equivalence_interval="boundary", n_boot_equivalence=99)
    assert r.equivalence["interval"] == "boundary"
    assert "equivalence (boundary)" in r.summary()
    assert r.profile()["marginal_difference"] in {"equivalent_within_tolerance", "inconclusive"}


def test_unknown_interval_raises():
    S, T = _pair(0.0, n=300)
    try:
        analyze_equivalence(S, T, ["a", "b"], 0.05, 50, np.random.default_rng(0), interval="nope")
    except ValueError:
        return
    raise AssertionError("expected ValueError for an unknown interval")
