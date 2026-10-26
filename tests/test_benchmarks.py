import numpy as np
from shiftcheck.benchmarks.runner import wilson, _rho_for_spearman_gap, _grid
from shiftcheck.simulate.generator import make_shift_pair


def test_wilson_basic():
    p, lo, hi = wilson(5, 100)
    assert p == 0.05 and lo < p < hi
    assert 0 <= wilson(0, 10)[1] and wilson(10, 10)[2] <= 1


def test_rho_for_spearman_gap_roundtrip():
    r0 = 0.3
    for gap in (0.0, 0.05, 0.1):
        rt = _rho_for_spearman_gap(r0, gap)
        back = 6 / np.pi * (np.arcsin(rt / 2) - np.arcsin(r0 / 2))
        assert abs(back - gap) < 1e-12


def test_grids_build_for_all_levels():
    for lvl in ("smoke", "quick", "full"):
        cfg = _grid(lvl)
        assert {c["exp"] for c in cfg} >= {"typeI", "power", "coverage", "transform", "equiv_marginal", "decomp"}


def test_location_kind_true_effect():
    S, T, tr = make_shift_pair(200000, 200000, d=2, kinds=["location"], strength=0.3, seed=2)
    shift = (np.median(T[:, 0]) - np.median(S[:, 0])) / tr["population_robust_scale"]
    assert abs(shift - 0.3) < 0.02
