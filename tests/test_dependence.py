import numpy as np
import pytest

from shiftcheck import diagnose
from shiftcheck.dependence.engine import _ecdf_at


def test_ecdf_at_matches_direct_calculation_with_empty_feature_grid():
    observations = np.empty((5, 0))
    grid = np.empty((3, 0))

    np.testing.assert_array_equal(_ecdf_at(observations, grid), np.ones(3))


def test_ecdf_at_matches_direct_empirical_copula_calculation():
    rng = np.random.default_rng(23)
    observations = rng.random((127, 7))
    grid = rng.random((139, 7))

    expected = np.array([
        np.all(observations <= point, axis=1).mean()
        for point in grid
    ])

    np.testing.assert_array_equal(_ecdf_at(observations, grid, chunk=13), expected)


def test_diagnose_exposes_dependence_grid_size():
    rng = np.random.default_rng(24)
    source = rng.normal(size=(30, 2))
    target = rng.normal(size=(30, 2))

    report = diagnose(
        source,
        target,
        n_boot=5,
        n_perm=5,
        n_boot_equivalence=5,
        max_n_decomp=50,
        dependence_n_grid=16,
        seed=24,
    )

    assert report.dependence["n_grid"] == 16


@pytest.mark.parametrize("value", [0, -1, 1.5, True])
def test_diagnose_rejects_invalid_dependence_grid_size(value):
    with pytest.raises(ValueError, match="dependence_n_grid must be a positive integer"):
        diagnose(np.zeros((2, 2)), np.zeros((2, 2)), dependence_n_grid=value)
