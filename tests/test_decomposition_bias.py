"""Regression test: the cross-fitted decomposition must be unbiased under no shift.

The first version (rank/(n+1) counterfactuals compared against the samples they were
built from) had a dependence contribution of about -3/n (-0.015 at n=200, d=3) under
no shift, found by the benchmark suite. The cross-fitted version is ~0.
"""
import numpy as np
from shiftcheck.core.scale import robust_scale
from shiftcheck.decomposition.engine import _components
from shiftcheck.simulate.generator import make_shift_pair


def test_null_contributions_are_unbiased():
    rng = np.random.default_rng(0)
    out = []
    for r in range(30):
        S, T, _ = make_shift_pair(200, 200, d=3, kinds=[], seed=900 + r)
        c = np.median(S, 0)
        sc = np.array([robust_scale(S[:, j]) for j in range(3)])
        A, B = (S - c) / sc, (T - c) / sc
        ia, ib = rng.permutation(200), rng.permutation(200)
        out.append(_components(A[ia[:100]], A[ia[100:]], B[ib[:100]], B[ib[100:]])[:4])
    mean, se = np.mean(out, 0), np.std(out, 0) / np.sqrt(len(out))
    # total, marginal, dependence, interaction: each within 3 SE of zero and far above the old -0.0148
    assert np.all(np.abs(mean) < 3 * se + 1e-3), (mean, se)
    assert mean[2] > -0.008
