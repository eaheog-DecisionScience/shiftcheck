"""Multiple-testing corrections. Layer 1."""
import numpy as np


def holm(p):
    p = np.asarray(p, dtype=float)
    m, order = len(p), np.argsort(p)
    adj = np.empty(m)
    run = 0.0
    for rank, i in enumerate(order):
        run = max(run, (m - rank) * p[i])
        adj[i] = min(1.0, run)
    return adj


def benjamini_hochberg(p):
    p = np.asarray(p, dtype=float)
    m, order = len(p), np.argsort(p)[::-1]
    adj = np.empty(m)
    run = 1.0
    for k, i in zip(range(m, 0, -1), order):
        run = min(run, p[i] * m / k)
        adj[i] = run
    return adj


def adjust(p, method="holm"):
    if method == "holm":
        return holm(p)
    if method in ("bh", "fdr_bh"):
        return benjamini_hochberg(p)
    raise ValueError(f"unknown multiplicity method: {method}")
