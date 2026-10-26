"""Bootstrap and permutation machinery. Layer 1 (depends on core)."""
import numpy as np
from shiftcheck.core.types import Estimate


def _resample_indices(n, rng, block_length=None):
    """Draw bootstrap indices: IID or Moving Block Bootstrap (MBB) for time series."""
    if block_length is None or block_length <= 1 or n <= block_length:
        return rng.integers(0, n, n)
    n_blocks = int(np.ceil(n / block_length))
    max_start = n - block_length + 1
    starts = rng.integers(0, max_start, n_blocks)
    idx = np.concatenate([np.arange(s, s + block_length) for s in starts])[:n]
    return idx


def bootstrap_estimates(stat, S, T, n_boot, alpha, rng, paired=False, block_length=None):
    """Percentile bootstrap supporting IID, paired data, and block bootstrap for time series.

    `stat(S, T)` returns a 1-D array; one Estimate is returned per output.
    When `paired=True`, S and T are resampled using identical row indices to preserve pairing.
    When `block_length` is provided, contiguous blocks are drawn to preserve temporal autocorrelation.
    """
    obs = np.atleast_1d(np.asarray(stat(S, T), dtype=float))
    draws = np.empty((n_boot, obs.size))
    n_s, n_t = len(S), len(T)
    for b in range(n_boot):
        if paired:
            idx = _resample_indices(n_s, rng, block_length)
            Sb, Tb = S[idx], T[idx]
        else:
            idx_s = _resample_indices(n_s, rng, block_length)
            idx_t = _resample_indices(n_t, rng, block_length)
            Sb, Tb = S[idx_s], T[idx_t]
        draws[b] = np.atleast_1d(stat(Sb, Tb))
    lo, hi = np.nanquantile(draws, [alpha / 2, 1 - alpha / 2], axis=0)
    method = ("paired block bootstrap" if (paired and block_length) else
              ("paired bootstrap" if paired else
               ("block bootstrap" if block_length else "percentile bootstrap")))
    return [Estimate(float(o), float(l), float(h), method, n_boot)
            for o, l, h in zip(obs, lo, hi)]


def permutation_null(stat, S, T, n_perm, rng, paired=False, block_length=None):
    """Permutation null distribution supporting independent, paired (sign-flip), and time-series shifts."""
    obs = np.atleast_1d(np.asarray(stat(S, T), dtype=float))
    n_s, n_t = len(S), len(T)
    null = np.empty((n_perm, obs.size))
    if paired:
        # Paired sign-flip permutation: randomly swap S_i and T_i with probability 0.5
        for b in range(n_perm):
            swap = rng.random(n_s) < 0.5
            if S.ndim == 1:
                Sp = np.where(swap, T, S)
                Tp = np.where(swap, S, T)
            else:
                Sp = np.where(swap[:, None], T, S)
                Tp = np.where(swap[:, None], S, T)
            null[b] = np.atleast_1d(stat(Sp, Tp))
    elif block_length and block_length > 1:
        # Circular shift permutation preserving internal temporal autocorrelation
        for b in range(n_perm):
            shift_s = rng.integers(0, n_s)
            shift_t = rng.integers(0, n_t)
            Sp = np.roll(S, shift_s, axis=0)
            Tp = np.roll(T, shift_t, axis=0)
            null[b] = np.atleast_1d(stat(Sp, Tp))
    else:
        pooled, n = np.vstack([S, T]), len(S)
        for b in range(n_perm):
            idx = rng.permutation(len(pooled))
            null[b] = np.atleast_1d(stat(pooled[idx[:n]], pooled[idx[n:]]))
    return obs, null


def p_value(obs, null):
    """One-sided Monte Carlo p-value with the +1 correction."""
    return float((1 + np.sum(null >= obs)) / (len(null) + 1))
