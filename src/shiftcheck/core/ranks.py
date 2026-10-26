"""Rank / pseudo-observation transforms. Layer 0."""
import numpy as np
from scipy.stats import rankdata


def pseudo_obs(X):
    """Column-wise ranks / (n+1): empirical-CDF coordinates in (0,1)."""
    n = len(X)
    return np.column_stack([rankdata(X[:, j]) for j in range(X.shape[1])]) / (n + 1.0)
