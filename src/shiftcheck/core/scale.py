"""Robust scale. Layer 0."""
import numpy as np


def robust_scale(x):
    """IQR/1.349, falling back to the SD, then to 1.0 for constant data."""
    q75, q25 = np.quantile(x, [0.75, 0.25])
    s = (q75 - q25) / 1.349
    if s > 0:
        return float(s)
    sd = float(np.std(x))
    return sd if sd > 0 else 1.0
