"""Input validation and preprocessing. Layer 0 (pure core, no dependencies)."""
import warnings
from dataclasses import dataclass, field
import numpy as np


@dataclass
class Prepared:
    source: np.ndarray
    target: np.ndarray
    names: list
    notes: list = field(default_factory=list)
    categorical_source: np.ndarray = None
    categorical_target: np.ndarray = None
    categorical_names: list = field(default_factory=list)
    paired: bool = False
    time_series: bool = False
    block_length: int = None


def _to_raw_columns(x, name):
    """Convert input (numpy array, pandas DataFrame, list) to (col_list, names)."""
    if hasattr(x, "columns"):  # pandas DataFrame
        names = list(x.columns)
        cols = [np.asarray(x[col]) for col in names]
        return cols, names
    a = np.asarray(x)
    if a.ndim == 1:
        a = a[:, None]
    if a.ndim != 2:
        raise ValueError(f"{name} must be 1-D or 2-D, got shape {a.shape}")
    d = a.shape[1]
    cols = [a[:, j] for j in range(d)]
    names = [f"x{j}" for j in range(d)]
    return cols, names


def prepare(source, target, feature_names=None, min_n=20, tie_warn=0.05,
            paired=False, time_series=False, block_length=None,
            categorical_features=None):
    """Validate and align source and target data across continuous, categorical, paired, and time-series modalities."""
    s_cols, s_names = _to_raw_columns(source, "source")
    t_cols, t_names = _to_raw_columns(target, "target")

    if len(s_cols) != len(t_cols):
        raise ValueError(f"feature count mismatch: {len(s_cols)} vs {len(t_cols)}")

    d = len(s_cols)
    names = list(feature_names) if feature_names is not None else s_names
    if len(names) != d:
        raise ValueError("feature_names length does not match number of features")

    notes = []
    n_s, n_t = len(s_cols[0]), len(t_cols[0])

    if paired and n_s != n_t:
        raise ValueError(f"paired=True requires equal sample sizes, got {n_s} in source vs {n_t} in target")

    # Determine which features are categorical vs continuous
    cat_mask = [False] * d
    if categorical_features is not None:
        cat_set = set(categorical_features)
        for j, nm in enumerate(names):
            if nm in cat_set or j in cat_set:
                cat_mask[j] = True
    else:
        for j in range(d):
            # Auto-detect string/object or non-numeric types
            s_dt = s_cols[j].dtype
            t_dt = t_cols[j].dtype
            if s_dt.kind in ("O", "U", "S", "b") or t_dt.kind in ("O", "U", "S", "b"):
                cat_mask[j] = True

    cont_indices = [j for j in range(d) if not cat_mask[j]]
    cat_indices = [j for j in range(d) if cat_mask[j]]

    cont_names = [names[j] for j in cont_indices]
    cat_names = [names[j] for j in cat_indices]

    # Process Continuous Features
    if cont_indices:
        S_cont = np.column_stack([np.asarray(s_cols[j], dtype=float) for j in cont_indices])
        T_cont = np.column_stack([np.asarray(t_cols[j], dtype=float) for j in cont_indices])
    else:
        S_cont = np.empty((n_s, 0), dtype=float)
        T_cont = np.empty((n_t, 0), dtype=float)

    # Process Categorical Features
    if cat_indices:
        S_cat = np.column_stack([np.asarray(s_cols[j], dtype=object) for j in cat_indices])
        T_cat = np.column_stack([np.asarray(t_cols[j], dtype=object) for j in cat_indices])
    else:
        S_cat = None
        T_cat = None

    # Handle Missing / Non-finite Values
    if cont_indices:
        if paired:
            # Synchronous pair-wise deletion
            s_fin = np.isfinite(S_cont).all(axis=1)
            t_fin = np.isfinite(T_cont).all(axis=1)
            valid = s_fin & t_fin
            dropped_s = n_s - int(np.sum(valid))
            dropped_t = dropped_s
            Sk = S_cont[valid]
            Tk = T_cont[valid]
            if S_cat is not None:
                S_cat = S_cat[valid]
                T_cat = T_cat[valid]
        else:
            s_fin = np.isfinite(S_cont).all(axis=1)
            t_fin = np.isfinite(T_cont).all(axis=1)
            dropped_s = n_s - int(np.sum(s_fin))
            dropped_t = n_t - int(np.sum(t_fin))
            Sk = S_cont[s_fin]
            Tk = T_cont[t_fin]
            if S_cat is not None:
                S_cat = S_cat[s_fin]
                T_cat = T_cat[t_fin]

        if dropped_s > 0 or dropped_t > 0:
            msg = (f"dropped rows with non-finite values (source={dropped_s}, target={dropped_t}); "
                   "listwise deletion assumes missingness is ignorable")
            warnings.warn(msg)
            notes.append(msg)
    else:
        Sk, Tk = S_cont, T_cont

    if len(Sk) < min_n or len(Tk) < min_n:
        if cont_indices or (S_cat is not None and (len(S_cat) < min_n or len(T_cat) < min_n)):
            raise ValueError(f"need at least {min_n} complete rows per sample (got {len(Sk)}, {len(Tk)})")

    # Tie warnings for continuous features
    for j_local, nm in enumerate(cont_names):
        for lab, A in (("source", Sk), ("target", Tk)):
            frac = 1.0 - len(np.unique(A[:, j_local])) / len(A)
            if frac > tie_warn:
                notes.append(f"feature '{nm}' has {frac:.0%} tied values in {lab}; "
                             "continuous-data methods may be unreliable")

    # Time series configuration
    if time_series:
        if block_length is None or block_length < 1:
            eff_n = min(len(Sk), len(Tk)) if cont_indices else min(len(S_cat), len(T_cat))
            block_length = max(2, int(np.cbrt(eff_n)))
        notes.append(f"time series mode enabled: moving block bootstrap with block_length={block_length}")

    if paired:
        notes.append(f"paired data mode enabled: intra-pair resampling and sign-flip permutation on {len(Sk)} pairs")

    if cat_indices:
        notes.append(f"categorical features ({len(cat_names)}): {cat_names}")

    return Prepared(
        source=Sk,
        target=Tk,
        names=cont_names,
        notes=notes,
        categorical_source=S_cat,
        categorical_target=T_cat,
        categorical_names=cat_names,
        paired=paired,
        time_series=time_series,
        block_length=block_length if time_series else None,
    )
