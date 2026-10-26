"""Tests for extended modalities: time series, paired data, and categorical data."""
import numpy as np
import pytest
from shiftcheck import diagnose

FAST = dict(n_boot=60, n_perm=60, max_n_decomp=200, n_boot_equivalence=80)


def test_time_series_modality():
    rng = np.random.default_rng(101)
    n = 200
    # Generate AR(1) process with autocorrelation rho=0.6
    eps_s = rng.normal(size=n)
    eps_t = rng.normal(size=n)
    s = np.zeros(n)
    t = np.zeros(n)
    for i in range(1, n):
        s[i] = 0.6 * s[i - 1] + eps_s[i]
        t[i] = 0.6 * t[i - 1] + eps_t[i] + 0.5  # Mean shift in target

    rep = diagnose(s[:, None], t[:, None], feature_names=["ts_feature"],
                   time_series=True, block_length=6, seed=42, **FAST)
    p = rep.profile()
    assert "time_series" in p["sampling_design"]
    f = rep.marginals["features"][0]
    assert "autocorrelation_lag1" in f
    assert f["autocorrelation_lag1"]["source"] > 0.3
    summary = rep.summary()
    assert "time series" in summary
    assert "autocorr lag-1" in summary


def test_paired_data_modality():
    rng = np.random.default_rng(202)
    n = 150
    # Within-subject pre vs post treatment
    pre = rng.normal(loc=10.0, scale=2.0, size=(n, 2))
    # Post treatment: feature 0 increases by +1.5, feature 1 has zero shift
    post = pre.copy()
    post[:, 0] += 1.5 + rng.normal(scale=0.3, size=n)
    post[:, 1] += rng.normal(scale=0.3, size=n)

    rep = diagnose(pre, post, feature_names=["bp_systolic", "control_biomarker"],
                   paired=True, seed=42, **FAST)
    p = rep.profile()
    assert p["sampling_design"] == "paired"
    f0 = rep.marginals["features"][0]
    assert "paired" in f0
    assert f0["paired"]["mean_difference"] > 1.0
    assert f0["paired"]["wilcoxon_p"] is not None
    assert f0["paired"]["wilcoxon_p"] < 0.001
    assert f0["paired"]["positive_shift_frac"] > 0.95
    summary = rep.summary()
    assert "paired" in summary
    assert "Wilcoxon p" in summary


def test_categorical_data_modality():
    rng = np.random.default_rng(303)
    n = 200
    # Source has balanced categories A, B, C (33% each)
    s_cat = rng.choice(["A", "B", "C"], size=(n, 1), p=[0.33, 0.33, 0.34])
    # Target has heavily shifted categories (70% A, 20% B, 10% C)
    t_cat = rng.choice(["A", "B", "C"], size=(n, 1), p=[0.70, 0.20, 0.10])

    rep = diagnose(s_cat, t_cat, feature_names=["patient_group"],
                   categorical_features=["patient_group"], seed=42, **FAST)
    p = rep.profile()
    assert p["categorical_difference"] == "detected"
    cf = rep.categorical_marginals["features"][0]
    assert cf["chi2_p_adjusted"] < 0.001
    assert cf["tvd"]["value"] > 0.25
    summary = rep.summary()
    assert "Categorical Marginals" in summary
    assert "TVD=" in summary


def test_mixed_continuous_and_categorical_modality():
    rng = np.random.default_rng(404)
    n = 150
    # Mixed data: 2 continuous columns, 1 categorical column
    s_cont = rng.normal(size=(n, 2))
    t_cont = rng.normal(size=(n, 2))  # no shift continuous

    s_cat = rng.choice(["Site_A", "Site_B"], size=(n, 1), p=[0.8, 0.2])
    t_cat = rng.choice(["Site_A", "Site_B"], size=(n, 1), p=[0.2, 0.8])  # shifted site

    S = np.column_stack([s_cont, s_cat])
    T = np.column_stack([t_cont, t_cat])

    rep = diagnose(S, T, feature_names=["biomarker_1", "biomarker_2", "site"],
                   categorical_features=["site"], seed=42, **FAST)
    p = rep.profile()
    assert p["categorical_difference"] == "detected"
    assert len(rep.marginals["features"]) == 2
    assert len(rep.categorical_marginals["features"]) == 1
    summary = rep.summary()
    assert "Continuous Marginals:" in summary
    assert "Categorical Marginals:" in summary
