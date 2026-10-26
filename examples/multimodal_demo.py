"""Multimodal Demonstration: Time Series, Paired Data, and Categorical Data.

Run:
  PYTHONPATH=src python examples/multimodal_demo.py
"""
import numpy as np
import pandas as pd
from shiftcheck import diagnose


def demo_multimodal():
    print("=" * 75)
    print("DEMO 1: TIME SERIES MODALITY (Moving Block Bootstrap & Serial Autocorrelation)")
    print("=" * 75)
    rng = np.random.default_rng(101)
    n = 300
    # AR(1) process with autocorrelation rho=0.65
    eps_s = rng.normal(size=n)
    eps_t = rng.normal(size=n)
    s = np.zeros(n)
    t = np.zeros(n)
    for i in range(1, n):
        s[i] = 0.65 * s[i - 1] + eps_s[i]
        t[i] = 0.65 * t[i - 1] + eps_t[i] + 0.4  # Location shift under autocorrelation

    rep_ts = diagnose(
        s[:, None], t[:, None],
        feature_names=["sensor_stream"],
        time_series=True,
        block_length=7,
        seed=42,
        n_boot=100,
        n_perm=100
    )
    print(rep_ts.summary())

    print("\n" + "=" * 75)
    print("DEMO 2: PAIRED DATA MODALITY (Within-Subject Pre-Post Longitudinal Trial)")
    print("=" * 75)
    n_patients = 200
    # Baseline systolic blood pressure and cholesterol
    pre_sbp = rng.normal(loc=140, scale=12, size=n_patients)
    pre_chol = rng.normal(loc=220, scale=20, size=n_patients)
    # Post antihypertensive treatment: SBP drops by ~12 mmHg within patient; cholesterol unchanged
    post_sbp = pre_sbp - 12.0 + rng.normal(loc=0, scale=3, size=n_patients)
    post_chol = pre_chol + rng.normal(loc=0, scale=4, size=n_patients)

    S_paired = np.column_stack([pre_sbp, pre_chol])
    T_paired = np.column_stack([post_sbp, post_chol])

    rep_paired = diagnose(
        S_paired, T_paired,
        feature_names=["systolic_bp", "serum_cholesterol"],
        paired=True,
        seed=42,
        n_boot=100,
        n_perm=100
    )
    print(rep_paired.summary())

    print("\n" + "=" * 75)
    print("DEMO 3: CATEGORICAL & MIXED MODALITY (Demographics, Survey & Biomarkers)")
    print("=" * 75)
    n_sample = 250
    # Continuous features: Age and Income
    age_s = rng.normal(loc=45, scale=10, size=n_sample)
    age_t = rng.normal(loc=46, scale=10, size=n_sample)  # stable age

    # Categorical feature: Clinical Trial Study Arm
    arm_s = rng.choice(["Standard_Care", "Experimental_Drug", "Placebo"], size=n_sample, p=[0.5, 0.25, 0.25])
    arm_t = rng.choice(["Standard_Care", "Experimental_Drug", "Placebo"], size=n_sample, p=[0.2, 0.60, 0.20])

    # Categorical feature: Disease Severity
    stage_s = rng.choice(["Mild", "Moderate", "Severe"], size=n_sample, p=[0.6, 0.3, 0.1])
    stage_t = rng.choice(["Mild", "Moderate", "Severe"], size=n_sample, p=[0.2, 0.4, 0.4])

    df_s = pd.DataFrame({"age": age_s, "trial_arm": arm_s, "severity": stage_s})
    df_t = pd.DataFrame({"age": age_t, "trial_arm": arm_t, "severity": stage_t})

    rep_cat = diagnose(
        df_s, df_t,
        categorical_features=["trial_arm", "severity"],
        seed=42,
        n_boot=100,
        n_perm=100
    )
    print(rep_cat.summary())


if __name__ == "__main__":
    demo_multimodal()
