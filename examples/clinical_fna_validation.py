"""Clinical Validation Case Study: Wisconsin Diagnostic Breast Cancer (WDBC)

This script validates `shiftcheck` on an authentic, peer-reviewed clinical benchmark:
fine-needle aspirate (FNA) biopsy cytology from 569 patients (Wolberg et al.,
University of Wisconsin Hospitals).

We compare two real-world clinical scenarios:
  1. The Biological Shift: Benign (n=357) vs Malignant (n=212) biopsies.
     Tests whether disease progression alters tumor scale (marginal) or disrupts
     cellular geometric coordination (copula dependence).
  2. The Empirical Null: Random 50/50 split of the benign cohort (n=178 vs 179).
     Validates Type-I error control and power-aware 'inconclusive' verdicts under H0.

Run:
  PYTHONPATH=src python examples/clinical_fna_validation.py
"""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from sklearn.datasets import load_breast_cancer

from shiftcheck import diagnose
from shiftcheck.viz.plots import (
    plot_quantile_map,
    plot_quantile_shift,
    plot_dependence_diff,
    plot_decomposition,
)


def run_clinical_validation():
    print("=" * 70)
    print("LOADING WISCONSIN DIAGNOSTIC BREAST CANCER (WDBC) CYTOLOGY BENCHMARK")
    print("=" * 70)

    data = load_breast_cancer(as_frame=True)
    df = data.data
    target = data.target

    # Focus on core morphological and geometric features of the cell nuclei
    features = [
        "mean radius",
        "mean texture",
        "mean perimeter",
        "mean area",
        "mean concavity",
    ]

    benign = df.loc[target == 1, features].to_numpy()
    malignant = df.loc[target == 0, features].to_numpy()

    print(f"\nCohort Overview:")
    print(f"  Source (Benign Biopsies):    n = {len(benign)}")
    print(f"  Target (Malignant Biopsies): n = {len(malignant)}")
    print(f"  Features: {features}")

    # -------------------------------------------------------------------------
    # Scenario 1: Biological Shift (Benign vs Malignant)
    # -------------------------------------------------------------------------
    print("\n" + "-" * 70)
    print("SCENARIO 1: BENIGN vs MALIGNANT COHORT COMPARISON")
    print("-" * 70)
    report_shift = diagnose(benign, malignant, feature_names=features, seed=42)
    print(report_shift.summary())
    report_shift.to_json("clinical_shift_report.json")

    # Generate multi-panel diagnostic figure
    fig, axes = plt.subplots(2, 2, figsize=(12, 9))
    plot_quantile_map(report_shift, feature=0, ax=axes[0, 0])
    plot_quantile_shift(report_shift, ax=axes[0, 1])
    plot_dependence_diff(report_shift, ax=axes[1, 0])
    plot_decomposition(report_shift, ax=axes[1, 1])
    fig.suptitle(
        "shiftcheck Clinical Cytology Diagnostic: Benign vs Malignant Biopsies",
        fontsize=13,
        weight="bold",
    )
    fig.tight_layout()
    fig.savefig("clinical_fna_shift.png", dpi=130)
    print("\nSaved diagnostic plot to 'clinical_fna_shift.png'")

    # -------------------------------------------------------------------------
    # Scenario 2: Empirical Null (Random Split of Benign Cohort)
    # -------------------------------------------------------------------------
    print("\n" + "-" * 70)
    print("SCENARIO 2: EMPIRICAL NULL (RANDOM 50/50 SPLIT OF BENIGN COHORT)")
    print("-" * 70)
    rng = np.random.default_rng(42)
    perm = rng.permutation(len(benign))
    s_null = benign[perm[: len(benign) // 2]]
    t_null = benign[perm[len(benign) // 2 :]]

    report_null = diagnose(s_null, t_null, feature_names=features, seed=42)
    print(report_null.summary())
    report_null.to_json("clinical_null_report.json")
    print("\nValidation completed successfully.")


if __name__ == "__main__":
    run_clinical_validation()
