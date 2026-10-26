"""Wine Chemical Composition Cultivar Shift Study
UCI Wine Recognition Dataset — Python shiftcheck validation.

Three scenarios:
  1. Cultivar 0 vs Cultivar 2  — two maximally different varietals
  2. Cultivar 0 vs Cultivar 1  — moderately different varietals
  3. Empirical null: random 50/50 split of Cultivar 1 (no true shift)

Run:
  PYTHONPATH=src python examples/wine_cultivar_shift.py
"""
import sys
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.datasets import load_wine
from shiftcheck import diagnose
from shiftcheck.viz.plots import (
    plot_quantile_map, plot_quantile_shift,
    plot_dependence_diff, plot_decomposition,
)


def run():
    print("=" * 70)
    print("UCI WINE RECOGNITION DATASET — CULTIVAR DISTRIBUTION SHIFT ANALYSIS")
    print("=" * 70)

    data = load_wine(as_frame=True)
    df, target = data.data, data.target

    # All 13 chemical features (fully continuous, no tie issues)
    features = list(data.feature_names)

    c0 = df[target == 0].to_numpy()   # Cultivar 0  (n=59)  — Barolo-type
    c1 = df[target == 1].to_numpy()   # Cultivar 1  (n=71)  — Grignolino-type
    c2 = df[target == 2].to_numpy()   # Cultivar 2  (n=48)  — Barbera-type

    print(f"\nDataset: {len(df)} wines, 3 cultivars, {len(features)} chemical features")
    print(f"  Cultivar 0 (Barolo-type):     n = {len(c0)}")
    print(f"  Cultivar 1 (Grignolino-type): n = {len(c1)}")
    print(f"  Cultivar 2 (Barbera-type):    n = {len(c2)}")

    FAST = dict(n_boot=300, n_perm=300, max_n_decomp=800, n_boot_equivalence=500)

    # -------------------------------------------------------------------------
    # Scenario 1: Cultivar 0 vs Cultivar 2 (maximal chemical contrast)
    # -------------------------------------------------------------------------
    print("\n" + "-" * 70)
    print("SCENARIO 1: CULTIVAR 0 vs CULTIVAR 2 (Maximum Varietal Contrast)")
    print("-" * 70)
    rep_02 = diagnose(c0, c2, feature_names=features, seed=42, **FAST)
    print(rep_02.summary())
    rep_02.to_json("wine_c0_vs_c2_report.json")

    # Four-panel diagnostic figure
    fig, axes = plt.subplots(2, 2, figsize=(13, 9))
    plot_quantile_map(rep_02, feature=0, ax=axes[0, 0])   # alcohol
    plot_quantile_shift(rep_02, ax=axes[0, 1])
    plot_dependence_diff(rep_02, ax=axes[1, 0])
    plot_decomposition(rep_02, ax=axes[1, 1])
    fig.suptitle(
        "shiftcheck: UCI Wine — Cultivar 0 vs Cultivar 2 Chemical Shift",
        fontsize=12, weight="bold"
    )
    fig.tight_layout()
    fig.savefig("wine_c0_vs_c2.png", dpi=130)
    plt.close(fig)
    print("\nSaved: wine_c0_vs_c2.png")

    # -------------------------------------------------------------------------
    # Scenario 2: Cultivar 0 vs Cultivar 1 (moderate contrast)
    # -------------------------------------------------------------------------
    print("\n" + "-" * 70)
    print("SCENARIO 2: CULTIVAR 0 vs CULTIVAR 1 (Moderate Varietal Contrast)")
    print("-" * 70)
    rep_01 = diagnose(c0, c1, feature_names=features, seed=42, **FAST)
    print(rep_01.summary())
    rep_01.to_json("wine_c0_vs_c1_report.json")

    fig2, axes2 = plt.subplots(2, 2, figsize=(13, 9))
    plot_quantile_map(rep_01, feature=0, ax=axes2[0, 0])
    plot_quantile_shift(rep_01, ax=axes2[0, 1])
    plot_dependence_diff(rep_01, ax=axes2[1, 0])
    plot_decomposition(rep_01, ax=axes2[1, 1])
    fig2.suptitle(
        "shiftcheck: UCI Wine — Cultivar 0 vs Cultivar 1 Chemical Shift",
        fontsize=12, weight="bold"
    )
    fig2.tight_layout()
    fig2.savefig("wine_c0_vs_c1.png", dpi=130)
    plt.close(fig2)
    print("\nSaved: wine_c0_vs_c1.png")

    # -------------------------------------------------------------------------
    # Scenario 3: Empirical null (random 50/50 split of Cultivar 1)
    # -------------------------------------------------------------------------
    print("\n" + "-" * 70)
    print("SCENARIO 3: EMPIRICAL NULL — Random 50/50 Split of Cultivar 1")
    print("-" * 70)
    rng = np.random.default_rng(99)
    perm = rng.permutation(len(c1))
    mid = len(c1) // 2
    null_s, null_t = c1[perm[:mid]], c1[perm[mid:]]
    rep_null = diagnose(null_s, null_t, feature_names=features, seed=42, **FAST)
    print(rep_null.summary())
    rep_null.to_json("wine_null_report.json")

    print("\nAll scenarios complete.")
    return rep_02, rep_01, rep_null


if __name__ == "__main__":
    run()
