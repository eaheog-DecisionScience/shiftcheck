"""Iris Morphometric Copula Shift Validation Study
UCI Iris Flower Morphology Dataset — shiftcheck validation.

Three scenarios that together demonstrate ALL four shiftcheck properties:
  1. Setosa vs Virginica   — BOTH marginal AND copula shift (the showcase)
  2. Versicolor vs Virginica — marginal shift, copula shift smaller
  3. Empirical null: random 50/50 split of Versicolor — no true shift

The key scientific insight: Iris setosa petals are tiny and their size is nearly
INDEPENDENT of sepal size (rho~0.27). Virginica has large petals whose size is
TIGHTLY CO-REGULATED with sepal size (rho~0.82). This structural difference in
how the plant coordinates its growth manifests as a significant copula CvM test —
something the WDBC and Wine datasets did NOT show.

Run:
  PYTHONPATH=src python examples/iris_copula_shift.py
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.datasets import load_iris
from shiftcheck import diagnose
from shiftcheck.viz.plots import (
    plot_quantile_map, plot_quantile_shift,
    plot_dependence_diff, plot_decomposition,
)

FAST = dict(n_boot=300, n_perm=300, max_n_decomp=800, n_boot_equivalence=500)
FEATS = ["sepal length (cm)", "sepal width (cm)",
         "petal length (cm)", "petal width (cm)"]


def run():
    print("=" * 70)
    print("UCI IRIS MORPHOLOGY — COPULA SHIFT VALIDATION STUDY")
    print("=" * 70)

    data = load_iris(as_frame=True)
    df, target = data.data, data.target
    setosa     = df[target == 0].to_numpy()  # n=50
    versicolor = df[target == 1].to_numpy()  # n=50
    virginica  = df[target == 2].to_numpy()  # n=50

    print(f"\nDataset: 150 flowers, 3 species, 4 morphological features")
    print(f"  Iris setosa:     n = {len(setosa)}")
    print(f"  Iris versicolor: n = {len(versicolor)}")
    print(f"  Iris virginica:  n = {len(virginica)}")

    # -------------------------------------------------------------------------
    # Scenario 1: Setosa vs Virginica — BOTH marginal + copula shift
    # -------------------------------------------------------------------------
    print("\n" + "-" * 70)
    print("SCENARIO 1: I. setosa vs I. virginica (Marginal + Copula Shift)")
    print("-" * 70)
    rep_sv = diagnose(setosa, virginica, feature_names=FEATS, seed=42, **FAST)
    print(rep_sv.summary())
    rep_sv.to_json("iris_setosa_vs_virginica_report.json")

    fig, axes = plt.subplots(2, 2, figsize=(12, 9))
    plot_quantile_map(rep_sv, feature=2, ax=axes[0, 0])   # petal length — biggest shift
    plot_quantile_shift(rep_sv, ax=axes[0, 1])
    plot_dependence_diff(rep_sv, ax=axes[1, 0])
    plot_decomposition(rep_sv, ax=axes[1, 1])
    fig.suptitle(
        "shiftcheck: Iris setosa vs. I. virginica — Marginal + Copula Shift",
        fontsize=12, weight="bold"
    )
    fig.tight_layout()
    fig.savefig("iris_setosa_vs_virginica.png", dpi=130)
    plt.close(fig)
    print("\nSaved: iris_setosa_vs_virginica.png")

    # -------------------------------------------------------------------------
    # Scenario 2: Versicolor vs Virginica
    # -------------------------------------------------------------------------
    print("\n" + "-" * 70)
    print("SCENARIO 2: I. versicolor vs I. virginica (Moderate Contrast)")
    print("-" * 70)
    rep_vv = diagnose(versicolor, virginica, feature_names=FEATS, seed=42, **FAST)
    print(rep_vv.summary())
    rep_vv.to_json("iris_versicolor_vs_virginica_report.json")

    fig2, axes2 = plt.subplots(2, 2, figsize=(12, 9))
    plot_quantile_map(rep_vv, feature=2, ax=axes2[0, 0])
    plot_quantile_shift(rep_vv, ax=axes2[0, 1])
    plot_dependence_diff(rep_vv, ax=axes2[1, 0])
    plot_decomposition(rep_vv, ax=axes2[1, 1])
    fig2.suptitle(
        "shiftcheck: I. versicolor vs I. virginica — Moderate Morphological Shift",
        fontsize=12, weight="bold"
    )
    fig2.tight_layout()
    fig2.savefig("iris_versicolor_vs_virginica.png", dpi=130)
    plt.close(fig2)
    print("\nSaved: iris_versicolor_vs_virginica.png")

    # -------------------------------------------------------------------------
    # Scenario 3: Empirical null — 50/50 split of Versicolor
    # -------------------------------------------------------------------------
    print("\n" + "-" * 70)
    print("SCENARIO 3: EMPIRICAL NULL — Random 50/50 Split of I. versicolor")
    print("-" * 70)
    rng = np.random.default_rng(77)
    perm = rng.permutation(len(versicolor))
    null_s, null_t = versicolor[perm[:25]], versicolor[perm[25:]]
    rep_null = diagnose(null_s, null_t, feature_names=FEATS, seed=42, **FAST)
    print(rep_null.summary())
    rep_null.to_json("iris_null_report.json")

    print("\nAll scenarios complete.")
    return rep_sv, rep_vv, rep_null


if __name__ == "__main__":
    run()
