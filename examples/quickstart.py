"""Run with `python examples/quickstart.py` from the repository root.

In Colab or Jupyter, run the diagnosis cell and display `report` (or call
`display(report)`) to see its interactive HTML summary.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from shiftcheck import diagnose, make_shift_pair
from shiftcheck.viz.plots import (plot_quantile_map, plot_quantile_shift, plot_dependence_diff,
                                  plot_decomposition, plot_power_curves)

# Known ground truth: nonlinear monotone marginal shift AND a copula shift together.
S, T, truth = make_shift_pair(n_source=1500, n_target=1500, d=3,
                              kinds=["monotone", "copula"], strength=0.8, seed=7)
report = diagnose(S, T, feature_names=["a", "b", "c"], seed=7, power=True)
print(report.summary())
report.to_json("shift_report.json")

fig, axes = plt.subplots(2, 3, figsize=(15, 8.4))
plot_quantile_map(report, feature=0, ax=axes[0, 0])
plot_quantile_shift(report, ax=axes[0, 1])
plot_dependence_diff(report, ax=axes[0, 2])
plot_decomposition(report, ax=axes[1, 0])
plot_power_curves(report, ax=axes[1, 1])
axes[1, 2].axis("off")
fig.tight_layout()
fig.savefig("shiftcheck_demo.png", dpi=130)
