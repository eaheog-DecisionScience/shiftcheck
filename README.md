# shiftcheck



> **Characterize, not just detect, distribution shift.**

`shiftcheck` is a statistical diagnostics library designed to explain **how** and **why** a target dataset $P_T(X)$ differs from a source dataset $P_S(X)$.

---

##  Motivation

Standard dataset drift tools return a single test statistic or a binary $p$-value indicating whether a shift occurred ($H_0: P_S = P_T$). However, a single alarm can be triggered by vastly different physical or empirical causes:

* A minor sensor calibration offset (location shift) across a few variables.
* A nonlinear monotone transformation in a measurement device.
* A breakdown in feature correlations (dependence shift).
* A combination of marginal and dependence changes.

```
       Source P_S(X)                 Target P_T(X)
   ┌──────────────────┐          ┌──────────────────┐
   │                  │          │  Marginal Shift  │  <-- Location/scale change
   │  Joint Feature   │  ──────> ├──────────────────┤
   │   Distribution   │          │ Transformation   │  <-- Affine vs Nonlinear
   │                  │          ├──────────────────┤
   └──────────────────┘          │ Dependence Shift │  <-- Correlation/Copula drift
                                 └──────────────────┘
```

These scenarios represent fundamentally different real-world conditions. `shiftcheck` decomposes dataset discrepancies into individual, mathematically rigorous, and interpretable components so you can take targeted action.

---

##  Key Features

### 1. Marginal Analysis & Marginal Form
For each individual feature $X_j$:
* **Statistical Divergence:** Kolmogorov–Smirnov (KS) tests with Holm multiple-testing correction, normalized Wasserstein distance, and quantile-shift profiles.
* **Form Classification:** Determines whether feature shifts are approximately affine ($X_T \approx a X_S + b$) or genuinely nonlinear by fitting robust quantile-quantile regressions against an empirical sampling noise floor.

Possible form outcomes:
* `affine_compatible` — Can be corrected via linear standardization.
* `nonlinear_monotone` — Requires rank-warping or quantile matching.
* `inconclusive` — Insufficient evidence to rule out noise.

### 2. Dependence Shift Analysis
Isolates cross-feature relationships from marginal changes by operating on rank-transformed pseudo-observations (empirical copulas):
* **Global Copula Shift:** Cramér–von Mises distance between source and target empirical copulas, calibrated via permutation tests.
* **Pairwise Correlation Discrepancies:** Pairwise Spearman rank correlation differences ($\Delta \rho$) with max-$|T|$ multiplicity adjustment to control family-wise error rates.

### 3. Cross-Fitted Shift Attribution
Uses two-player Shapley decomposition built on squared energy distance $\mathcal{E}^2(P, Q)$ to quantify the relative contribution of each driver:
$$\text{Total Discrepancy} = \text{Marginal Share} + \text{Dependence Share} + \text{Interaction}$$

To prevent coordinate-sharing bias, counterfactual samples are evaluated using out-of-fold cross-fitting, yielding stable percentages such as:
> **73% Marginal Contribution | 19% Dependence Contribution | 8% Interaction**

### 4. Equivalence Certification (TOST)
A failure to reject $H_0$ ($p > 0.05$) does **not** mean datasets are equivalent—it often means the test is underpowered. `shiftcheck` uses Two One-Sided Tests (TOST) to certify true equivalence:
* `equivalent_within_tolerance` — Certified only when the entire $(1 - \alpha)$ confidence interval falls strictly below a user-defined equivalence tolerance $\delta$.
* `inconclusive` — Reported alongside the **smallest certifiable effect size** ($\delta^*$) attainable with the available sample sizes.

---

##  Installation

Install `shiftcheck` via `pip`:

```bash
pip install shiftcheck
```

Or install with optional dependencies for visualization and development:

```bash
pip install "shiftcheck[plots,dev]"
```

**Requirements:** `numpy`, `scipy`

**Optional:** `matplotlib` (for diagnostic plots)

---

##  Quick Start

```python
import numpy as np
from shiftcheck import diagnose

# Generate synthetic source and target data
rng = np.random.default_rng(42)
X_source = rng.normal(loc=0.0, scale=1.0, size=(1000, 4))
X_target = rng.normal(loc=0.5, scale=1.0, size=(1000, 4)) # Location shift added

# Run diagnostic shift analysis
report = diagnose(
    source=X_source,
    target=X_target,
    seed=0
)

# Print high-level summary
print(report.summary())

# Export detailed diagnostic breakdown to JSON
report.to_json("shift_report.json")
```

---

##  Extended Modalities

`shiftcheck` adapts its underlying resampling techniques to match your data structure.

### Time Series Data
Accounts for temporal autocorrelation using Moving Block Bootstrap (MBB):

```python
report = diagnose(
    source=stream_train,
    target=stream_deploy,
    time_series=True,
    block_length=8
)
```

### Paired / Longitudinal Samples
Evaluates within-subject changes via sign-flip permutation and Wilcoxon signed-rank testing:

```python
report = diagnose(
    source=pre_treatment,
    target=post_treatment,
    paired=True
)
```

### Categorical and Mixed Data
Processes discrete coordinates using Total Variation Distance (TVD), Chi-square tests of homogeneity, and Cramér's $V$:

```python
report = diagnose(
    source=source_df,
    target=target_df,
    categorical_features=["site_id", "treatment_group"]
)
```

---

##  Scope & Boundaries

To ensure statistical validity, `shiftcheck` explicitly bounds its scope:

| Included in `shiftcheck` | Outside `shiftcheck` Scope |
| :--- | :--- |
| **Unsupervised Feature Shift $P(X)$** | **Concept Drift $P(Y \mid X)$** (Requires labels $Y$) |
| Marginal & Dependence Decomposition | **Label Shift $P(Y)$** (Requires ground-truth labels) |
| Equivalence & Power Certification | Black-box performance estimation |

---

##  Architecture & Design Principles

`shiftcheck` is built around a layered, acyclic architecture:

```text
L0 Core         (validation, rank transforms, scale estimators, types)
  │
L1 Uncertainty  (bootstrap resampling, block bootstrap, multiplicity adjustments)
  │
L2 Engines      (marginal, transform, dependence, attribution, equivalence)
  │
L3 Reporting    (unified diagnose API, JSON export, diagnostic visualization)
```

**Design Invariants:**
1. **Engine Isolation:** Diagnostic engines in $L2$ operate independently and never import each other.
2. **Representation Decoupling:** Visualization routines consume plain report dictionaries, separating compute from render logic.
3. **No Uninformative Warnings:** Every alert is backed by an effect size, confidence interval, or power diagnostic.

---

##  Benchmarking & Validation

The repository includes a benchmark suite to verify empirical false-positive control, CI coverage, and calibration across various sample sizes and tail distributions:

```bash
# Run quick benchmark validation
python -m shiftcheck.benchmarks --level quick

# Run full Monte Carlo validation suite
python -m shiftcheck.benchmarks --level full
```

---

##  Citation & License

If you use `shiftcheck` in academic research, please cite:

```bibtex
@software{pandey2026shiftcheck,
  author = {Deepak Pandey},
  title = {shiftcheck: A Python Package for Characterizing, Decomposing, and Certifying Distribution Shift},
  year = {2026},
  publisher = {GitHub},
  url = {https://github.com/eaheog-DecisionScience/shiftcheck}
}
```

This project is licensed under the terms of the **MIT License**.