---
title: 'shiftcheck: A Python Package for Characterizing, Decomposing, and Certifying Distribution Shift'
tags:
  - Python
  - distribution shift
  - dataset drift
  - copula
  - equivalence testing
  - energy distance
  - Shapley attribution
  - biostatistics
  - machine learning safety
authors:
  - name: Deepak Pandey
    affiliation: 1
affiliations:
  - name: Independent Researcher
    index: 1
date: 24 October 2026
bibliography: paper.bib
---

# Summary

In empirical science and machine learning, researchers routinely compare pairs of datasets: evaluating multi-site clinical trial measures, assessing experimental batch-effect normalizations in genomics, auditing demographic survey reweighting, or monitoring production machine learning models against training baselines. Standard distribution shift detectors—such as Kolmogorov–Smirnov (KS) screens, Maximum Mean Discrepancy (MMD), and Classifier Two-Sample Tests (C2ST)—typically treat distribution shift as a binary hypothesis test ($H_0: P_S = P_T$). When sample sizes are large, these tests almost invariably reject $H_0$ with extreme significance ($p < 10^{-15}$), leaving practitioners blind to whether the shift is scientifically negligible, driven by simple affine sensor recalibration, caused by a breakdown in cross-feature dependencies, or restricted to categorical demographics [@rabanser2019failing]. Conversely, at small sample sizes, a non-significant test ($p > 0.05$) is routinely misconstrued as "proven equivalence," falling victim to the classic equivalence fallacy.

`shiftcheck` is a Python package engineered to **characterize, decompose, and certify** distribution shifts rather than merely detect them. Operating under an acyclic layered software architecture, `shiftcheck` evaluates independent properties of the joint distribution: per-feature marginal divergence with bootstrap confidence intervals, robust affine versus nonlinear monotone transformation compatibility, empirical-copula dependence divergence, two-order Shapley attribution on squared energy distance with cross-fitting, and Two One-Sided Tests (TOST) equivalence certification. Furthermore, `shiftcheck` supports time series (via moving block bootstrap), paired longitudinal observations (via sign-flip permutation), and mixed continuous-categorical data.

# Statement of Need

Existing open-source drift detection tools—including `scipy.stats` [@scipy], `alibi-detect`, and `evidently`—suffer from several critical methodological limitations in scientific applications:

1. **Conflation of Marginal and Dependence Shifts:** A shift in joint distribution $P(X)$ can stem from changes in marginal distributions $F_j(x_j)$, shifts in the underlying copula $C(u_1, \dots, u_d)$, or both [@nelsen2006introduction]. Existing omnibus tests sound an uninformative alarm that conflates these phenomena, obscuring whether features can be corrected via simple marginal standard-scaling or require re-learning multivariate relationships.
2. **Absence of Transformation Form Classification:** When marginal drift is detected, practitioners need to know *how* the feature shifted. An affine shift ($x \mapsto a x + b$) can be resolved by linear standardization, whereas a nonlinear monotonic shift requires rank warping or quantile matching.
3. **Finite-Sample Bias in Counterfactual Attribution:** Quantile replacement counterfactuals that share coordinate values with comparison sets suffer from an $O(d/n)$ negative bias in dependence attribution under $H_0$. Most ad-hoc implementations ignore this structural bias.
4. **The Equivalence Fallacy:** Standard software reports $p$-values against the null of strict equality. When underpowered, practitioners erroneously conclude that populations are equivalent. Legitimate scientific certification requires establishing that the effect size lies strictly within an equivalence tolerance $(-\delta, \delta)$ at nominal confidence [@lakens2017equivalence].

`shiftcheck` addresses each of these gaps in a unified, mathematically defensible framework.

# Mathematical Formulation & Architecture

`shiftcheck` decomposes the evaluation into four mutually independent statistical engines:

## 1. Marginal Difference and Transformation Form
For continuous features, `shiftcheck` computes the two-sample Kolmogorov–Smirnov statistic with Holm family-wise error rate control, paired with normalized Wasserstein distances scaled by the source's robust scale $\hat{\sigma}_{\text{IQR}} = (Q_{0.75} - Q_{0.25}) / 1.349$:
$$\mathcal{W}_{\text{norm}}(S_j, T_j) = \frac{1}{\hat{\sigma}_{\text{IQR}}} \int_{-\infty}^{\infty} |F_{S, j}(x) - F_{T, j}(x)| \, dx$$
Percentile bootstrap confidence intervals accompany each distance. To determine the functional form of the marginal shift, `shiftcheck` fits a robust Theil–Sen estimator [@theil1950rank] to empirical quantile-quantile curves $q_T(p)$ vs. $q_S(p)$ on $p \in [0.10, 0.90]$. It measures the root-mean-square deviation (RMSD) from linearity and compares the bootstrap confidence bounds against an empirical sampling noise floor plus user tolerance $\delta_{\text{tol}}$:
$$\text{RMSD} = \frac{1}{\hat{\sigma}_T} \sqrt{\frac{1}{K} \sum_{k=1}^K \left(q_T(p_k) - (\hat{\alpha} + \hat{\beta} q_S(p_k))\right)^2}$$
Classifications commit to `affine_compatible`, `nonlinear_monotone`, or `inconclusive`.

## 2. Dependence Diagnostics via Empirical Copula
To isolate multivariate dependence from marginal shifts, observations are transformed to within-sample pseudo-observations $U_{i, j} = R_{i, j} / (n + 1)$, where $R_{i, j}$ is the rank of observation $i$ along coordinate $j$ [@genest2009goodness]. `shiftcheck` evaluates:
* **Global Copula Shift:** Cramér–von Mises (CvM) distance between empirical copulas $\mathbb{C}_S$ and $\mathbb{C}_T$ over a fixed quasi-random grid in $[0, 1]^d$, calibrated via permutation testing.
* **Pairwise Spearman Differences:** Pairwise rank correlation discrepancies $\Delta \rho_{jk} = \rho_T(j, k) - \rho_S(j, k)$, adjusted using the max-$|t|$ permutation distribution to guarantee strong family-wise error control across all $\binom{d}{2}$ pairs.

## 3. Cross-Fitted Shapley Discrepancy Attribution
`shiftcheck` quantifies total distribution divergence using squared energy distance [@szekely2013energy]:
$$\mathcal{E}^2(P, Q) = 2 \mathbb{E}\|X - Y\| - \mathbb{E}\|X - X'\| - \mathbb{E}\|Y - Y'\|$$
Using Hazen quantile functions $Q_S$ and $Q_T$, counterfactual distributions are constructed: $M = Q_T(U_S)$ (target marginals with source copula) and $C = Q_S(U_T)$ (source marginals with target copula). To prevent coordinate-sharing bias, `shiftcheck` cross-fits the attribution: splitting samples into disjoint folds, constructing counterfactuals on one fold, and evaluating energy distances exclusively against the opposite fold. The two attribution paths yield a two-player Shapley decomposition [@shapley1953value]:
$$\text{Marginal Share} = \frac{\mathcal{E}^2(S, M) + \mathcal{E}^2(C, T)}{2}, \quad \text{Dependence Share} = \frac{\mathcal{E}^2(S, C) + \mathcal{E}^2(M, T)}{2}$$
with the non-additive interaction discrepancy reported explicitly.

## 4. Equivalence Testing and Certifiable Tolerance
Rather than testing the null of equality, `shiftcheck` tests the equivalence null $H_0: \text{Effect} \ge \delta$ using Two One-Sided Tests (TOST) [@lakens2017equivalence]. Equivalence (`equivalent_within_tolerance`) is certified only when the $(1 - \alpha)$ bootstrap upper bound on effect size lies strictly below $\delta$. `shiftcheck` also reports the **smallest certifiable tolerance** $\delta^*$, providing transparent power calibration for negative findings.

## 5. Extended Modalities
* **Time Series (`time_series=True`):** Replaces I.I.D. resampling with a Moving Block Bootstrap (MBB) [@kunsch1989jackknife] with block length $L \approx \lfloor \sqrt[3]{n} \rfloor$ and circular shift permutations, preserving serial autocorrelation up to lag $L$.
* **Paired Data (`paired=True`):** Performs intra-pair joint resampling and exact sign-flip permutation tests ($S_i \leftrightarrow T_i$ with $p = 0.5$), reporting Wilcoxon signed-rank statistics and within-subject difference profiles.
* **Categorical Data (`categorical_features=[...]`):** Computes Total Variation Distance (TVD) with bootstrap confidence intervals, Chi-squared tests of homogeneity with Holm adjustment, and pairwise Cramér's V association shifts.

# Empirical Validation

`shiftcheck` has undergone extensive empirical calibration:
1. **Benchmark Suite:** Validated across 2,340 Monte Carlo replicates covering gamma, log-normal, student-$t_3$, and Gaussian distributions ($n \in [100, 3000]$, $d \in [2, 5]$). The suite established:
   * Empirical false-detection rates strictly $\le \alpha = 0.05$ under all no-shift null scenarios.
   * Zero false certifications of equivalence when true effect sizes equaled or exceeded $\delta = 0.1$.
   * Elimination of the $O(d/n)$ dependence bias through cross-fitting.
2. **Clinical Biopsy Benchmark:** Applied to the Wisconsin Diagnostic Breast Cancer (WDBC) fine-needle aspirate cohort ($n_{\text{benign}} = 357, n_{\text{malignant}} = 212$) [@wolberg1995breast]. `shiftcheck` revealed that the difference between benign and malignant tissue is **$100\%$ marginal-dominated** (marginal share: $1.00$ $[0.98, 1.00]$, dependence CvM $p = 0.096$). Physical coordinate laws ($\text{radius} \sim \text{perimeter} \sim \text{area}$) were preserved with rank difference $\Delta \rho = -0.000, p_{\text{adj}} = 1.00$, and discrete measurement ties were flagged automatically.

# Software Architecture

The package strictly enforces an acyclic dependency structure verified by automated AST analysis:
```
L0 core         (validation, ranks, scale, types)
L1 uncertainty  (resampling, multiplicity) -> core
L2 engines      (marginal, transform, dependence, decomposition, power) -> core, uncertainty
L3 report       (unified diagnose API and diagnostic reports) -> core, uncertainty, L2
L4 benchmarks   (Monte Carlo validation suite) -> all layers
```
The codebase has 100% test passing status across 24 test suites, zero cyclical dependencies, and is released under the permissive MIT License.

# References
