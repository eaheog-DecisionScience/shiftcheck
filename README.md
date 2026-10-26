shiftcheck

Characterize, not just detect, distribution shift.

shiftcheck is a statistical diagnostics library for understanding how a target dataset differs from a source dataset.

Unlike many drift-detection tools that return a single test statistic or p-value, shiftcheck decomposes dataset shift into multiple interpretable components:

Are feature marginals different?
Is the difference approximately affine or genuinely nonlinear?
Has the dependence structure changed?
How much of the total discrepancy is attributable to marginal versus dependence shifts?
Is the observed shift practically negligible for a user-defined tolerance?
Given the available sample sizes, what shifts would actually be detectable?

The library focuses on distribution comparison, not model performance estimation.

Motivation

A "dataset shift detected" warning provides limited information.

The same alarm may be triggered by:

a small location shift in a few variables,
a nonlinear monotone transformation,
a change in feature correlations,
or a combination of marginal and dependence changes.

These represent fundamentally different deployment scenarios.

shiftcheck aims to produce a statistically defensible characterization of shift rather than a binary detection outcome.

Core Concepts

The source distribution Ps(X)P_s(X)Ps​(X) and target distribution Pt(X)P_t(X)Pt​(X) are analyzed across several independent properties.

A dataset may simultaneously exhibit:

marginal shift,
nonlinear transformation,
dependence shift,
or any combination thereof.

The package intentionally does not collapse these into a single verdict.

Features
Marginal Analysis

For each feature:

Kolmogorov-Smirnov test
Holm multiple-testing correction
Normalized Wasserstein distance
Bootstrap confidence intervals
Quantile-shift profiles

Example questions:

Which variables changed?

How large are the changes?

Are the observed differences statistically distinguishable from sampling noise?

Marginal Form Analysis

Tests whether changes are approximately affine or require a nonlinear explanation.

Methodology:

Construct feature-wise Q-Q curves.
Fit robust affine transformations.
Measure residual nonlinearity.
Compare against a sampling-noise floor estimated from resampling.

Possible outcomes:

Plain Text
affine_compatible
nonlinear_monotone
inconclusive

This distinguishes:

Plain Text
X_target ≈ aX_source + b

from genuinely nonlinear transformations.

Dependence Shift Analysis

Changes in dependence are assessed independently of marginal changes.

Methods:

empirical copula comparison
Cramér-von Mises statistic
permutation testing

Additional diagnostics:

pairwise Spearman correlation shifts
max-T multiplicity adjustment

Example questions:

Did relationships between variables change?

Are marginals stable while dependency structure drifted?

Shift Attribution

Total discrepancy is decomposed into:

marginal contribution
dependence contribution
interaction component

Methodology:

counterfactual redistribution
two-order Shapley decomposition
squared energy distance
bootstrap uncertainty intervals

This allows interpretation of shifts such as:

Plain Text
73% marginal
19% dependence
8% interaction

rather than a single undifferentiated discrepancy metric.

Equivalence Assessment

Failure to detect a shift does not imply equivalence.

shiftcheck explicitly distinguishes

Plain Text
not detected

from

Plain Text
small enough to certify as negligible

A property is reported as:

Plain Text
equivalent_within_tolerance

only when the entire confidence interval lies below a user-specified tolerance.

Otherwise:

Plain Text
inconclusive

is reported together with the smallest effect size that current data can certify.

Detectability Analysis (Optional)

With

Python
diagnose(..., power=True)

the package estimates:

power curves
minimal detectable location shifts
minimal detectable dependence changes

for the observed sample sizes.

This helps distinguish:

"No evidence of shift"

from

"The experiment is underpowered."

What shiftcheck Does Not Infer

Some properties cannot be determined from features alone.

Specifically:

P(Y∣X)P(Y|X)P(Y∣X)

cannot be assessed without outcome labels.

Therefore:

concept drift
label shift
outcome mechanism changes

are outside the scope of the library unless labels are explicitly incorporated.

Quick Start
Python
from shiftcheck import diagnose
 
report = diagnose(
X_source,
X_target,
seed=0
)
 
print(report.summary())
 
report.to_json("shift_report.json")
Extended Modalities
Time Series

Supports serially dependent observations through Moving Block Bootstrap.

Python
report = diagnose(
stream_train,
stream_deploy,
time_series=True,
block_length=8
)

Features:

dependence-aware uncertainty estimation
autocorrelation diagnostics
Paired Samples

For within-subject comparisons:

Python
report = diagnose(
pre_treatment,
post_treatment,
paired=True
)

Methods:

sign-flip permutation tests
Wilcoxon signed-rank tests
Categorical and Mixed Data
Python
report = diagnose(
source_df,
target_df,
categorical_features=[
"site",
"treatment_arm"
]
)

Methods:

chi-square testing
total variation distance
Cramér's V
Architecture

The package follows a layered architecture.

Plain Text
L0 core
│
├── validation
├── ranking utilities
└── Estimate types
 
L1 uncertainty | simulation | visualization
│
└── depends on core
 
L2 engines
│
├── marginal
├── transform
├── dependence
├── decomposition
└── power
│
└── depends on core + uncertainty
 
L3 reporting
│
└── depends on engines

Design invariants:

engines never import each other
visualization consumes report dictionaries only
layering verified by automated tests
Assumptions

Current implementation is designed primarily for:

continuous tabular data
source/target comparison
moderate dimensions
independent samples

Missing values are removed via listwise deletion.

Statistical Philosophy

The library prioritizes:

effect sizes over p-values,
uncertainty intervals over point estimates,
practical significance over statistical significance,
characterization over detection.

A statistically significant result may still be practically negligible.

A non-significant result may simply indicate insufficient power.

Benchmarking

Run validation experiments:

Shell
python -m shiftcheck.benchmarks --level smoke
python -m shiftcheck.benchmarks --level quick
python -m shiftcheck.benchmarks --level full

Outputs:

Plain Text
results.json
report.md
figures/

Benchmarks evaluate:

type-I error
power
confidence interval coverage
decomposition calibration
transformation classification
Current Limitations

The following capabilities are planned but experimental or not yet implemented:

classifier two-sample tests (C2ST)
tail dependence diagnostics
high-dimensional copula approximations
expanded power simulations
broader time-series support
Installation
Shell
pip install shiftcheck

or

Shell
pip install "shiftcheck[plots,dev]"

Requirements:

numpy
scipy

Optional:

matplotlib
Citation

If you use shiftcheck in research, please cite the repository and specify the version used in analysis.

scikit-learn Integration

The core API is the function `diagnose(source, target)`. An optional fitted
monitor is also available for applications using scikit-learn estimator
conventions:

```python
from shiftcheck.sklearn import ShiftCheck

monitor = ShiftCheck(diagnose_kwargs={"alpha": 0.05, "seed": 0}).fit(X_train)
report = monitor.compare(X_deploy)
```

Install the integration with `pip install "shiftcheck[sklearn]"`. `ShiftCheck`
implements `BaseEstimator` parameter handling and `fit`; `compare` returns a
`DiagnosticReport`. Because it compares datasets rather than transforming or
predicting samples, it is not a `Pipeline` transformer.

License

MIT License.