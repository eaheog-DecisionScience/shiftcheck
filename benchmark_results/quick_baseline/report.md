# shiftcheck benchmark report

level `quick` · 2340 replicates · 1083s · alpha=0.05 · equivalence tolerance=0.1

Rates are shown as `rate [Wilson 95% interval]`. With few replicates intervals are wide; an unflagged cell means *no significant departure was found*, not that the property is proven.

## Flags (rates significantly off expectation)

- ⚠ coverage [location 0.3 SD n=1000] W_covers_truth = 0.87 [0.70,0.95] below 0.95

## typeI

No shift is present, so every detection rate should be <= alpha (flag: Wilson lower bound above alpha).

| scenario | reps | results |
|---|---|---|
| gamma d=2 n=100 | 30 | marginal_detected=0.00 [0.00,0.11]; copula_detected=0.03 [0.01,0.17]; any_pair_detected=0.03 [0.01,0.17] |
| gamma d=5 n=100 | 30 | marginal_detected=0.07 [0.02,0.21]; copula_detected=0.03 [0.01,0.17]; any_pair_detected=0.03 [0.01,0.17] |
| gamma d=2 n=300 | 30 | marginal_detected=0.00 [0.00,0.11]; copula_detected=0.10 [0.03,0.26]; any_pair_detected=0.07 [0.02,0.21] |
| gamma d=5 n=300 | 30 | marginal_detected=0.00 [0.00,0.11]; copula_detected=0.03 [0.01,0.17]; any_pair_detected=0.03 [0.01,0.17] |
| gamma d=2 n=1000 | 30 | marginal_detected=0.00 [0.00,0.11]; copula_detected=0.00 [0.00,0.11]; any_pair_detected=0.00 [0.00,0.11] |
| gamma d=5 n=1000 | 30 | marginal_detected=0.00 [0.00,0.11]; copula_detected=0.07 [0.02,0.21]; any_pair_detected=0.00 [0.00,0.11] |
| normal d=3 n=300 | 30 | marginal_detected=0.10 [0.03,0.26]; copula_detected=0.03 [0.01,0.17]; any_pair_detected=0.07 [0.02,0.21] |
| t3 d=3 n=300 | 30 | marginal_detected=0.10 [0.03,0.26]; copula_detected=0.07 [0.02,0.21]; any_pair_detected=0.00 [0.00,0.11] |
| lognormal d=3 n=300 | 30 | marginal_detected=0.03 [0.01,0.17]; copula_detected=0.07 [0.02,0.21]; any_pair_detected=0.00 [0.00,0.11] |

## power

Informational: detection rates should rise with effect size and sample size.

| scenario | reps | results |
|---|---|---|
| location 0.05 SD n=100 | 30 | marginal_detected=0.10 [0.03,0.26] |
| location 0.1 SD n=100 | 30 | marginal_detected=0.20 [0.10,0.37] |
| location 0.2 SD n=100 | 30 | marginal_detected=0.30 [0.17,0.48] |
| location 0.3 SD n=100 | 30 | marginal_detected=0.57 [0.39,0.73] |
| copula rho 0.3->0.5 n=100 | 30 | marginal_detected=0.00 [0.00,0.11]; copula_detected=0.33 [0.19,0.51]; any_pair_detected=0.30 [0.17,0.48] |
| copula rho 0.3->0.7 n=100 | 30 | marginal_detected=0.00 [0.00,0.11]; copula_detected=0.97 [0.83,0.99]; any_pair_detected=1.00 [0.89,1.00] |
| copula rho 0.3->0.9 n=100 | 30 | marginal_detected=0.00 [0.00,0.11]; copula_detected=1.00 [0.89,1.00]; any_pair_detected=1.00 [0.89,1.00] |
| location 0.05 SD n=300 | 30 | marginal_detected=0.07 [0.02,0.21] |
| location 0.1 SD n=300 | 30 | marginal_detected=0.20 [0.10,0.37] |
| location 0.2 SD n=300 | 30 | marginal_detected=0.70 [0.52,0.83] |
| location 0.3 SD n=300 | 30 | marginal_detected=1.00 [0.89,1.00] |
| copula rho 0.3->0.5 n=300 | 30 | marginal_detected=0.07 [0.02,0.21]; copula_detected=0.90 [0.74,0.97]; any_pair_detected=0.87 [0.70,0.95] |
| copula rho 0.3->0.7 n=300 | 30 | marginal_detected=0.00 [0.00,0.11]; copula_detected=1.00 [0.89,1.00]; any_pair_detected=1.00 [0.89,1.00] |
| copula rho 0.3->0.9 n=300 | 30 | marginal_detected=0.03 [0.01,0.17]; copula_detected=1.00 [0.89,1.00]; any_pair_detected=1.00 [0.89,1.00] |
| location 0.05 SD n=1000 | 30 | marginal_detected=0.13 [0.05,0.30] |
| location 0.1 SD n=1000 | 30 | marginal_detected=0.57 [0.39,0.73] |
| location 0.2 SD n=1000 | 30 | marginal_detected=1.00 [0.89,1.00] |
| location 0.3 SD n=1000 | 30 | marginal_detected=1.00 [0.89,1.00] |
| copula rho 0.3->0.5 n=1000 | 30 | marginal_detected=0.07 [0.02,0.21]; copula_detected=1.00 [0.89,1.00]; any_pair_detected=1.00 [0.89,1.00] |
| copula rho 0.3->0.7 n=1000 | 30 | marginal_detected=0.03 [0.01,0.17]; copula_detected=1.00 [0.89,1.00]; any_pair_detected=1.00 [0.89,1.00] |
| copula rho 0.3->0.9 n=1000 | 30 | marginal_detected=0.10 [0.03,0.26]; copula_detected=1.00 [0.89,1.00]; any_pair_detected=1.00 [0.89,1.00] |

## coverage

Nominal 95% intervals (flag: Wilson upper bound below 0.95). The normalized-Wasserstein interval is expected to behave worst near a true effect of 0 because the statistic is non-negative and noise-inflated.

| scenario | reps | results |
|---|---|---|
| location 0.1 SD n=1000 | 30 | W_covers_truth=0.93 [0.79,0.98]; W_estimate_bias=0.00932 (sd 0.041) |
| location 0.3 SD n=1000 | 30 | W_covers_truth=0.87 [0.70,0.95]; W_estimate_bias=-0.000433 (sd 0.049) |
| affine s=0.2 n=1000 | 30 | slope_covers_truth=0.93 [0.79,0.98]; intercept_covers_truth=0.97 [0.83,0.99] |
| affine s=0.4 n=1000 | 30 | slope_covers_truth=0.90 [0.74,0.97]; intercept_covers_truth=0.97 [0.83,0.99] |

## transform

Under a true affine map `nonlinear_monotone` should be rare (flag: lower bound above alpha). Under monotone maps the rate of `nonlinear_monotone` is the power.

| scenario | reps | results |
|---|---|---|
| affine s=0.4 n=300 | 30 | verdict: affine_compatible=0.07 [0.02,0.21], inconclusive=0.93 [0.79,0.98] |
| monotone s=0.1 n=300 | 30 | verdict: affine_compatible=0.13 [0.05,0.30], inconclusive=0.87 [0.70,0.95] |
| monotone s=0.25 n=300 | 30 | verdict: inconclusive=1.00 [0.89,1.00] |
| monotone s=0.5 n=300 | 30 | verdict: inconclusive=0.97 [0.83,0.99], nonlinear_monotone=0.03 [0.01,0.17] |
| affine s=0.4 n=1000 | 30 | verdict: affine_compatible=0.67 [0.49,0.81], inconclusive=0.33 [0.19,0.51] |
| monotone s=0.1 n=1000 | 30 | verdict: affine_compatible=0.60 [0.42,0.75], inconclusive=0.40 [0.25,0.58] |
| monotone s=0.25 n=1000 | 30 | verdict: affine_compatible=0.33 [0.19,0.51], inconclusive=0.67 [0.49,0.81] |
| monotone s=0.5 n=1000 | 30 | verdict: affine_compatible=0.03 [0.01,0.17], inconclusive=0.97 [0.83,0.99] |
| affine s=0.4 n=3000 | 30 | verdict: affine_compatible=1.00 [0.89,1.00] |
| monotone s=0.1 n=3000 | 30 | verdict: affine_compatible=0.97 [0.83,0.99], inconclusive=0.03 [0.01,0.17] |
| monotone s=0.25 n=3000 | 30 | verdict: affine_compatible=0.60 [0.42,0.75], inconclusive=0.40 [0.25,0.58] |
| monotone s=0.5 n=3000 | 30 | verdict: inconclusive=0.80 [0.63,0.90], nonlinear_monotone=0.20 [0.10,0.37] |

## equiv_marginal

Once the true effect >= the tolerance, certifying equivalence is an error and should occur <= alpha (flag: lower bound above alpha). Below the tolerance the rate is power to certify.

| scenario | reps | results |
|---|---|---|
| true W=0.0 n=300 | 30 | certified_equivalent=0.00 [0.00,0.11]; smallest_certifiable_delta=0.3 (sd 0.043) |
| true W=0.05 n=300 | 30 | certified_equivalent=0.00 [0.00,0.11]; smallest_certifiable_delta=0.325 (sd 0.057) |
| true W=0.1 n=300 | 30 | certified_equivalent=0.00 [0.00,0.11]; smallest_certifiable_delta=0.346 (sd 0.072) |
| true W=0.15 n=300 | 30 | certified_equivalent=0.00 [0.00,0.11]; smallest_certifiable_delta=0.388 (sd 0.074) |
| true W=0.0 n=1000 | 30 | certified_equivalent=0.00 [0.00,0.11]; smallest_certifiable_delta=0.162 (sd 0.027) |
| true W=0.05 n=1000 | 30 | certified_equivalent=0.00 [0.00,0.11]; smallest_certifiable_delta=0.185 (sd 0.036) |
| true W=0.1 n=1000 | 30 | certified_equivalent=0.00 [0.00,0.11]; smallest_certifiable_delta=0.232 (sd 0.042) |
| true W=0.15 n=1000 | 30 | certified_equivalent=0.00 [0.00,0.11]; smallest_certifiable_delta=0.284 (sd 0.032) |
| true W=0.0 n=3000 | 30 | certified_equivalent=0.50 [0.33,0.67]; smallest_certifiable_delta=0.0978 (sd 0.015) |
| true W=0.05 n=3000 | 30 | certified_equivalent=0.13 [0.05,0.30]; smallest_certifiable_delta=0.132 (sd 0.024) |
| true W=0.1 n=3000 | 30 | certified_equivalent=0.00 [0.00,0.11]; smallest_certifiable_delta=0.179 (sd 0.028) |
| true W=0.15 n=3000 | 30 | certified_equivalent=0.00 [0.00,0.11]; smallest_certifiable_delta=0.223 (sd 0.023) |

## equiv_dependence

As for the marginal case, with the effect measured as max |dSpearman|.

| scenario | reps | results |
|---|---|---|
| true dSpearman=0.0 n=300 | 30 | certified_equivalent=0.00 [0.00,0.11]; smallest_certifiable_delta=0.259 (sd 0.04) |
| true dSpearman=0.05 n=300 | 30 | certified_equivalent=0.00 [0.00,0.11]; smallest_certifiable_delta=0.278 (sd 0.052) |
| true dSpearman=0.1 n=300 | 30 | certified_equivalent=0.00 [0.00,0.11]; smallest_certifiable_delta=0.302 (sd 0.053) |
| true dSpearman=0.15 n=300 | 30 | certified_equivalent=0.00 [0.00,0.11]; smallest_certifiable_delta=0.336 (sd 0.047) |
| true dSpearman=0.0 n=1000 | 30 | certified_equivalent=0.00 [0.00,0.11]; smallest_certifiable_delta=0.135 (sd 0.016) |
| true dSpearman=0.05 n=1000 | 30 | certified_equivalent=0.00 [0.00,0.11]; smallest_certifiable_delta=0.152 (sd 0.035) |
| true dSpearman=0.1 n=1000 | 30 | certified_equivalent=0.00 [0.00,0.11]; smallest_certifiable_delta=0.217 (sd 0.029) |
| true dSpearman=0.15 n=1000 | 30 | certified_equivalent=0.00 [0.00,0.11]; smallest_certifiable_delta=0.26 (sd 0.033) |
| true dSpearman=0.0 n=3000 | 30 | certified_equivalent=0.93 [0.79,0.98]; smallest_certifiable_delta=0.0826 (sd 0.012) |
| true dSpearman=0.05 n=3000 | 30 | certified_equivalent=0.20 [0.10,0.37]; smallest_certifiable_delta=0.114 (sd 0.019) |
| true dSpearman=0.1 n=3000 | 30 | certified_equivalent=0.00 [0.00,0.11]; smallest_certifiable_delta=0.161 (sd 0.02) |
| true dSpearman=0.15 n=3000 | 30 | certified_equivalent=0.00 [0.00,0.11]; smallest_certifiable_delta=0.212 (sd 0.023) |

## decomp

No-shift row: total-discrepancy test is a type-I check. Other rows: marginal share should be high for marginal-only, low for copula-only, and in between for both; `interaction` shows non-additivity.

| scenario | reps | results |
|---|---|---|
| no shift n=200 | 30 | total_detected=0.03 [0.01,0.17]; verdict: marginal_dominated=0.03 [0.01,0.17], not_detected=0.97 [0.83,0.99]; marginal_share=1 (sd 0); marginal=-0.00341 (sd 0.011); dependence=-0.0149 (sd 0.0013); interaction=0.0193 (sd 0.0023) |
| marginal only (monotone 0.4) n=200 | 30 | total_detected=1.00 [0.89,1.00]; verdict: marginal_dominated=1.00 [0.89,1.00]; marginal_share=1 (sd 0); marginal=0.95 (sd 0.2); dependence=-0.024 (sd 0.0035); interaction=0.0172 (sd 0.031) |
| copula only (0.3->0.8) n=200 | 30 | total_detected=0.90 [0.74,0.97]; verdict: dependence_dominated=0.10 [0.03,0.26], mixed_or_unclear=0.80 [0.63,0.90], not_detected=0.10 [0.03,0.26]; marginal_share=0.162 (sd 0.25); marginal=0.00117 (sd 0.019); dependence=0.0246 (sd 0.0098); interaction=0.0173 (sd 0.0064) |
| both n=200 | 30 | total_detected=1.00 [0.89,1.00]; verdict: marginal_dominated=1.00 [0.89,1.00]; marginal_share=0.957 (sd 0.017); marginal=0.899 (sd 0.22); dependence=0.0396 (sd 0.015); interaction=-0.14 (sd 0.037) |
| no shift n=400 | 30 | total_detected=0.03 [0.01,0.17]; verdict: marginal_dominated=0.03 [0.01,0.17], not_detected=0.97 [0.83,0.99]; marginal_share=1 (sd 0); marginal=-0.00234 (sd 0.005); dependence=-0.00777 (sd 0.00058); interaction=0.0105 (sd 0.00096) |
| marginal only (monotone 0.4) n=400 | 30 | total_detected=1.00 [0.89,1.00]; verdict: marginal_dominated=1.00 [0.89,1.00]; marginal_share=1 (sd 0); marginal=0.981 (sd 0.12); dependence=-0.0124 (sd 0.0011); interaction=0.00915 (sd 0.019) |
| copula only (0.3->0.8) n=400 | 30 | total_detected=1.00 [0.89,1.00]; verdict: dependence_dominated=0.73 [0.56,0.86], mixed_or_unclear=0.27 [0.14,0.44]; marginal_share=0.0316 (sd 0.073); marginal=-0.00225 (sd 0.0051); dependence=0.0313 (sd 0.0062); interaction=0.0105 (sd 0.0044) |
| both n=400 | 30 | total_detected=1.00 [0.89,1.00]; verdict: marginal_dominated=1.00 [0.89,1.00]; marginal_share=0.944 (sd 0.012); marginal=0.917 (sd 0.14); dependence=0.0532 (sd 0.011); interaction=-0.154 (sd 0.029) |
