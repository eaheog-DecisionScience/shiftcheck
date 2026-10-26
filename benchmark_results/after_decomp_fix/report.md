# shiftcheck benchmark report

level `quick` · 690 replicates · 296s · alpha=0.05 · equivalence tolerance=0.1

Rates are shown as `rate [Wilson 95% interval]`. With few replicates intervals are wide; an unflagged cell means *no significant departure was found*, not that the property is proven.

## Flags (rates significantly off expectation)

- none

## transform

Each scenario states the true (noise-free) gap of the simulated map. If the true gap is below the 0.1 tolerance, `nonlinear_monotone` is an error; if it is at or above, `affine_compatible` is an error (flag: Wilson lower bound above alpha). Beyond tolerance, the `nonlinear_monotone` rate is power.

| scenario | reps | results |
|---|---|---|
| affine s=0.4 (true gap 0) n=300 | 30 | verdict: affine_compatible=0.03 [0.01,0.17], inconclusive=0.97 [0.83,0.99] |
| monotone s=0.1 (true gap 0.030) n=300 | 30 | verdict: affine_compatible=0.10 [0.03,0.26], inconclusive=0.90 [0.74,0.97] |
| monotone s=0.25 (true gap 0.075) n=300 | 30 | verdict: affine_compatible=0.03 [0.01,0.17], inconclusive=0.97 [0.83,0.99] |
| monotone s=0.5 (true gap 0.147) n=300 | 30 | verdict: inconclusive=1.00 [0.89,1.00] |
| monotone s=1.0 (true gap 0.285) n=300 | 30 | verdict: inconclusive=0.83 [0.66,0.93], nonlinear_monotone=0.17 [0.07,0.34] |
| affine s=0.4 (true gap 0) n=1000 | 30 | verdict: affine_compatible=0.77 [0.59,0.88], inconclusive=0.23 [0.12,0.41] |
| monotone s=0.1 (true gap 0.030) n=1000 | 30 | verdict: affine_compatible=0.53 [0.36,0.70], inconclusive=0.47 [0.30,0.64] |
| monotone s=0.25 (true gap 0.075) n=1000 | 30 | verdict: affine_compatible=0.37 [0.22,0.54], inconclusive=0.63 [0.46,0.78] |
| monotone s=0.5 (true gap 0.147) n=1000 | 30 | verdict: inconclusive=1.00 [0.89,1.00] |
| monotone s=1.0 (true gap 0.285) n=1000 | 30 | verdict: inconclusive=0.20 [0.10,0.37], nonlinear_monotone=0.80 [0.63,0.90] |
| affine s=0.4 (true gap 0) n=3000 | 30 | verdict: affine_compatible=1.00 [0.89,1.00] |
| monotone s=0.1 (true gap 0.030) n=3000 | 30 | verdict: affine_compatible=0.90 [0.74,0.97], inconclusive=0.10 [0.03,0.26] |
| monotone s=0.25 (true gap 0.075) n=3000 | 30 | verdict: affine_compatible=0.50 [0.33,0.67], inconclusive=0.50 [0.33,0.67] |
| monotone s=0.5 (true gap 0.147) n=3000 | 30 | verdict: inconclusive=0.80 [0.63,0.90], nonlinear_monotone=0.20 [0.10,0.37] |
| monotone s=1.0 (true gap 0.285) n=3000 | 30 | verdict: nonlinear_monotone=1.00 [0.89,1.00] |

## decomp

No-shift row: total-discrepancy test is a type-I check. Other rows: marginal share should be high for marginal-only, low for copula-only, and in between for both; `interaction` shows non-additivity.

| scenario | reps | results |
|---|---|---|
| no shift n=200 | 30 | total_detected=0.07 [0.02,0.21]; verdict: mixed_or_unclear=0.07 [0.02,0.21], not_detected=0.93 [0.79,0.98]; marginal_share=0.466 (sd 0.38); marginal=0.000379 (sd 0.015); dependence=-0.000807 (sd 0.014); interaction=0.00147 (sd 0.017) |
| marginal only (monotone 0.4) n=200 | 30 | total_detected=1.00 [0.89,1.00]; verdict: marginal_dominated=1.00 [0.89,1.00]; marginal_share=0.99 (sd 0.02); marginal=1.01 (sd 0.19); dependence=-0.00172 (sd 0.028); interaction=-0.00338 (sd 0.044) |
| copula only (0.3->0.8) n=200 | 30 | total_detected=0.70 [0.52,0.83]; verdict: dependence_dominated=0.17 [0.07,0.34], mixed_or_unclear=0.53 [0.36,0.70], not_detected=0.30 [0.17,0.48]; marginal_share=0.0979 (sd 0.15); marginal=-0.000975 (sd 0.017); dependence=0.0375 (sd 0.017); interaction=0.00149 (sd 0.022) |
| both n=200 | 30 | total_detected=1.00 [0.89,1.00]; verdict: marginal_dominated=1.00 [0.89,1.00]; marginal_share=0.928 (sd 0.044); marginal=0.936 (sd 0.23); dependence=0.0722 (sd 0.045); interaction=-0.174 (sd 0.068) |
| no shift n=400 | 30 | total_detected=0.03 [0.01,0.17]; verdict: mixed_or_unclear=0.03 [0.01,0.17], not_detected=0.97 [0.83,0.99]; marginal_share=0.306 (sd 0.41); marginal=-0.00166 (sd 0.0068); dependence=-0.00168 (sd 0.0054); interaction=0.00237 (sd 0.0063) |
| marginal only (monotone 0.4) n=400 | 30 | total_detected=1.00 [0.89,1.00]; verdict: marginal_dominated=1.00 [0.89,1.00]; marginal_share=0.993 (sd 0.012); marginal=0.945 (sd 0.15); dependence=0.00315 (sd 0.015); interaction=-0.000916 (sd 0.022) |
| copula only (0.3->0.8) n=400 | 30 | total_detected=1.00 [0.89,1.00]; verdict: dependence_dominated=0.87 [0.70,0.95], mixed_or_unclear=0.13 [0.05,0.30]; marginal_share=0.0513 (sd 0.082); marginal=-0.000807 (sd 0.0068); dependence=0.0413 (sd 0.01); interaction=0.000619 (sd 0.012) |
| both n=400 | 30 | total_detected=1.00 [0.89,1.00]; verdict: marginal_dominated=1.00 [0.89,1.00]; marginal_share=0.936 (sd 0.014); marginal=0.906 (sd 0.14); dependence=0.0617 (sd 0.016); interaction=-0.16 (sd 0.026) |
