# Boundary validity and certification power of equivalence constructions

alpha=0.05, tolerance=0.1, bootstrap draws=500, boundary simulations=499, scale=1.0. Rates are certification rates `rate [Wilson 95%]`.

**validity** rows: the true effect equals the tolerance, so a valid method must certify <= alpha.
**power** rows: the true effect is below the tolerance, so higher is better (among valid methods).

| type | scenario | n | reps | conservative | tost | boundary | verdict |
|---|---|---|---|---|---|---|---|
| validity | marginal location (single feature, true W = 0.1) | 1000 | 300 | 0.000 [0.000,0.013] | 0.000 [0.000,0.013] | 0.040 [0.023,0.069] | no method significantly > alpha |
| validity | dependence Gaussian (single pair, true dSpearman = 0.1) | 1000 | 300 | 0.013 [0.005,0.034] | 0.047 [0.028,0.077] | 0.050 [0.031,0.081] | no method significantly > alpha |
| validity | marginal location (single feature, true W = 0.1) | 3000 | 300 | 0.017 [0.007,0.038] | 0.040 [0.023,0.069] | 0.067 [0.044,0.101] | no method significantly > alpha |
| validity | dependence Gaussian (single pair, true dSpearman = 0.1) | 3000 | 300 | 0.027 [0.014,0.052] | 0.033 [0.018,0.060] | 0.037 [0.021,0.064] | no method significantly > alpha |
| validity | marginal SCALE shift (single feature, true W = 0.1) | 1000 | 150 | 0.000 [0.000,0.025] | 0.000 [0.000,0.025] | 0.007 [0.001,0.037] | no method significantly > alpha |
| validity | dependence t-copula nu=4 (single pair, true dSpearman = 0.1) | 1000 | 150 | 0.020 [0.007,0.057] | 0.040 [0.018,0.085] | 0.060 [0.032,0.110] | no method significantly > alpha |
| validity | marginal SCALE shift (single feature, true W = 0.1) | 3000 | 150 | 0.007 [0.001,0.037] | 0.020 [0.007,0.057] | 0.020 [0.007,0.057] | no method significantly > alpha |
| validity | dependence t-copula nu=4 (single pair, true dSpearman = 0.1) | 3000 | 150 | 0.020 [0.007,0.057] | 0.040 [0.018,0.085] | 0.060 [0.032,0.110] | no method significantly > alpha |
| power | marginal location true W = 0.0 | 1000 | 100 | 0.000 [0.000,0.037] | 0.010 [0.002,0.054] | 0.630 [0.532,0.718] |  |
| power | dependence Gaussian true dSpearman = 0.0 | 1000 | 100 | 0.220 [0.150,0.311] | 0.460 [0.366,0.557] | 0.430 [0.337,0.528] |  |
| power | marginal location true W = 0.0 | 3000 | 100 | 0.900 [0.826,0.945] | 0.940 [0.875,0.972] | 0.980 [0.930,0.994] |  |
| power | dependence Gaussian true dSpearman = 0.0 | 3000 | 100 | 0.980 [0.930,0.994] | 0.990 [0.946,0.998] | 0.990 [0.946,0.998] |  |
| power | marginal location true W = 0.05 | 1000 | 100 | 0.000 [0.000,0.037] | 0.020 [0.006,0.070] | 0.390 [0.300,0.488] |  |
| power | dependence Gaussian true dSpearman = 0.05 | 1000 | 100 | 0.100 [0.055,0.174] | 0.260 [0.184,0.354] | 0.280 [0.201,0.375] |  |
| power | marginal location true W = 0.05 | 3000 | 100 | 0.350 [0.264,0.447] | 0.480 [0.385,0.577] | 0.580 [0.482,0.672] |  |
| power | dependence Gaussian true dSpearman = 0.05 | 3000 | 100 | 0.570 [0.472,0.663] | 0.710 [0.615,0.790] | 0.700 [0.604,0.781] |  |

Note: the final power cell (dependence, true dSpearman = 0.05, n=3000) was re-run separately with identical seeds after the original run was interrupted.
