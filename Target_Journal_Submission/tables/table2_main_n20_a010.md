| method | label | n_systems | coverage_mean | coverage_ci_lo | coverage_ci_hi | coverage_median | coverage_p10 | share_severe_under | share_below_nominal | width_mean_min | width_median_min | width_rel_median | winkler_mean | winkler_median | share_infinite |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| P | Proposed (meta-learned σ) | 228 | 0.911 | 0.908 | 0.913 | 0.911 | 0.888 | 0.000 | 0.219 | 6.213 | 4.677 | 0.350 | 9.070 | 6.738 | 0.000 |
| B1_gauss_rmse | Gaussian ±z·RMSE | 228 | 0.885 | 0.883 | 0.888 | 0.887 | 0.859 | 0.044 | 0.754 | 5.640 | 4.248 | 0.302 | 9.649 | 7.309 | 0.000 |
| B5_jackknife_plus | Jackknife+ | 228 | 0.909 | 0.907 | 0.911 | 0.910 | 0.886 | 0.013 | 0.263 | 6.470 | 5.058 | 0.351 | 9.865 | 7.392 | 0.000 |
| B4_split_conformal | Split conformal | 228 | 0.909 | 0.907 | 0.912 | 0.909 | 0.883 | 0.004 | 0.298 | 8.706 | 7.037 | 0.453 | 11.967 | 9.758 | 0.000 |
| B3_gp_PI | GP projection PI | 228 | 0.822 | 0.818 | 0.826 | 0.823 | 0.783 | 0.825 | 1.000 | 4.182 | 3.401 | 0.224 | 9.511 | 7.444 | 0.000 |
| B7_gp_meta | GP projection, meta-learned prior | 228 | 0.973 | 0.968 | 0.978 | 0.990 | 0.937 | 0.026 | 0.061 | 12.447 | 12.001 | 0.848 | 14.163 | 13.156 | 0.000 |
| B2_pspline_PI | P-spline PI | 228 | 0.851 | 0.847 | 0.854 | 0.852 | 0.811 | 0.491 | 0.956 | 4.693 | 3.726 | 0.256 | 9.925 | 7.684 | 0.000 |
| B6_jk+_ensemble | Jackknife+ (ensemble σ) | 228 | 0.918 | 0.916 | 0.920 | 0.919 | 0.900 | 0.000 | 0.101 | 26.764 | 9.967 | 0.824 | 28.987 | 11.438 | 0.000 |
