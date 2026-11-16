| scenario | method | label | n_systems | coverage_mean | coverage_ci_lo | coverage_ci_hi | coverage_median | coverage_p10 | share_severe_under | share_below_nominal | width_mean_min | width_median_min | width_rel_median | winkler_mean | winkler_median | share_infinite |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| random calibrants (main) | P | Proposed (meta-learned σ) | 228 | 0.911 | 0.908 | 0.913 | 0.911 | 0.888 | 0.000 | 0.219 | 6.213 | 4.677 | 0.350 | 9.070 | 6.738 | 0.000 |
| random calibrants (main) | B1_gauss_rmse | Gaussian ±z·RMSE | 228 | 0.885 | 0.883 | 0.888 | 0.887 | 0.859 | 0.044 | 0.754 | 5.640 | 4.248 | 0.302 | 9.649 | 7.309 | 0.000 |
| random calibrants (main) | B5_jackknife_plus | Jackknife+ | 228 | 0.909 | 0.907 | 0.911 | 0.910 | 0.886 | 0.013 | 0.263 | 6.470 | 5.058 | 0.351 | 9.865 | 7.392 | 0.000 |
| random calibrants (main) | B4_split_conformal | Split conformal | 228 | 0.909 | 0.907 | 0.912 | 0.909 | 0.883 | 0.004 | 0.298 | 8.706 | 7.037 | 0.453 | 11.967 | 9.758 | 0.000 |
| random calibrants (main) | B3_gp_PI | GP projection PI | 228 | 0.822 | 0.818 | 0.826 | 0.823 | 0.783 | 0.825 | 1.000 | 4.182 | 3.401 | 0.224 | 9.511 | 7.444 | 0.000 |
| random calibrants (main) | B7_gp_meta | GP projection, meta-learned prior | 228 | 0.973 | 0.968 | 0.978 | 0.990 | 0.937 | 0.026 | 0.061 | 12.447 | 12.001 | 0.848 | 14.163 | 13.156 | 0.000 |
| random calibrants (main) | B2_pspline_PI | P-spline PI | 228 | 0.851 | 0.847 | 0.854 | 0.852 | 0.811 | 0.491 | 0.956 | 4.693 | 3.726 | 0.256 | 9.925 | 7.684 | 0.000 |
| random calibrants (main) | B6_jk+_ensemble | Jackknife+ (ensemble σ) | 228 | 0.918 | 0.916 | 0.920 | 0.919 | 0.900 | 0.000 | 0.101 | 26.764 | 9.967 | 0.824 | 28.987 | 11.438 | 0.000 |
| commonly measured calibrants | P | Proposed (meta-learned σ) | 228 | 0.815 | 0.794 | 0.833 | 0.840 | 0.605 | 0.535 | 0.662 | 4.946 | 3.863 | 0.203 | 10.446 | 7.931 | 0.000 |
| commonly measured calibrants | B1_gauss_rmse | Gaussian ±z·RMSE | 228 | 0.724 | 0.703 | 0.745 | 0.763 | 0.464 | 0.772 | 0.877 | 3.477 | 2.628 | 0.165 | 12.527 | 9.859 | 0.000 |
| commonly measured calibrants | B5_jackknife_plus | Jackknife+ | 228 | 0.746 | 0.725 | 0.767 | 0.791 | 0.477 | 0.706 | 0.829 | 3.926 | 2.788 | 0.191 | 12.001 | 9.889 | 0.000 |
| commonly measured calibrants | B4_split_conformal | Split conformal | 228 | 0.801 | 0.786 | 0.815 | 0.833 | 0.619 | 0.579 | 0.807 | 6.685 | 4.529 | 0.355 | 13.505 | 10.664 | 0.000 |
| commonly measured calibrants | B3_gp_PI | GP projection PI | 228 | 0.668 | 0.646 | 0.689 | 0.688 | 0.416 | 0.803 | 0.921 | 3.316 | 2.459 | 0.173 | 12.493 | 8.767 | 0.000 |
| commonly measured calibrants | B2_pspline_PI | P-spline PI | 228 | 0.670 | 0.647 | 0.691 | 0.725 | 0.413 | 0.873 | 0.961 | 3.234 | 2.269 | 0.146 | 13.675 | 10.797 | 0.000 |
| commonly measured calibrants | B6_jk+_ensemble | Jackknife+ (ensemble σ) | 228 | 0.902 | 0.891 | 0.912 | 0.933 | 0.781 | 0.224 | 0.368 | 42.034 | 10.888 | 0.636 | 44.360 | 13.087 | 0.000 |
| calibrant RT noise (2% of span) | P | Proposed (meta-learned σ) | 228 | 0.928 | 0.924 | 0.931 | 0.928 | 0.898 | 0.009 | 0.114 | 7.003 | 5.435 | 0.382 | 9.374 | 6.796 | 0.000 |
| calibrant RT noise (2% of span) | B1_gauss_rmse | Gaussian ±z·RMSE | 228 | 0.892 | 0.888 | 0.896 | 0.896 | 0.858 | 0.057 | 0.557 | 5.867 | 4.431 | 0.311 | 9.669 | 7.365 | 0.000 |
| calibrant RT noise (2% of span) | B5_jackknife_plus | Jackknife+ | 228 | 0.913 | 0.910 | 0.916 | 0.915 | 0.883 | 0.018 | 0.263 | 6.642 | 5.100 | 0.360 | 9.882 | 7.421 | 0.000 |
| calibrant RT noise (2% of span) | B4_split_conformal | Split conformal | 228 | 0.912 | 0.908 | 0.915 | 0.912 | 0.875 | 0.018 | 0.307 | 8.878 | 6.738 | 0.463 | 12.047 | 8.948 | 0.000 |
| calibrant RT noise (2% of span) | B3_gp_PI | GP projection PI | 228 | 0.844 | 0.839 | 0.848 | 0.843 | 0.799 | 0.575 | 0.956 | 4.427 | 3.530 | 0.241 | 9.312 | 7.000 | 0.000 |
| calibrant RT noise (2% of span) | B2_pspline_PI | P-spline PI | 228 | 0.868 | 0.864 | 0.872 | 0.870 | 0.831 | 0.259 | 0.846 | 4.946 | 3.862 | 0.268 | 9.680 | 7.554 | 0.000 |
| calibrant RT noise (2% of span) | B6_jk+_ensemble | Jackknife+ (ensemble σ) | 228 | 0.940 | 0.937 | 0.943 | 0.939 | 0.911 | 0.000 | 0.061 | 33.385 | 13.002 | 1.059 | 34.920 | 14.458 | 0.000 |
| one outlying calibrant (+10% span) | P | Proposed (meta-learned σ) | 228 | 0.926 | 0.922 | 0.929 | 0.927 | 0.892 | 0.004 | 0.127 | 6.879 | 5.199 | 0.388 | 9.259 | 6.821 | 0.000 |
| one outlying calibrant (+10% span) | B1_gauss_rmse | Gaussian ±z·RMSE | 228 | 0.895 | 0.892 | 0.899 | 0.899 | 0.862 | 0.053 | 0.504 | 5.906 | 4.593 | 0.313 | 9.619 | 7.398 | 0.000 |
| one outlying calibrant (+10% span) | B5_jackknife_plus | Jackknife+ | 228 | 0.918 | 0.915 | 0.921 | 0.920 | 0.889 | 0.009 | 0.206 | 6.738 | 5.218 | 0.363 | 9.855 | 7.549 | 0.000 |
| one outlying calibrant (+10% span) | B4_split_conformal | Split conformal | 228 | 0.914 | 0.910 | 0.918 | 0.917 | 0.876 | 0.018 | 0.263 | 8.835 | 6.719 | 0.464 | 11.902 | 9.188 | 0.000 |
| one outlying calibrant (+10% span) | B3_gp_PI | GP projection PI | 228 | 0.842 | 0.837 | 0.846 | 0.844 | 0.791 | 0.544 | 0.952 | 4.451 | 3.566 | 0.242 | 9.417 | 7.118 | 0.000 |
| one outlying calibrant (+10% span) | B2_pspline_PI | P-spline PI | 228 | 0.870 | 0.866 | 0.874 | 0.872 | 0.831 | 0.232 | 0.829 | 5.000 | 3.898 | 0.273 | 9.683 | 7.354 | 0.000 |
| one outlying calibrant (+10% span) | B6_jk+_ensemble | Jackknife+ (ensemble σ) | 228 | 0.932 | 0.929 | 0.935 | 0.933 | 0.904 | 0.000 | 0.079 | 31.268 | 12.547 | 1.001 | 33.091 | 14.052 | 0.000 |
| RT-spanning calibrants | P | Proposed (meta-learned σ) | 228 | 0.928 | 0.925 | 0.930 | 0.930 | 0.902 | 0.000 | 0.096 | 6.267 | 4.600 | 0.337 | 8.648 | 6.135 | 0.000 |
| RT-spanning calibrants | B1_gauss_rmse | Gaussian ±z·RMSE | 228 | 0.914 | 0.910 | 0.917 | 0.913 | 0.878 | 0.009 | 0.320 | 6.011 | 4.503 | 0.318 | 9.156 | 6.822 | 0.000 |
| RT-spanning calibrants | B5_jackknife_plus | Jackknife+ | 228 | 0.929 | 0.927 | 0.932 | 0.931 | 0.899 | 0.000 | 0.110 | 6.779 | 5.064 | 0.361 | 9.452 | 7.052 | 0.000 |
| RT-spanning calibrants | B4_split_conformal | Split conformal | 228 | 0.927 | 0.923 | 0.930 | 0.928 | 0.895 | 0.004 | 0.140 | 9.185 | 7.303 | 0.471 | 11.736 | 9.070 | 0.000 |
| RT-spanning calibrants | B3_gp_PI | GP projection PI | 228 | 0.858 | 0.853 | 0.862 | 0.861 | 0.811 | 0.373 | 0.912 | 4.366 | 3.443 | 0.241 | 8.597 | 6.269 | 0.000 |
| RT-spanning calibrants | B2_pspline_PI | P-spline PI | 228 | 0.885 | 0.881 | 0.889 | 0.887 | 0.843 | 0.127 | 0.675 | 4.922 | 3.766 | 0.262 | 8.854 | 6.867 | 0.000 |
| RT-spanning calibrants | B6_jk+_ensemble | Jackknife+ (ensemble σ) | 228 | 0.936 | 0.934 | 0.939 | 0.939 | 0.910 | 0.000 | 0.061 | 29.708 | 9.708 | 0.852 | 31.479 | 10.774 | 0.000 |
