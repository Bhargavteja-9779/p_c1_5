---
title: "Distribution-free retention-time windows for candidate filtering in LC–HRMS non-target screening on new chromatographic systems"
---

[First Author]^a^, [Second Author]^a,\*^

^a^ [Department, Institution, Address, Country]

\* Corresponding author: [name, e-mail address]

# Highlights {.unnumbered}
- Parametric and fixed RT windows lose correct structures on new LC systems
- Calibrant-based conformal windows keep a user-chosen retention rate
- At least 1/α − 1 calibrants are needed for a finite window
- Calibrants must resemble the analytes; known selection rules can be corrected
- Tested on {{n_target_systems}} RepoRT systems and {{norman.n_systems}} interlaboratory NTS methods

# Abstract {.unnumbered}
**Background.** Retention-time (RT) filters discard candidate structures in liquid chromatography–high-resolution mass spectrometry (LC–HRMS) non-target screening, but their tolerance windows are set heuristically, and the probability that the correct structure survives on a laboratory's own chromatographic system is unknown. We aimed to construct RT windows that retain the correct structure at a user-chosen rate using only the calibrants that a laboratory measures.

**Results.** A retention-order model trained on public RepoRT data was mapped to each target system through its calibrants, and jackknife+ conformal windows were computed from the calibrants, optionally normalised by a difficulty model meta-learned on other systems. On {{n_target_systems}} reversed-phase systems in leave-source-out validation (20 calibrants, nominal 90 %), Gaussian-process, smoothing-spline and ±z·RMSE windows retained the correct structure in {{main.B3_gp_PI.coverage_mean|pct}}, {{main.B2_pspline_PI.coverage_mean|pct}} and {{main.B1_gauss_rmse.coverage_mean|pct}} of cases and fixed ±1 min windows in {{fw.FIX_1min.n20.coverage_mean|pct}}, whereas conformal windows retained {{main.B5_jackknife_plus.coverage_mean|pct}}–{{main.P.coverage_mean|pct}}; this held for HILIC and for {{norman.n_systems}} interlaboratory NTS methods. Finite windows require at least $1/\alpha-1$ calibrants. Normalisation improved efficiency for longer-than-median runs (median change in interval score versus jackknife+ {{span.median_rel_change_span_above_median|pct0}}) but not for short runs, HILIC or the external systems. Calibrant mixes unlike the analytes broke validity ({{nor.fixed_mix.mix.a0.10.B5_jackknife_plus.coverage_mean|pct}} for an interlaboratory calibrant mix); replaying a known selection rule on public data restored {{sasc.common.P+SASC.coverage_mean|pct}} for common-compound calibrants.

**Significance.** Calibrant-based conformal RT windows turn the probability of keeping the correct structure from an unknown into a design choice, provided that the calibrants represent the analytes; conventional windows silently discard many correct structures on new systems. An open implementation needs only calibrant RTs and candidate SMILES.

**Keywords:** non-target screening; retention-time prediction; conformal prediction; uncertainty quantification; compound identification; LC–HRMS
