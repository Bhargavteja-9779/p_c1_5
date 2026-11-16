---
title: "Distribution-free retention-time windows for candidate filtering in LC–HRMS non-target screening on new chromatographic systems"
---

[First Author]^a^, [Second Author]^a,\*^

^a^ [Department, Institution, Address, Country]

\* Corresponding author: [name, e-mail address]

# Highlights {.unnumbered}
- Parametric and fixed RT windows lose correct structures on new LC systems
- Conformal windows from representative calibrants keep a user-chosen retention rate
- At least 1/α − 1 calibrants are needed for a finite window
- Calibrants must resemble analytes; known selection rules are corrected on average
- Tested on {{n_target_systems}} RepoRT systems and {{norman.n_systems}} interlaboratory NTS methods

# Abstract {.unnumbered}
**Background.** Retention-time (RT) filters discard candidate structures in liquid chromatography–high-resolution mass spectrometry non-target screening (NTS), but their tolerance windows are heuristic and the probability that the correct structure survives on a laboratory's own system is unknown. We aimed at RT windows that keep the correct structure at a user-chosen rate, using only a laboratory's calibrants.

**Results.** A retention-order model trained on public RepoRT data was mapped to each target system through its calibrants, which also defined jackknife+ conformal windows, optionally normalised by a meta-learned difficulty model. On {{n_target_systems}} reversed-phase systems in leave-source-out validation (20 calibrants, nominal 90 %), Gaussian-process, smoothing-spline and ±z·RMSE windows retained the correct structure in {{main.B3_gp_PI.coverage_mean|pct}}, {{main.B2_pspline_PI.coverage_mean|pct}} and {{main.B1_gauss_rmse.coverage_mean|pct}} and fixed ±1 min windows in {{fw.FIX_1min.n20.coverage_mean|pct}}, whereas conformal windows retained {{confrange.main.min|pct}}–{{confrange.main.max|pct}} (HILIC: {{confrange.hilic.min|pct}}–{{confrange.hilic.max|pct}}; {{norman.n_systems}} interlaboratory NTS methods with random calibrants: {{confrange.norrand.min|pct}}–{{confrange.norrand.max|pct}}). Finite windows require at least $1/\alpha-1$ calibrants. Normalisation improved efficiency for longer-than-median runs (median interval-score change versus jackknife+ {{span.median_rel_change_span_above_median|pct0}}) but not for shorter runs or the external methods, and reduced it for HILIC. Calibrant mixes unlike the analytes broke validity ({{nor.fixed_mix.mix.a0.10.P.coverage_mean|pct}}–{{nor.fixed_mix.mix.a0.10.B5_jackknife_plus.coverage_mean|pct}} for an interlaboratory mix); replaying a known selection rule on public data restored {{sasc.common.P+SASC.coverage_mean|pct}} for common-compound calibrants, but not for the interlaboratory mix.

**Significance.** Conformal RT windows make the probability of keeping the correct structure a design choice, provided the calibrants represent the analytes; conventional windows discard many correct structures on new systems. The implementation needs only calibrant structures and RTs, routine chromatographic metadata and candidate structures.

**Keywords:** non-target screening; retention-time prediction; conformal prediction; uncertainty quantification; compound identification; LC–HRMS
