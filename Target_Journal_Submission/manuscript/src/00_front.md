---
title: "Distribution-free retention-time windows for candidate filtering in LC–HRMS non-target screening on new chromatographic systems"
---

[First Author]^a^, [Second Author]^a,\*^

^a^ [Department, Institution, Address, Country]

\* Corresponding author: [name, e-mail address]

# Highlights

- Parametric retention-time windows lose correct structures on new LC systems
- Calibrant-based conformal windows keep a user-chosen retention rate
- At least 1/α − 1 calibrants are needed for a finite window
- A difficulty model learned on other systems narrows windows where it matters
- Validated on {{n_target_systems}} RepoRT systems and {{norman.n_systems}} interlaboratory NTS methods

# Abstract

**Background.** Retention-time (RT) filters discard candidate structures in liquid chromatography–high-resolution mass spectrometry (LC–HRMS) non-target screening, but their tolerance windows are set heuristically, and the probability that the correct structure survives on a laboratory's own chromatographic system is unknown. We aimed to construct RT windows that retain the correct structure at a user-chosen rate, using only the few calibrants that a laboratory measures.

**Results.** A retention-order model trained on public RepoRT data was mapped to the target system through the calibrants, and jackknife+ conformal windows were computed with nonconformity scores normalised by a difficulty model meta-learned on other systems. Across {{n_target_systems}} reversed-phase systems in leave-source-out validation (20 calibrants, nominal 90 %), Gaussian-process, smoothing-spline and ±z·RMSE windows retained the correct structure in {{main.B3_gp_PI.coverage_mean|pct}}, {{main.B2_pspline_PI.coverage_mean|pct}} and {{main.B1_gauss_rmse.coverage_mean|pct}} of cases, whereas the proposed windows retained {{main.P.coverage_mean|pct}} with the lowest interval score; the median gain over plain jackknife+ reached {{het.B5_jackknife_plus.Folberth.median_rel_change|pct0}} for long gradients but was negligible for short ones. Finite windows require at least $1/\alpha-1$ calibrants. For compounds outside the calibrant range, coverage rose from {{strata.region.outside_cal_range.B5_jackknife_plus|pct0}} to {{strata.region.outside_cal_range.P|pct0}}. {{abstract_external_sentence}} At matched retention, all windows eliminated similar fractions of isomeric decoys.

**Significance.** Calibrant-based conformal RT windows let laboratories fix, before the analysis, the probability of keeping the correct structure on their own system, independently of the RT model used; an open implementation requires only calibrant RTs and candidate SMILES.

**Keywords:** non-target screening; retention-time prediction; conformal prediction; uncertainty quantification; compound identification; LC–HRMS
