---
title: "Laboratory-level guarantees for retention-time filtering in non-target screening by cross-system conformal calibration"
---

P N Bhargav Teja^a,\*^, [Second Author]^b^

^a^ [Department, Institution, City, Country]

^b^ [Department, Institution, City, Country]

\* Corresponding author: P N Bhargav Teja, [e-mail address]

# Highlights {.unnumbered}
- Conformal RT windows keep 90 % only on average; a third of laboratories fall below
- Replaying a laboratory's calibrant protocol on public LC systems calibrates its window
- Guarantee: ≥90 % of correct structures retained in ≥90 % of laboratories
- Finite guaranteed windows from three calibrants and for non-random calibrant mixes
- Validated on 228 systems, PubChemLite candidates and 37 NORMAN laboratories

# Abstract {.unnumbered}
**Background.** Retention-time (RT) filters discard candidate structures in LC–HRMS non-target screening (NTS), but their tolerance windows are heuristic. Even conformal windows, which guarantee a retention rate, do so only on average over laboratories: an individual laboratory cannot know whether its own filter keeps the correct structures. We aimed at RT windows with a guarantee for the laboratory itself, from a handful of calibrants.

**Results.** We introduce cross-system conformal calibration (XSC). The laboratory's calibrants place and scale the window, with the scale shrunk towards a meta-learned difficulty model, and the window width is calibrated by replaying the laboratory's calibrant protocol on hundreds of public chromatographic systems. On {{x.n_systems}} reversed-phase systems in leave-source-out validation with 20 calibrants, jackknife+ windows met the average target but {{x.random.20.B5_jackknife_plus.fail|pct0}} of laboratories retained fewer than 90 % of correct structures ({{x.random.20.B1_gauss_rmse.fail|pct0}} with ±z·RMSE windows). With XSC, {{x.random.20.XSC_sigs_lab.fail|pct}} did (target ≤ 10 %), with windows {{x.ratio_XSC_sigs_lab_vs_JP_lab.median|.2f}} times as wide as equally guaranteed jackknife+ windows. The control held with three calibrants, for RT-spanning, common or polar calibrant mixes, and in 37 NORMAN interlaboratory methods using their own calibration chemicals, where XSC retained {{xk.XSC_sigs_lab.kept|pct0}} of correct suspect structures and eliminated {{xk.XSC_sigs_lab.elim|pct0}} of PubChemLite isomers.

**Significance.** XSC turns RT filtering into a procedure with a stated, laboratory-level error rate: at least 90 % of correct structures retained, with 90 % confidence, for the laboratory's actual calibrants. It needs only calibrant structures and RTs, routine metadata and public data, and it shows before measurement which calibrant designs are informative.

**Keywords:** non-target screening; retention-time prediction; conformal prediction; uncertainty quantification; compound identification; LC–HRMS

