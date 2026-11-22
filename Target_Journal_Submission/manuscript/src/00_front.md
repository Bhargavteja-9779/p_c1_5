---
title: "Laboratory-level guarantees for retention-time filtering in non-target screening by cross-system conformal calibration"
---

P N Bhargav Teja^a,\*^, [Second Author]^b^

^a^ [Department, Institution, City, Country]

^b^ [Department, Institution, City, Country]

\* Corresponding author: P N Bhargav Teja, [e-mail address]

# Highlights {.unnumbered}
- Conformal RT windows hold 90 % on average, yet a third of laboratories fall short
- Replaying a lab's calibrant protocol on public LC systems calibrates its window
- {{x.random.20.XSC_sigs_lab.fail|pct0}} of labs below 90 % retention in unseen source groups (target ≤ 10 %)
- Finite guaranteed windows from three calibrants and for non-random calibrant mixes
- In 37 NORMAN laboratories, current RT filters silently lose correct structures

# Abstract {.unnumbered}
**Background.** Retention-time (RT) filters discard candidate structures in LC–HRMS non-target screening (NTS), but their tolerance windows are heuristic. Even conformal windows, which guarantee a retention rate, do so only on average over laboratories: an individual laboratory cannot know whether its own filter keeps the correct structures. We aimed at RT windows with a guarantee for the laboratory itself, from a handful of calibrants.

**Results.** We introduce cross-system conformal calibration (XSC). The laboratory's calibrants place and scale the window, with the scale shrunk towards a meta-learned difficulty model. The window width is calibrated by replaying the laboratory's calibrant protocol on hundreds of public chromatographic systems. On 228 reversed-phase systems in leave-source-out validation with 20 calibrants, jackknife+ met the average target, but {{x.random.20.B5_jackknife_plus.fail|pct0}} of laboratories retained fewer than 90 % of correct structures. With XSC, {{x.random.20.XSC_sigs_lab.fail|pct}} did (target ≤ 10 %), with windows {{x.ratio_XSC_sigs_lab_vs_JP_lab.median|.2f}} times as wide as equally guaranteed jackknife+ windows. Control held with three calibrants, for non-random calibrant mixes and in 37 external NORMAN laboratories ({{xn.random_all.20.XSC_sigs_lab.fail|pct0}} failing). In a retrospective case study, these laboratories screened suspects against PubChemLite isomers using their own calibration mixes. With current windows, {{xk.RMSE_x2.labfail|pct0}}–{{xk.B1_gauss_rmse.labfail|pct0}} of laboratories lost more than 10 % of correct structures. XSC without the meta-learned scale retained {{xk.XSC_flat_lab.kept|pct0}} ({{xk.XSC_flat_lab.labfail|pct0}} of laboratories below 90 %); with it, {{xk.XSC_sigs_lab.labfail|pct0}} fell short.

**Significance.** XSC turns RT filtering into a procedure with a stated, laboratory-level error rate: at least 90 % of correct structures retained, with 90 % confidence, for the laboratory's actual calibrants. It needs only calibrant structures and RTs, routine metadata and public data, and it shows before measurement which calibrant designs are informative.

**Keywords:** non-target screening; retention-time prediction; conformal prediction; uncertainty quantification; compound identification; LC–HRMS

