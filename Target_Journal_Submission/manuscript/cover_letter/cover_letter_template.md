[Date]

The Editors
*Analytica Chimica Acta*

Dear Editor,

We submit the manuscript **"Distribution-free retention-time windows for candidate filtering in LC–HRMS non-target screening on new chromatographic systems"** for consideration as a Research Article in *Analytica Chimica Acta*.

**The problem.** In suspect and non-target screening, candidate structures are filtered by comparing predicted with observed retention times. The tolerance windows used for this filter are set heuristically (fixed tolerances, multiples of an RMSE, or parametric prediction intervals), so a laboratory does not know how many correct structures its filter removes on its own chromatographic system.

**What we did.** We evaluated the window constructions in current use on {{n_target_systems}} reversed-phase systems from the public RepoRT repository in leave-source-out validation, with HILIC systems, a second base model and {{norman.n_systems}} routine NTS methods from the NORMAN interlaboratory comparison as independent checks. We propose calibrant-based conformal windows (jackknife+) for a new system, optionally normalised by a difficulty model that is meta-learned on other systems, and a selection-aware calibration for calibrant mixes chosen by a known rule. The comparators, metrics and statistical tests were specified in an analysis plan committed before the full run; later additions are labelled post hoc.

**Main findings.**

- Current windows are miscalibrated on new systems: at a nominal 90 %, Gaussian-process projection, smoothing-spline and ±z·RMSE windows retained the correct structure in {{main.B3_gp_PI.coverage_mean|pct}}, {{main.B2_pspline_PI.coverage_mean|pct}} and {{main.B1_gauss_rmse.coverage_mean|pct}} of cases, and fixed ±1 min windows in {{fw.FIX_1min.n20.coverage_mean|pct}}.
- Conformal windows, which carry a distribution-free finite-sample guarantee, retained {{main.B5_jackknife_plus.coverage_mean|pct}}–{{main.P.coverage_mean|pct}}, also for HILIC and, with randomly drawn calibrants, for the interlaboratory data. At least $1/\alpha-1$ calibrants are required, which gives laboratories a simple design rule.
- The guarantee requires calibrants that resemble the analytes. We show where this fails (an interlaboratory calibrant mix) and how a known selection rule can be corrected (from {{rob.common.P.coverage_mean|pct}} to {{sasc.common.P+SASC.coverage_mean|pct}} for commonly measured calibrants).
- The learned normalisation improves efficiency only for longer gradients and is less efficient than plain jackknife+ for HILIC and the external data; we report this as a delimitation of the method rather than a general gain.

**Fit to the journal.** *Analytica Chimica Acta* publishes data-treatment methods with quantifiable advantages for analytical practice. Our work concerns a decision made in every LC–HRMS screening workflow, quantifies its consequences at a scale not previously used for retention-time windows, and comes with code that regenerates every number, table and figure from the public data.

**Declarations.** The manuscript is original, has not been published and is not under consideration elsewhere. All authors have approved the submission. The authors declare no competing interests. The use of generative AI is declared in the manuscript. All data are public; code, configuration, raw results and the pre-specified analysis plan are provided as Supplementary Material.

**Suggested reviewers.** [To be completed by the authors: names, affiliations and e-mail addresses of independent experts without conflicts of interest.]

Yours sincerely,

[Corresponding author]
[Affiliation, address, e-mail]
