[Date]

The Editors
*Analytica Chimica Acta*

Dear Editor,

We submit the manuscript **"Laboratory-level guarantees for retention-time filtering in non-target screening by cross-system conformal calibration"** for consideration as a Research Article in *Analytica Chimica Acta*.

**The problem.** In suspect and non-target screening, candidate structures are filtered by comparing predicted with observed retention times. The tolerance windows used for this filter are heuristic, so a laboratory does not know how many correct structures its filter removes. Recent conformal windows guarantee a retention rate only on average over laboratories. We show that with 20 calibrants, {{x.random.20.B5_jackknife_plus.fail|pct0}} of laboratories using such windows retain fewer than 90 % of the correct structures, and that they cannot tell whether they are among them.

**What is new.** We introduce cross-system conformal calibration (XSC). The laboratory's own calibrant protocol is replayed on hundreds of public chromatographic systems, and systems, rather than compounds, serve as the calibration units. This gives a guarantee for the individual laboratory, which to our knowledge has not been available for retention-based filtering: at least 90 % of the correct structures are retained with 90 % confidence, for the laboratory's actual calibrant set. The guarantee holds with as few as three calibrants and for calibrant mixes that do not resemble the analytes. A proof is given in the Supplementary Material.

**Evidence.**

- In leave-source-out validation on {{x.n_systems}} reversed-phase systems, {{x.random.20.XSC_sigs_lab.fail|pct}} of laboratories fell below 90 % retention (target ≤ 10 %). The windows were narrower than equally guaranteed alternatives (median ratio {{x.ratio_XSC_sigs_lab_vs_JP_lab.median|.2f}}).
- Control was at or close to the target with 3–40 calibrants, with RT-spanning, commonly measured and polar calibrant mixes, and for each held-out source group ({{x.src.min_fail_proposed|pct}}–{{x.src.max_fail_proposed|pct}} failing).
- External validation on the 37 laboratories of the NORMAN interlaboratory comparison, with no NORMAN data used for calibration: {{xn.random_all.20.XSC_sigs_lab.fail|pct0}} of laboratories failed with random calibrants.
- Retrospective case study in the same 37 laboratories, each using its own calibration mix and screening its suspects against all PubChemLite isomers. Current RT windows silently discarded the correct structure in many laboratories ({{xk.RMSE_x2.labfail|pct0}}–{{xk.B1_gauss_rmse.labfail|pct0}} retained fewer than 90 %). XSC (σ ≡ 1) retained {{xk.XSC_flat_lab.kept|pct}} and gave the most correct top-ranked identifications. The case study also shows how much filtering such a mix can safely support. The meta-learned variant did not transfer to these methods, and we report this.

**Why Analytica Chimica Acta.** The journal welcomes data-treatment methods that demonstrate new principles with clear advantages over existing approaches and proper statistical treatment. The method gives a quantifiable advantage: a stated laboratory-level error rate where current practice has none. It uses only data a laboratory already has (calibrant structures and RTs, routine metadata), and it lets laboratories compare calibrant designs before measuring them.

The manuscript has not been published and is not under consideration elsewhere. A related article of ours on conformal calibration of near-infrared spectra under deployment shift (Talanta, 2026) is cited; it does not address chromatographic transfer, laboratory-level guarantees or candidate filtering. All data used are public; the code is available on request. The authors declare no competing interests, and the work received no specific funding.

Yours sincerely,

P N Bhargav Teja (corresponding author), on behalf of all authors
