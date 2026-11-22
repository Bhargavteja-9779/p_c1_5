# PROJECT_STATE.md — persistent research checkpoint

Last updated: 2026-11-15 22:25 UTC (session 2, after container restart at 22:18 UTC)

## Selected target journal (Phase 0 — DONE)
Analytica Chimica Acta (Elsevier; ISSN 0003-2670/1873-4324). JCR 2025 (released 2026): JIF 6.1, Q1, Chemistry,
Analytical, rank 23/111. Verified received→accepted (PubMed history, articles published Nov 2025–Oct 2026): median 60 d
(all; 90.5 % ≤100 d), 72 d (ML/chemometrics/RT/uncertainty papers; 77.5 % ≤100 d). Details: `Target_Journal_Submission/journal_selection/`.

## Selected topic (Phases 1–2 — DONE)
Coverage-guaranteed retention-time (RT) windows for LC-HRMS non-target screening on new chromatographic systems from a
handful of calibrants: monotone calibrant mapping of a multi-system retention-order score + normalised jackknife+ with a
difficulty model meta-learned across systems. Record: `Target_Journal_Submission/supplementary/topic_selection.md`.
NIR conformal idea abandoned (pre-empted by the author's own Talanta 2026 paper, doi 10.1016/j.talanta.2026.130941).

## Data
RepoRT v1.0 (Zenodo 10.5281/zenodo.16267768, sha256 4eeeb46a…6091, CC BY-SA 4.0) → `data_raw/` (gitignored).
Processed: `data_processed/` (rt_long.parquet 155,485 rows; 368 RP+HILIC systems; 93,275 compounds; RDKit features).
Rebuild: `python code/scripts/build_dataset.py --repo <RepoRT root> --out data_processed`.

## Code
`Target_Journal_Submission/code/` (package `conformrt`): data loader, features, base models (XGBoost pairwise ranking;
SMRT LightGBM), monotone/P-spline/GP mappings, jackknife+ variants, difficulty model, protocol, driver, statistics,
figures. Tests: 9 passing (`python -m pytest -q tests`). Pre-registered plan: `code/EXPERIMENT_PLAN.md` (+ Addendum A).

## Experiments
* Session 1 run died at 18:07 UTC (container suspended at usage limit). Cached: fold-0 base models, fold-0 σ rows, 3 ensemble reps.
* Session 2: per-fold checkpointing added; `code/experiments/run_pipeline.sh` relaunched 22:2x UTC (main → robust → smrt → hilic → sens).
  Resume after any interruption: `cd Target_Journal_Submission/code && nohup setsid experiments/run_pipeline.sh &` (skips finished folds).
* Status file: `Target_Journal_Submission/results/experiment_logs/pipeline_status.txt` ("PIPELINE_DONE" when complete).

## Literature / benchmark
* ACA journal-quality benchmark (43 coded ACA papers) — DONE by agent: `journal_selection/aca_benchmark_papers.csv`,
  `journal_selection/aca_journal_quality_benchmark.md`.
* Literature matrix + verified BibTeX + prior-art novelty report — NOT DONE (agent died at usage limit) → re-run.

## Next steps (in order)
1. Literature matrix / prior-art workflow (agents) while pipeline computes.
2. When main results exist: `code/experiments/analyze.py` (to write) → tables, figures, statistics.json, metrics.json.
3. Manuscript (ACA format: structured abstract Background/Results/Significance; highlights; graphical abstract).
4. 10-agent × 10-cycle review/revision (workflows), fixes + reruns.
5. DOCX/PDF, cover letter, final audit, risk report, decision; commit + push.

## Update 2026-11-15 22:32 UTC
* Pipeline (main → robust/smrt/hilic/sens → fixed-window baselines (Addendum B) → small-n) running; fold-0 evaluation started 22:20.
* Literature/novelty workflow `wf_12511b4b-242` running (8 agents, concurrency cap 2). Resume: Workflow({scriptPath: <session>/workflows/scripts/rt-conformal-literature-and-novelty-wf_12511b4b-242.js, resumeFromRunId: "wf_12511b4b-242"}).
* Scheduled self-resume message at 2026-11-16 03:10 UTC (trigger trig_015c48QwnsDts8kYv2Fvjciw).
* Manuscript sections drafted: 01_introduction (to verify citations), 02_theory, 03_experimental, 90_declarations. Build chain ready
  (`code/scripts/build_manuscript.py` + `manuscript_facts.py`, pandoc + reference_aca.docx + elsevier-with-titles.csl).
* Claims to verify in the citation audit: typical RT tolerance windows (±0.5–2 min, ±2–3 SD) need sources.

## Update 2026-11-16 03:29 UTC (session 3, after usage-limit pause 22:58→03:10)
* Container was reclaimed during the pause; pipeline relaunched 03:11 (resumed at fold 2). Folds 0–2 DONE; fold 3 running.
* Literature workflow wf_12511b4b-242: 2/8 agents completed (rt_models, nts_practice sweeps; 40 verified refs) — 6 failed at the
  usage limit and are NOT re-run on purpose (token budget protects the compute job); their role (adversarial novelty search,
  citation checks) is covered by reviewers R1/R8 in every review cycle.
* Done this session: verified bibliography (manuscript/references.bib; supplementary/citation_verification.csv, 55 Crossref/arXiv
  + Holm manual), claim-source verification (supplementary/claim_source_verification.md), literature matrix (44), prior-art map,
  NORMAN external validation script (Addendum C) appended to run_pipeline.sh, analysis extended (fixed windows, small n, NORMAN),
  Introduction revised (Ferrando Plo & Kruve 2026 CCS conformal acknowledged).
* Scheduled resumes: trig_015c48QwnsDts8kYv2Fvjciw (fired 03:10), trig_018aBvaNf3AJEmPCc4NykmHv (08:15 UTC).
* NEXT: when pipeline_status.txt shows PIPELINE_DONE → `cd code && python3 experiments/analyze.py --raw ../results/raw --out ..`
  → write 04_results.md (tokens only) + highlights/abstract → `python3 scripts/build_manuscript.py` → review cycle 1 via
  Workflow({scriptPath: REVIEW/workflow/review_cycle.js, args: {cycle: 1}}).

## Update 2026-11-16 04:20 UTC
* Main experiment DONE (228 targets). Results analysed; Results/front matter/conclusions drafted (token-based, manuscript/src).
  Key: conformal windows valid (P 0.911), parametric under-cover (GP 0.822); P best interval score but gain heterogeneous
  (long gradients −18 %, short ≈0/worse vs Gaussian); conditional coverage gains outside calibrant range; decoy trade-off equal at matched retention.
* Robustness: span → conformal conservative (0.93); COMMON calibrants → all under-cover (P 0.815, J+ 0.746) [mechanism: 27 % of tests
  outside calibrant range + easy calibrants] → Phase-19 improvement SASC (EXPERIMENT_PLAN Addendum D) implemented + queued.
* Pipeline order now: robust(noise, outlier) → smrt → hilic → sens → fixed windows → small-n → NORMAN → SASC → NORMAN(+SASC) → compute benchmark.
* NEXT after PIPELINE_DONE: analyze.py → write robustness/generalisation/SASC/NORMAN/compute subsections + limitations →
  build_manuscript.py → review cycle 1 (Workflow REVIEW/workflow/review_cycle.js, args {cycle: 1}).
* TODO after PIPELINE_DONE: remove the TEMPORARY raw-results block from .gitignore and git add -f results/raw/*_rows.parquet, *_decoys.parquet, *_timing.parquet, sasc_*.csv/json

## Update 2026-11-16 05:32 UTC — NORMAN external validation done (37 systems)
* Exchangeable calibrants (random_all n=20): conformal valid (P 0.912, J+ 0.909); GP 0.839, Gaussian 0.881 — replicates main finding.
  Efficiency: P NOT better than J+ (interval score 11.15 vs 10.70) — learned normalisation gives no gain externally.
* Fixed lab mix (calibrants = NORMAN chemicals, test = suspects): all under-cover (P 0.851, J+ 0.869, Gauss 0.851, GP 0.777).
* FRAMING DECISION PENDING SASC: core contributions = failure of current windows (large scale) + calibrant-based conformal
  windows + n >= 1/alpha - 1 rule + non-exchangeability diagnosis/correction (SASC); learned normalisation = context-dependent
  enhancement (long RP gradients). Check: NORMAN systems have missing temp/flow/t0 descriptors (never missing in training).

## Update 2026-11-16 06:35 UTC — manuscript v1 + Supplementary Material built; audits running
* PIPELINE_DONE 06:05:56Z; analysis complete (re-run 06:17 after NORMAN-strata fix). All results/tables/figures committed.
* Manuscript v1: `manuscript/Final_Manuscript.{md,docx,pdf}` (32 pp; numbered sections, line + page numbers, portable math).
  Supplementary: `supplementary/Supplementary_Material.{md,docx,pdf}` (S1–S6, Tables S1–S26, Figs. S1–S2), built by
  `code/scripts/build_supplementary.py` (table_index.csv maps S-numbers to analysis outputs). Both builders are token-based.
* Phase 22–23 audits: Workflow run wf_0e2f7cad-fff (claim–evidence auditor + citation auditor) — results → REVIEW/claim_evidence_audit.md, REVIEW/citation_audit.md.
* NEXT: fix audit findings → rebuild → review cycle 1: Workflow({scriptPath: "Target_Journal_Submission/REVIEW/workflow/review_cycle.js" (absolute path), args: {cycle: 1}})
  → write REVIEW/review_01.md (reviews + consolidated actions + responses) → implement → cycles 2–10 → final_audit.md, cover letter, risk report, decision.
* Scheduled self-resume: trig_018aBvaNf3AJEmPCc4NykmHv at 08:15 UTC (re-schedule if work continues past the next usage-limit pause).

## Update 2026-11-16 13:55 UTC — audits resolved; comparator B7 added; review cycle 1 running
* Phase 22–23 audits (wf_0e2f7cad-fff, rerun after the 06:38 usage-limit failure): 106 claims audited (51 corrected), 57→80
  references (80/80 verified). Records: REVIEW/claim_evidence_audit.md, REVIEW/citation_audit.md (with resolution tables).
* Analysis now follows the pre-specified plan (|gap| test, one-sided coverage test); Table 1 RT span fixed; fold-adjusted
  gain–span correlation; SASC width statistic aligned; NORMAN strata per scenario.
* Addendum E (post hoc, pre-run commit 34ea26f, bound amendment 38de8b5): comparator B7 = GP projection with meta-learned
  hyper-parameters → coverage 0.973 at nominal 0.90, 2.6× wider windows, 3 % of systems < 0.85 (identical draws verified).
* Clean-clone quick reproduction PASSED (1238 s). Graphical abstract and cover letter (token-based) added.
* Review cycle 1: Workflow run wf_42a67b1d-ef3 (started 13:31; 10 reviewers + consolidation + adversarial verification).
  When done: `python REVIEW/workflow/journal_to_md.py review <transcript dir of wf_42a67b1d-ef3> 1 REVIEW/review_01.md`,
  implement the actions, rebuild, write responses into review_01.md, commit, then cycle 2 (args {cycle: 2}).
* Resume after a usage-limit pause: re-run any failed workflow with resumeFromRunId (completed agents are cached).

## Novelty revision (started 2026-11-21 20:50 UTC, session 4) — user: desk rejections everywhere → more novelty
* Journal-level scan + prior art: supplementary/journal_level_and_prior_art_2026-11-21.md. Accepted RT/NTS papers show
  identification gains on real candidate lists and/or real samples; borrowing calibration across systems (B2) has no prior art
  in chemistry; weighted conformal (B1) and conformal FDR for annotation (B3) also open.
* NEW CORE CONTRIBUTION (Addendum G in code/EXPERIMENT_PLAN.md): cross-system conformal calibration (XSC) →
  laboratory-level (ε, δ) retention guarantee + finite windows from 3–5 calibrants. Code: src/conformrt/evaluation/crosssys.py,
  experiments/novelty/run_xsc.py. Cache (base models) re-created by experiments/novelty/prewarm.py (cache is not in git).
* Running: `nohup setsid code/experiments/novelty/run_novelty_pipeline.sh &` (prewarm → XSC; resumable). Status:
  results/experiment_logs/novelty_status.txt (XSC_STAGE_DONE), logs prewarm.stdout, run_xsc_full.stdout.
* NEXT: (1) when XSC done: analyse (lab failure rates, widths, n=3/5, rules) → decide framing; (2) NORMAN external XSC
  (calibrate on all RepoRT, apply to 37 NORMAN methods, random + lab mix); (3) real candidate lists (PubChemLite isomers)
  → candidates removed / true-structure rank at guaranteed retention vs ±k·RMSE; (4) rewrite manuscript around XSC
  (new title/abstract/highlights), demote meta-learned J+ to comparator; implement remaining review_01 actions;
  (5) rebuild DOCX/PDF, cover letter; commit + push. Self-resume: scheduled task trig_012DjbqH1ygJ8njjao5anxKe (every 2 h).

### Update 2026-11-21 22:45 UTC
* Fold 0 XSC result: random n=20 → J+ 33 % of labs < 0.90, XSC_sigs_lab 5 % (target ≤ 10 %), width 3.1 vs 3.9 min (lab-calibrated J+);
  n=3 finite & valid with shrinkage; polar calibrants: J+ fails 91 %, XSC_sigs_lab 0 %.
* Done: new 01_introduction, 02_theory (XSC section), 03_experimental (XSC/candidates/NORMAN), supplementary/S7 proof,
  manuscript_facts.xsc_facts (tokens x.*, xa*.*, xc.*, xn.*; tables xsc_table_*), analyze_xsc.py, refs (PubChemLite, park2022pac).
* Background chain: run_novelty_pipeline.sh (XSC) → run_stage2.sh (candidates → NORMAN XSC). If killed: re-launch both with
  `nohup setsid code/experiments/novelty/run_novelty_pipeline.sh &` and `... run_stage2.sh &` (both resumable).
* NEXT after NORMAN_STAGE_DONE: `cd code && python3 experiments/novelty/analyze_xsc.py` → write 04_results.md (new structure:
  4.1 marginal guarantees leave labs below target; 4.2 XSC lab-level per source; 4.3 few calibrants; 4.4 selection rules;
  4.5 PubChemLite candidates; 4.6 NORMAN; 4.7 old comparisons condensed → SI) + 00_front (title/abstract/highlights),
  05_conclusions, fig1 workflow redraw, cover letter; build_manuscript.py; commit/push.

## RESUME HERE (2026-11-22 03:05 UTC) — handover to a new session/account
User instructions (keep): first author **P N Bhargav Teja**, second author placeholder; **no funding, no competing interests**;
**code available on request**; maximise novelty; real-world framing welcome; work fully autonomously.
### Done
* XSC (cross-system conformal calibration, laboratory-level guarantee) run on all 5 folds → results/raw/xsc/fold*_eval.parquet,
  fold*_multipliers.csv, fold*_sim.pkl (replays); analysis: `cd code && python3 experiments/novelty/analyze_xsc.py`
  → tables/xsc_*.md, results/xsc_metrics.json, figures/figX1_lab_level, figX2_selection_rules.
  KEY: random n=20: lab failure (<0.90) J+ 34.9 %, meta-σ J+ 33.6 %, Gaussian 48.0 %; XSC lab-level 8.0 % [7.1–9.0], width 5.70 min
  vs lab-calibrated J+ 7.01 (9.2 %) and flat XSC 6.81; ≤10 % in nearly all rule×n cells (polar n=20 11.4 %); n=3 finite (7.8 %).
* Manuscript src rewritten: 01_introduction, 02_theory (XSC section), 03_experimental, 90_declarations (authors/funding/code).
  STILL OLD: 00_front (title/abstract/highlights), 04_results, 05_conclusions, fig1 workflow → must be rewritten.
* supplementary/S7_laboratory_level_guarantee.md (proof); supplementary/journal_level_and_prior_art_2026-11-21.md.
### Remaining (in order)
1. Environment in a NEW container: `cd Target_Journal_Submission/code && pip install -r requirements.txt && mkdir -p data &&
   ln -sfn /home/user/p_c1_5/data_processed data/processed && python3 scripts/download_data.py` (RepoRT+NORMAN), then
   PubChemLite: `mkdir -p data/raw/pubchemlite && curl -L -o data/raw/pubchemlite/PubChemLite_exposomics_20260925.csv
   "https://zenodo.org/records/22953038/files/PubChemLite_exposomics_20260925.csv?download=1"` (sha256 6a258459…89cd).
2. Rebuild model cache (~1 h): `nohup setsid python3 experiments/novelty/prewarm.py config/full.yaml > ../results/experiment_logs/prewarm.stdout 2>&1 &`
3. Then (each resumable, skip finished folds): `python3 experiments/novelty/run_xsc_pivots.py config/full.yaml` (fold0 done),
   `run_candidates.py` (PubChemLite candidate lists), `run_norman_xsc.py` (NORMAN external), `run_norman_case.py`
   (REAL-WORLD case study: 37 NORMAN labs, own calibration mix, suspects, PubChemLite candidates). Write "XSC_STAGE_DONE" into
   results/experiment_logs/novelty_status.txt first if using run_stage2.sh/run_stage3.sh chains.
4. `python3 experiments/novelty/analyze_xsc.py` and `analyze_xsc_offline.py` (lab-calibration curve figX3, localized XSC K=60).
5. Write 04_results (4.1 average guarantees leave a third of labs below target; 4.2 XSC lab-level per source; 4.3 few calibrants;
   4.4 selection rules; 4.5 PubChemLite candidates; 4.6 NORMAN real-world case study; 4.7 localized XSC/limitations), new title
   ("Laboratory-level guarantees for retention-time filtering in non-target screening by cross-system conformal calibration"),
   abstract/highlights, conclusions, redraw fig1; move old main comparison to SI; `python3 scripts/build_manuscript.py`;
   rebuild supplementary + cover letter (no funding/competing interests; code on request); commit + push.
6. Scheduled task trig_012DjbqH1ygJ8njjao5anxKe belongs to the old account/session — create a new one if needed.

### Update 2026-11-22 06:57 UTC
* Worker restarts kill background jobs (~every few hours). Single resumable chain: `cd Target_Journal_Submission/code &&
  nohup setsid experiments/novelty/run_chain.sh >/dev/null 2>&1 &` (NORMAN XSC → NORMAN case → candidates → pivots; status
  in results/experiment_logs/novelty_status.txt). Relaunch after every restart.
* DONE since last update: 04_results (new), 00_front (title/highlights/abstract), 05_conclusions, Fig. 1 redrawn, Supplementary
  (S7 proof, S8 tables S27–S31, Fig. S3), cover letter template, offline analysis (calibration curve: observed ≤ nominal δ).
* Pending tokens (need chain outputs): xc.* (PubChemLite), xk.* (NORMAN case), xn.* (NORMAN XSC). After chain: run
  analyze_xsc.py + analyze_xsc_offline.py, VERIFY claims in 04_results §"Filtering real candidate lists" and §"Real-world case"
  against numbers (text was written before the data; rewrite if not supported), then build_manuscript.py, build_supplementary.py,
  build_cover_letter.py; adversarial review; commit/push.

## Status 2026-11-22 08:30 UTC — novelty revision COMPLETE (submission-ready apart from author details)
* Manuscript (34 pp), Supplementary Material (31 tables, S7 proof, Fig. S3), cover letter and graphical abstract rebuilt from
  machine-generated results; all claims checked against tables after an internal hostile review (fixes applied).
* Author-side items left: second author name, affiliations, e-mail; date in cover letter; optionally trim body (8.3k words).
* Honest limitations stated in the text: shrunk-σ XSC failed on NORMAN own mixes (43 %; masked-metadata replay did not fix it
  → σ ≡ 1 recommended for fixed mixes); polar rule n=20 11.4 % and one source group 10.8 % slightly above δ.
