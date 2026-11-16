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

## Update 2026-11-16 07:00 UTC — manuscript v1 + Supplementary Material built; audits running
* PIPELINE_DONE 06:05:56Z; analysis complete (re-run 06:40 after NORMAN-strata fix). All results/tables/figures committed.
* Manuscript v1: `manuscript/Final_Manuscript.{md,docx,pdf}` (32 pp; numbered sections, line + page numbers, portable math).
  Supplementary: `supplementary/Supplementary_Material.{md,docx,pdf}` (S1–S6, Tables S1–S26, Figs. S1–S2), built by
  `code/scripts/build_supplementary.py` (table_index.csv maps S-numbers to analysis outputs). Both builders are token-based.
* Phase 22–23 audits: Workflow run wf_0e2f7cad-fff (claim–evidence auditor + citation auditor) — results → REVIEW/claim_evidence_audit.md, REVIEW/citation_audit.md.
* NEXT: fix audit findings → rebuild → review cycle 1: Workflow({scriptPath: "Target_Journal_Submission/REVIEW/workflow/review_cycle.js" (absolute path), args: {cycle: 1}})
  → write REVIEW/review_01.md (reviews + consolidated actions + responses) → implement → cycles 2–10 → final_audit.md, cover letter, risk report, decision.
* Scheduled self-resume: trig_018aBvaNf3AJEmPCc4NykmHv at 08:15 UTC (re-schedule if work continues past the next usage-limit pause).
