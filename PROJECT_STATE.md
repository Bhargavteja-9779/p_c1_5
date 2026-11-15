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
