# REPRODUCIBILITY_STATUS

Created 2026-11-15.

## Status 2026-11-16 07:00 UTC
* One command: `cd Target_Journal_Submission/code && python run_all.py --mode full` (data → tests → experiments → analysis →
  manuscript + supplementary). Quick check: `--mode quick` (~5 min, separate output dir).
* Data: RepoRT v1.0 (Zenodo 10.5281/zenodo.16267768) and NORMAN ILC SI archive downloaded with SHA-256 verification
  (`code/scripts/download_data.py`).
* Determinism: fixed seeds from config; per-system seeds via zlib.crc32 (not Python hash); deterministic figure PDFs.
  Raw per-draw results committed (`results/raw/*_rows.parquet` etc.); `results/raw_results.csv` exported.
* Unit tests: 9/9 pass (`python -m pytest -q tests`), incl. Monte-Carlo coverage checks of jackknife+ and normalised jackknife+.
* Environment: Python 3.13.18, pinned requirements.txt / environment.yml; pandoc 3.1.3; LibreOffice 24.2 (+ Math) for PDF.
* Not yet done: clean-room re-run of `run_all.py --mode quick` from a fresh clone (planned for the final reproducibility audit).
