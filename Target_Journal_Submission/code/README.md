# conformrt — coverage-guaranteed retention-time windows for LC–HRMS non-target screening

Code, configuration and tests for the study *"Distribution-free retention-time windows for candidate filtering on new
chromatographic systems from a handful of calibrants"* (submitted to *Analytica Chimica Acta*).

The method transfers a multi-system retention-order model to a new chromatographic system using only n measured
calibrants, and returns, for any candidate structure, an RT window that contains the compound's true RT with probability
≥ 1 − 2α (≈ 1 − α in practice) when calibrants are randomly drawn from the analysed compounds (normalised jackknife+;
proof in `../supplementary/S2_coverage_proof.md`). Window widths adapt to the compound and to the calibrant geometry via a
difficulty model meta-learned on other systems.

## 1. Installation
```bash
python -m venv .venv && source .venv/bin/activate      # Python 3.13 tested (3.10+ should work)
pip install -r requirements.txt                         # or: conda env create -f environment.yml
```
Hardware: CPU only (tested on 4 cores, 15 GB RAM). `conformrt.utils.detect_hardware()` reports CUDA/MPS availability;
XGBoost can use CUDA if present, everything else runs on CPU (Apple Silicon falls back to CPU).

## 2. Data
RepoRT v1.0 (Kretschmer et al., *Nat. Methods* 2024; Zenodo 10.5281/zenodo.16267768; CC BY-SA 4.0).
```bash
python scripts/download_data.py --dest data/raw          # downloads 176.7 MB, verifies SHA-256 4eeeb46a…6091
python scripts/build_dataset.py --repo data/raw/repoRT/f-kretschmer-RepoRT-f4a2171 --out data/processed
```
`build_dataset.py` standardises compound identity (InChIKey connectivity block), removes inconsistent stereoisomer
entries, derives mobile-phase pH class / organic modifier, and computes RDKit descriptors, Morgan fingerprints and
ionisable-group counts for all 93,275 compounds (~4 min on 4 cores).

## 3. Reproduce everything
```bash
python run_all.py --mode quick    # pipeline check: tiny models, 2 systems per fold, ~5 min -> ../results_quick
python run_all.py --mode full     # full pre-registered study -> ../results, ../tables, ../figures (several CPU-hours)
```
Stages (`--stages data,tests,experiments,analysis,manuscript`): download/build → unit tests → experiments
(`experiments/run_main.py`, `experiments/run_followups.py`; resumable per fold) → `experiments/analyze.py` (all tables,
figures, `results/metrics.json`, `results/statistics.json`, `results/raw_results.csv`) → manuscript and Supplementary
Material build (`scripts/build_manuscript.py`, `scripts/build_supplementary.py`; every number is inserted from the result
files and unresolved tokens abort the build). The document stage additionally needs pandoc (≥ 3.1; tested 3.1.3) and,
for the PDF export, LibreOffice with its Math component (tested 24.2; on Debian/Ubuntu: `apt-get install
libreoffice-writer libreoffice-math`); without LibreOffice the DOCX files are still produced.
If a run is interrupted, re-run the same command: trained base models (`cache/`) and finished folds
(`../results/raw/checkpoints/`) are reused.

## 4. Configuration
`config/full.yaml` — the pre-registered settings (see `EXPERIMENT_PLAN.md`): target systems = RP systems with ≥60
compounds (SMRT excluded), leave-source-out folds, n ∈ {10, 15, 20, 30, 40} calibrants, 20 random draws, α ∈ {0.05, 0.10, 0.20},
base model (XGBoost rank:pairwise, 600 trees), difficulty model (LightGBM, 400 trees), seeds fixed.
`config/quick.yaml` — reduced settings for pipeline verification only.

## 5. Layout
```
src/conformrt/data/repo_rt.py            RepoRT loading and cleaning
src/conformrt/preprocessing/features.py  RDKit features, ionisable-group SMARTS, Tanimoto similarity
src/conformrt/models/base_models.py      retention-order ranking model; SMRT regressor
src/conformrt/models/mapping.py          monotone (isotonic+PCHIP), P-spline (+PI), GP (+PI), linear maps
src/conformrt/models/difficulty.py       meta-learned difficulty model and calibrant-geometry features
src/conformrt/evaluation/intervals.py    jackknife+, normalised jackknife+, split conformal, Gaussian, Winkler score
src/conformrt/evaluation/protocol.py     leave-source-out folds, inner cross-fitting, calibrant-protocol simulation
src/conformrt/evaluation/system_eval.py  evaluation of all methods on one target system (coverage, width, strata, decoys)
src/conformrt/evaluation/driver.py       experiment orchestration with caching/checkpointing
src/conformrt/statistics/analysis.py     per-system aggregation, Wilcoxon/Holm, rank-biserial, bootstrap CIs
src/conformrt/visualization/figures.py   publication figures
experiments/                             run_main.py, run_followups.py, analyze.py, smoke_test.py, run_pipeline.sh
tests/test_core.py                       unit tests incl. Monte-Carlo coverage checks
```

## 6. Using the method on your own system
Provide (i) SMILES and measured RTs of ≥ 1/α − 1 calibrants (e.g. ≥ 19 for α = 0.05), chosen *at random* from the kind
of compounds you analyse if you want the formal guarantee, and (ii) candidate SMILES. The window for each candidate is
computed with `jackknife_plus(...)` using leave-one-out monotone maps and the difficulty model trained on RepoRT
(see `evaluation/system_eval.py` for the exact sequence of calls).

## 7. License
Code: MIT. Data derived from RepoRT remain under CC BY-SA 4.0.
