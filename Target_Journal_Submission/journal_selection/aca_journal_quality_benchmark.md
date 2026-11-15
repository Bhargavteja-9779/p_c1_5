# Analytica Chimica Acta: journal-quality benchmark for data-treatment papers (2023–2026)

Built 2026-11-15 from PubMed E-utilities and Crossref (no ScienceDirect access). Per-paper coding: `aca_benchmark_papers.csv`; tallies: `data/aca_benchmark_raw/summary_stats.json`; pipeline: `scripts/aca_benchmark.py` (stages search → fetch → crossref → rank → coverage → context → pmc → build); manual coding sheet: `data/aca_benchmark_manual_coding.json`.

## Evidence base

* **Pool:** 469 ACA research articles (2023–2026, PubMed `journal article`, reviews/editorials excluded) matched by six theme queries (RT/QSRR 30, NTS/annotation 99, ML/DL 211, chemometrics 159, uncertainty/validation 37, transfer/domain shift 47; overlapping). Every record has a Crossref citation count (snapshot 2026-11-15) and PubMed received/accepted dates.
* **Completeness:** PubMed indexes all of ACA (yearly PubMed totals ≥ Crossref totals). Of Crossref's 400 most-cited ACA articles (2023–2026, n = 4,791), no on-theme *research article* was missing from the pool; the only on-theme misses were tutorials (listed at the end).
* **Coded set:** 43 papers (7 from 2023, 14 from 2024, 13 from 2025, 9 from 2026), covering every theme and including the nine most-cited on-theme papers. Coding is **abstract-based** (only one paper, PMID 38438240, has PMC full text), so "not determinable" means *not stated in the abstract*, not *absent from the paper*.

## 1. Novelty level

ACA's accepted novelty in these themes is mostly **"new to analytical chemistry"**, not new statistical theory. In the coded set it falls into five patterns:

1. Incremental algorithm upgrades: a filter step added to SPA, generalized Procrustes cross-validation, XYOnion data splitting.
2. New combinations or adaptations of existing methods: mid-level fusion strategies, Monte Carlo dropout, Gramian-angular-field encoding, semi-supervised domain adaptation, LLM agents.
3. Existing ML put to new analytical use: an RPLC-subspace classifier, QSRR-predicted scheduled-MRM windows, CCS prediction-interval filtering.
4. Statistical re-framing: nonparametric ANOVA, asymmetric coverage intervals, Bayesian hierarchical calibration, delta-method plus bootstrap variance for MLP calibration.
5. Pure benchmarks with no new algorithm: metabolomics processing modules, deep learning vs chemometrics, data-partitioning strategies.

Nearly every abstract claims "novel" or "first". What these papers have in common is a measurable analytical benefit, which matches the scope wording recorded in `selected_journal.md` ("quantifiable advantages over existing approaches").

## 2. Sophistication of accepted methods

Workhorse models dominate: MLR, PLS/iPLS, random forest, gradient boosting, SVM, logistic regression and shallow ANNs. Deep architectures do appear (GNN/GIN, Faster R-CNN, SE-ResNet50, LLM agents), but sophistication is not what gets a paper accepted. Two papers make this point directly:

* 42887958 found that descriptor-based gradient boosting matched a graph network within experimental reproducibility.
* 40024653 found no universally best model in low-data spectroscopy.

Statistical depth is concentrated in the metrology and uncertainty papers.

## 3. Experimental depth

* **Number of datasets:** 79% use two or more datasets or case studies, usually 2–3. Only LUMIR (41330688) uses six.
* **RT papers:** data are single-laboratory, with 535–2,003 compounds (41167895, 42887958), 604 standards (40288861) or 1,217 compounds × 51 conditions (41526108).
* **Real samples:** 63% apply the method to real samples. Data-only papers are accepted when they combine simulation with real or published data (40221201, 38677839, 42951562).
* **Sample size:** small numbers pass in applied papers. The most-cited paper (39788666) used eight steel samples.

## 4. Common baselines

Only 42% of abstracts name a competing method. The baselines that do appear are:

* the PLS family (PLS, iPLS, SPA-MLR);
* established open tools (XCMS, MZmine, MS-DIAL, Compound Discoverer, Peakonly);
* conventional practice (VIP > 1, LOESS drift correction, classical ANOVA, symmetric root-sum-of-squares uncertainty);
* single-technique models vs fused models;
* PDS/DOP for calibration transfer;
* same-instrument performance as the upper reference.

The RT papers compare only their own model variants (PLS+NN vs GNN; descriptors vs fingerprints vs GIN; condition-specific vs multi-condition models). **None benchmarks against the RT-window practice used in screening.**

## 5. How results are presented

* **Figures of merit:** R², RMSE/RMSEP/RMSECV, MAE in minutes, accuracy, MCC, recall/F1, top-k accuracy, FDR, RSD%, empirical coverage and interval width.
* **Abstracts:** quantitative, but almost always as point estimates.
* **Structure:** 67% of the coded set (64% of the pool) use structured Background / Results / Significance (and Novelty) abstracts.
* **Interpretation devices:** SHAP, Grad-CAM, stability and diagnostic plots, and method × dataset comparison tables.

## 6. Statistical evidence: typical vs missing

* **Formal inference is rare.** Only 21% (9/43) of abstracts report any (paired error test, p-values, bootstrap or permutation inference, coverage, FDR). For 60%, the statistical validation is "not determinable".
* **Reported metrics carry no uncertainty.** Confidence intervals on RMSE or accuracy are essentially absent. For example, 41167895 reports cross-validated R² 0.58 against test R² 0.85 without comment in the abstract.
* **Uncertainty without guarantees is accepted.** Monte Carlo dropout (39580162) and CCS prediction-interval filtering (37709461) were published without any coverage calibration. Coverage studies appear only in the statistics/metrology papers (40221201 and 42951562, both against a nominal 95%).
* **Full-text check (38438240):** only descriptive dispersion statistics, no tests or confidence intervals, and no data/code statement.
* **Code and data:** only 12% of abstracts mention them; computational cost is reported in 9%.

## 7. Generalization testing

* **Overall:** 44% (19/43) evaluate on external or independent data; 9% use an internal split or CV only; 47% do not say.
* **Transfer papers** test across instruments, laboratories or seasons by design.
* **RT/QSRR papers** are mostly internal. The exceptions are 38609268 (external validation set plus unseen standards) and 40480691 (a second ionization system).
* **The gap:** no ACA paper evaluates RT prediction on a *new chromatographic system calibrated with a few compounds*, with intervals whose coverage is checked. That is the gap our paper targets.

## 8. What distinguishes the most-cited papers

Rigour does not drive citations. Compared with the other 35 coded papers, the eight most-cited:

| Feature | Top 8 | Other 35 |
|---|---|---|
| External validation | 12.5% | 51% |
| Formal statistics in the abstract | 12.5% | 23% |
| Code or data statement | 0% | 14% |

Three things do drive citations:

* **Large, fast-growing application communities:** LIBS data fusion, smartphone and wearable sensors, Raman deep learning. Paper 39788666 ranks 6th of 4,791 ACA articles.
* **Practice-guiding messages or reusable tools:** the PLS-DA prefiltration rule (91 citations), Procrustes cross-validation, Raman spike removal, the modular comparison of metabolomics processing.
* **Age:** the 2023 pool has a median of 12 citations against 1 for 2026.

The two most-cited on-theme items overall are **validation-strategy tutorials** (124 and 38 citations), which shows the community wants validation guidance.

## Implications for designing our study

1. **Novelty framing.** Present the work as the first distribution-free, coverage-guaranteed RT-window method for NTS on a new system with few calibrants. Prove a quantifiable advantage over current practice, not over a straw man.
2. **Named baselines.** Compare against:
   * fixed ±Δt and ±% windows;
   * calibrant-based RT projection (linear, LOESS, GAM) with ad hoc tolerances;
   * QSRR point models (gradient boosting or random forest; optionally a GNN);
   * bootstrap, quantile-regression and Gaussian-process/Bayesian intervals;
   * MC dropout or ensembles.
3. **Datasets.** Use at least 3 public multi-system datasets (e.g., RepoRT, already in `data_raw/`) with leave-one-system-out evaluation. Add one real NTS case study with confirmed identities, since 63% of accepted papers include real samples.
4. **Leakage-proof splits.** Split by compound or scaffold and by chromatographic system, and cite the validation tutorials (37524478, 41513337) and 42838683.
5. **Core metrics.** Report empirical vs nominal coverage (80/90/95%), interval width (minutes and % of gradient), candidate-list reduction at fixed true-ID retention, and conditional coverage by chemical class and RT region.
6. **Formal statistics.** Give bootstrap CIs across systems, paired Wilcoxon tests vs baselines, and binomial checks of coverage. This puts the paper in the top fifth for statistical evidence.
7. **Ablations.** Vary the number of calibrants (e.g., 3–20), calibrant-selection strategy, base model, conformal variant (split, CV+/jackknife+, normalized, Mondrian) and calibration-set size.
8. **Robustness.** Test column-chemistry and gradient shifts, erroneous reference RTs (cf. 38438243) and exchangeability violations.
9. **Practitioner rule.** Include a "how many calibrants for a target width" guideline. Practice-guiding results earn citations in ACA.
10. **Reporting.** State the computational cost, make code and data openly available, and say so in the abstract.
11. **Limitations.** Report them explicitly: marginal vs conditional coverage, calibrant representativeness and the exchangeability assumption.
12. **Positioning.** Position against 37709461 (prediction-interval filtering without guarantees), 39580162 (uncertainty without calibration), 40221201 and 42951562 (coverage-evaluation templates), 38462330 (one-reference-sample transfer), 42702424 (cross-platform LC-MS with significance tests), 40288861 (RT reproducibility of 9–15%) and 38609268 (QSRR RT windows for suspect screening).
13. **Timeline.** The coded set has a median received→accepted time of 92 d (IQR 77–138; 25/43 ≤ 100 d). The whole pool fell to a median of 68 d for 2026 publications, against 97–104 d for 2023–2025. Statistics-heavy papers can run long (40903110 took 265 d).

## Most-cited papers in these themes (Crossref, 2026-11-15)

| # | PMID | Year | Paper (short title) | DOI | Citations (per year) | Rank in ACA 2023–26 |
|---|---|---|---|---|---|---|
| 1 | 39788666 | 2025 | LIBS + plasma acoustic emission fusion for steel classification | 10.1016/j.aca.2024.343496 | 106 (56.7) | 6 |
| 2 | 38182346 | 2024 | Univariate filtration before vs after PLS-DA in untargeted LC-MS | 10.1016/j.aca.2023.342103 | 91 (31.7) | 11 |
| 3 | 38834276 | 2024 | Wearable sweat-glucose sensor with machine learning | 10.1016/j.aca.2024.342761 | 71 (29.9) | 26 |
| 4 | 36657888 | 2023 | Closed-loop LC gradient design by Bayesian optimization | 10.1016/j.aca.2023.340789 | 47 (12.4) | 93 |
| 5 | 36641144 | 2023 | Mid-level Raman–LIBS fusion with CNN for ores | 10.1016/j.aca.2022.340772 | 45 (11.9) | 115 |
| 6 | 37455078 | 2023 | Improved successive projections algorithm (fSPA-MLR) | 10.1016/j.aca.2023.341560 | 42 (13.1) | 133 |
| 7 | 38182389 | 2024 | LC-HRMS + ML screening of unknown food contaminants | 10.1016/j.aca.2023.342116 | 37 (12.9) | 187 |
| 8 | 37032062 | 2023 | Procrustes cross-validation of multivariate regression | 10.1016/j.aca.2023.341096 | 37 (10.4) | 189 |

**Context (tutorials, PubMed type "Review"; excluded from coding):**

* 37524478, validation strategy in predictive models, Part 1 (10.1016/j.aca.2023.341532; 124 citations, ACA rank 3).
* 41513337, Part 2, recipes for avoiding overfitting (10.1016/j.aca.2025.344838; 38 citations).
* 37311606, SIMCA tutorial (10.1016/j.aca.2023.341304; 75 citations).
* 36641149, LC-MS calibration tutorial (10.1016/j.aca.2022.340711; 74 citations).

**Caveats:**

* Coding is abstract-based.
* Citation counts are an age-biased snapshot.
* Only accepted papers are visible (survivorship), so nothing here measures rejection risk.
* Theme assignment involved manual triage of broad queries.
