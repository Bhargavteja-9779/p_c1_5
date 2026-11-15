# Phases 1–2: candidate research directions, stress test and topic selection

Constraints that shaped the search: (i) target journal fixed first = Analytica Chimica Acta (chemometrics / data-treatment
scope, verified Q1, verified fast acceptance); (ii) compute = 4 CPU cores, 15 GB RAM, no GPU; (iii) only public data;
(iv) contribution must be a new principle in data treatment with "obvious and quantifiable advantages" (ACA scope).

## Candidate directions considered (19)

| # | Direction | Data (public, CPU-feasible?) | Novelty status after search | Main risk | Verdict |
|---|---|---|---|---|---|
| 1 | **Coverage-guaranteed retention-time (RT) windows for LC-HRMS non-target screening (NTS) on new chromatographic systems from few calibrants** | RepoRT v1.0 (373 systems, 155k RTs, CC BY-SA) – yes | No conformal/coverage-guaranteed RT work found in PubMed (query "conformal AND retention time" = 0 hits, 2026-11-15) or preprint searches; closest: probabilistic GP projections (García 2022), accuracy of projection vs prediction (Malm & Kruve 2025), uncertainty-aware multi-condition RT (Beck 2026) | Could be viewed as "applying conformal prediction"; must show quantifiable advantage | **SELECTED** |
| 2 | Conformal prediction intervals for NIR multivariate calibration under season/instrument shift | Corn, IDRC-2002 tablets, mango DMC, OSSL – yes | **Pre-empted**: "Group-conformal calibration of NIR prediction intervals for new seasons, instruments and regions using PLS diagnostics", Talanta 2026 (doi 10.1016/j.talanta.2026.130941) — by the same author as this project; duplicating it would be salami publication | — | Abandoned |
| 3 | Conformal calibration transfer with few transfer standards (NIR) | small sets (80–655 samples) | Overlaps strongly with #2 | small data | Abandoned |
| 4 | Distribution-free control limits for PCA/PLS process monitoring (conformal T²/Q) | Tennessee Eastman (simulated), metal-etch (129 wafers) | Conformal SPC literature exists; classic small data | incremental | Rejected |
| 5 | Conformal FDR control for MS/MS spectral-library annotation | GNPS/MassBank – yes | Conformal retrieval sets already published (Rakhshaninejad et al., JCIM 2026); target–decoy FDR established | pre-empted | Rejected |
| 6 | Coverage-guaranteed collision-cross-section (CCS) windows across IMS platforms | METLIN-CCS/AllCCS – partly | open but CCS errors are small (~2 %), weaker practical need | weak significance | Rejected |
| 7 | Semi-quantification (ionisation-efficiency) prediction intervals for NTS | limited public multi-lab response data | partly addressed by MS2Quant-type work | data scarcity | Rejected |
| 8 | Electronic-nose drift compensation (UCI gas-sensor drift) | yes | saturated (dozens of domain-adaptation papers) | incremental | Rejected |
| 9 | Open-set Raman pathogen identification | Bacteria-ID data availability uncertain | open-set/conformal Raman/IR papers 2025–2027 exist | pre-empted/data | Rejected |
| 10 | Microplastic IR identification with abstention | yes | 2026–2027 conformal/dual-evidence microplastic papers exist | pre-empted | Rejected |
| 11 | Spectral pre-processing selection / ensembles | yes | saturated (SPORT/PROSAC/ensemble pre-processing) | incremental | Rejected |
| 12 | Wavelength-selection algorithms | yes | saturated | incremental | Rejected |
| 13 | Deep learning for NIR calibration | yes | saturated | incremental, CPU-heavy | Rejected |
| 14 | LIBS calibration transfer (ChemCam) | yes | niche | weak fit | Rejected |
| 15 | Retention-order-based joint annotation | yes | LC-MS2Struct (Bach et al. 2022) | pre-empted | Rejected |
| 16 | Optimal calibrant selection for RT transfer (design of calibrant mixes) | RepoRT | FastRet (JCIM 2026) "selective measuring" partly covers; kept as an analysis inside #1 | partial pre-emption | Merged into #1 |
| 17 | New point-prediction RT model (GNN/transformer) | RepoRT/SMRT | saturated (Uni-RT 2026, 2-step 2026, GNN-RT, RT-Transformer, FastRet …) | marginal accuracy gains | Rejected |
| 18 | GC retention-index prediction intervals | open RI data limited (NIST not open) | open | data licensing | Rejected |
| 19 | Soil MIR global-to-local transfer | OSSL | overlaps with #2 (soil regions) | self-overlap | Rejected |

## Stress test of the selected topic (skeptical ACA editor)
* *Too incremental / just combining techniques?* Risk acknowledged. The contribution is defined as an analytical
  data-treatment principle with measurable consequences: (a) quantify how often current RT-window practice discards the
  correct structure on new systems; (b) a calibrant-only procedure with a finite-sample guarantee; (c) a difficulty model
  meta-learned across hundreds of systems that adapts window width to compound and calibrant geometry; (d) the
  calibrant-number requirement (n ≥ 1/α − 1 for finite windows) and filtering-power trade-off. Each claim has an experiment.
* *Established solution exists?* Current practice uses fixed or ±k·RMSE windows, smoothing-spline (PredRet) or GP
  intervals; none guarantees coverage, and coverage across many systems has not been measured.
* *Simpler method achieves the same?* Plain jackknife+ is included as a baseline; if the meta-learned normalisation
  does not beat it, this is reported (negative result) and claims are restricted accordingly.
* *Enough data / leakage?* 239 RP target systems (≥60 compounds) from 5 source groups; leave-source-out folds; inner
  cross-fitting for difficulty-model training; compounds-seen analysis; strict novel-compound strata.
* *Can improvement disappear under proper statistics?* Unit = system; paired Wilcoxon with Holm correction; bootstrap CIs;
  per-fold consistency.
* *Useful beyond one dataset?* Second base model (SMRT regressor), HILIC systems (different separation mechanism),
  calibrant-selection strategies and perturbations.
* *Journal fit?* ACA scope explicitly invites "data treatment that truly demonstrate new principles and/or important
  applications"; ACA published 11 RT-prediction and 39 NTS papers in 2023–2026.
* *Lack of wet-lab validation?* Unavoidable here (no laboratory); mitigated by the scale and diversity of real
  inter-laboratory data. Stated as a limitation.

Decision: proceed with #1.
