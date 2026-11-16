# Experimental

## Data

All RT data were taken from RepoRT v1.0 [@kretschmer2024report], a curated repository of small-molecule RTs with standardised structures and chromatographic metadata (Zenodo record 10.5281/zenodo.16267768, CC BY-SA 4.0; archive SHA-256 4eeeb46a…6091). Compounds were identified by the connectivity block of the standard InChIKey; stereoisomer entries within a system were averaged when their RTs differed by ≤0.2 min and removed otherwise. Mobile-phase pH class and organic modifier were derived from the reported additives and solvents. Target systems were all reversed-phase (RP) systems with at least 60 compounds after cleaning, excluding METLIN SMRT, which served only as an auxiliary training library: {{n_target_systems}} systems with {{n_rt_values_targets}} retention times (median {{median_compounds}} compounds per system; Table 1, Fig. S1). HILIC systems with ≥60 compounds were analysed separately as an out-of-mechanism test ({{n_hilic_systems}} systems).

{{table:table1}}

Table: Table 1. Reversed-phase target systems (≥60 compounds) by source group, which defines the leave-source-out folds.

## Evaluation protocol and leakage control

Systems were partitioned by data source into five leave-source-out folds ({{fold_list}}). For each fold, the base model, the bootstrap ensemble (B6) and the difficulty model were trained exclusively on systems from the other four folds (and SMRT; numbers of training systems per fold in Table S1); the only information taken from a target system was its calibrants. Within each target system, calibrant sets of size $n\in\{10,15,20,30,40\}$ were drawn uniformly at random (20 draws per $n$, fixed seeds), and all remaining compounds of the system formed the test set; windows were computed for $\alpha\in\{0.05,0.10,0.20\}$. The primary setting ($n=20$, $\alpha=0.10$), all comparators, metrics and statistical tests were specified in a plan written before any target-system result was inspected (code/EXPERIMENT_PLAN.md). Because many compounds of a target system also occur in public data from other laboratories ({{share_seen_test}} of test compounds at $n=20$ were present in the base model's training data), test compounds were additionally stratified into those seen and not seen by the base model, and by their maximum Tanimoto similarity to its training set.

Pre-specified robustness experiments used $n=20$ and 10 draws per system with (i) calibrants spread over the RT range (one random compound per RT quantile stratum), (ii) the 20 compounds of each system most frequently measured across RepoRT, (iii) Gaussian noise with a standard deviation of 2 % of the RT span added to the calibrant RTs and (iv) one calibrant shifted by 10 % of the span. Generalisation was assessed with a second base model (a LightGBM regressor trained on SMRT only) and on HILIC systems (base and difficulty models trained on HILIC systems from other sources), and sensitivity to the difficulty-model hyper-parameters with 15 or 63 leaves and 200 or 800 trees. Analyses added after the plan was written but before the corresponding results were computed are labelled post hoc: fixed-tolerance windows on the same calibrant draws, very small calibrant sets ($n=5$ and 8), the external validation on NORMAN data and the selection-aware calibration (plan addenda B–D).

## Metrics

For each system and draw we computed empirical coverage (fraction of test compounds whose RT lies in their window), mean and median window width (min), and the Winkler interval score [@winkler1972decision; @gneiting2007strictly]
$$\mathrm{IS}_\alpha(l,u;t)=(u-l)+\frac{2}{\alpha}(l-t)\,\mathbb{1}\{t<l\}+\frac{2}{\alpha}(t-u)\,\mathbb{1}\{t>u\},$$
a proper score that rewards narrow windows but penalises misses, so that an under-covering method cannot appear efficient merely by being narrow. Candidate-filtering performance was evaluated with *isomeric decoys*: for each test compound, up to 30 RepoRT compounds with the same molecular formula but a different connectivity block were treated as incorrect candidates; a decoy was eliminated when the observed RT of the true compound fell outside the decoy's window (draws 0–4 at $n=20$).

## Statistical analysis

The unit of analysis is the target system; metrics were averaged over draws within each system. Mean coverage is reported with 95 % bootstrap confidence intervals (10,000 resamples of systems) and, as a sensitivity analysis, with a cluster bootstrap over source folds. The proposed method was compared with each comparator by two-sided Wilcoxon signed-rank tests on per-system interval scores, widths and absolute coverage gaps, with Holm correction within each metric family [@holm1979simple] and the matched-pairs rank-biserial correlation as effect size [@kerby2014simple]. Robustness experiments used 10 draws per system.

## Implementation

All computations ran on a 4-core CPU without GPU. Code (Python 3.13; RDKit, XGBoost, LightGBM, scikit-learn, SciPy), configuration files, unit tests (including a simulation check of the jackknife+ coverage guarantee), the pre-registered plan and all raw outputs are provided (see Data availability); every number, table and figure in this article is regenerated by `python run_all.py --mode full`.
