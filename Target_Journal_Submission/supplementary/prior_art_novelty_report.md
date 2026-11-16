# Prior-art investigation and novelty map (Phases 3–4)

Sources: PubMed/Crossref/arXiv/Europe PMC queries run by the author (session 1–2, 2026-11-15/16) and two verified literature
sweeps by agents (`literature_raw/sweep_rt_models.json`, `literature_raw/sweep_nts_practice.json`; 40 references, all
Crossref/PubMed-verified). Four adversarial "refute-the-claim" agents and two further sweeps were terminated by a usage
limit before returning; their role is taken over by the novelty reviewer in every review cycle (REVIEW/).
Key negative search: PubMed `"conformal"[tiab] AND "retention time"[tiab]` → 0 records (2026-11-15).

## What exists (closest works)

| Existing work | What it contributes | What it does not do (relative to this study) |
|---|---|---|
| García et al. 2022, J Cheminform (meta-learned GP projections) | Projects SMRT-model predictions to a new chromatographic method with a GP whose prior is meta-learned on other methods; probabilistic annotation from ~10 calibrants; interval-score evaluation | Bayesian/Gaussian intervals without finite-sample coverage guarantee; evaluated on a limited set of PredRet methods; no candidate-filter coverage analysis across hundreds of systems; no calibrant-number requirement |
| Ferrando Plo & Kruve 2026, JCIM (conformal CCS) | Conformal predictive system for **CCS**-based candidate prioritisation in NTS with guaranteed fraction of true structures prioritised; MCES applicability domain | CCS needs no per-system transfer: calibration uses a large held-out set. Does not address RT, the few-calibrant regime, re-use of calibrants for mapping and calibration, or meta-learned difficulty across systems |
| Kretschmer et al. 2026, Nat Methods (2-step) | Condition-aware retention-order index learned by pairwise ranking on 171 RepoRT RP datasets + monotone polynomial anchor mapping (~15 anchors); public code | Point predictions; fixed 0.40 min acceptance window; no prediction intervals or coverage analysis |
| Zhang et al. 2024, Commun Chem (post-projection calibration) | 35 calibrants correct cross-method projections; reports 6–20 % of true identities filtered at a 2×RMSE threshold | No coverage control; few methods |
| Malm & Kruve 2025, Analyst | Projection vs prediction accuracy on 37 NORMAN systems; recommends ±2 RMSE windows (~95 % expected) | Accuracy only; coverage of windows not evaluated; no guarantee |
| Domingo-Almenara et al. 2019, Nat Commun (SMRT) | Projection with ROC-optimised thresholds; 21–41 % of correct identities filtered on four methods | Heuristic thresholds; 4 systems |
| PredRet (Stanstrup et al. 2015) | Pairwise projection models with prediction intervals | Requires overlapping measured compounds; intervals not coverage-guaranteed |
| Bach et al. 2018 (retention-order ranking) | Pairwise ranking of retention order from multi-system data | No absolute RT, no windows |
| Matyushin et al. 2024 | Flags unreliable RT/RI/CCS predictions from ensemble spread, similarity, clusters | Classification of unreliable predictions; no calibrated windows |
| Beck et al. 2026 (MC-Retention) | Heteroscedastic, uncertainty-aware GNN across 8 conditions for method development | Calibrated uncertainty within known conditions; no new-system transfer guarantee |
| Bhargav Teja & S 2026, Talanta (same author group) | Group-conformal NIR calibration intervals under season/instrument/region shift | Spectroscopic calibration, not chromatographic transfer or candidate filtering |

## Novelty map → genuinely different elements of this study
1. **Conformal RT windows for the few-calibrant transfer setting.** The window is built only from the calibrants measured on
   the target system, with the calibrants used both to fit the monotone map and to calibrate (jackknife+), giving a finite-sample
   guarantee (≥1−2α; ≈1−α empirically) under exchangeability. *Closest:* Ferrando Plo & Kruve 2026 (CCS, no transfer), García
   2022 (RT, Bayesian, no guarantee). **Claim wording:** "to our knowledge, the first conformal (distribution-free) RT windows for
   transfer to a new chromatographic system from a handful of calibrants" — not "first conformal method in NTS".
2. **Difficulty normalisation meta-learned across chromatographic systems by simulating the calibrant protocol**, with
   leave-one-out-consistent calibrant-geometry features, keeping the guarantee (Supplementary S2). *Related:* meta-learned GP
   priors (García 2022), few-shot conformal with auxiliary tasks in ML (Fisch et al. 2021), uncertainty-flagging features
   (Matyushin 2024). **Claim wording:** "a meta-learned normalisation" (new combination), not "first meta-learned conformal".
3. **Scale of the coverage evaluation**: 239 RP target systems in leave-source-out folds + 37 NORMAN interlaboratory systems with
   a real fixed calibrant mix, + HILIC; previous true-identity-loss measurements covered ≤ a few dozen methods with heuristic windows.
4. **Calibrant-number requirement** n ≥ 1/α − 1 for a finite window (a direct consequence of conformal theory, stated and verified
   for analytical practice) with width-vs-n guidance. Modest novelty — presented as practical guidance, not as a theoretical result.

## Elements that are NOT new (stated as such in the manuscript)
* Pairwise retention-order learning (Bach 2018; Kretschmer 2026); monotone calibrant projection (PredRet; 2-step anchors).
* Jackknife+ and normalised conformal prediction as general techniques (Barber 2021; Lei 2018; Papadopoulos 2002).
* Using similarity-to-training-set as an uncertainty signal (Matyushin 2024; Ferrando Plo & Kruve 2026 MCES).
* Conformal guarantees for NTS candidate prioritisation as a concept (Ferrando Plo & Kruve 2026, CCS).

## Residual novelty risk
Medium: an editor may view the work as "applying jackknife+ to RT". Mitigation: the quantified failure of current practice at scale,
the few-calibrant transfer formulation, the LOO-consistent meta-learned normalisation with proof, external NORMAN validation and
practical calibrant guidance; honest statement of what is not new.
