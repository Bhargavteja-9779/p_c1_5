# EXPERIMENT_LOG

Created 2026-11-15.

## 2026-11-15
* 17:41 Built processed RepoRT tables and RDKit feature cache (234 s; 0 featurisation failures).
* 17:52 Quick-mode smoke test (tiny models; 10 systems) — pipeline verification only; found B6 LOO-consistency bug → fixed
  (EXPERIMENT_PLAN Addendum A). Unit tests 9/9 pass (incl. synthetic coverage of jackknife+ and normalised jackknife+).
* 17:58 Full main run launched; died 18:07 when the container was suspended (usage limit). Fold-0 caches survived.
* 22:2x Relaunched with per-fold checkpointing as a single sequential pipeline (main → follow-ups).

## 2026-11-16
* 03:11 Pipeline relaunched after container reclaim (resumed at fold 2 from per-fold checkpoints).
* 03:53 Main experiment done (5 folds, 228 RP target systems; 2541 s for the resumed part).
* 04:xx–05:23 Follow-ups done: robustness (span, common, noise, outlier), SMRT base model, HILIC, sensitivity (4 settings),
  fixed-tolerance windows on identical draws (Addendum B), small n (5, 8).
* 05:58 SASC (Addendum D): c per fold/rule from training systems; NORMAN-mix factor from 36 RepoRT systems (c_P 0.92, c_J 0.64).
* 06:01 NORMAN external validation (Addendum C): 37 systems, 82 compounds (95 % seen by base model); incl. P+SASC / J+SASC rows.
* 06:05 Compute benchmark (1 thread) → PIPELINE_DONE 2026-11-16T06:05:56Z.
* 06:08 analyze.py (B = 10,000) → 57 table files, 18 figure files, metrics.json, statistics.json, raw_results.csv.
* 06:17 analyze.py re-run after fixing the NORMAN per-stratum table (strata were pooled over the three calibration scenarios;
  now reported per scenario). No other output changed in substance.
* 06:28 Diagnostic (pending since 05:32): NORMAN systems lack temperature/flow/dead time (NaN). In RepoRT these are missing for
  116/33/37 of 297 RP systems and encoded as 0. LightGBM converts NaN to 0 for splits trained without NaN (empirical check,
  lightgbm 4.7.0: max |pred(NaN) − pred(0)| = 0.0), so NORMAN systems are encoded exactly like RepoRT systems with missing
  metadata — no encoding mismatch explains the missing efficiency gain on NORMAN.
