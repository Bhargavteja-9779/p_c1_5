# EXPERIMENT_LOG

Created 2026-11-15.

## 2026-11-15
* 17:41 Built processed RepoRT tables and RDKit feature cache (234 s; 0 featurisation failures).
* 17:52 Quick-mode smoke test (tiny models; 10 systems) — pipeline verification only; found B6 LOO-consistency bug → fixed
  (EXPERIMENT_PLAN Addendum A). Unit tests 9/9 pass (incl. synthetic coverage of jackknife+ and normalised jackknife+).
* 17:58 Full main run launched; died 18:07 when the container was suspended (usage limit). Fold-0 caches survived.
* 22:2x Relaunched with per-fold checkpointing as a single sequential pipeline (main → follow-ups).
