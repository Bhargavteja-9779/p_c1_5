# S7. Laboratory-level guarantee of cross-system calibration (XSC)

**Setting.** A *laboratory* is a pair (T, D): a chromatographic system T with its population of analytes (the compounds whose
structures are to be filtered on T) and a draw D of n calibrants by the laboratory's selection rule ρ. The laboratory fits the
monotone map ĝ and the scale τ̂ (Section 2.6) on D; for a multiplier c its windows are
C_c(x) = ĝ(s(x)) ± c·τ̂·σ(x), and its *own retention* is

  cov_T(c) = P_{x ~ analytes of T}{ t(x) ∈ C_c(x) | D } ,

a random quantity (it depends on T and D). Because C_c(x) grows monotonically in c, the smallest multiplier that gives the
laboratory a retention of at least 1 − ε,

  V_T(ε) = inf{ c : cov_T(c) ≥ 1 − ε } ,

is well defined, and cov_T(c) ≥ 1 − ε ⇔ V_T(ε) ≤ c.

**Assumption (system exchangeability).** The target laboratory (T, D) and S training laboratories (s_1, D_1), …, (s_S, D_S)
are exchangeable, where D_k is drawn on s_k by the same rule ρ and the same n, and the base score and the difficulty model used
on s_k did not use data from s_k (cross-fitting by source group; Section 3.4). No assumption is made on the distribution of
compounds within a system, on the shape of the error distribution or on the quality of the base model.

**Proposition 1 (exact level, one draw per system).** Let V_k = V_{s_k}(ε) computed on the held-out compounds of s_k and
ĉ = the ⌈(1 − δ)(S + 1)⌉-th smallest of V_1, …, V_S (ĉ = ∞ if this index exceeds S). Then

  P{ cov_T(ĉ) ≥ 1 − ε } ≥ 1 − δ.

*Proof.* V_T, V_1, …, V_S are exchangeable, so the rank of V_T among the S + 1 values is uniform (ties broken at random), and
P{V_T ≤ ĉ} ≥ ⌈(1 − δ)(S + 1)⌉/(S + 1) ≥ 1 − δ — the standard split-conformal argument [Vovk et al. 2005; Lei et al. 2018]
applied to laboratories instead of compounds. The claim follows from cov_T(ĉ) ≥ 1 − ε ⇔ V_T ≤ ĉ. ∎

The guarantee is *training-conditional with respect to the laboratory's own calibrants*: it holds for the laboratory's actual
calibrant draw (with probability 1 − δ over laboratories and draws), not on average over hypothetical repetitions, which is the
only kind of guarantee that jackknife+ or split conformal windows give with n ≈ 20 calibrants (Section 4.2). It needs no
minimum n: n enters only through the replayed protocol. For n ≥ 1/α − 1 this is the cross-system analogue of PAC/
training-conditional conformal prediction [Vovk 2012; Park et al. 2022, meta-learning PAC prediction sets].

**Draw-averaged version (used in the main text).** With several replay draws per system (8 here), we use the system-level CDF
F_k(c) = (1/R)Σ_r 1{V_{k,r} ≤ c} and choose the smallest c with (1/S)Σ_k F_k(c) ≥ (1 − δ)(1 + 1/S). This reduces the
Monte-Carlo noise of a single draw per system; its validity is approximate (it is exact as R → ∞ for the population version
of the pooled CDF, cf. the pooling estimators of Dunn et al. 2023). Both versions are reported (Table S-XSC); they differ by
at most a few percent in c.

**Average-retention version.** Choosing the smallest c with (1/S)Σ_k cov_{s_k}(c) ≥ 1 − α controls the expected retention
over laboratories, E cov_T(c) ≥ 1 − α, up to an O(1/S) term (selection-aware calibration of Section 2.7 is the special case
in which the window shape is that of jackknife+).

**What the guarantee does not cover.** (i) Laboratories whose system or analyte population is systematically unlike the
public systems (e.g. a new column chemistry or a mobile-phase pH absent from RepoRT) violate exchangeability; leave-source-out
evaluation and the external NORMAN methods test this empirically (Sections 4.3, 4.6). (ii) The analytes of a laboratory are
represented by the measured compounds of public systems; analytes that are systematically harder than these (e.g. outside the
applicability domain of the base model) are covered only to the extent that the difficulty model σ widens their windows.
