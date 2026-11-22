# S7. Laboratory-level guarantee of cross-system calibration (XSC)

**Setting.** A *laboratory* is a pair (T, D): a chromatographic system T with its population of analytes (the compounds whose
structures are to be filtered on T) and a draw D of n calibrants by the laboratory's selection rule ρ. The laboratory fits the
monotone map ĝ and the scale τ̂ (Section 2.6) on D; for a multiplier c its windows are
C_c(x) = ĝ(s(x)) ± c·τ̂·σ(x), and its *own retention* is

  cov_T(c) = P_{x ~ analytes of T}{ t(x) ∈ C_c(x) | D } ,

a random quantity (it depends on T and D). In the evaluation, cov_T is the empirical retention on the analytes measured on T (all compounds other than the calibrants), and V_k below is computed in the same way; the guarantee therefore concerns the empirical retention on the measured analytes. Because C_c(x) grows monotonically in c, the smallest multiplier that gives the
laboratory a retention of at least 1 − ε,

  V_T(ε) = inf{ c : cov_T(c) ≥ 1 − ε } ,

is well defined, and cov_T(c) ≥ 1 − ε ⇔ V_T(ε) ≤ c.

**Assumption (system exchangeability).** The target laboratory (T, D) and S training laboratories (s_1, D_1), …, (s_S, D_S)
are exchangeable, where D_k is drawn on s_k by the same rule ρ and the same n, and the base score and the difficulty model used
on s_k did not use data from s_k (cross-fitting by source group; main text, Section 3.3). No assumption is made on the distribution of
compounds within a system, on the shape of the error distribution or on the quality of the base model.

**Proposition 1 (exact level, one draw per system).** Let V_k = V_{s_k}(ε) computed on the held-out compounds of s_k and
ĉ = the ⌈(1 − δ)(S + 1)⌉-th smallest of V_1, …, V_S (ĉ = ∞ if this index exceeds S). Then

  P{ cov_T(ĉ) ≥ 1 − ε } ≥ 1 − δ.

*Proof.* V_T, V_1, …, V_S are exchangeable, so the rank of V_T among the S + 1 values is uniform (ties broken at random), and
P{V_T ≤ ĉ} ≥ ⌈(1 − δ)(S + 1)⌉/(S + 1) ≥ 1 − δ — the standard split-conformal argument [Vovk et al. 2005; Lei et al. 2018]
applied to laboratories instead of compounds. The claim follows from cov_T(ĉ) ≥ 1 − ε ⇔ V_T ≤ ĉ. ∎

The probability in Proposition 1 is over the laboratory (T, D) *and* the S public systems. The guarantee is therefore of PAC type
with respect to the laboratory's own calibrant draw: it bounds the probability that the retention achieved with the actual
calibrants is inadequate, not merely the average retention, which is all that jackknife+ or split conformal windows provide with
n ≈ 20 calibrants (main text, Section 4.1). It is marginal with respect to the public systems. A guarantee conditional on the
public set as well would require a binomial correction of the order index (Vovk 2012; Park et al. 2022). With S ≈ 200 it would
raise the index by a few ranks, and we did not apply it. It needs no
minimum n: n enters only through the replayed protocol. For n ≥ 1/α − 1 this is the cross-system analogue of PAC/
training-conditional conformal prediction [Vovk 2012; Park et al. 2022, meta-learning PAC prediction sets].

**Draw-averaged version (used in the main text).** With several replay draws per system (R = 8 here), we use the system-level
CDF F_k(c) = (1/R)Σ_r 1{V_{k,r} ≤ c} and choose the smallest c with (1/S)Σ_k F_k(c) ≥ (1 − δ)(1 + 1/S). This reduces the
Monte-Carlo noise of a single draw per system, but its finite-sample validity is only approximate. For averages of
within-group CDFs the worst-case level is 1 − 2δ (Dunn et al. 2023). The exact rule is therefore evaluated alongside it. For 90 % of the
protocols the two multipliers differ by less than {{x.mult_diff_q90|pct0}}, and their failure rates on held-out systems are reported in
the main text (Section 4.2) and Fig. S3. The grid of multipliers ends at c = 40: an entry ∞ in Table S28 means that the required
multiplier exceeds 40, so that the windows are effectively unbounded (this occurs only for σ ≡ 1 with the 'most polar' rule).

**Average-retention version.** Choosing the smallest c with (1/S)Σ_k cov_{s_k}(c) ≥ 1 − α controls the expected retention
over laboratories, E cov_T(c) ≥ 1 − α, up to an O(1/S) term (selection-aware calibration of Section S6 is the special case
in which the window shape is that of jackknife+).

**What the guarantee does not cover.** (i) Laboratories whose system or analyte population is systematically unlike the
public systems (e.g. a new column chemistry or a mobile-phase pH absent from RepoRT) violate exchangeability; leave-source-out
evaluation and the external NORMAN methods test this empirically (main text, Sections 4.2 and 4.6). (ii) The analytes of a laboratory are
represented by the measured compounds of public systems; analytes that are systematically harder than these (e.g. outside the
applicability domain of the base model) are covered only to the extent that the difficulty model σ widens their windows.
