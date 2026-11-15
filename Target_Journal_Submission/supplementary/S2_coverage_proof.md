# Supplementary Section S2 — Coverage of the normalised jackknife+ window

**Setting.** Let $Z_i=(X_i,Y_i)$, $i=1,\dots,n+1$, be exchangeable, where $Z_1,\dots,Z_n$ are the calibrants and
$Z_{n+1}$ is the compound whose retention time $Y_{n+1}$ is to be covered. Let $\mathcal{A}$ be any algorithm that maps a
training multiset to a pair $(\hat\mu,\hat\sigma)$ with $\hat\sigma(\cdot)>0$ and that is *symmetric* (invariant to the
order of its training points). In our method $\hat\mu=g$ is the monotone calibrant map and
$\hat\sigma(x)=\exp\{h(\phi(x,\mathcal{D},z_T))\}+10^{-3}$, where $h$ is fixed (trained on other chromatographic systems
only) and $\phi$ depends on the training set $\mathcal{D}$ only through symmetric summaries (score range, nearest-calibrant
gap, RT span, set size).

For $i\neq j$ let $(\tilde\mu_{-(i,j)},\tilde\sigma_{-(i,j)})=\mathcal{A}(\{Z_k:k\neq i,j\})$ and define the leave-two-out
scaled residuals
$$\tilde R_{ij}=\frac{|Y_i-\tilde\mu_{-(i,j)}(X_i)|}{\tilde\sigma_{-(i,j)}(X_i)} .$$
The quantities used by the method are $\hat\mu_{-i}=\tilde\mu_{-(i,n+1)}$, $\hat\sigma_{-i}=\tilde\sigma_{-(i,n+1)}$ and
$R_i=\tilde R_{i,n+1}$, and the window is
$$C_\alpha(x)=\big[\hat q^-_\alpha\{\hat\mu_{-i}(x)-\hat\sigma_{-i}(x)R_i\},\ \hat q^+_\alpha\{\hat\mu_{-i}(x)+\hat\sigma_{-i}(x)R_i\}\big].$$

**Proposition.** $\mathbb{P}\{Y_{n+1}\in C_\alpha(X_{n+1})\}\ge 1-2\alpha$.

**Proof.** *(i) Deterministic step.* Suppose $Y_{n+1}>\hat q^+_\alpha\{v_i\}$ with $v_i=\hat\mu_{-i}(X_{n+1})+\hat\sigma_{-i}(X_{n+1})R_i$.
Since $\hat q^+_\alpha$ is the $\lceil(1-\alpha)(n+1)\rceil$-th smallest $v_i$, $Y_{n+1}>v_i$ for at least
$\lceil(1-\alpha)(n+1)\rceil$ indices $i$; for each of them $|Y_{n+1}-\hat\mu_{-i}(X_{n+1})|/\hat\sigma_{-i}(X_{n+1})>R_i$,
i.e. $\tilde R_{n+1,i}>\tilde R_{i,n+1}$. If instead $Y_{n+1}<\hat q^-_\alpha\{u_i\}$ with
$u_i=\hat\mu_{-i}(X_{n+1})-\hat\sigma_{-i}(X_{n+1})R_i$ (the $\lfloor\alpha(n+1)\rfloor$-th smallest $u_i$), then
$Y_{n+1}<u_i$ for at least $n+1-\lfloor\alpha(n+1)\rfloor\ge(1-\alpha)(n+1)$ indices, and again
$\tilde R_{n+1,i}>\tilde R_{i,n+1}$ for each. Hence miscoverage implies $\sum_{i\le n}A_{n+1,i}\ge(1-\alpha)(n+1)$ for the
tournament matrix $A_{ij}=\mathbb{1}\{\tilde R_{ij}>\tilde R_{ji}\}$.

*(ii) Tournament step.* For any tournament matrix, the set $S(A)=\{i:\sum_jA_{ij}\ge(1-\alpha)(n+1)\}$ satisfies
$|S(A)|\le 2\alpha(n+1)$ (Barber et al., Ann. Stat. 49 (2021), proof of Theorem 1). Because $\mathcal{A}$ is symmetric and
the $Z_i$ are exchangeable, the distribution of $A$ is invariant under simultaneous permutations of its rows and
columns, so $\mathbb{P}\{n+1\in S(A)\}=\mathbb{E}|S(A)|/(n+1)\le 2\alpha$. Combining (i) and (ii) gives the claim. $\square$

**Remarks.** (1) Nothing in the argument uses the form of $\hat\sigma$ except positivity and symmetry, so any difficulty
model trained independently of the target data may be used; this is why the difficulty model is learned exclusively on
other systems and its calibrant-dependent inputs are recomputed for every leave-one-out set. (2) Plugging in the
full-calibrant-set $\hat\sigma$ for the calibrants' own residuals (the "no leave-one-out refit" ablation) breaks the
symmetry and is not covered by the proposition; empirically it under-covers. (3) The guarantee is marginal over the
random draw of calibrants and compound from the compounds analysed on the target system; it does not hold when
calibrants are selected non-randomly (e.g. a fixed standard mixture), which is why such selections are evaluated
empirically (robustness experiments). (4) The window is finite iff $\lceil(1-\alpha)(n+1)\rceil\le n$, i.e. $n\ge 1/\alpha-1$.
(5) A Monte-Carlo check with a known heteroscedastic scale and $n=20$, $\alpha=0.10$ is included in the unit tests
(`code/tests/test_core.py::test_normalised_jackknife_plus_coverage_exchangeable`).
