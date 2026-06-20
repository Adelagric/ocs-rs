# On the support size of the genomic OCS optimum

> A self-contained note — a proof, a counterexample, and the empirical regime — for the
> number of individuals an optimum-contribution-selection optimum activates. Companion to
> the support-first solver (`MANUSCRIPT.md`). Every empirical claim is reproducible; the
> scripts are listed under *Numerical evidence*.

## 1. The question

Genomic optimum contribution selection (OCS) solves, for n candidates with breeding values
**b** ∈ ℝⁿ and relationship matrix **G = ZZᵀ/s + εI** (with **Z** the n×m centred genotype
matrix, s a scale, ε a small ridge for positive-definiteness),

  maximise **bᵀc**  s.t.  **Ac = d**,  **0 ≤ c ≤ u**,  **cᵀGc ≤ k**,    (OCS)

where **A** is the q×n budget matrix (q = 1 for the simplex 𝟙ᵀc = 1; q = 2 for the sexed
form Σ_males = Σ_females = ½). Empirically the optimum **c\*** activates very few
candidates: its support S = {i : c\*ᵢ > 0} is ≈ 15–30 on real panels and stays bounded as n
grows at a fixed cap. What bounds |S|?

The answer is a **bracket**. In the no-ridge limit there is a clean n-independent bound
(Theorem 1). Once the ridge is present there is provably **no** universal bound — the
support can be the whole population (Theorem 2). The practical regime between them is set by
the joint geometry of (spectrum, **b**, cap k); it admits no single-scalar law, and we map
it empirically (§5). A KKT identity (§4) bridges the two theorems.

## 2. Theorem 1 — the ε = 0 bound

**Theorem 1.** *For ε = 0 (so* **G** = **G₀** = **ZZᵀ**/s *of rank r ≤ m), (OCS) has an
optimal* **c\*** *with* **|S| ≤ q + r + 1**, *independent of n.*

*Proof.* By strong duality (Slater holds at any interior-feasible instance) there is a
multiplier λ\* ≥ 0 such that **c\*** maximises the Lagrangian L(**c**) = **bᵀc** −
λ\* **cᵀG₀c** over the polytope {**Ac = d**, **c ≥ 0**}. Since **cᵀG₀c** = ‖**Zᵀc**‖²/s, the
value L(**c**) depends on **c** only through the pair (**Zᵀc**, **bᵀc**) ∈ ℝ^{r+1}. Hence
the slice

  P = { **c ≥ 0** : **Ac = d**, **Zᵀc = Zᵀc\***, **bᵀc = bᵀc\*** }

is **entirely optimal**: every **c** ∈ P is feasible (**cᵀG₀c** = ‖**Zᵀc\***‖²/s ≤ k) and
attains the optimal value **bᵀc\***. P is a nonempty (it contains **c\***), bounded (the
budget plus **c ≥ 0**) polytope cut by q + r + 1 equality rows, so it has a vertex, and a
vertex of {**c ≥ 0** : **Mc = e**} with q + r + 1 rows has at most q + r + 1 nonzeros. That
vertex is an optimum with the stated support. ∎

This is the Carathéodory / Barvinok–Pataki idea applied correctly: one does *not* count on
the curved boundary of the ellipsoid (where, **G** being positive definite, every point is
already an extreme point of the feasible set and the support is unconstrained) — instead one
**fixes the low-rank image Zᵀc and the objective**, which makes the optimal slice affine, and
counts LP vertices there. Active-set solvers (support-first; the critical line algorithm)
return such a vertex, whereas interior-point and ADMM solvers return non-sparse interior
points and threshold post hoc — so the bound is a statement about precisely what an
active-set solver computes.

## 3. Theorem 2 — no universal bound for ε > 0

**Theorem 2.** *For ε > 0 no bound on |S| of the form f(q, r) independent of n exists: the
support can equal n.*

*Proof (counterexample).* Take **G** = εI — the degenerate case m = 0 (no markers, every pair
equally unrelated; equivalently rank(**G₀**) = 0). Then (OCS), simplex form, is

  maximise **bᵀc**  s.t.  𝟙ᵀ**c** = 1, **c ≥ 0**, ε‖**c**‖² ≤ k,

and ‖**c**‖² = Σ c²ᵢ is exactly the rate-of-inbreeding proxy. The minimum of ‖**c**‖² on the
simplex is 1/n, attained only at the uniform plan **c** = 𝟙/n. For k slightly above ε/n the
feasible set is the simplex intersected with a ball that shrinks onto 𝟙/n, forcing |S| = n.
Numerically (`bound_validation.py`, block 7, n = 300): |S| = 300, 295, 256, 151, 52, 12 as k
loosens from 1.05× to 50× the minimum. So |S| ranges over the whole interval [1, n] with no
n-independent ceiling. ∎

Genetically, **G** = εI is the *no-structure* limit; with no relatedness pattern the
diversity cap can be met only by spreading contributions over the whole population. The
small support seen on real data is therefore a property of the **structure** of **G** (and of
**b**), not a worst-case guarantee.

## 4. The bridge — a KKT identity

**Proposition.** *At an optimum with the kinship constraint active, the contributions on the
support are an affine function of each candidate's augmented feature* (bᵢ, **zᵢ**) ∈ ℝ^{m+1}
*(* **zᵢ** *the i-th row of* **Z** *):*  **c\*ᵢ = α bᵢ + wᵀzᵢ + β_{sex(i)}**  *for i ∈ S.*

*Proof.* KKT gives μ ∈ ℝ^q, λ ≥ 0, **s ≥ 0** with **b** = **Aᵀμ** + 2λ**Gc\*** − **s** and
sᵢc\*ᵢ = 0. On S (sᵢ = 0): bᵢ = (**Aᵀμ**)ᵢ + 2λ(**Gc\***)ᵢ. Writing **y** := **Zᵀc\*** and using
**Gc\*** = **Z**y/s + ε**c\*** gives, for i ∈ S,
ε c\*ᵢ = (bᵢ − (**Aᵀμ**)ᵢ)/(2λ) − (**zᵢ**ᵀ**y**)/s, i.e. **c\*ᵢ = α bᵢ + wᵀzᵢ + β_{sex(i)}**
with α = 1/(2λε), **w** = −**y**/(sε), β_{sex(i)} = −(**Aᵀμ**)ᵢ/(2λε). ∎

The support's contributions thus live on an (m + 2)-parameter family. As ε → 0 the εc\*ᵢ term
vanishes and the identity forces **b** − **Aᵀμ** into the rank-r row space of **Z**,
recovering Theorem 1 and exhibiting m (the marker count) as the relevant dimension. For
ε > 0 the same identity does **not** cap |S| — consistent with Theorem 2.

## 5. The empirical regime between the theorems

Between the two clean statements lies the regime breeders actually use, which we mapped
across synthetic spectra, b-alignments and caps and on the real panels. The findings:

- **What drives |S| is b's alignment with the spectrum, inversely.** With **G** and k held
  fixed, putting **b** on the dominant (coancestry-expensive) eigendirection forces a large
  support (mean |S| ≈ 146 of n = 800 for **b** = top eigenvector); spreading **b** across many
  cheap directions collapses it to ≈ 4. Gain sought in an expensive direction must be diluted
  over many candidates to respect the cap (`bound_balign.py`).
- **No single scalar predicts |S| across regimes.** Neither rank(**G₀**), the effective rank
  / participation ratio of the spectrum, the directional cost **bᵀGb/bᵀb**, nor the support of
  the unconstrained direction (**G⁻¹b**)₊ tracks |S| across both synthetic and real instances;
  the best dimensionless combination, |S| ∼ (**bᵀGb**/**bᵀb** · 1/k)^{0.8}, reaches only
  R² ≈ 0.58 and mispredicts the real panels' relative support (`bound_predictor.py`,
  `bound_lawfit.py`).
- **The curve |S|(k) is a per-instance power law of non-universal exponent.** |S| ∼ k^{−α}
  fits each instance well (R² ≈ 0.9–0.98 where the support is well resolved), with α ≈ 1 on
  both real panels but ranging 0.5–1.7 across synthetic spectra and alignments
  (`bound_curve.py`).
- **Real panels.** Wheat (n = 599) and mouse (n = 1814): |S| = 25 / 26 and 17 / 26
  (simplex / sexed-optiSel) at the working cap, both of the order of tens — small and stable
  in n, as the spectrum-with-decay regime predicts qualitatively, though (per the above) no
  formula pins the value (`bound_real.py`).

Theorem 2 explains *why* no scalar law exists: there is no universal bound to predict, so the
practical value is genuinely a function of the full (spectrum, **b**, k) geometry.

## 6. A conditional bound: the obstruction and the conjecture

The universal question is closed negatively (Theorem 2); the useful one is **conditional**. Two
clean bounds are provable and bracket the rest: at a **loose cap** (kinship slack, cᵀGc\* < k) the
active set is the polytope {Ac = d, c ≥ 0} alone and c\* is an LP vertex, so |S| ≤ q; at **exact
rank** d, Theorem 1 gives |S| ≤ q + d + 1.

A bound in terms of *effective* rank — the realistic ask, since real spectra decay — does **not**
follow from the same argument, and the obstruction is precise. Theorem 1 linearises the optimal
slice by pinning Zᵀc, i.e. the **entire** range of G₀; for a full-rank G one would have to pin all
n projections uᵢᵀc, giving only the trivial |S| ≤ q + n + 1. One cannot pin just the top-d
projections and drop the tail, because the tail contributes cᵀEc = Σ_{i>d} λᵢ(uᵢᵀc)² — a *genuine
quadratic*, not a low-rank term — so fixing it is not a linear constraint and the slice is not a
polytope. The tail's curvature is exactly what lets the support spread; effective rank cannot enter
a vertex-counting argument. A conditional bound therefore needs a different (spectral-perturbation)
argument.

It also needs two assumptions, **both necessary**:
- **(A1) No large degenerate cheap plateau.** A block of many near-equal small eigenvalues is a
  cheap subspace the optimum spreads into; G = εI is the extreme case, giving |S| = n (Theorem 2).
  A smoothly decaying spectrum has no such plateau.
- **(A2) b avoids the dominant directions.** Even with a gap, b on the top (expensive) eigenvector
  forces dilution and a large support (mean |S| ≈ 146 of 800 in the alignment sweep).

Both counterexamples are established, so neither assumption can be dropped. Under (A1)+(A2) — the
regime of real genomic matrices, where |S| ≈ 15–30 — we conjecture |S| is small and n-independent,
of the order of the dominant eigendirections b actually excites. The growth half (|S| rising as
k → 0) is governed by the classical contributions ↔ ΔF ↔ Nₑ link (Wray & Thompson 1990; Woolliams &
Bijma 2000). Proving the bound — coupling the spectral gap with b's projection onto the top
eigenspace, robust to the tail's curvature — is the open problem, the meeting point of optimisation
(the perturbed-cone face) and quantitative genetics (effective lineages).

**An exact handle on the effective support.** Write the *effective* number of contributors as the
participation ratio PR(c\*) = 1/‖c\*‖² (it equals |S| for uniform contributions, less when they are
uneven). At the active cap with Σc\* = 1 an identity holds: **PR(c\*) = R(c\*)/k**, where
R(c\*) = c\*ᵀGc\*/c\*ᵀc\* is the optimum's own Rayleigh quotient (since R·‖c\*‖² = c\*ᵀGc\* = k).
Hence **PR(c\*) ≤ λ₁/k** unconditionally, and the conditional bound *reduces to bounding R(c\*)* —
the directional coancestry cost the optimum is forced to incur. This is exactly where (A2) bites:
when b avoids the dominant directions, c\* lives in the cheap part of the spectrum and R(c\*) is
small (measured: R ≈ 0.02–0.4 on the structured cases versus λ₁ ≈ 45–143, giving PR ≈ 2–12), whereas
G = εI forces R = ε = λ₁ and PR = n. The residual open piece is now sharp and scalar — bound R(c\*)
a priori under (A1)+(A2). (The raw cardinality |S| runs ≈ 1.5–2.5× PR here; PR is the quantity with
the clean law.) Verified in `bound_kernel.py`.

## Numerical evidence (reproducible)

- `bound_validation.py` — solver vs SciPy; the ε = 0 bound; n-independence; the effective-rank
  behaviour; the ridge sweep; a spectral-decay sweep; and (block 7) the Theorem-2
  counterexample.
- `bound_real.py` — wheat and mouse kinship matrices (after `research/repro/*_export.R`).
- `bound_balign.py` — the b-alignment sweep at fixed **G** and cap.
- `bound_predictor.py` — the cross-regime scalar-predictor search.
- `bound_lawfit.py` — the dimensionless law fit.
- `bound_curve.py` — the per-instance |S|(k) power-law fit.

## References

- Pataki, G. (1998). On the rank of extreme matrices in semidefinite programs…
  *Math. Oper. Res.* 23(2):339–358. DOI 10.1287/moor.23.2.339. — extreme-point rank counting.
- Markowitz, H. (1956). *Naval Res. Logist. Q.* 3:111–133. DOI 10.1002/nav.3800030110. — critical-line / corner-portfolio structure.
- Wray, N.R. & Thompson, R. (1990). *Genet. Res.* 55(1):41–54. DOI 10.1017/S0016672300025180. — ΔF ∝ Σ(contribution²).
- Woolliams, J.A. & Bijma, P. (2000). *Genetics* 154(4):1851–1864. DOI 10.1093/genetics/154.4.1851. — contributions ↔ ΔF ↔ Nₑ.
- Yamashita, M., Mullin, T.J. & Safarina, S. (2018). *Optim. Lett.* 12(7):1683–1697. DOI 10.1007/s11590-018-1229-y. — OCS as a single second-order-cone program.
- Waldmann, P. (2025). *Bioinform. Adv.* 5(1):vbaf259. DOI 10.1093/bioadv/vbaf259. — independent empirical support counts; post-hoc truncation.
