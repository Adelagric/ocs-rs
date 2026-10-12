# On the support size of the genomic OCS optimum

> A self-contained note — a proof, a counterexample, and the empirical regime — for the
> number of individuals an optimum-contribution-selection optimum activates. Companion to
> the support-first solver (`MANUSCRIPT.tex`). Historical experiment results are retained with context; the
> scripts and their prerequisites are listed under *Numerical evidence*.

## 1. The question

Genomic optimum contribution selection (OCS) solves, for n candidates with breeding values
**b** ∈ ℝⁿ and relationship matrix **G = ZZᵀ/s + εI** (with **Z** the n×m centred genotype
matrix, s a scale, ε a small ridge for positive-definiteness),

  maximise **bᵀc**  s.t.  **Ac = d**,  **0 ≤ c ≤ u**,  **cᵀGc ≤ k**,    (OCS)

where **A** is the q×n budget matrix (q = 1 for the simplex 𝟙ᵀc = 1; q = 2 for the sexed
form Σ_males = Σ_females = ½). Empirically the optimum **c\*** activates very few
candidates: its total support S = {i : c\*ᵢ > 0} was ≈ 15–30 in the historical loose-cap runs. At the operational contribution-proxy points in the manuscript it is tens to hundreds (mouse: 47/89/182 at ΔF = 2/1/0.5%). A structured n-sweep at a fixed proxy point gives roughly 80–103. These observations do not establish a universal bound. What bounds |S|?

The answer is a **bracket**. In the uncapped no-ridge case there is an existential n-independent bound
(Theorem 1). Once the ridge is present there is provably **no** universal bound — the
support can be the whole population (Theorem 2). The practical regime between them is set by
the joint geometry of (spectrum, **b**, cap k); the tested scalar predictors are not reliable across regimes, and we map
it empirically (§5). A conditional KKT identity (§4) describes free contributions.

## 2. Theorem 1 — the ε = 0 bound

**Theorem 1 (uncapped contributions).** *Assume ε = 0, s > 0, a nonempty feasible set, and the simplex or sex-budget equalities so that {Ac = d, c ≥ 0} is bounded. With no effective upper bounds, an optimum exists with*
**|S| ≤ rank([A; Zᵀ; bᵀ]) ≤ q + r + 1**, *where r = rank(Z) ≤ m.*

*Proof.* Start with any attained optimum **c\***. Fix its low-rank image and objective in the slice

  P = { **c ≥ 0** : **Ac = d**, **Zᵀc = Zᵀc\***, **bᵀc = bᵀc\*** }.

Every point of P has the same kinship ‖**Zᵀc\***‖²/s and the same objective as **c\***, and is therefore optimal. The slice is nonempty and bounded. Let M = [A; Zᵀ; bᵀ], whose rank is at most q + r + 1. At a vertex, more than rank(M) positive coordinates would give a nonzero null direction supported on those coordinates; sufficiently small perturbations in both directions preserve nonnegativity and all equalities, contradicting extremality. A vertex of P proves the bound. No duality or constraint qualification is needed. ∎

**Finite-cap extension.** Intersect the slice with **c ≤ u** and choose a vertex there. Only its *free* coordinates F = {i : 0 < cᵢ < uᵢ} satisfy **|F| ≤ rank(M) ≤ q + r + 1**: a null direction on more free coordinates again permits perturbations both ways. Positive upper-bound coordinates U₊ = {i : cᵢ = uᵢ > 0} also contribute to the support, so **|S| = |F| + |U₊| ≤ rank(M) + |U₊|**. There is no n-independent total-support bound with arbitrary finite caps. For example, q = 1, r = 0, G = 0, k = 0 and uᵢ = 2/n force every feasible plan to have at least ceil(n/2) positive coordinates, even though q + r + 1 = 2. A strictly decreasing b attains that minimum at an optimal capped LP vertex.

This is an existence result about vertices of an *optimal slice*. It does not guarantee that the implemented active-set method returns such a vertex, nor does it imply that an interior-point or ADMM method cannot return sparse solutions. A separate vertex extraction would be needed to enforce this cardinality guarantee on an arbitrary solver output.

## 3. Theorem 2 — no universal bound for ε > 0

**Theorem 2.** *For ε > 0 no bound on total support of the form f(q, r) independent of n exists, even with ineffective upper bounds: the support can equal n.*

*Proof (counterexample).* Take the abstract OCS instance q = 1, G = εI, r = 0 and uᵢ = 1. This algebraic counterexample is not a usable marker panel for a public entry point. The constraint is ε‖c‖² ≤ k on the simplex. Cauchy–Schwarz gives ‖c‖² ≥ 1/n. A plan with any zero coordinate has ‖c‖² ≥ 1/(n−1). Consequently, for

  **ε/n < k < ε/(n−1)**,

the uniform plan is strictly kinship-feasible and every feasible plan, hence every optimum for any b, has full support. This proves the claim for arbitrary n ≥ 2. ∎

The older numerical sweep (`bound_validation.py`, block 7, n = 300) reported |S| = 300, 295, 256, 151, 52, 12 as k relaxed from 1.05× to 50× the minimum. These are observations for one b; 1.05× is outside the analytic full-support interval and is not a universal witness. Block 8 checks the rigorous interval and the finite-cap obstruction directly.

Here Σcᵢ² measures concentration of contributions. With Nₑ = 1/Σcᵢ², the quantity ΔF = 1/(2Nₑ) is a contribution-based proxy, not a measured rate of realised genomic inbreeding. Small observed support depends on the joint instance geometry and is not a worst-case guarantee.

## 4. A conditional KKT identity

**Proposition.** *Suppose ε > 0 and an optimum has KKT multipliers with λ > 0 for kinship. On the free coordinates F = {i : 0 < c\*ᵢ < uᵢ}, contributions are affine in (bᵢ, zᵢ, aᵢ):*

  **c\*ᵢ = α bᵢ + wᵀzᵢ + vᵀaᵢ**,  *i ∈ F,*

*where aᵢ is column i of A. For sex-incidence budgets the final term is β_{sex(i)}.*

*Proof.* Write lower and upper multipliers ℓ,h ≥ 0. KKT stationarity is **b = Aᵀμ + 2λGc\* − ℓ + h**, with ℓᵢc\*ᵢ = 0, hᵢ(uᵢ−c\*ᵢ) = 0 and λ(c\*ᵀGc\*−k) = 0. On F both bound multipliers vanish. Let y = Zᵀc\*. Since Gc\* = Zy/s + εc\*, divide free-coordinate stationarity by 2λε to obtain α = 1/(2λε), w = −y/(sε) and v = −μ/(2λε). ∎

This identity uses m + q + 1 coefficients (m + 2 for a simplex, m + 3 for two sex budgets). Active kinship alone does not imply λ > 0. KKT multipliers require an appropriate constraint qualification, or must be supplied as a certificate; the proposition assumes their existence. It does not apply to coordinates at an upper cap, does not bound their number, and proves no cardinality bound on F. Taking ε → 0 cannot recover Theorem 1 from this formula: multipliers and supports may change or diverge, and stationarity holds only on the current free coordinates. The optimal-slice proof stands independently.

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
  (simplex / sexed-optiSel) at the historical loose cap, historical loose-cap observations of order tens; these two fixed panels alone do not establish stability in n, and the operational points in the manuscript have larger supports (`bound_real.py`).

Theorem 2 rules out a universal bound depending only on q and r. The empirical sweeps separately show that the scalar predictors tested here do not reliably predict support across regimes.

## 6. A conditional bound: the obstruction and the conjecture

A conditional bound remains an open question. If the kinship constraint is *redundant on the whole budget polytope*, an uncapped LP has an optimal vertex with |S| ≤ rank(A) ≤ q; with finite caps, the same argument bounds only |F| ≤ rank(A). Slack kinship at one optimum is insufficient: for constant b and G = I, the uniform optimum is slack when 1/n < k < 1/(n−1), but every feasible optimum has full support. At exact rank r and ε = 0, Theorem 1 gives the uncapped existential bound; its finite-cap extension bounds free coordinates.

A bound in terms of *effective* rank — the realistic ask, since real spectra decay — does **not**
follow from the same argument, and the obstruction is precise. Theorem 1 linearises the optimal
slice by pinning Zᵀc, i.e. the **entire** range of G₀; for a full-rank G one would have to pin all
n projections uᵢᵀc, giving only the trivial |S| ≤ q + n + 1. One cannot pin just the top-d
projections and drop the tail, because the tail contributes cᵀEc = Σ_{i>d} λᵢ(uᵢᵀc)² — a *genuine
quadratic*, not a low-rank term — so fixing it is not a linear constraint and the slice is not a
polytope. The tail's curvature is exactly what lets the support spread; effective rank alone cannot replace exact rank in this vertex-counting argument. A conditional bound therefore needs a different (spectral-perturbation)
argument.

Two qualitative conditions motivate a possible conditional theorem; they have not been specified quantitatively or proved sufficient:
- **(A1) No large degenerate cheap plateau.** A block of many near-equal small eigenvalues is a
  cheap subspace the optimum spreads into; G = εI is the extreme case, giving |S| = n (Theorem 2).
  A smoothly decaying spectrum has no such plateau.
- **(A2) b avoids the dominant directions.** Even with a gap, b on the top (expensive) eigenvector
  forces dilution and a large support (mean |S| ≈ 146 of 800 in the alignment sweep).

The pure-ridge counterexample and alignment sweep motivate (A1) and (A2), but do not prove necessity of these informal conditions. Under quantitative versions of them, at a specified operating point, we conjecture a small support; the historical 15–30 values describe loose-cap runs, not all real genomic instances,
of the order of the dominant eigendirections b actually excites. The observed growth as k approaches its feasible minimum is indexed here by the contribution-based ΔF = 1/(2Nₑ) proxy; this is distinct from realised genomic inbreeding. The genetic contributions literature provides context (Wray & Thompson 1990; Woolliams &
Bijma 2000). Proving the bound — coupling the spectral gap with b's projection onto the top
eigenspace, robust to the tail's curvature — is the open problem, the meeting point of optimisation
(the perturbed-cone face) and quantitative genetics (effective lineages).

**An exact handle on the effective support.** Write the *effective* number of contributors as the
participation ratio PR(c\*) = 1/‖c\*‖² (it equals |S| for uniform contributions, less when they are
uneven). At a positive active cap k > 0 with Σc\* = 1 an identity holds: **PR(c\*) = R(c\*)/k**, where
R(c\*) = c\*ᵀGc\*/c\*ᵀc\* is the optimum's own Rayleigh quotient (since R·‖c\*‖² = c\*ᵀGc\* = k).
Hence **PR(c\*) ≤ λ₁/k** for every such active positive cap, and a bound on *effective* support reduces to bounding R(c\*) —
the directional coancestry cost the optimum is forced to incur. This is exactly where (A2) bites:
the measured cases where b avoids the dominant directions have small R(c\*) (observed: R ≈ 0.02–0.4 on the structured cases versus λ₁ ≈ 45–143, giving PR ≈ 2–12), whereas
G = εI forces R = ε = λ₁. At an active cap, PR = ε/k lies in (n−1, n) throughout
the full-support interval ε/n < k < ε/(n−1), and approaches n as k decreases to ε/n.
The total support is already n throughout that interval; PR = n holds only for uniform
contributions. For example, n = 3 and k = 5ε/12 give |S| = 3 but PR = 12/5.
The residual open piece is now sharp and scalar — bound R(c\*)
a priori under (A1)+(A2). (The historical raw cardinality |S| ran ≈ 1.5–2.5× PR in these cases. In general PR ≤ |S|, so an upper bound on PR does not upper-bound |S| without another assumption on contribution sizes.) Verified in `bound_kernel.py`.

**Attacking this residual.** For a positive-definite G, the simplex relaxation without lower or upper bounds has, in the active-kinship case with λ > 0, the Rayleigh-quotient expression
R_relax = (βᵀG⁻¹β)/(βᵀG⁻²β) with β = b − μ𝟙 — a λ-average weighted by βᵢ²/λᵢ², so small under
(A2). But it predicts R(c\*) only when c ≥ 0 is near-inactive: on the real panels the
unconstrained optimum carries large negative entries (min c_relax ≈ −1 to −2.6) and R_relax
falls 60–300× below R(c\*). So **nonnegativity is the irreducible core** — fittingly, |S| is the
cardinality of a nonnegativity-constrained QP optimum, hard in general — and a proof must be
genuinely c ≥ 0-aware, not a relaxation shortcut. The clean a priori handle (the relaxation
formula) and the irreducible obstacle (c ≥ 0) are now both identified.

## Numerical evidence (reproducible)

- `bound_validation.py` — solver vs SciPy; ridged low-rank observations (not a test of an ε = 0 vertex theorem); empirical n-sweep; the effective-rank
  behaviour; the ridge sweep; a spectral-decay sweep; and (block 7) the Theorem-2
  counterexample; and analytic finite-cap / full-support witnesses in block 8.
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
