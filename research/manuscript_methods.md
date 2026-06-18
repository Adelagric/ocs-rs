# Methods (manuscript draft)

> Companion to `manuscript_intro.md`. Math matches the Rust implementation
> (`src/support_first.rs`, `solve` and `solve_sexed`). Display equations are in
> readable plain math, to be set in LaTeX for the manuscript.

## Optimum contribution selection

For n selection candidates with estimated breeding values **b** ∈ ℝⁿ and a
genomic relationship matrix **G** ∈ ℝⁿˣⁿ, optimum contribution selection chooses
proportional genetic contributions **c** that solve

  maximise **bᵀc**  subject to  **Ac = d**,  **c ≥ 0**,  **cᵀGc ≤ k**.

The affine constraints **Ac = d** encode the contribution budget. In the simplex
form a single row **A = 𝟙ᵀ**, **d = 1**, imposes Σcᵢ = 1. In the *sexed* form —
the true OCS, since each mating draws one parent of each sex — **A** is the 2×n
sex-incidence matrix (row 1 the male indicator, row 2 the female indicator) and
**d** = (½, ½)ᵀ, imposing Σ_{males}cᵢ = Σ_{females}cᵢ = ½. The kinship bound
**cᵀGc ≤ k** caps the mean coancestry of the offspring. We take **G** as the
VanRaden genomic relationship matrix, **G = ZZᵀ/s + εI**, where **Z** is the n×m
matrix of genotypes centred by twice the allele frequencies, s = 2Σⱼpⱼ(1−pⱼ), and
ε a small ridge ensuring positive definiteness (the same ridged matrix the conic
baselines constrain, so all solvers face an identical problem). The program is a
convex quadratically constrained quadratic program, equivalently a second-order
cone program.

## Two structural facts

Support-first rests on two properties of the optimum **c\***. First, at a binding
kinship cap **c\*** is supported on a small set S = {i : c\*ᵢ > 0}, and S is
empirically bounded as n grows (Results). Second, on the face where the c ≥ 0
constraints inactive on S are dropped, the problem restricted to S — maximise a
linear form over an ellipsoid intersected with the affine constraints — has a
closed-form solution. The whole difficulty is therefore identifying S; everything
else is a small direct solve.

## Closed form on a fixed support

Fix a support S of size |S|, and write **G_S**, **b_S**, **A_S** for the
restrictions to S (**A_S** is q×|S| with q = 1 or 2 equality rows). At the optimum
of the restricted problem the kinship constraint is active, and the stationarity
conditions read

  **b_S = A_Sᵀ μ + 2λ G_S c_S**,   **A_S c_S = d**,   **c_Sᵀ G_S c_S = k**,

with multipliers **μ** ∈ ℝ^q and λ ≥ 0. Solving stationarity for **c_S** gives
**c_S = (1/2λ) G_S⁻¹(b_S − A_Sᵀμ)**. Imposing the affine constraints eliminates
**μ**: with the q×q matrix **P = A_S G_S⁻¹ A_Sᵀ** and vector **q_v = A_S G_S⁻¹ b_S**,

  **μ = P⁻¹(q_v − 2λ d)**,   so   **c_S = g/(2λ) + h**,

where **g = G_S⁻¹(b_S − A_Sᵀ P⁻¹ q_v)** and **h = G_S⁻¹ A_Sᵀ P⁻¹ d** satisfy
**A_S g = 0** and **A_S h = d**, so **A_S c_S = d** for every λ. Substituting into
the active ellipsoid **c_Sᵀ G_S c_S = k** and writing α = **gᵀG_S g**,
β = **gᵀG_S h**, γ = **hᵀG_S h** yields a single scalar quadratic in the
multiplier:

  **4(γ − k) λ² + 4β λ + α = 0.**

We take the positive root maximising **b_Sᵀ c_S**. Because **G_S g = b_S − A_SᵀP⁻¹q_v**
and **G_S h = A_SᵀP⁻¹d** are already available, α, β, γ cost only inner products —
no second factorisation. A restricted solve is thus: one Cholesky factorisation of
**G_S**, back-substitutions for the q+1 right-hand sides **[A_Sᵀ | b_S]**, a q×q
inverse (q ≤ 2), and a scalar quadratic — no inner iteration. The simplex case
(q = 1, **A_S = 𝟙ᵀ**) is identical and reduces to a scalar quadratic in μ. When the
equalities fully determine **c_S** (|S| = q + 1) the contributions are forced —
(½, ½) in the sexed case — and we set λ = 0 (the binding subcase is measure-zero,
as for a singleton support). The solve returns "infeasible on S" when S lacks a
sex or the ellipsoid does not meet the affine hull, the signal to enlarge S.

## The support-first algorithm

The support is found by active-set column generation. It is seeded with the best
candidate of each sex, S = {argmax_{males} b, argmax_{females} b}. Each iteration
solves the closed form on the current S and branches:

- **Infeasible on S** (the support is too related to meet the cap k): add the
  least-related candidates, those with the smallest **(Gc)ⱼ**, in a chunk that
  doubles |S|. This feasibility phase, when done one candidate at a time, was the
  dominant cost; chunking cut total iterations roughly thirty-fold with an
  identical optimum.
- **A contribution is negative**: drop every i with c_Sᵢ ≤ 0, marking them taboo
  to re-entry until the next genuine improvement, and re-solve; if dropping
  empties a sex, its best candidate is re-seeded.
- **Otherwise**: price each candidate j ∉ S by its reduced cost
  **rⱼ = bⱼ − μ_{sex(j)} − 2λ (Gc)ⱼ**, and add the maximiser if rⱼ > tol;
  if none is positive, stop. A genuine addition clears the taboo set.

Because the problem is convex, a feasible point at which no candidate has a
positive reduced cost satisfies the Karush–Kuhn–Tucker conditions and is the
global optimum; correctness is therefore not heuristic. The taboo set together
with a degeneracy threshold scaled to the coefficient magnitudes guarantees finite
termination. Reduced-cost pricing toward a *single* fixed cap k distinguishes the
support update from the critical line algorithm, which changes the active set by
±1 along a swept multiplier to trace the whole efficient frontier.

## Matrix-free kinship products

Every full kinship product is formed without materialising **G**:

  **Gc = ε c + Z(Zᵀc) / s**,

two matrix–vector products against the n×m genotype matrix, at O(n·m) time and
O(n·m) memory for **Z** versus O(n²) for a dense **G**. The restricted **G_S** is
assembled only on the support. This is the enabler at large n, where the dense
n×n **G** is infeasible to store or costly (O(n²m)) to build; it is not an
inner-loop speed-up when m > n, and the order-of-magnitude advantage over the
conic baselines is algorithmic — the tiny active set — independent of this choice.

## Implementation and reproducibility

The solver is implemented in Rust over a pure-Rust stack: dense linear algebra
(the support Cholesky and back-substitutions) via `faer`, with no system
BLAS/LAPACK dependency, no `unsafe`, and a seeded reproducible RNG for the
synthetic generators. A conic interior-point solver (`clarabel`) is included only
as an independent oracle for cross-checking. Correctness is asserted by
KKT-certificate unit tests across a range of caps k — a `Solved` result is always
feasible, on the budget, and sex-split — and by cross-language agreement: on a
binding instance the Rust `solve_sexed` reproduces a NumPy reference optimum to a
gain difference of 1.5×10⁻¹⁴. The build is gated on `cargo fmt`,
`cargo clippy -D warnings`, and the test suite.

## Data and baselines

The method is evaluated on three public marker panels: a CIMMYT wheat panel
(`BGLR`, n = 599), a PIC pig panel (n = 3534, 52k SNP, with real estimated
breeding values), and a heterogeneous-stock mouse panel (`BGLR`, n = 1814, with
recorded sex, 934 males and 880 females, using body-mass index as the selection
criterion). For each, **G** is the VanRaden matrix with ridge ε = 10⁻⁵. Baselines
are optiSel (R, the `cccp` Nesterov–Todd interior-point solver; the exact domain
reference), Clarabel (the conic cross-check), and AlphaMate (Fortran differential
evolution; run from its Linux binary under emulation, as no macOS build exists and
the source is locked to an Intel toolchain). Honest caveats carried into the
Discussion: recorded sex is available only for the mouse panel (an arbitrary
balanced split is used elsewhere — wheat is autogamous and the PIC panel ships no
usable sex); the selection criterion **b** is a recorded phenotype or EBV standing
in for a true genomic breeding value on these public panels; and the prototype
timings compare a NumPy support-first against R/optiSel, an algorithmic rather
than a language difference, with the single-language Rust timing reported
separately.
