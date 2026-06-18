# Exact optimum contribution selection in milliseconds: a matrix-free, support-first solver validated against optiSel and AlphaMate on real genomic data

**Adel Kaleche**

*Affiliation to be completed.*


## Abstract

Optimum contribution selection (OCS) maximises expected genetic gain subject to a cap
on the average coancestry of the next generation, the central tool for balancing gain
against the loss of diversity in breeding programmes. At genomic scale the relationship
matrix G is dense and n×n, and the standard exact tool (the optiSel R package) solves a
quadratically-constrained program whose cost grows steeply with the number of candidates,
while the most widely used alternative (AlphaMate) is a stochastic heuristic. We present
**support-first**, an exact solver that exploits two structural facts of OCS: the optimal
contribution vector is supported on a tiny subset of candidates (typically 2–50 of n),
and the only constraint that is hard to satisfy is non-negativity. Support-first grows the
support by an active-set / column-generation rule and solves each fixed-support subproblem
— maximise a linear objective over an ellipsoid intersected with the affine sum (and sex)
constraints — in closed form, with no inner iterative solver. It is **matrix-free**: it
never forms or stores the dense n×n G, computing G·c from the raw centred genotype matrix Z
— the memory and large-n enabler — while the speed advantage over existing solvers comes
from the **tiny active set** (a handful of cheap iterations instead of a full conic
interior-point solve), not from the matrix-free product itself. On real data — a CIMMYT
wheat panel (n=599), a PIC pig (n=3534, 52k SNP), and a heterogeneous-stock mouse panel
(n=1814, with real sex) — support-first reaches the exact optimum (agreeing with a conic
interior-point solver to 1e-8, and saturating the kinship bound those solvers leave slightly
slack) while running 90×–2280× faster (0.008 s vs 6.96 s on the sexed mouse instance), and
~37000× faster than a general conic solver (Clarabel) at n=10000. Against AlphaMate the exact
optimum dominates at every matched coancestry (mouse: Δgain +0.004 at the 45° tradeoff,
larger elsewhere), at a small fraction of the run time (882 s CPU for the frontier vs ≤1.1 s
per point).
Support-first makes exact, reproducible OCS practical at genomic scale on a laptop.

## 1. Introduction

Genetic improvement of livestock and crops must hold two objectives in tension:
maximising the genetic merit of the next generation, and conserving the genetic
diversity whose erosion drives inbreeding depression and forecloses future
response to selection. Optimum contribution selection (OCS; Meuwissen 1997) makes
that trade-off explicit. Given each candidate's estimated breeding
value and the matrix of relationships among candidates, OCS chooses the
proportional genetic contributions **c** that maximise expected gain **bᵀc**
subject to a cap on the mean coancestry of the offspring, **cᵀGc ≤ k**, with the
contributions non-negative and summing to one — and, in any real mating scheme,
split so that sires and dams each supply half. OCS is the de facto framework for
managing diversity in both commercial breeding and conservation programmes, and
the genomic relationship matrix (VanRaden 2008) has largely replaced its
pedigree-based predecessor as the **G** it constrains.

With genomic data this is a convex quadratically constrained quadratic program,
equivalently a second-order cone program, and two costs come to dominate as
candidate numbers grow. The first is the relationship matrix itself: dense and
n×n, it costs O(n²) to store and O(n²m) to build from m markers, and at tens of
thousands of candidates it no longer fits in a workstation's memory. The second
is the conic solve. Yet the optimal contribution vector is almost always
supported on a handful of individuals — a sparsity that, as we show, persists as
n grows — and general-purpose solvers exploit neither this nor the cheap
structure of the kinship constraint.

The tools in use span the exact and the heuristic, but share a template: assemble
the full relationship matrix and hand the whole problem to a generic optimiser
over every candidate. Meuwissen's original method (1997), and its successor
Gencont2 (Dagnachew & Meuwissen 2016), iterate the Lagrangian
conditions over the full candidate set. The optiSel package (Wellmann 2019), the
standard exact tool, casts OCS as a cone program for a primal-dual interior-point
solver. Semidefinite formulations (Pong-Wong & Woolliams 2007) and recent
ADMM/JuMP solvers (Waldmann 2025) follow suit, the latter forming the dense **G**
explicitly and truncating small contributions only after the fact. AlphaMate
(Gorjanc & Hickey 2018) trades exactness for a differential-evolution heuristic
that also allocates matings. The lone attempt to exploit sparsity, Yamashita,
Mullin & Safarina (2018), exploits the sparsity of the *pedigree inverse* **A⁻¹**
— a property of the data matrix — inside a full-candidate interior-point solve;
it is neither the sparsity of the *solution* nor a genotype-based, matrix-free
method. None of these exploits the two facts that make genomic OCS cheap: that
the optimal support is tiny and bounded in n, and that within a fixed support the
kinship-constrained subproblem is solvable in closed form. None avoids
materialising **G**, which becomes the binding constraint — in memory and in
build time — exactly at the population sizes where OCS matters.

We present **support-first**, an exact OCS solver built on those two facts. It
maintains a small working support, seeded with the best candidate of each sex and
grown by an active-set, reduced-cost rule: candidates are priced against the
current multipliers and the one most worth adding is brought in, candidates that
turn negative are dropped, all toward a single fixed coancestry cap k. Each
fixed-support subproblem — maximise a linear form over an ellipsoid intersected
with the affine sum (and sex) constraints — is solved in closed form by
eliminating the equality constraints and reducing to a scalar quadratic in the
binding multiplier, with no inner iterative solver. The kinship products **G·c**
are formed matrix-free from the centred genotype matrix **Z** as ε**c** +
**Z**(**Zᵀc**)/s, so the solver never builds or stores the dense n×n **G**.

The individual ingredients are classical, and we claim none of them. Active-set
methods that track a small working set of held assets descend from Markowitz's
critical line algorithm (1956) for long-only mean–variance selection, of which
OCS is the kinship-constrained analogue; the closed form for a linear objective
over an ellipsoid intersected with linear constraints is a constrained-eigenvalue
/ secular-equation result (Gander, Golub & von Matt 1989); and the matrix-free
product **Z**(**Zᵀc**) is standard in genomic prediction (VanRaden 2008;
Legarra & Misztal 2008). Our contribution is their **synthesis, specialised to
genomic OCS**: reduced-cost column generation that grows a tiny support toward a
*single* coancestry cap — as opposed to the critical line algorithm's parametric
sweep of the *entire* efficient frontier — with the kinship product evaluated
matrix-free, so that the cost of a solve follows the support size and the marker
count rather than the dense n×n matrix. To our knowledge no prior OCS method
exploits solution-support sparsity, and none is genotype-matrix-free; the
benefit, as the matrix-free product is what makes a support-bounded solve pay off
at genomic n, is the combination rather than any one part.

On real data — a CIMMYT wheat panel (n = 599), a PIC pig (n = 3534, 52k SNP) and a
heterogeneous-stock mouse panel (n = 1814, with real sex) — support-first reaches
the exact optimum, agreeing with the conic optimum to 1e-8 and saturating the
kinship bound that the interior-point methods leave slightly slack, while running
90×–2280× faster, and ~37000× faster than a general conic interior-point solver at
n = 10000. Against AlphaMate it attains strictly higher gain at every matched
coancestry, at a small fraction of the run time. Across synthetic populations the
optimal support stays 14–19 as n grows from
1000 to 40000, while the dense **G** the alternatives must form reaches 11.9 GiB
— a 40× larger footprint than **Z**, past the working memory of an 8–16 GB laptop
— in a regime where support-first still solves in under 0.1 s. Support-first makes
exact, reproducible optimum contribution
selection practical at genomic scale on a laptop.

## 2. Methods

### Optimum contribution selection

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

### Two structural facts

Support-first rests on two properties of the optimum **c\***. First, at a binding
kinship cap **c\*** is supported on a small set S = {i : c\*ᵢ > 0}, and S is
empirically bounded as n grows (Results). Second, on the face where the c ≥ 0
constraints inactive on S are dropped, the problem restricted to S — maximise a
linear form over an ellipsoid intersected with the affine constraints — has a
closed-form solution. The whole difficulty is therefore identifying S; everything
else is a small direct solve.

### Closed form on a fixed support

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

### The support-first algorithm

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

### Matrix-free kinship products

Every full kinship product is formed without materialising **G**:

  **Gc = ε c + Z(Zᵀc) / s**,

two matrix–vector products against the n×m genotype matrix, at O(n·m) time and
O(n·m) memory for **Z** versus O(n²) for a dense **G**. The restricted **G_S** is
assembled only on the support. This is the enabler at large n, where the dense
n×n **G** is infeasible to store or costly (O(n²m)) to build; it is not an
inner-loop speed-up when m > n, and the order-of-magnitude advantage over the
conic baselines is algorithmic — the tiny active set — independent of this choice.

### Implementation and reproducibility

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

### Data and baselines

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

## 3. Results

### Exactness

Support-first is exact by construction: for a convex program a feasible point at
which no candidate has a positive reduced cost satisfies the KKT conditions, and
this is what the active set certifies on termination (and what the unit tests
assert across a range of kinship caps). Empirically, on the optiSel `Cattle`
example (Angler cattle, n = 268) with the package's own real recorded sex,
support-first and optiSel select the *same* 36-individual support and agree on the
contributions to within 3×10⁻⁴ (maximum contribution 0.0664 vs 0.0661). Where they
differ is instructive: support-first reaches the exact optimum — gain 1.8315 with
the kinship constraint saturated at its bound — while optiSel's interior-point
solver stops about 0.4 % short (gain 1.8246, group kinship 0.0576 against a bound
of 0.0578), leaving feasible gain on the table. Against an independent conic
interior-point solver (Clarabel) on synthetic data the two agree to a gain
difference below 10⁻⁸ across kinship caps, and the Rust implementation of the
sexed solver reproduces a NumPy reference optimum to 1.5×10⁻¹⁴ on a binding
instance. Support-first is therefore at least as accurate as the domain tool, and
exact where the interior-point methods are merely close.

### Speed against a generic conic solver

On synthetic genomic instances support-first is faster than the Clarabel conic
solver by a factor that grows steeply with n (Table 1): from ~130× at n = 1000 to
~37000× at n = 10000, where Clarabel takes 26 minutes and support-first 43 ms. The
scaling is structural — Clarabel factors a dense O(n³) KKT system at every
interior-point iteration, whereas support-first performs a near-constant number of
matrix–vector products (3–6 here), each O(n·m).

**Table 1.** support-first vs Clarabel (synthetic, identical optimum).

| n | Clarabel | support-first | speed-up |
|---|---|---|---|
| 1000 | 1.40 s | 0.011 s | 126× |
| 2000 | 10.6 s | 0.022 s | 472× |
| 5000 | 160 s | 0.036 s | 4474× |
| 10000 | 1579 s | 0.043 s | **37090×** |

### Speed against the domain tool optiSel

A generic solver is a soft target; the informative comparison is against optiSel,
the standard exact OCS tool, on its own formulation. After extending support-first
to the two sex-equality constraints, it returns the same optimum as optiSel while
running 90–2280× faster on three real marker panels (Table 2), and 22–97× faster
on structured synthetic populations. The largest margin is on the PIC pig panel
(n = 3534), where optiSel takes 55 s and support-first 24 ms.

**Table 2.** support-first vs optiSel (R/`cccp`), same optimum throughout.

| dataset | n | support-first | optiSel | speed-up |
|---|---|---|---|---|
| synthetic (structured) | 1000 | 0.09 s | 2.0 s | 22× |
| synthetic (structured) | 2000 | 0.29 s | 11.9 s | 41× |
| synthetic (structured) | 5000 | 1.47 s | 143 s | 97× |
| CIMMYT wheat (real GRM) | 599 | 0.007 s | 0.63 s | 90× |
| PIC pig (real GRM, 52k SNP) | 3534 | 0.024 s | 54.8 s | **2280×** |
| HS mouse (real GRM, real sex) | 1814 | 0.008 s | 6.96 s | 870× |

The mouse row is the true sexed OCS, on a genuine recorded sex (934 males, 880
females). Caveats carried to the Discussion: the support-first timings are the
NumPy prototype and optiSel is R, so the factor reflects the algorithm (the active
set on a tiny support) rather than the language; sex is real only for the mouse
panel; and the selection criterion is a recorded phenotype or EBV used as a proxy
for a genomic breeding value.

### Comparison with the heuristic AlphaMate

AlphaMate optimises a related but distinct problem — discrete mate allocation by a
stochastic evolutionary algorithm — so the fair comparison evaluates *its*
contribution vector in our metric and pits it against support-first at the same
group coancestry. On the mouse panel support-first attains strictly higher genetic
gain at every point of AlphaMate's frontier (Table 3): a small margin at the
angle-45° trade-off optimum (Δgain +0.004) and larger margins at the
diversity-control and gain-maximising corners. AlphaMate, being a discrete
heuristic, leaves gain on the table everywhere. It was also markedly fragile on
real genomic data: a successful run required six configurations and three distinct
work-arounds — capping matings below n; restoring the full parent set, to avoid a
setup segmentation fault that the reduced parent count triggered; and positively
shifting the selection criterion, to undo a value-over-maximum sign inversion that
made the heuristic maximise in the wrong direction on the centred, negative EBVs —
whereas support-first and optiSel ran unmodified. AlphaMate computed its whole
frontier in 882 s of CPU time (an emulated x86 binary, no native build existing);
support-first traces the exact frontier at ≤ 1.1 s per point.

**Table 3.** Genetic gain at matched group coancestry, mouse panel (scored in the
same metric from each method's contribution vector).

| group coancestry cᵀKc | AlphaMate | support-first | Δgain |
|---|---|---|---|
| 0.000272 (AlphaMate min-coancestry) | −0.45885 | **−0.37843** | +0.080 |
| 0.001317 (AlphaMate 45° optimum) | −0.35748 | **−0.35318** | +0.004 |
| 0.007574 (AlphaMate max-criterion) | −0.34066 | **−0.32240** | +0.018 |

### Scaling and the matrix-free advantage (Figure 1)

The cost advantage holds — and grows — at scale. Sweeping the candidate count from
1000 to 40000 at a fixed marker panel (m = 1000) under a binding kinship cap, the
optimal support stays between 14 and 19 (Figure 1A) and the matrix-free solve
stays under 0.1 s. The dense relationship matrix that every other solver must
materialise tells the opposite story. Merely *building* it costs O(n²m): 3.6 s at
n = 30000, already 63× the entire support-first solve at that size, and rising
quadratically (Figure 1A). *Storing* it costs O(n²): 11.9 GiB at n = 40000, where
the matrix-free Z footprint is 0.30 GiB and the dense matrix no longer fits in a
laptop's working memory (Figure 1B). The dense matrix becomes the binding
constraint, in setup time and in memory, exactly in the regime where support-first
remains cheap — because the solver's cost follows the support and the marker count,
never the n×n matrix.

### Support behaviour

The advantage rests on the small, bounded support. At the mouse operating
coancestry (0.0346) the optimum places weight on 19 of the 1814 candidates; the
support enlarges only as the cap is driven toward zero — about 1163 individuals to
force group coancestry near 0, where the solution must spread over much of the
population to minimise relatedness. Across the synthetic sweep the support stays in
the low tens as n grows forty-fold (Figure 1A), which is what makes the per-solve
cost scale with the support rather than with n.

## 4. Discussion

Support-first makes exact optimum contribution selection cheap at genomic scale.
On the panels tested it reaches the same optimum as the domain's exact tool while
running two to three orders of magnitude faster, and it stays cheap precisely
where the dense relationship matrix that every other solver forms becomes
infeasible — the large-candidate regime that motivates genomic OCS in the first
place. Because the method is exact and deterministic — a Karush–Kuhn–Tucker–
certified active set rather than a stochastic search — its output is reproducible
to the last digit, unlike the heuristic mate-selection tools it is measured
against.

The method's ingredients are individually classical, and we claim only their
synthesis. Active-set tracking of a small working set descends from the critical
line algorithm for long-only mean–variance selection; the per-support closed form
is a constrained-eigenvalue result; the matrix-free product is standard in genomic
prediction. The contribution is their combination, specialised to OCS:
reduced-cost column generation that grows a tiny support toward a single coancestry
cap, with the kinship product evaluated matrix-free so that cost follows the
support and the marker count rather than the n×n matrix. Framed against the closest
prior work, this is neither the pedigree-matrix sparsity of an interior-point OCS
solver nor the full-population dense solves of the standard tools.

Several limitations bound these results, and we state them plainly. First, the
public panels carry a recorded phenotype or EBV that we use as a stand-in for a
true genomic breeding value; the speed and exactness results are unaffected, but
the contribution vectors themselves are illustrative, not breeding
recommendations. Second, a genuine recorded sex is available only for the mouse
panel; elsewhere we impose an arbitrary balanced split, so only the mouse result
exercises the true sexed constraints on real data. Third, the head-to-head
timings compare a NumPy prototype against R/optiSel: the order-of-magnitude gap is
algorithmic — both realise the same active set — but a single-language comparison
would place it beyond doubt, and our own measurements are explicit that the
matrix-free product is *not* an inner-loop speed-up when markers outnumber
candidates (m > n), where streaming the genotype matrix costs more than a resident
dense product. The matrix-free route is the memory and large-n enabler; the speed
advantage over the conic solvers is the small active set. Fourth, the
boundedness of the optimal support as n grows is, here, an empirical observation
across a synthetic sweep and the real panels, not a theorem; how the bound depends
on population structure and linkage disequilibrium is open. Finally, we solve a
single quadratic kinship constraint and return continuous contributions, whereas
breeding programmes may impose several constraints at once — multiple relationship
matrices, own-relationship caps, group-specific limits — and ultimately require an
integer mate allocation, which AlphaMate provides and support-first does not.

Each limitation points to an extension. The closed-form-per-support core already
admits more than one equality constraint through the same elimination
P = A G⁻¹ Aᵀ; additional active quadratic caps turn the scalar root-finding into a
small low-degree system rather than changing the structure, so multiple kinship
constraints are within reach. A rounding or branch-and-bound layer over the
continuous optimum would add mate allocation while keeping the exact relaxation as
a tight bound — combining support-first's exactness with the discrete plan that
heuristics target directly. The matrix-free crossover in m versus n deserves an
explicit characterisation, as does the support bound itself. And the method should
be validated on datasets carrying true genomic breeding values and on populations
an order of magnitude larger, where the memory gap of Figure 1B turns from visible
into decisive.

Support-first does one thing — exact, single-constraint, continuous optimum
contribution selection — at a scale and speed that bring genomic OCS within reach
of a laptop and a reproducible script. That narrow, verified claim, and the open
extensions it invites, are the contribution.

## References

### OCS methods and tools

1. **Meuwissen, T.H.E. (1997).** Maximizing the response of selection with a
   predefined rate of inbreeding. *Journal of Animal Science* 75(4):934–940.
   DOI 10.2527/1997.754934x.

2. **Dagnachew, B.S. & Meuwissen, T.H.E. (2016).** A fast Newton–Raphson based
   iterative algorithm for large scale optimal contribution selection. *Genetics
   Selection Evolution* 48(1):70. DOI 10.1186/s12711-016-0249-2.

3. **Pong-Wong, R. & Woolliams, J.A. (2007).** Optimisation of contribution of
   candidate parents to maximise genetic gain and restricting inbreeding using
   semidefinite programming. *Genetics Selection Evolution* 39(1):3–25.
   DOI 10.1186/1297-9686-39-1-3.

4. **Wellmann, R. (2019).** Optimum contribution selection for animal breeding and
   conservation: the R package optiSel. *BMC Bioinformatics* 20(1):25.
   DOI 10.1186/s12859-018-2450-5.

5. **Gorjanc, G. & Hickey, J.M. (2018).** AlphaMate: a program for optimizing
   selection, maintenance of diversity and mate allocation in breeding programs.
   *Bioinformatics* 34(19):3408–3411. DOI 10.1093/bioinformatics/bty375.

6. **Waldmann, P. (2025).** Genomic optimum contribution selection and mate
   allocation using JuMP. *Bioinformatics Advances* 5(1):vbaf259.
   DOI 10.1093/bioadv/vbaf259.

7. **Yamashita, M., Mullin, T.J. & Safarina, S. (2018).** An efficient
   second-order cone programming approach for optimal selection in tree breeding.
   *Optimization Letters* 12(7):1683–1697. DOI 10.1007/s11590-018-1229-y;
   arXiv:1506.04487. — Exploits *pedigree-inverse* sparsity, not solution support.

### Genomic relationships and matrix-free products

8. **VanRaden, P.M. (2008).** Efficient methods to compute genomic predictions.
   *Journal of Dairy Science* 91(11):4414–4423. DOI 10.3168/jds.2007-0980.

9. **Legarra, A. & Misztal, I. (2008).** Technical note: Computing strategies in
   genome-wide selection. *Journal of Dairy Science* 91(1):360–366.
   DOI 10.3168/jds.2007-0403.

### Active-set and closed-form lineage (conceded prior art)

10. **Markowitz, H. (1956).** The optimization of a quadratic function subject to
    linear constraints. *Naval Research Logistics Quarterly* 3(1–2):111–133.
    DOI 10.1002/nav.3800030110. — The critical line algorithm.

11. **Gander, W., Golub, G.H. & von Matt, U. (1989).** A constrained eigenvalue
    problem. *Linear Algebra and its Applications* 114–115:815–839.
    DOI 10.1016/0024-3795(89)90494-1. — The per-support closed form (secular equation).

### Software and data

12. **Goulart, P.J. & Chen, Y. (2024).** Clarabel: an interior-point solver for
    conic programs with quadratic objectives. arXiv:2405.12762. (Peer-reviewed
    version: *Mathematical Programming Computation*, 2026,
    DOI 10.1007/s12532-026-00320-7.) — Used only as an independent cross-check oracle.

13. **Quiñones, S.** faer: a linear-algebra library for the Rust programming
    language. Software, version 0.24.0 (pinned in `Cargo.toml`);
    <https://github.com/sarah-quinones/faer-rs>. (JOSS submission under review; no
    published DOI at time of writing.)

14. **Pérez, P. & de los Campos, G. (2014).** Genome-wide regression and
    prediction with the BGLR statistical package. *Genetics* 198(2):483–495.
    DOI 10.1534/genetics.114.164442. — Source of the wheat (CIMMYT) and mouse panels.

15. **Cleveland, M.A., Hickey, J.M. & Forni, S. (2012).** A common dataset for
    genomic analysis of livestock populations. *G3: Genes|Genomes|Genetics*
    2(4):429–435. DOI 10.1534/g3.111.001453. — The PIC pig panel.

## Figure

**Figure 1.** Matrix-free vs dense-G scaling (`research/fig_scaling.pdf`); described in Results, section *Scaling and the matrix-free advantage*.
