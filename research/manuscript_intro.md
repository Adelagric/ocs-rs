# Introduction (manuscript draft)

> Draft prose for the bioRxiv preprint. English, manuscript-ready. Numbers are
> locked from the benchmark tables (see PREPRINT.md ledger). Citations are in
> (Author Year) form; a few bodies were paywalled at drafting time and are
> flagged in PREPRINT.md — verify the exact solver each used before camera-ready.

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
the exact optimum, agreeing with the conic optimum to 1e-8; at matched realised
coancestry it agrees with the interior-point methods, and where they halt just
inside the constraint support-first reaches the boundary, so its small edge is the
diversity budget they leave unspent rather than a different optimum — all while
running 90×–2280× faster, and ~37000× faster than a general conic interior-point
solver at n = 10000. Against AlphaMate, a heuristic for the distinct problem of
discrete mate allocation, the exact optimum is no worse at matched coancestry on
the continuous relaxation the two share, at a small fraction of the run time. Across synthetic populations the
optimal support stays 14–19 as n grows from
1000 to 40000, while the dense **G** the alternatives must form reaches 11.9 GiB
— a 40× larger footprint than **Z**, past the working memory of an 8–16 GB laptop
— in a regime where support-first still solves in under 0.1 s. Support-first makes
exact, reproducible optimum contribution
selection practical at genomic scale on a laptop.
