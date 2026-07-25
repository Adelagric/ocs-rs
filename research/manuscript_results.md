# Results (manuscript draft)

> Companion to `manuscript_intro.md` / `manuscript_methods.md`. Numbers locked
> from the benchmark tables (`research/repro/table2_numbers.md`; Figure 1 =
> `research/fig_scaling.pdf`). Table 2 reports two support-first timings — the active
> set given a dense `G` (the regime optiSel runs in), and the shipped matrix-free
> solver end to end — both against optiSel (R/`cccp`).

## Exactness

Support-first is exact by construction: for a convex program a feasible point at
which no candidate has a positive reduced cost satisfies the KKT conditions, and
this is what the active set certifies on termination (and what the unit tests
assert across a range of kinship caps). Empirically, on the optiSel `Cattle`
example (Angler cattle, n = 268) with the package's own real recorded sex,
support-first and optiSel select the *same* 36-individual support and agree on the
contributions to within 3×10⁻⁴ (maximum contribution 0.0664 vs 0.0661). They differ
only at the kinship boundary: the OCS optimum sits on the cap (gain is monotone in
it), and support-first reaches it (cᵀGc = k), whereas optiSel's interior-point
solver halts at its convergence tolerance just inside the feasible region (group
kinship 0.0576 against a bound of 0.0578). At matched realised coancestry the two
agree — on the mouse panel both reach gain −0.29754 at coancestry 0.0346 — so
support-first's small edge at its own operating point is the diversity budget
optiSel leaves unspent on the boundary, not a different optimum: the difference is
boundary-saturation versus an interior stop. Against an independent conic
interior-point solver (Clarabel) on synthetic data the two agree to a gain
difference below 10⁻⁸ across kinship caps, and the Rust implementation of the
sexed solver reproduces a NumPy reference optimum to 1.5×10⁻¹⁴ on a binding
instance. Support-first is therefore at least as accurate as the domain tool, and
exact where the interior-point methods are merely close.

## Speed against a generic conic solver

This comparison also settles a question the optiSel head-to-head cannot: both
solvers here are Rust, called in the same process on the same generated data, so no
language difference can be mistaken for an algorithmic one. On synthetic genomic
instances support-first is faster than the Clarabel conic solver by a factor that
grows steeply with n (Table 1): from ~290× at n = 1000 to ~76000× at n = 10000, where
the conic route takes half an hour and support-first 24 ms. The scaling is
structural — Clarabel factors a dense O(n³) KKT system at every interior-point
iteration, whereas support-first performs a near-constant number of matrix–vector
products (3–6 here), each O(n·m).

**Table 1.** support-first vs Clarabel, identical optimum (agreeing to 1.2×10⁻⁹–
4.8×10⁻⁹ in gain), one serial sweep on an idle machine. Both solvers are Rust in one
process on the same synthetic instance (m = 20000, cap k = 0.6·mean diag **G**), so
the comparison carries no language confound.

| n | m | Clarabel¹ | support-first | speed-up |
|---|---|---|---|---|
| 1000 | 20000 | 1.70 s | 0.006 s | 290× |
| 2000 | 20000 | 12.66 s | 0.012 s | 1093× |
| 5000 | 20000 | 182.9 s | 0.018 s | 9917× |
| 10000 | 20000 | 1792 s | 0.024 s | **76179×** |

¹ the whole conic route: forming **G** (7.0 s at n = 10000), its Cholesky, cone
assembly and the interior-point solve. The solve alone is 1.60 / 12.31 / 180.9 /
1784.2 s, i.e. 274×–75831×. support-first's column is its entire run, since it builds
nothing.

Two honest qualifications. The support-first timings are deterministic at this
precision (three runs at n = 1000 and at n = 2000 reproduced 0.006 s and 0.012 s to
the millisecond) while Clarabel varies a few percent; the large-n cells are single
runs, as Clarabel needs half an hour at n = 10000. And the cap here is loose, so the
optimal support is 2–5 — the easy end of the spectrum, which is why so few kinship
products suffice. The tight-cap regime, where the support runs to 59–168, is the
synthetic half of Table 2, and the factor shrinks there accordingly; we report both
rather than only the favourable one.

## Speed against the domain tool optiSel

Table 1 isolates the algorithm; what matters to a breeder is the comparison against
optiSel, the standard exact OCS tool, on its own formulation. After extending support-first
to the two sex-equality constraints, it returns the same optimum as optiSel. The
shipped matrix-free solver — forming no `G` at all — runs **13×–182×** faster across
real and synthetic panels (Table 2), the factor widening with n as optiSel's
interior-point cost climbs. Given `G` precomputed, as the interior-point tools
require, the active set alone reaches the same optimum **~90×–2500×** faster: that is
the algorithmic advantage this section dissects, and the matrix-free solver trades
part of it for never building the matrix at all.

**Table 2.** support-first vs optiSel, same optimum throughout, measured serially on
one idle machine (Apple M4 Max): the matrix-free solver and optiSel in one R session,
the algorithm prototype in Python on the same instances. Each time is a median
(matrix-free and algorithm 5×, optiSel 3×). Two support-first columns: the active set
given a dense `G` (the regime optiSel runs in), and the shipped matrix-free solver,
which forms no `G`, timed end to end from genotypes including the R-binding copy.
optiSel constrains `cᵀ sKin c ≤ ub` on `sKin = G/2 + 10⁻⁵ I`, so support-first is
handed the identical problem as `k = 2·ub` on `G + 2×10⁻⁵ I`.

| dataset | n | m | algorithm¹ | matrix-free² | optiSel | speed-up³ |
|---|---|---|---|---|---|---|
| synthetic (structured) | 1000 | 500 | 0.023 s | 0.117 s | 2.06 s | 18× |
| synthetic (structured) | 2000 | 500 | 0.102 s | 0.280 s | 13.16 s | 47× |
| synthetic (structured) | 5000 | 500 | 0.506 s | 0.902 s | 163.8 s | 182× |
| CIMMYT wheat (real GRM) | 599 | 1279 | 0.003 s | 0.008 s | 0.595 s | 74× |
| PIC pig (real GRM, 52k SNP) | 3534 | 52843 | 0.020 s | 0.622 s | 49.9 s | 80× |
| HS mouse (real GRM, real sex) | 1814 | 10346 | 0.008 s | 0.057 s | 6.93 s | 122× |
| PIC pig, sexed sub-panel (real EBV **and** real sex) | 1194 | 52843 | 0.004 s | 0.227 s | 2.85 s | 13× |

¹ active set given a precomputed dense `G`, solve only — the regime optiSel runs in
(NumPy prototype); ² matrix-free, forming no `G`, end to end including the R-binding
genotype copy (shipped Rust); ³ optiSel / matrix-free.

The speed-up column is the shipped solver, which forms nothing, against optiSel:
13×–182×. The algorithm column — the active set handed `G` as the interior-point
tools are — reaches the same optimum ~90×–2500× faster than optiSel (the pig's
~2500× the largest, and the origin of this work's earlier figure), but that column
pays the `O(n²m)` to build and `O(n²)` to store `G` that the matrix-free solver never
does; handed `G`, the solve avoids streaming `Z`, so it is the faster of the two on
every panel, and the matrix-free column is what that costs to keep nothing dense in
memory. On the pig much of the shipped 0.622 s is the R binding copying the 1.5 GB
genotype matrix, so at m = 52 843 data movement, not the solve, sets the end-to-end
number. Both columns return the same optimum on every row, reaching the boundary
optiSel stops just inside.

The last row is the one instance carrying a real breeding value **and** a real sex at
once, and so the one closest to what a breeding programme actually runs. Sex is not
assumed there but read off the pedigree the PIC download ships — sire ⇒ male, dam ⇒
female — which resolves 1194 of the 3534 genotyped animals (390 sires, 804 dams, none
appearing as both), every one of them carrying a real EBV from the accompanying
evaluation at a mean accuracy of 0.811. On it support-first reaches gain 1.844614
against optiSel's 1.842975, on the constraint boundary where optiSel stops at 99.53%
of it, with the sex budget split exactly in half and 6 males and 21 females selected
out of 1194 candidates. It is also the least favourable row for the matrix-free
route: at m/n ≈ 44 the marker matrix dominates, which is why its speed-up (13×) is
the smallest in the table while the same solve given `G` takes 4 ms. That is the
trade-off the Methods state, visible in the one row where it bites hardest — and the
mouse row (934 males, 880 females) remains the other genuinely sexed instance.
Caveats to the Discussion: on wheat and mouse the selection criterion is a recorded
phenotype rather than a breeding value, and in no panel is it a *genomic* breeding
value; and the sexed pig sub-panel, being the animals that became parents, is a
post-selection subset.

## Comparison with the heuristic AlphaMate

AlphaMate optimises a related but distinct problem — discrete mate allocation by a
stochastic evolutionary algorithm — so this is not a like-for-like contest, and we
read it accordingly. We score *its* contribution vector in our metric and compare,
at matched group coancestry, only on the continuous-contribution relaxation the two
methods share. On that relaxation support-first's exact optimum has higher gain than
AlphaMate's vector at every point of its frontier (Table 3): a small margin at the
angle-45° trade-off (Δgain +0.004) and larger ones at the corners — as an exact
convex optimum should against a stochastic heuristic that additionally carries
integer mating constraints. We read this as a consistency check (support-first is
exact on the shared relaxation), not as AlphaMate being beaten at its own task,
discrete mate allocation, which support-first does not do and which is out of scope
here. AlphaMate was also markedly fragile on
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

## Scaling and the matrix-free advantage (Figure 1)

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

## Support behaviour

The advantage rests on the support, whose size is a property of the operating point
rather than a universal constant — so we characterise it along the whole frontier,
not at a single cap. On the mouse panel the support grows smoothly and monotonically
as the coancestry cap is tightened: 19 of 1814 candidates at the working coancestry
(0.0346), 61 at 0.010, 133 at 0.003, 189 at 0.0015, 473 at 0.0002, and about 1163 as
the cap approaches zero — where minimising relatedness forces the solution to spread
over much of the population — and it collapses to two at a loose cap. The claim is
therefore two-fold, and we keep the axes separate: at a *fixed* operating cap the
support is bounded as n grows (14–19 up to n = 40000, Figure 1A), while *along* the
frontier it grows as diversity is pushed. The per-solve cost tracks the support,
which is small in the regime breeders actually use, and which we report across the
frontier rather than at one point. Whether this is provable — n-independent at a
fixed cap, growing as the cap tightens — is taken up in the Discussion.
