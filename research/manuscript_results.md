# Results (manuscript draft)

> Companion to `manuscript_intro.md` / `manuscript_methods.md`. Numbers locked
> from the benchmark tables (PREPRINT.md ledger). Figure 1 = `research/fig_scaling.pdf`.
> Timings: support-first is a NumPy prototype unless noted "Rust"; optiSel is
> R/`cccp`; the gap is algorithmic, not a language artefact (see Methods).

## Exactness

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

## Speed against a generic conic solver

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

## Speed against the domain tool optiSel

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

## Comparison with the heuristic AlphaMate

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

The advantage rests on the small, bounded support. At the mouse operating
coancestry (0.0346) the optimum places weight on 19 of the 1814 candidates; the
support enlarges only as the cap is driven toward zero — about 1163 individuals to
force group coancestry near 0, where the solution must spread over much of the
population to minimise relatedness. Across the synthetic sweep the support stays in
the low tens as n grows forty-fold (Figure 1A), which is what makes the per-solve
cost scale with the support rather than with n.
