# Discussion (manuscript draft)

> Companion to the intro/methods/results drafts. Limitations are stated plainly
> and up front — under-promising is the credible posture for a methods paper.

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
boundedness of the optimal support is, in this paper, an empirical observation
across a synthetic sweep and the real panels, not yet a theorem, and the
route is subtler than a constraint count. With no ridge (ε = 0, G = ZZᵀ/s of rank
r) a clean argument bounds it: the optimum maximises a Lagrangian depending on c
only through (Zᵀc, bᵀc), so the optimal slice {Ac = d, c ≥ 0, Zᵀc = Zᵀc\*,
bᵀc = bᵀc\*} is an LP polytope whose vertices carry ≤ q + r + 1 nonzeros — an
extreme-point / Carathéodory bound, independent of n, on exactly the vertex an
active-set solver such as support-first returns (interior-point and ADMM solvers
instead return non-sparse interior points and threshold post hoc). The operative
ridge ε > 0 makes G full rank, however, and we find numerically that the realised
support is then not governed by rank(G₀): it stays small (a few dozen at most) and
flat in n on the panels here, but does not reduce to a single clean spectral
quantity — the effective rank is suggestive yet, across a sweep of spectra, not a
reliable predictor — so a tight bound for the ridged problem remains open. The genetics accounts for the growth half: ΔF ∝ Σcᵢ²
(Wray & Thompson 1990; Woolliams & Bijma 2000) with Ne ≈ 1/(2ΔF), so tightening the
cap spreads contributions and grows the support, as observed. We develop this
characterisation in follow-on work. Finally, the solver handles a single
quadratic kinship constraint, per-candidate contribution caps (0 ≤ c ≤ u), and the
sex equalities, but not yet several quadratic constraints at once — multiple
relationship matrices, group-specific coancestry limits — nor the integer mate
allocation that AlphaMate provides; it returns continuous contributions.

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
