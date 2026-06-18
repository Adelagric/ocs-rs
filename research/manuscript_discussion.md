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
