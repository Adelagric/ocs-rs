# Discussion (manuscript draft)

> Companion to the intro/methods/results drafts. Limitations are stated plainly
> and up front — under-promising is the credible posture for a methods paper.

Support-first makes exact optimum contribution selection cheap at genomic scale.
On the panels tested it reaches the same optimum as the domain's exact tool while
running one to two orders of magnitude faster end to end, forming no relationship
matrix at all — and, handed that matrix as the interior-point tools require it, its
active set alone reaches the same optimum up to three orders of magnitude faster. It
stays cheap precisely where the dense relationship matrix that every other solver
forms becomes infeasible — the large-candidate regime that motivates genomic OCS in
the first place. Because the method is exact and deterministic — a Karush–Kuhn–Tucker–
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
selection criterion is a real estimated breeding value on the pig panels and a
recorded phenotype standing in for one on wheat and mouse. On the sexed pig sub-panel
we additionally ran a genuine GEBV, fitted by GBLUP on that panel's own phenotypes
and relationship matrix, and the result was unchanged, so the conclusions do not turn
on which breeding value serves as the criterion; on wheat and mouse, however, the
criterion remains a phenotype, and the contribution vectors are in every case
illustrative rather than breeding recommendations.
Second, a genuine recorded sex is available for the mouse panel and for the sexed pig
sub-panel — the single instance combining a real breeding value with a real sex,
which is therefore the closest thing here to an operational OCS problem — while wheat
and the full pig panel take an arbitrary balanced split. That sub-panel carries its
own restriction: sex is recoverable from the pedigree precisely for the animals that
became parents, so it is a post-selection subset of the candidates rather than a
random sample of them. Third, Table 2 separates two things the speed claim used to conflate. The shipped
matrix-free solver, timed end to end against optiSel including the genotype copy
across the R binding, is 12×–180× faster; the active set given a dense `G` — the
regime optiSel runs in, but forming and storing the `O(n²)` matrix — reaches the same
optimum ~90×–2500× faster, and it is this algorithmic figure the pig's ~2500×
represents. The two coincide in mechanism (the same tiny active set) and differ only
in whether `G` is materialised, which is precisely the cost the matrix-free route
avoids. That route is *not* an inner-loop speed-up when markers outnumber candidates
(m > n): on the pig its end-to-end cost is dominated by streaming and copying the
52k-marker matrix, and its 80× understates the algorithm — a true zero-copy binding,
of the kind the Python path uses, would lift it. The matrix-free route is the memory
and large-n enabler; the speed advantage over the conic solvers is the small active
set. Fourth, the
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

One boundary of the study should be stated plainly, because it decides which claims
are tested and which are not. This work was carried out without institutional data
agreements, so every panel used is one any reader can download and every number in
the paper is regenerable from public sources — a deliberate gain in reproducibility,
and a real ceiling on scale. The largest real populations in animal breeding sit in
national evaluation databases holding hundreds of thousands of genotyped animals, and
those are access-controlled; the largest openly available cattle sequence resource
holds a few thousand individuals, fewer than the pig panel used here. The scaling
results at n = 10⁴–4×10⁴ are therefore synthetic by necessity rather than by
preference. What they establish — that the dense relationship matrix becomes
infeasible before the solver does — is a property of that matrix and does not depend
on the genotypes being real. What they do not establish is the method running inside
a production evaluation at that scale, and no public dataset can establish it. That
comparison needs a breeding organisation or a national evaluation centre, and we
would welcome it.

Each limitation points to an extension. The closed-form-per-support core already
admits more than one equality constraint through the same elimination
P = A G⁻¹ Aᵀ; additional active quadratic caps turn the scalar root-finding into a
small low-degree system rather than changing the structure, so multiple kinship
constraints are within reach. A rounding or branch-and-bound layer over the
continuous optimum would add mate allocation while keeping the exact relaxation as
a tight bound — combining support-first's exactness with the discrete plan that
heuristics target directly. The matrix-free crossover in m versus n deserves an
explicit characterisation, as does the support bound itself. The criterion question is
now closed on one instance — a GBLUP-fitted GEBV on the sexed pig sub-panel changes
nothing — so what remains is validation on populations an order of magnitude larger,
where the memory gap of Figure 1B turns from visible into decisive.

Support-first does one thing — exact, single-constraint, continuous optimum
contribution selection — at a scale and speed that bring genomic OCS within reach
of a laptop and a reproducible script. That narrow, verified claim, and the open
extensions it invites, are the contribution.
