# Toward a bound on the OCS support size

> Working note — turning the empirical observation (|S| small, bounded in n) into
> theory, with the subtleties stated honestly and **checked numerically**
> (`research/bound_validation.py`). Not a finished paper. The headline correction:
> the clean rank bound is an **ε = 0** statement; at the operative ridge the realised
> support is governed by the **effective (spectral) rank**, not by rank(G₀). Seed for
> the theory follow-on and the INRAE / optimisation (Alliot) collaboration.

## The object

Sexed OCS is the second-order cone program

  maximise **bᵀc**  s.t.  **A c = d**  (q rows: q = 1 simplex, q = 2 sexed),
  **c ≥ 0**,  **cᵀGc ≤ k**,   with **G = G₀ + εI**, **G₀ = ZZᵀ/s**, **Z** ∈ ℝ^{n×m}.

Let S = supp(c\*) at the optimum. We observe (Results; Waldmann 2025): |S| is small
(≈15 at a working cap), bounded as n grows (14–19 up to n = 40000 at a fixed cap),
and grows monotonically as k is tightened toward zero coancestry.

## Loose cap — clean

If cᵀGc\* < k the active set near c\* is just {A c = d, c ≥ 0}: a polytope, c\* an LP
vertex with q equality rows, so **|S| ≤ q**. The unconstrained sexed optimum puts all
mass on the best male and the best female (|S| = 2). The loose end of the frontier.

## Binding cap — the trap

When cᵀGc\* = k, c\* lies on the *curved* ellipsoid boundary. A naive count (|S| free,
minus q equalities, minus 1 active quadratic ⇒ |S| = q + 1) is **wrong**: that is the
count for *linear* active constraints. The quadratic surface is strictly curved, so
**every** point on it (∩ a coordinate face) is already extreme — curvature, not a
0-dimensional vertex, makes it extreme. So "extreme point of F" alone gives **no**
support bound. Any honest statement must dodge this.

## The ε = 0 theorem (this part is proved)

Take ε = 0, so G = G₀ = ZZᵀ/s with rank r ≤ m.

**Claim.** *There exists an optimal c\* with* **|S| ≤ q + r + 1**, *independent of n.*

**Proof.** By strong duality (Slater) there is λ\* ≥ 0 with c\* maximising the
Lagrangian L(c) = bᵀc − λ\*·cᵀG₀c over the polytope {A c = d, c ≥ 0}. Since
cᵀG₀c = ‖Zᵀc‖²/s, **L depends on c only through (Zᵀc, bᵀc) ∈ ℝ^{r+1}.** Define

  P = { c ≥ 0 : A c = d,  Zᵀc = Zᵀc\*,  bᵀc = bᵀc\* }.

Every c ∈ P is feasible (cᵀG₀c = ‖Zᵀc\*‖²/s ≤ k) and attains the optimum (bᵀc = bᵀc\*),
so **all of P is optimal.** P is a bounded, nonempty polytope cut by q + r + 1 equality
rows (A, then Zᵀ contributing r independent rows, then bᵀ), so it has a vertex with
≤ q + r + 1 nonzeros. That vertex is an optimum with the stated support. ∎

This is the Carathéodory / Barvinok–Pataki idea done **correctly**: don't count on the
curved boundary — **fix the low-rank image Zᵀc and the objective**, which linearises the
optimal slice, then count LP vertices. And it is a statement about *exactly* what an
active-set solver returns (a vertex), unlike interior-point / ADMM solvers, which return
non-sparse interior points and threshold post hoc (Waldmann 2025 truncates at 1e-4).

## The ε > 0 reality — the ridge rewrites the support (numerically established)

The theorem is for ε = 0. The solver — and any PD-requiring method — uses ε > 0, which
makes **G full rank n**, so rank(G₀) no longer bounds the *solved* problem. What governs
the realised support is then **spectral**, and the experiments (`bound_validation.py`)
are unambiguous:

- **Exact low rank + tiny ridge ⇒ large, n-growing support.** With G₀ of exact rank
  r = 10 and ε = 1e-5, |S| = 124, 265, 398, 616, 629 for n = 500 … 10000 — it climbs
  with n, far above q + r + 1 = 12. The cheap, near-degenerate ε-floor subspace gives a
  flat optimal face the solution spreads across.
- **The ε = 0 bound is the ε → 0 limit, reached only for ε ≪ the spectral floor.**
  Fixed r = 10, n = 1500, sweeping the ridge: |S| = 435 (ε=1e-2), 422 (1e-4), 66 (1e-6),
  **11 (1e-8)** — collapsing to ≤ q+r+1 = 12 only deep below the default ε. Gain is
  identical across these (~1.574): same flat optimal face, ε selects which representative.
- **Decaying full spectrum with a floor ≫ ε ⇒ small support — but no clean law.**
  Full-rank G₀ (rank 799), unit-strength factors: participation ratio
  3.5 / 6.4 / 12.3 / 24.5 → |S| = 5 / 8 / 12 / 20 (|S| of the order of the effective
  rank). **But a decay sweep breaks the clean version:** across spectra with
  participation ratio 1.2–55 at a fixed binding cap, |S| stayed **6–10** (not 1–55) — so
  the effective rank is *suggestive, not a reliable predictor*. What is robust is only
  that |S| is **small** and is **not** rank(G₀) or n; no single spectral scalar tracks it.

This is the operative regime for **real** GRMs: a smoothly decaying spectrum with no cheap
degenerate floor to spread into. There the support is robustly **small and flat in n** (the
paper's |S| ≈ 15), of the order of the effective rank / number of effective lineages — but,
per the decay sweep, that is a heuristic, not a proven scalar law. The clean rank theorem
was an ε = 0 truth oversold as governing the practical problem; the ridge term ε‖c‖² is not
a harmless regulariser for the support.

**Real GRMs — the missing datapoint** (`bound_real.py`): on the exported wheat (n = 599)
and mouse (n = 1814) kinship matrices, the support at the working cap is small and of the
*order* of the participation ratio, but not exactly it. Wheat — PR 28.5, |S| = 25 (simplex)
/ 26 (optiSel, sexed): a near-match. Mouse — PR 92.3, |S| = 17 (simplex) / 26 (sexed): same
order, PR over-predicting ~3–5×. Stable rank (1.5, 3.8, dominated by λ₁ ≈ 45–50) is far too
small to be the proxy. The cap moves |S| strongly at fixed spectrum (wheat 67 → 25 → 13 → 6
as it loosens). So on real data too the effective rank sets the *scale* (tens) while the cap
and b's alignment set the rest — no single scalar pins |S|, consistent with the synthetic
sweep. (A cross-check falls out: simplex |S| = 25 vs optiSel's sexed 26 on the real wheat K
— the support is robust across the formulation.)

**Isolating b's alignment** (`bound_balign.py`): with G *and* the cap held fixed, sweeping
only how the objective b aligns with G's eigenbasis settles which factor moves |S| — and the
relation is **inverse**. Putting b on the single dominant eigenvector (the most
coancestry-expensive direction, λ₁) forces a large support (mean |S| ≈ 146 of n = 800);
spreading b across many eigendirections, including cheap low-λ ones, collapses it to ≈ 4–5.
Mechanism: gain sought in an expensive direction must be diluted over many candidates to
respect the cap, whereas an objective weighted toward cheap directions concentrates. So the
support is driven by **b's projection onto the spectrum — alignment with the top (expensive)
eigenspace inflates it** — not by G's effective rank alone. This is exactly what the
spectrum-only proxies (rank, participation ratio) miss, and it explains how two panels of
very different effective rank (wheat 28, mouse 92) can share |S| ≈ 26: what matters is where
each b sits in its own spectrum. A predictive bound must therefore couple b and Λ — the
cost-adjusted object G^{-1/2}b is the natural candidate.

## The growth half (cap → 0)

As k → 0 the optimum approaches minimum coancestry and must spread mass to drive cᵀGc
down. The genetics quantifies it: ΔF ∝ E[Σcᵢ²] (Wray & Thompson 1990; Woolliams & Bijma
2000), Ne ≈ 1/(2ΔF). Tightening k forces Σcᵢ² down, i.e. spreads contributions, so |S|
grows monotonically, as observed (mouse: 19 → 61 → 133 → 189 → 473 → ~1163 as k → 0).

## Status and the open prize

- **Proved:** the ε = 0 bound |S| ≤ q + r + 1, n-independent (linearise-on-low-rank-image
  + LP vertex). Confirmed numerically as the ε → 0 limit.
- **Open (the prize):** a predictive bound on the **ridged** support. Effective rank /
  spectral gap is the natural suspect — small support coincides with a decaying spectrum and
  no cheap degenerate floor — but the decay sweep above shows that no single spectral scalar
  (participation ratio included) tracks |S| across regimes, so the right quantity is not yet
  identified. It plausibly couples the spectral gap, b's alignment with the top eigenspace,
  and the cap k. A spectral-perturbation question — the natural meeting point of optimisation
  (Alliot: face dimension of the perturbed cone) and genetics (Bouchet: effective number of
  contributing lineages).
- **Marrying the halves:** |S| as a function of the target ΔF / Ne, capped by the (still to
  be identified) spectral quantity above — the publishable characterisation.

## References (verified)

- Pataki, G. (1998). On the rank of extreme matrices in semidefinite programs…
  *Math. Oper. Res.* 23(2):339–358. DOI 10.1287/moor.23.2.339. — extreme-point rank counting.
- Markowitz, H. (1956). *Naval Res. Logist. Q.* 3:111–133. DOI 10.1002/nav.3800030110. — critical-line / corner-portfolio structure (the "grows along the frontier" analogue).
- Wray, N.R. & Thompson, R. (1990). *Genet. Res.* 55(1):41–54. DOI 10.1017/S0016672300025180. — ΔF ∝ Σ(contribution²).
- Woolliams, J.A. & Bijma, P. (2000). *Genetics* 154(4):1851–1864. DOI 10.1093/genetics/154.4.1851. — contributions ↔ ΔF ↔ Ne.
- Yamashita, M., Mullin, T.J. & Safarina, S. (2018). *Optim. Lett.* 12(7):1683–1697. DOI 10.1007/s11590-018-1229-y. — OCS as a single-SOC program (the cone setup).
- Waldmann, P. (2025). *Bioinform. Adv.* 5(1):vbaf259. DOI 10.1093/bioadv/vbaf259. — independent empirical support counts; post-hoc truncation.

*Numerical evidence: `research/bound_validation.py` (blocks 1–6: solver vs scipy, the
ε = 0 bound sweep, n-independence, the effective-rank behaviour, the ridge sweep, and a
spectral-decay sweep), `research/bound_real.py` (wheat and mouse GRMs, after the
`research/repro/*_export.R` exports), and `research/bound_balign.py` (the b-alignment sweep).*
