# Revision roadmap — response to the severe review

Truth-first. Claims are scoped strictly to what is demonstrated; headline numbers
are re-established at breeder-defensible operating points **before** any manuscript
edit. Whatever the re-benchmarks show is what gets reported, even if the headline
speed-ups shrink. This file tracks each review point → verdict → fix → status.

## Load-bearing defects (confirmed against code/data)

| # | Point | Verdict | Fix | Status |
|---|-------|---------|-----|--------|
| 1 | Operating point never translated to ΔF/Ne; headline timings sit at the loosest end of our own frontier | **Confirmed.** "Ne" appears 6× in the whole manuscript. Pig working point (½/½ over ~6+21) ⇒ Σc²≈0.054 ⇒ **Ne≈19**, below the FAO floor (50) and commercial ΔF targets (0.5–1%). §3.6 already shows the operating regime is \|S\|=61–133 (coancestry 0.003–0.010), not \|S\|=19 | Report Ne_c=1/Σc² and ΔF≈1/(2Ne_c) at every operating point; re-benchmark at ΔF ∈ {0.5,1,2%} (Ne ∈ {100,50,25}); report time + \|S\| + speed-up **along** the frontier, secondary Ne axis | in progress |
| 3 | "Support bounded as n grows" is an artifact of a fixed **absolute** cap while the feasible floor collapses | **Confirmed in code**: `examples/scaling_matrixfree.rs` uses `k = 0.1*mean_diag` ≈ const, while `k_min ≈ 1/n` falls from 1e-3 to 2.5e-5 across n=1000→40000. Deeper cause: the synthetic generator is **HWE with independent loci** (no structure/LD/founders), so growing n at fixed m just dilutes relatedness. Fixing Ne bounds \|S\| by construction, so the honest test needs a **structured** population (fixed base Ne) via AlphaSimR | (a) annotate current Fig 1A points with Ne_c/ΔF/realised f; (b) structured-population scaling at a fixed operating point; report the result whatever it is | in progress |
| 4 | Theorem 1 (\|S\| ≤ q+r+1) is vacuous on the panels: m≫n ⇒ r=n−1 ⇒ bound = n | Correct | State it explicitly; frame the appendix as "no useful bound in the genomic regime, here is why both obvious directions fail" | todo |
| 5b | "None exploits solution-support sparsity" is false — Meuwissen (1997) prunes negative contributions (a decreasing dense active set) | Correct (intro l.46, l.77) | Reframe: increasing active set by column generation that never forms the n×n system, vs decreasing pruning on the full dense system | todo |
| 6b | Table 2 puts NumPy and Rust in adjacent, cross-readable columns — the language confound we claim to avoid in Table 1 | Correct | Port the prototype to Rust, or rename the column to forbid the cross-read | todo |
| 6c | Table 2 gives no \|S\| per row — the variable that explains the times | Correct | Add \|S\|, effective cap, ΔF per row | todo |
| 7a | Abstract "a solve costs the support and the marker count rather than the n×n matrix" elides the n | Correct: KKT pricing is a full G·c over n candidates = **O(nm) per iteration** (the code labels `products` "the dominant O(n·m) cost") | Write O(iters·nm + \|S\|³): linear in n, never quadratic, never formed | todo |
| 7b | No numerical-robustness story: Cholesky on G_S with near-duplicate genotypes | Partially correct: the ridge εI floors the spectrum (κ(G_S) ≤ (λmax+ε)/ε), so precision degrades, it does not "explode". But no conditioning is reported and finite-termination is asserted, not proved | Spike: panel with exact-duplicate genotypes, measure κ(G_S) and solution accuracy; soften "guarantees finite termination" to observed-no-cycling unless proved | **done** (Exp 5): with 5 exact clones the solve is `Solved`+feasible; G_S min eigenvalue is −3.4e-16 without the ridge (singular) but floored at ε=1e-5 with it, so **κ(G_S)=2.7e5, bounded** — no explosion; clones get exactly equal contributions. `research/repro/robustness_duplicates.R`. Soften "guarantees finite termination" → "terminated at a feasible KKT optimum on every instance tested, incl. exact-duplicate degeneracy" |

## Findings so far

**Experiment 1 — `examples/scaling_ne.rs`, synthetic HWE panels, m=1000** (artifact `artifacts/scaling_ne.csv`):

- **Figure 1A sits at Ne ≈ 10, ΔF ≈ 5 %.** At the current cap (k = 0.1·mean_diag) the realised
  operating point is Ne_c 9–10 and ΔF ≈ 5 % per generation across n = 1000→40000 — *worse* than
  the reviewer's estimate (Ne ≈ 15), and far below the FAO floor (50) and commercial targets
  (0.5–1 %). The "support ~15" headline is a ΔF = 5 % regime. Point 1 confirmed and quantified.

- **At a defensible operating point the core survives.** Re-measured at fixed Ne (n = 1000→40000):

  | target | ΔF | \|S\| range | solve s @ n=40000 |
  |---|---|---|---|
  | Ne = 25 | 2 % | 41–49 | 0.12 |
  | Ne = 50 | 1 % | 86–99 | 0.31 |
  | Ne = 100 | 0.5 % | 171–203 | 0.86 |

  The support is ~90 at Ne = 50, **not 15** — roughly 6× the headline. But the solve stays
  sub-second at n = 40000 with **no G formed** (the dense G is 11.9 GiB, unbuildable), so the
  matrix-free / memory argument holds at every operating point, and the speed almost certainly
  survives too. The honest headline moves to defensible ground and the core gets *stronger*, not
  weaker.

- **\|S\| stays bounded in n even at fixed Ne — but on THIS generator that is largely because
  G ≈ I.** With independent HWE loci, cᵀGc ≈ Σcᵢ² = 1/Ne_c, so at fixed coancestry the parent
  count is set by the cap, not by n. Realised f (0.007 / 0.016 / 0.037) tracks Σc² (0.010 / 0.020 /
  0.040), confirming the near-identity regime. The genuine "bounded as n grows" test therefore
  needs a STRUCTURED population where G carries real off-diagonal mass — next step (AlphaSimR).

**Experiment 2 — `research/repro/r_operational_ne.R`, real panels, optiSel head-to-head at
operational Ne** (artifact `artifacts/operational_ne_realpanels.csv`; R BLAS = OpenBLAS 0.3.33,
so the baseline is not a slow reference BLAS — the speed-up is conservative):

Mouse (HS, real sex, n=1814) — the clean comparison, gain gaps ≤ 2.6e-3 (same optimum):

| Ne | ΔF | \|S\| | optiSel | ocsrs | speed-up |
|---|---|---|---|---|---|
| 25 | 2 % | 47 | 9.4 s | 0.28 s | **34×** |
| 50 | 1 % | 89 | 8.1 s | 1.46 s | **6×** |
| 100 | 0.5 % | 182 | 7.3 s | 8.4 s | **~1× (parity)** |

Wheat (n=599): 9× / 3× / 0.6×, but at tight caps optiSel converges ("cccp2 optimal", no warning)
to a point using only **42 % of the budget** with lower gain (0.87 vs 1.12). The two formulations
diverge there, so wheat is **not** a clean same-optimum comparison at operational caps — mouse is.
(Open item: why optiSel leaves budget on wheat — likely cccp tolerance on an ill-conditioned
autogamous panel with an artificial sex split.)

**Conclusion (load-bearing).** The "12–180× faster on real panels" headline was a **loose-cap
artifact**. At operational ΔF on modest real panels support-first is a few × faster (Ne=50) and at
**parity or slower** at the tightest defensible point (Ne=100). The durable contribution is **not
speed on small panels**; it is (a) **exactness** (KKT-optimal; optiSel halts inside the constraint
on wheat) and (b) **memory at scale** — at n=40000 the dense G is 11.9 GiB (unbuildable), which
optiSel cannot run **at all**, while matrix-free solves in <1 s. The thesis must move from "orders
of magnitude faster" to "exact, and enables population-scale OCS the dense-G tools cannot open."
The real-panel speed claim is now: 6× at Ne=50, parity at Ne=100 (mouse), reported honestly.

**Experiment 3 — `research/repro/scaling_structured.R`, AlphaSimR structured panels (fixed
100-parent pool, coalescent LD, m=1000), operating point Ne=50** (artifact
`artifacts/scaling_structured.csv`):

| n | \|S\| | solve s | dense G | optiSel @ Ne=50 |
|---|---|---|---|---|
| 1000 | 86 | 0.053 | 0.01 GiB | 2.1 s → **40×** |
| 2000 | 103 | 0.062 | 0.03 GiB | 11.5 s → **186×** |
| 5000 | 84 | 0.103 | 0.19 GiB | — |
| 10000 | 81 | 0.105 | 0.75 GiB | — |
| 20000 | 94 | 0.207 | 2.98 GiB | — |
| 40000 | 88 | 0.646 | 11.92 GiB | — |

- **Enablement holds on realistic structured data**: matrix-free solves n=40000 in 0.65 s at
  ΔF=1 % while the dense G is 11.9 GiB (unbuildable). This is the durable contribution, now shown
  on AlphaSimR data, not just the HWE toy.
- **\|S\| stays bounded (~80–103) at fixed Ne as n grows** → reframe "bounded support" as **bounded
  at a fixed operating point**, honest and demonstrated — not "intrinsically bounded", which the
  fixed-absolute-cap Figure 1A conflated.
- **Metric caveat (self-corrected)**: the uniform-allocation coancestry floor is ≈0 by VanRaden
  centring (collapses to 3.8e-11 at n=40000) — a centring artifact, NOT a structure indicator. Do
  not report it. The true intrinsic-boundedness test needs the actual min-coancestry (min c'Gc), a
  matrix-free QP — deeper open item.
- **Speed-up is strongly m-dependent** (matrix-free cost O(iters·nm)): at Ne=50, m=1000 gives
  40–186×, whereas mouse (m=10346) gave 6×. The operational-ΔF speed-up is **not uniformly dead** —
  it shrinks with marker density. **Pig (m=52843) is the acid test** and is the next run.
- optiSel again leaves gain at the operational cap (gaingap ~1e-2 at n=1000–2000, like wheat, unlike
  mouse) — the inexactness recurs, strengthening the exactness claim (pending a clean "why").

**Experiment 4 — pig sexed sub-panel (real EBV + real sex), n=1194, m=52843, the high-marker acid
test** (`research/repro/r_operational_ne_pig.R`; probe numbers, single timed rep):

| Ne | ΔF | \|S\| | optiSel | ocsrs matrix-free | speed-up |
|---|---|---|---|---|---|
| 25 | 2 % | 56 | 2.89 s | 1.42 s | 2.0× |
| 50 | 1 % | 106 | 2.77 s | **10.06 s** | **0.3× (slower)** |
| 100 | 0.5 % | 199 | 2.77 s | **51.3 s** | **0.1× (18× slower)** |

At the manuscript's loose pig cap the operating point is again **Ne=10.3, ΔF=4.87 %** (point 1
reconfirmed on the flagship panel). At operational ΔF the **matrix-free path is SLOWER than
optiSel**, badly so at tight caps: optiSel forms G once (n=1194 is small, so the wall does not bite)
and its IPM is a fixed ~2.8 s reusing that G, while matrix-free re-derives kinship from Z at
O(iters·nm) with m=52843 and iters growing with \|S\| (56→106→199).

**The unified speed law (load-bearing).** At a fixed operating point the matrix-free speed-up is
governed by the marker count m (its per-iteration cost is O(nm)):

| panel | m | speed-up @ Ne=50 |
|---|---|---|
| AlphaSimR synthetic | 1000 | 40–186× |
| wheat | 1279 | 3× (messy comparison) |
| mouse | 10346 | 6× |
| **pig (flagship)** | **52843** | **0.3× (slower)** |

So "orders of magnitude faster on real panels" is **false at operational ΔF on the flagship**. The
matrix-free path wins at low-to-moderate m and, decisively, at **large n where G cannot be formed at
all** (n=40000 → 11.9 GiB); it loses at small-n/high-m/tight-cap, where forming G once and reusing
it (given-G) is better. This is a real, clean characterization — *when* matrix-free wins — and it is
a stronger scientific contribution than a raw speed claim.

**What survives, firmly:** (1) memory-at-scale enablement (n-driven, m-independent) — the durable
result; (2) exactness (ocsrs at the boundary; optiSel leaves gain). **Not yet measured:** the
"algorithm given G" active-set column at operational Ne on the pig — given G, both avoid the O(nm)
recompute, so the active-set may still beat optiSel's IPM even at high m. And the matrix-free path
may be optimisable (cache/update Zᵀc across pricing instead of recomputing) — an engineering
opportunity, not claimed.

**Experiment 6 — given-G active-set at operational Ne on the pig** (`research/repro/sf_operational_ne_pig.py`,
K precomputed; the algorithmic column, m-independent because pricing is K·c = O(n²)):

| Ne | ΔF | \|S\| | iters | algo time | vs optiSel IPM (2.77 s) |
|---|---|---|---|---|---|
| 25 | 2 % | 57 | 79 | 0.016 s | **173×** |
| 50 | 1 % | 110 | 124 | 0.078 s | **36×** |
| 100 | 0.5 % | 194 | 76 | 0.238 s | **12×** |

**Given G, the active-set is 12–173× faster than the IPM at operational ΔF even at m=52843** — while
the matrix-free path at the same points was *slower* than optiSel (10 s at Ne=50). The ~128× gap is
entirely the matrix-free O(nm) reconstruction penalty, **not** the algorithm. (Prototype in NumPy, so
the absolute times are an upper bound; a Rust given-G path would be faster — but the algorithmic win is
unambiguous from the iteration counts × O(n²) vs the IPM's O(n³).)

**The honest two-part thesis (now established):**
1. **Algorithm** — the support-first active set exploits solution-support sparsity and beats a general
   IPM at operational ΔF (12–173× on the pig, given G). n-bound, m-independent. This is the speed result,
   correctly scoped to the given-G comparison.
2. **Matrix-free representation** — never forming G *enables* population scale (n=40000, dense G=11.9 GiB,
   unbuildable) at O(nm) per iteration; the right representation at large n / moderate m, and deliberately
   *not* the fast choice at small n / very large m, where forming G once wins.

Table 2 must separate the two axes cleanly (algorithm vs representation) and stop inviting a cross-read of
NumPy-given-G against Rust-matrix-free absolute times (review point 6b).

## Strategic (venue-dependent)

| # | Point | Note |
|---|-------|------|
| 5a | Gencont2 not benchmarked | It **is** cited (intro l.35); the benchmark is missing. Nuance: Gencont2 exploits pedigree **A⁻¹** sparsity, not the dense genomic **G** — not a clean genomic head-to-head. State this |
| 5c | Waldmann 2025 dismissed in one line | Deserves a benchmark or a paragraph on why not comparable. Verify DOI vbaf259 before freezing |
| 8 | No multi-generation validation | AlphaSimR over 10–20 generations (support-first vs optiSel vs AlphaMate vs truncation) is the domain standard. Required to retain any breeder-regime claim |

## Writing (mechanical, correct regardless)

- Abstract 450→250 words, structured, not one paragraph with three-level nesting.
- Drop the honesty tic ("honest caveats", "we state them plainly", …, 6+ occurrences) — keep the caveats, cut the announcements.
- The 12–180× / ~2500× pair appears 5×; once suffices.
- Unify notation: Table 3 uses c'Kc, the rest G.
- Add two figures: gain/coancestry frontier with the three methods; \|S\| + time along that frontier (§3.6 is prose today).
- Tone down essayistic lines ("doing what is otherwise infeasible…", "the speed is the lesser result") for a defensive readership.
- Declare R's BLAS + full sessionInfo. Cut the AlphaMate Rosetta timing; keep one factual sentence on fragility. Single-run-at-large-n → table note.

## What survives everything

The dense-G memory wall is **ΔF-independent**: at any operating point the dense G is
11.9 GiB at n=40000 and O(n²m) to build; an exact solver that never forms it opens a
problem the existing tools cannot load. That is the durable contribution.

## Sequencing

1. Quantify Ne/ΔF at current operating points (Rust, matrix-free) — makes point 1 concrete. ← **now**
2. Operational re-benchmark at ΔF ∈ {0.5,1,2%} on the real panels.
3. Structured-population scaling (AlphaSimR panels) for the honest point-3 test.
4. Robustness spike (duplicate genotypes).
5. Framing corrections (thm 1, Meuwissen, O(nm), Table 2 columns/|S|).
6. Writing pass.
7. Baselines (Gencont2/Waldmann) + AlphaSimR multi-generation validation — per venue.

Manuscript edits happen **after** 1–4, never before the numbers exist.
