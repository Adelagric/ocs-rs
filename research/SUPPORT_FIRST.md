# Support-first OCS — research note

Status: **validated and written up — superseded by**
[`MANUSCRIPT.md`](MANUSCRIPT.md), which carries the final numbers, the `c ≤ u`
extension, and the Rust implementation. This note is kept as the development
log: the lead, the reduction, and how it was confirmed (synthetic → real wheat,
pig, and mouse panels). Numbers below are reproducible with
`python3 research/support_first.py`.

## Why

The Clarabel spike (see top-level `REPORT.md`) returns a clean **GO**: the conic
interior-point method solves genomic OCS reliably and exactly. But it is slow at
scale — n=10000 takes ~30 min and 8.8 GB, because each of its ~17 iterations
factors a dense `n×n` KKT system (O(n³)). Meanwhile the *answer* is tiny: at a
moderate kinship bound the optimal contribution vector activates only **2–5**
candidates out of n. The IPM pays O(n³) to describe a solution that lives in a
~5-dimensional face. This note asks whether exploiting that sparsity breaks the
wall — without giving up exactness (unlike the evolutionary heuristic AlphaMate).

## The reduction

OCS is `max bᵀc s.t. 1ᵀc=1, c≥0, cᵀGc≤k`. Drop `c≥0` for a moment and fix a
support `S`. Eliminate the equality with `c = c₀ + N y` (`c₀=1/|S|·𝟙`, `N` an
orthonormal basis of `ker(𝟙ᵀ)`):

```
max b̃ᵀy   s.t.   yᵀG̃y + 2g̃ᵀy + (q₀−k) ≤ 0
```

with `G̃=NᵀGN`, `b̃=Nᵀb`, `g̃=NᵀGc₀`, `q₀=c₀ᵀGc₀`. This is **maximising a
linear form over an ellipsoid** — closed form: centre `y_c=−G̃⁻¹g̃`, radius
`ρ²=g̃ᵀG̃⁻¹g̃−(q₀−k)`, optimum `y* = y_c + √ρ²·G̃⁻¹b̃ / √(b̃ᵀG̃⁻¹b̃)`. Two
linear solves with `G̃`, no iteration.

Because the OCS objective is **linear**, this is *not* a hard generalized
trust-region subproblem; the eigenvalue machinery (Adachi–Nakatsukasa) is only
needed if the objective becomes quadratic (robust / variance-aware / multi-trait
OCS — a separate lead). For standard OCS the binding difficulty is entirely
`c≥0`: dropping it changes the optimum drastically (e.g. n=50, gain jumps
1.75 → 4.55 with a contribution of −0.24).

## The algorithm

`c≥0` is handled by an **active set / column generation** outer loop; each inner
solve is the closed form above. `G` is never formed — every `Gc` is
`ridge·c + Z(Zᵀc)/s` in `O(n·m)`.

```
S ← {argmax b}
loop:
    solve closed form on S
    if infeasible on S (ellipsoid ∩ affine empty):   # support too related
        add the candidate least related to current c   (argmin (Gc)_j)
    elif some c_i < 0:    drop those i from S
    else:
        reduced costs r_j = b_j − μ − 2λ (Gc)_j      # μ,λ from KKT on S
        if max_j r_j ≤ tol:   return c                # KKT optimal + feasible
        else: add argmax_j r_j to S
```

Cost ≈ `(#products) × O(n·m)`, with `#products ≈ |S|`, versus `~17 × O(n³)` for
the IPM.

## Measured (synthetic, M4 Max)

**Exactness** — matches reference solvers:
- vs SciPy SLSQP (iid, n≤1000): `Δgain ≈ 1e-12`
- vs Clarabel (crate dump, n=2000): `Δgain ≈ 1e-8`
- vs SciPy (pedigree + correlated EBV, n=400, |S|=7): `Δgain ≈ 1e-11`

**Aligned against optiSel** (the domain reference, R `cccp` IPM). Reproduced
optiSel's marker-based OCS on its `Cattle` example (Angler, n=268) with the *real*
reproduction constraints `Σmales = Σfemales = 0.5` — which required extending
support-first to two equality constraints (kernel of `[𝟙ₘ; 𝟙f]`). Result:
identical active support (36/36), near-identical contributions (max 0.0664 vs
0.0661), and support-first reaches the **exact** optimum (kinship saturated at
the bound, gain 1.8315) where optiSel's IPM stops ~0.4% short (gain 1.8246,
kinship 0.0576 < ub 0.0578). Same problem, same solution; support-first is at
least as good. Caveat: optiSel's `Cattle` data is itself *simulated* and small
(268 individuals, 800 markers) — a real dataset and a larger-scale *timing*
comparison are the open items.

An independent review verified the closed-form algebra (Vieta + Cauchy–Schwarz:
exactly one root has λ>0 and it is the gain-maximiser) and `Solved ⇒ feasible &
optimal` on 3400 brute-force instances (zero counter-examples). It also found one
real defect — an active-set cycle on *degenerate* inputs (random `G`, tiny `k`,
no non-negative optimum) that looped to `MaxIter` with an infeasible point. Fixed
with anti-cycling (taboo dropped indices until the next progress step) plus an
explicit "`c` valid only when `Solved`" contract; the VanRaden regime never
triggered it (400/400 `Solved`).

**Cost, matrix-free, iid** (`m=20000`, numpy prototype):

| n | \|S\| | G·c products | wall (python) | vs Clarabel |
|---|---|---|---|---|
| 1000 | 3 | 4 | 0.009 s | 1.66 s |
| 5000 | 3 | 4 | 0.037 s | 191 s |
| 10000 | 3 | 4 | **0.072 s** | **1815 s** |
| 20000 | 2 | 3 | 0.111 s | (IPM impractical) |

**Rust head-to-head, apples-to-apples** (`cargo run --release -- compare --n N`,
same data, same machine, both in Rust; gain agrees to ~1e-9 and the active
support is identical at every size):

| n | Clarabel solve | support-first solve | **speedup** |
|---|---|---|---|
| 1000 | 1.40 s | 0.011 s | 126× |
| 2000 | 10.6 s | 0.022 s | 472× |
| 5000 | 160 s | 0.036 s | 4474× |
| 10000 | 1579 s | 0.043 s | **37090×** |

The speedup grows steeply (Clarabel is O(n³) per iteration, support-first is
roughly constant in products: 3–6). **Caveat — this is vs Clarabel, a *generic*
conic solver, on *synthetic* data.** A specialised solver is expected to beat a
generic one; the comparison that matters for the field is vs the domain tools
(optiSel / AlphaMate / Gencont2) on *real* data, which is not done here.

**`|S|` stays bounded** — the key scaling result. Worst regime tested
(pedigree, EBV correlated to genotype, tight k), `|S|(n)`:

| k fraction | n=500 | n=2000 | n=5000 | n=10000 |
|---|---|---|---|---|
| 0.05 (very tight) | 32 | 40 | 44 | 31 |
| 0.15 (tight) | 11 | 12 | 12 | 11 |
| 0.30 (moderate) | 5 | 4 | 7 | 5 |

`|S|` is flat in `n` (it tracks the number of *effective lineages*, set by the
population structure and the diversity demanded, not by candidate count). The
product count is bounded likewise (~17 at k=0.15, ~50–60 at the tightest),
independent of n.

## Benchmark vs optiSel at scale — the test that matters

The Clarabel comparison pits us against a *generic* solver. The real test is vs
the **domain tool, optiSel** (R, `cccp` IPM), on its own formulation: the true
sex-constrained OCS `max bᵀc s.t. Σmales=Σfemales=0.5, c≥0, cᵀGc≤k`. This
required extending support-first to **two equality constraints** (eliminate
`[𝟙ₘ; 𝟙f]` by its kernel; reduced costs become sex-dependent), and exposed — then
fixed — a real efficiency bug: the feasibility phase added candidates **one at a
time** (675 of 768 iterations at n=1000). Adding them **in chunks** (double the
support until feasible) cuts iterations ~30× with the *identical* optimum.

After that fix, support-first solves the same problem as optiSel, reaching the
**exact** optimum on the kinship boundary where optiSel's IPM halts just inside it
(equal gain at matched realised coancestry; the edge is unspent diversity budget),
at:

| dataset | n | support-first | optiSel | speedup |
|---|---|---|---|---|
| synthetic (structured) | 1000 | 0.09 s | 2.0 s | 22× |
| synthetic (structured) | 2000 | 0.29 s | 11.9 s | 41× |
| synthetic (structured) | 5000 | 1.47 s | 143 s | 97× |
| **wheat** (real CIMMYT GRM) | 599 | 0.007 s | 0.63 s | 90× |
| **PIC pig** (real GRM, 52k SNP) | **3534** | **0.024 s** | **54.8 s** | **~2280×** |
| **HS mice** (real GRM, **REAL sex** 934♂/880♀) | **1814** | **0.008 s** | **6.96 s** | **~870×** |

The **mouse row is the *true* sexed OCS** — a genuine recorded sex (BGLR `mice`,
GENDER complete, 0 missing), not an arbitrary partition. Caveats kept explicit:
support-first is the numpy prototype, optiSel is R/`cccp` (the gap is algorithmic,
not language); the synthetic / wheat / pig rows use an **arbitrary 2-group sex
partition** (wheat is autogamous; PIC ships no usable sex — chromosomes removed,
only 390/3534 identifiable sires) so those benchmark the *solvers* on a real GRM
rather than the true sexed OCS; `b` is a real phenotype/EBV (mouse BMI, pig
trait-3 EBV); and the speedup grows with how *small* the active support is (a very
tight `k` enlarges it and shrinks the factor). The first revision (proto "naively
slower than optiSel") was entirely the one-at-a-time feasibility phase —
corrected, the advantage is real, exact, and grows with scale.

## Where it sits vs the state of the art

| method | exact? | regime | what it exploits | limit |
|---|---|---|---|---|
| Meuwissen 1997 | yes | dense | Lagrangian on λ | iterative, dense |
| Gencont2 (Pong-Wong) | yes | **pedigree sparse** | Newton-λ + Gauss-Seidel, matrix-free | not the dense genomic `G`; first-order on c |
| optiSel `cccp` / this crate's Clarabel | yes | genomic dense | SOCP + IPM | O(n³) |
| AlphaMate | **no** (heuristic) | any | differential evolution | no optimality |
| **support-first (here)** | yes | **genomic dense** | **solution sparsity + factored Z + closed form** | single quadratic constraint, continuous c |

The novelty is not "use a cone" (optiSel does), "matrix-free Lagrangian"
(Gencont2 does on pedigree), or "scale via heuristic" (AlphaMate does). It is the
specific assembly — **exact column generation on `c≥0` + closed form per support
+ factored `Z` operator** — applied to the **dense genomic** regime, with cost
governed by `|S|` rather than `n³`.

## Open — must be cleared before any claim

1. **Exactness & wall-clock — cleared.** Head-to-head vs Clarabel (pool-unique,
   identical data, n=5000: Δgain 1.9e-9, same support) and vs optiSel (sex
   constraint, real wheat and pig GRMs): same optimum every time, measured
   speedups from 22× to ~2280×. Not extrapolation.
2. **True sexed OCS — cleared.** HS mice (BGLR `mice`, real GENDER 934♂/880♀,
   1814 individuals, real 10k-SNP GRM): support-first matches optiSel's optimum
   and solves in 0.008 s vs 6.96 s (~870×). The genuine reproduction constraint,
   not an arbitrary partition.
3. **`|S|` on broad real populations / tight k.** On real pig data `|S|`≈28 (huge
   speedup); a conservation-grade tight `k` would enlarge `|S|` and shrink the
   factor — still bounded by structure, but to be mapped.
4. **Rust port of the sexed version — cleared.** [`src/support_first.rs`](../src/support_first.rs)
   carries the sexed solver (`solve_sexed`, two equality rows) and the `c ≤ u`
   variants (`solve_capped` / `solve_sexed_capped`), all KKT-certified and
   cross-checked against Clarabel.
5. **AlphaMate — cleared.** AlphaMate targets a *distinct* problem (discrete mate
   allocation); on the continuous relaxation the two share, scored at matched
   coancestry, support-first's exact optimum is no worse — a consistency check, not
   a duel — at a fraction of the run time (see [`MANUSCRIPT.md`](MANUSCRIPT.md)).

## References

- Meuwissen 1997, *Maximizing the response of selection with a predefined rate of
  inbreeding*, J. Anim. Sci.
- Pong-Wong & Woolliams 2014, *A fast Newton–Raphson based iterative algorithm for
  large scale optimal contribution selection* (Gencont2).
- Wellmann 2019, *optiSel*, BMC Bioinformatics — `cccp` SOCP path.
- Gorjanc & Hickey 2018, *AlphaMate*, Bioinformatics.
- Goulart & Chen, *Clarabel*, Mathematical Programming Computation.
- Adachi, Iwata, Nakatsukasa 2017, *Solving the trust-region subproblem by a
  generalized eigenvalue problem*, SIAM J. Optim. (for the quadratic-objective
  variant).
