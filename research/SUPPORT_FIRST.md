# Support-first OCS — research note

Status: **prototype, synthetic data only.** A promising lead, not a validated
result. Numbers below are reproducible with `python3 research/support_first.py`.

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

**Cost, matrix-free, iid** (`m=20000`):

| n | \|S\| | G·c products | wall (python) | vs Clarabel |
|---|---|---|---|---|
| 1000 | 3 | 4 | 0.009 s | 1.66 s |
| 5000 | 3 | 4 | 0.037 s | 191 s |
| 10000 | 3 | 4 | **0.072 s** | **1815 s** |
| 20000 | 2 | 3 | 0.111 s | (IPM impractical) |

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

## Where it sits vs the state of the art

| method | exact? | regime | what it exploits | limit |
|---|---|---|---|---|
| Meuwissen 1997 | yes | dense | Lagrangian on λ | iterative, dense |
| Gencont2 (Pong-Wong) | yes | **pedigree sparse** | Newton-λ + Gauss-Seidel, matrix-free | not the dense genomic `G`; first-order on c |
| optiSel `cccp` / this crate's Clarabel | yes | genomic dense | SOCP + IPM | O(n³) |
| AlphaMate | **no** (heuristic) | any | differential evolution | no optimality |
| **support-first (here)** | yes | **genomic dense** | **solution sparsity + factored Z + closed form** | validated on synthetic only |

The novelty is not "use a cone" (optiSel does), "matrix-free Lagrangian"
(Gencont2 does on pedigree), or "scale via heuristic" (AlphaMate does). It is the
specific assembly — **exact column generation on `c≥0` + closed form per support
+ factored `Z` operator** — applied to the **dense genomic** regime, with cost
governed by `|S|` rather than `n³`.

## Open — must be cleared before any claim

1. **Synthetic only.** The pedigree sim has 24 founders, so effective lineages
   (hence `|S|` at very tight k) are mechanically capped. A real, broader genetic
   base could give larger `|S|` at very tight k — still bounded by structure, not
   by n, but to be measured on real genotypes (e.g. optiSel's cattle data).
2. **Head-to-head exactness at n=10000** is KKT-certified, not yet checked
   against Clarabel on identical data at that scale.
3. **Wall-clock claim.** The ~10³–10⁴× gap is in operation counts; a Rust port
   (closed form + `Z` products via `faer`) is needed for a defensible time figure
   apples-to-apples with the crate's Clarabel path.

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
