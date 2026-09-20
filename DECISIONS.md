# Design decisions

Dated log of non-trivial choices for the OCS / Clarabel go-no-go spike. Newest
last within each section. Re-read before revisiting a related decision.

## Scope

This crate is a **decision spike**, not the OCS tool. It exists to answer one
question with measured numbers: *does Clarabel's sparse interior-point method
solve genomic-scale OCS reliably and fast, given that the second-order-cone
block it must factor is effectively dense?* Everything not needed to answer that
is out of scope (real data IO, multi-trait, mate allocation, a polished CLI).

## Dependencies (pinned exactly, 2026-06-17)

| crate | version | why | notes |
|-------|---------|-----|-------|
| `clarabel` | =0.11.1 | the solver under evaluation | pure Rust; default `qdldl` direct solver (also pure Rust) |
| `faer` | =0.24.0 | dense GEMM (GRM build) + Cholesky | `default-features = false`, only `linalg`+`std`+`rayon`. Drops `sparse`, `npy`, and faer's internal `rand` (which pulled a second `rand` 0.9 major). |
| `rand` | =0.10.1 | seeded RNG (`StdRng`) | |
| `rand_distr` | =0.6.0 | Uniform / Binomial / Normal samplers | resolves to the same `rand_core` 0.10 as `rand`, so samplers interop with our `StdRng` |

No system BLAS/LAPACK, no `unsafe`, no FFI. `criterion` deliberately *not*
used: the operations timed here are large one-shot pipeline stages (seconds),
not micro-kernels, so wall-clock with explicit warm-up/repeat is the honest tool;
criterion's resampling would re-run multi-second solves dozens of times for no
added signal. Timing methodology lives in BENCHES.md.

### Dual-`rand` coexistence — resolved, not worked around

`cargo add` initially produced `rand` 0.9 *and* 0.10 in the tree. Tracing showed
0.9 came **only** from faer's optional `rand` feature, which this crate never
touches. Disabling faer's default features removed it entirely (`cargo tree -i
rand@0.9.4` now errors — gone), so there is a single `rand`/`rand_core` 0.10.

## Numerics

- **VanRaden GRM**: `Z = M - 2p` (column-centred dosages), `G = Z Zᵀ / s` with
  `s = 2 Σ pⱼ(1-pⱼ)`. `G` is symmetric PSD but only positive *semi*-definite
  when `m < n` (rank ≤ m), so a ridge `G + εI` is added before Cholesky. Default
  `ε = 1e-5`. The code reports and auto-escalates ε if `llt` returns
  `NonPositivePivot`, and records the final ε used.
- **f64 throughout.** Genomic relationship values are O(1); the conditioning
  risk is from near-duplicate individuals, addressed by the ridge, not by
  precision.

## Two SOC formulations (the comparison the spike is built around)

The kinship bound `cᵀGc ≤ k` becomes a second-order cone `‖Fᵀc‖₂ ≤ r`:

- **Route A (Cholesky)**: `G = LLᵀ`, `F = L`, `r = √k`. SOC dimension `n+1`.
  The conic constraint matrix is `Lᵀ` — dense lower-triangular transposed, i.e.
  `n×n` with `n(n+1)/2` nonzeros.
- **Route B (raw factor)**: `cᵀGc = ‖Zᵀc‖²/s`, `F = Z`, `r = √(k·s)`. SOC
  dimension `m+1`. No Cholesky. The conic matrix is `Zᵀ` — fully dense `m×n`.

Default route = whichever of `{n+1, m+1}` is smaller. With `m ≥ 20000` and
`n ≤ 10000`, that is always Route A in the scaling sweep; Route B is still run
at two sizes for the head-to-head the brief requires.

## Clarabel cone mapping (verified against a hand-solved optimum)

Solver form: `min ½xᵀPx + qᵀx s.t. Ax + s = b, s ∈ K`. With `P = 0`, `q = -b`:

| constraint | rows | A block | b block | cone |
|------------|------|---------|---------|------|
| `Σcᵢ = 1` | 1 | `1ᵀ` | `1` | `ZeroConeT(1)` |
| `c ≥ 0` | n | `-I` | `0` | `NonnegativeConeT(n)` |
| `‖Fᵀc‖ ≤ r` | d+1 | `[0ᵀ; -Fᵀ]` | `[r, 0…]` | `SecondOrderConeT(d+1)` |
| `c ≤ u` (opt) | n | `I` | `u` | `NonnegativeConeT(n)` |

`s = b − Ax` makes the SOC slack `(r, Fᵀc)`, so `‖Fᵀc‖ ≤ r`. Validated in
`examples/spike.rs`: n=2, G=I, b=[1,2], k=0.6 → Clarabel returns
x=[0.276393, 0.723607], matching the closed-form c₂=(1+√0.2)/2 to 1e-6.

### CSC assembly

`A` is assembled **directly in compressed-sparse-column form** (not triplets):
for each candidate column j the entries land in strictly increasing row order by
construction (`0 < 1+j < 2+n+i < …`), which is exactly what `CscMatrix::new`
requires and avoids the O(nnz log nnz) sort and the ~3× triplet memory that
`new_from_triplets` would cost at Route-B scale (nnz ≈ n·m ≈ 2×10⁸). A unit test
cross-checks the direct builder against `new_from_triplets` on small random
instances.

## Feasibility is checked from ORIGINAL data, never solver internals

`Σc`, `c ≥ 0`, and `cᵀGc` are recomputed from the returned `c` and the original
`G` (not from Clarabel's `s`/`z`). This is what makes a "Solved" status
trustworthy: a wrong CSC assembly would still report Solved but fail this check.

## Peak memory

`/usr/bin/time -l` wraps each scaling run (the `maximum resident set size`
line). Reading RSS in-process would need `getrusage` via libc → `unsafe`, which
is banned, so the external wrapper is the authoritative source. The binary also
prints an *analytical* lower bound (G + Z + L + A nnz) labelled as an estimate.

## The support-Gram cache keeps a contiguous copy of the dense support rows (2026-09-20)

`Z` is column-major, so a row of it strides by `n`; building `G_S` one entry at a
time through `MatRef` indexing misses the cache on every marker. The `Kinship`
refactor (`e10e2f0`) dropped the row buffer `GramCache::zs` for a representation-
agnostic `gram(i, j)` and cost 4.6× on the n=1000 / m=20000 sweep (BENCHES.md).
Restored through `Kinship::dense_z() -> Option<(MatRef, f64)>`: the cache gathers
each support row once and takes contiguous dot products when the operator exposes a
dense `Z`, and keeps the entry-wise path otherwise (the packed columns are already
contiguous, and `PackedGeno::gram_entry` is the right primitive there). Alternatives
rejected: interior mutability inside `DenseZ` (a `RefCell` row cache behind a `&self`
trait method — hides state, adds a borrow on the hot path); `&mut self` trait methods
(forces a wrapper for the shared `&PackedGeno`). Summation order over the markers
is ascending on both paths, so the Gram entries — and the optimum — are
bit-identical across representations (`packed_solver_matches_dense`). Invariant to
preserve: **whatever the operator, `G_S` must never be built by striding through a
column-major `Z`.**
