# ocs-rs — support-first optimum contribution selection

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20746987.svg)](https://doi.org/10.5281/zenodo.20746987)

An **exact, matrix-free** solver for genomic Optimum Contribution Selection (OCS)
that exploits the sparsity of the *solution*: the optimal contribution vector
activates only a handful of candidates, so support-first finds that small support
by active-set column generation, solves each fixed support in closed form, and
never forms the dense `n×n` relationship matrix. It reaches the **same optimum**
as the domain's exact tool (optiSel) **orders of magnitude faster**, and is
validated on real genomic panels.

> Full write-up — derivation, tables, figure, references —
> [`research/MANUSCRIPT.md`](research/MANUSCRIPT.md).
> Reproduce every number in one command:
> [`bash research/repro/repro.sh`](research/repro/REPRODUCE.md).

## Results

On real marker panels (CIMMYT wheat n=599, PIC pig n=3534, heterogeneous-stock
mouse n=1814 with real sex):

- **Exact.** Agrees with a conic interior-point optimum to `1e-8`, and *saturates*
  the kinship bound that optiSel's IPM leaves slightly slack — at least as
  accurate as the domain tool, exact where it is merely close. Cross-language
  agreement with a NumPy reference: `1.5e-14`.
- **vs optiSel** (the standard exact tool): **90×–2280×** faster, same optimum
  (mouse: 0.008 s vs 6.96 s).
- **vs Clarabel** (a generic conic solver): up to **37090×** at n=10000
  (26 minutes → 43 ms).
- **vs AlphaMate** (the heuristic): strictly higher gain at *every* matched
  coancestry, at a small fraction of the run time.
- **Scales.** The optimal support stays ~15 as n grows to 40000, while the dense
  `G` every other solver forms reaches **11.9 GiB** (past laptop memory) — and
  merely *building* it costs more than the whole support-first solve — whereas
  the matrix-free `Z` footprint stays 0.30 GiB and the solve stays under 0.1 s
  (Figure 1: [`research/fig_scaling.pdf`](research/fig_scaling.pdf)).

## The model

OCS picks how much each of `n` candidates contributes to the next generation:
maximise genetic gain subject to a cap on average coancestry (to preserve
diversity).

```
maximize    bᵀc                  b = genomic estimated breeding values (GEBV)
subject to  A c = d              budget: Σcᵢ = 1, or sexed Σ_males = Σ_females = ½
            c ≥ 0
            cᵀ G c ≤ k           G = VanRaden genomic relationship matrix, k = kinship bound
```

`G = ZZᵀ/s + εI` is symmetric positive definite (ridge `ε`); `cᵀGc ≤ k` is a
convex quadratic constraint (a second-order cone). The *sexed* form — the true
OCS, each mating drawing one parent of each sex — replaces the simplex with two
equality rows.

## How it works

- **Active-set column generation.** Seed the support with the best male and
  female; price candidates by reduced cost and add the best, toward a *single*
  fixed coancestry cap `k`; drop negatives. A feasible point with no positive
  reduced cost is KKT-optimal — exact, not heuristic.
- **Closed form per support.** Eliminate the equality constraints through a
  `q×q` reduction `P = A_S G_S⁻¹ A_Sᵀ` (q ≤ 2) and read the multiplier off a
  scalar quadratic — one Cholesky, no inner iteration.
- **Matrix-free.** `Gc = εc + Z(Zᵀc)/s`; the dense `n×n` `G` is never formed.

Derivation in [`research/MANUSCRIPT.md`](research/MANUSCRIPT.md) (Methods);
implementation in [`src/support_first.rs`](src/support_first.rs)
(`solve` simplex, `solve_sexed` sexed).

## Running

```sh
# support-first head-to-head and benchmarks
cargo run --release -- compare --n 5000          # Clarabel vs support-first, same optimum
cargo run --release --example bench_sexed         # sexed solver: release timing + cross-language export
cargo run --release --example scaling_matrixfree  # Figure 1 data: support/time bounded vs dense-G blow-up

# the original Clarabel evaluation harness
cargo run --release                               # = `all`: correctness + frontier + scaling -> REPORT.md
cargo run --release -- scale --n 5000 --route a   # one timing point (CSV row on stdout)
```

`stdout` is data; progress goes to `stderr`. Flags:
`--n --m --seed --k-frac --route a|b --max-iter --time-limit --max-n`.

## Reproduction

```sh
bash research/repro/repro.sh
```

Runs everything the local toolchain allows (Rust timings, Figure 1, the R GRM
exports via `BGLR`, support-first vs optiSel) and skips — with a message —
anything whose dependency or dataset is absent. The PIC pig panel needs a
one-time manual download (URL + layout in
[`research/repro/REPRODUCE.md`](research/repro/REPRODUCE.md)); AlphaMate runs from
a Linux binary under Colima + Rosetta (recipe in the same file).

## Verification

```sh
cargo test                                 # KKT certificates, feasibility, monotonicity, CSC properties
cargo clippy --all-targets -- -D warnings
cargo fmt --check
```

A `Solved` result is always feasible and on the budget (sex-split for the sexed
solver), certified across a range of caps `k`. The optimum is cross-checked
independently of this crate's arithmetic: against Clarabel (conic IPM), a NumPy
reference (`1.5e-14`), SciPy SLSQP, and optiSel on real data.

## Origin — the Clarabel spike

This began as a go/no-go on whether the [Clarabel](https://clarabel.org) conic
interior-point solver handles genomic-scale OCS (verdict: **GO**, see
[`REPORT.md`](REPORT.md), [`BENCHES.md`](BENCHES.md), [`DECISIONS.md`](DECISIONS.md)).
That spike exposed the opening: the conic IPM pays `O(n³)` per iteration to
describe a solution that activates a handful of candidates. support-first
exploits exactly that. Clarabel is kept as an independent cross-check oracle.

## Honest caveats

`b` is a recorded phenotype or EBV standing in for a true genomic breeding value
on the public panels; a genuine recorded sex exists only for the mouse panel
(arbitrary balanced split elsewhere); the optiSel head-to-head times a NumPy
prototype against R (the gap is algorithmic, not language); the support bound is
empirical, not proven; and the solver handles a single quadratic constraint and
continuous contributions (no integer mate allocation). These are stated in the
manuscript's Discussion.

## Constraints honoured

Pure-Rust dependencies only (`clarabel`, `faer`, `rand`, `rand_distr`) — no
system BLAS/LAPACK, no FFI, no `unsafe`. Versions pinned exactly. Single-author.
