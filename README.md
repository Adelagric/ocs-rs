# ocs-rs

A go/no-go **spike**, not a product. It answers one question with measured
numbers: *can the [Clarabel](https://clarabel.org) conic interior-point solver
handle genomic-scale Optimum Contribution Selection (OCS) reliably and fast,
given that the second-order-cone block the problem produces is effectively
dense?*

The verdict, with the actual numbers from this machine, is written to
[`REPORT.md`](REPORT.md) by the `all` command. Timing methodology and hardware
are in [`BENCHES.md`](BENCHES.md); design rationale in [`DECISIONS.md`](DECISIONS.md).

## The model

OCS picks how much each of `n` candidates contributes to the next generation:
maximise genetic gain subject to a cap on average coancestry (to preserve
diversity).

```
maximize    bᵀc                  b = genomic estimated breeding values (GEBV)
subject to  Σ cᵢ = 1
            c ≥ 0
            cᵀ G c ≤ k            G = genomic relationship matrix (VanRaden), k = kinship bound
            c ≤ u                 optional per-candidate cap (off by default)
```

`G` is symmetric positive semi-definite. `cᵀGc ≤ k` is a convex quadratic
constraint, recast here as a **second-order cone** `‖Fᵀc‖₂ ≤ r` in two ways:

| | factor `F` | radius `r` | SOC dimension | conic matrix |
|---|---|---|---|---|
| **Route A** (Cholesky) | `L` where `G = LLᵀ` | `√k` | `n+1` | `Lᵀ` (dense lower-tri, transposed) |
| **Route B** (raw) | `Z` (centred genotypes), `G = ZZᵀ/s` | `√(k·s)` | `m+1` | `Zᵀ` (fully dense) |

Route B skips the Cholesky entirely (`cᵀGc = ‖Zᵀc‖²/s`), at the cost of an
`m+1` cone. With `m` (markers) ≫ `n` (individuals) — the genomic regime — Route
A's `n+1` cone is far smaller, so it is the default; Route B is benchmarked at
two sizes for the head-to-head.

### Clarabel mapping

Clarabel solves `min ½xᵀPx + qᵀx s.t. Ax + s = b, s ∈ K`. OCS maps to it with
`P = 0`, `q = -b`, and a stacked `A`/`b`/cone list:

| block | rows | `A` | `b` | cone |
|---|---|---|---|---|
| `Σcᵢ = 1` | 1 | `1ᵀ` | `1` | `ZeroConeT(1)` |
| `c ≥ 0` | `n` | `-I` | `0` | `NonnegativeConeT(n)` |
| `‖Fᵀc‖ ≤ r` | `d+1` | `[0ᵀ; -Fᵀ]` | `[r, 0…]` | `SecondOrderConeT(d+1)` |
| `c ≤ u` (opt) | `n` | `I` | `u` | `NonnegativeConeT(n)` |

The SOC slack is `s = b − Ax = (r, Fᵀc)`, hence `‖Fᵀc‖ ≤ r`. `A` is built
directly in compressed-sparse-column form. **The point of the spike**: that
`-Fᵀ` block (`Lᵀ` or `Zᵀ`) is dense, so the conic part of `A` is near-dense —
whether Clarabel's sparse IPM copes with that at scale is exactly what is
measured.

### Synthetic data (seeded, reproducible)

- allele frequencies `pⱼ ~ Uniform(0.05, 0.5)`
- dosages `Mᵢⱼ ~ Binomial(2, pⱼ) ∈ {0,1,2}`
- VanRaden: `Z = M − 2p`, `G = ZZᵀ / (2 Σⱼ pⱼ(1−pⱼ))`, ridged `G + εI`
  (`ε = 1e-5` default; auto-escalated and reported if Cholesky needs more)
- `b ~ N(0, 1)`

## Running

```sh
cargo run --release                       # = `all`: correctness + frontier + scaling sweep + REPORT.md
cargo run --release -- all --max-n 2000   # cap the sweep (faster smoke run)

cargo run --release -- correctness        # n=50, dump artifacts, assert feasibility (exit 0/1)
cargo run --release -- frontier           # sweep k, assert gain monotone in k -> artifacts/frontier.csv
cargo run --release -- scale --n 5000 --route a   # one point; prints a CSV row to stdout
cargo run --release -- report             # rebuild REPORT.md from existing artifacts/ (no re-solving)
```

Flags: `--n --m --seed --k-frac --route a|b --max-iter --time-limit --max-n`.
`stdout` is data (the `scale` CSV row); progress goes to `stderr`.

Artifacts land in `artifacts/`:
- `correctness/` — `grm_true.csv`, `b.csv`, `c.csv`, `meta.csv` for independent
  cross-check (the dumped problem is reproducible in cvxpy/optiSel/numpy)
- `frontier.csv` — `(k, diversity cᵀGc, gain, feasible)`
- `scaling.csv` — per-(n, route) timing split, status, iterations, peak RSS

Peak RSS is measured by re-spawning each scaling point under `/usr/bin/time -l`
(macOS), so each number is that process's own high-water mark.

## Verification

```sh
cargo test                                # feasibility + monotonicity + CSC-assembler properties
cargo clippy --all-targets -- -D warnings
cargo fmt --check
cargo run --release --example spike       # hand-checked 2-variable optimum (kept API regression test)
```

The n=50 solution has been cross-checked two ways independent of this crate's
own arithmetic: feasibility recomputed in numpy from the dumped CSVs, and the
optimum re-solved with SciPy SLSQP (a different solver) — both agree to ~1e-9.
See [`REPORT.md`](REPORT.md).

## Constraints honoured

Pure-Rust dependencies only (`clarabel`, `faer`, `rand`, `rand_distr`) — no
system BLAS/LAPACK, no FFI, no `unsafe`. Versions pinned exactly. Single-author.
