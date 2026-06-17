# Benchmark methodology and results

## Hardware / build

- **Machine**: Apple M4 Max, 14 cores, 36 GB unified memory (macOS, arm64).
- **Build profile**: `release` with `opt-level = 3`, `lto = "thin"` (see
  `Cargo.toml`). All numbers below are release numbers; debug numbers are not
  representative and are never quoted. A maximally-tuned `lto = "fat"` +
  `codegen-units = 1` build was *not* used — `thin` keeps compile times sane and
  is documented here so the figures are not mistaken for a peak-tuned build.
- **Toolchain**: stable `rustc` (edition 2024).

## Run hygiene (important — measurements at scale are environment-sensitive)

A first sweep was run with **LM Studio holding a ~9 GB local model resident**,
which drove macOS into heavy memory compression (~10 GB of compressed pages,
several GB of swap in use). Re-running on a quiet machine (LM Studio and Ollama
quit; compressor back to ~1.5 GB, no swap) isolated what the contamination
actually distorted — and it was **memory, not time**:

| n=10000 (Route A) | contaminated | quiet machine |
|---|---|---|
| solve time | 1851 s | 1815 s (within noise) |
| growth exponent (top two) | 3.31 | 3.25 |
| **peak RSS** | **6.1 GB** | **8.8 GB** |

Counter-intuitively the solve *time* barely moved, but peak RSS jumped by ~45%:
under memory pressure macOS compressed the benchmark's own pages, so the
contaminated run **under-reported the true working set**. The headline figures
below are from the quiet-machine run; its 8.8 GB peak is the honest memory
number. Lesson recorded so figures are never compared across unequal machine
states: *a memory measurement is only valid on a machine that is not already
under pressure.* (RSS is reported in GiB throughout — bytes ÷ 1024³ — matching
`REPORT.md`.)

## What is timed, and how

Wall-clock via `std::time::Instant` around each pipeline stage, reported per
`(n, route)` in `artifacts/scaling.csv`:

- `t_datagen_s` — synthetic genotype + GEBV generation.
- `t_grm_s` — VanRaden GRM build (`ZZᵀ/s` GEMM + ridge). Route A only.
- `t_factor_s` — dense Cholesky (Route A); includes any ridge-escalation
  re-factorings. ~0 for Route B (no factorization).
- `t_assemble_s` — building the Clarabel CSC `A` matrix.
- `t_solve_s` — `DefaultSolver::new` + `solve()` wall-clock.
- `clarabel_solve_s` — Clarabel's own internal solve-time counter (cross-check).

No `criterion`: these stages are multi-second one-shots, not micro-kernels, so a
single timed run per point is the honest measurement. `criterion`'s statistical
resampling would re-run multi-second/multi-minute solves dozens of times for no
extra signal and enormous wall cost. Each point is run **once**; the question is
order-of-magnitude scaling, not nanosecond variance.

### Peak memory

Each scaling point is re-spawned as a child process under `/usr/bin/time -l`;
the `maximum resident set size` line (bytes, macOS) is parsed into
`peak_rss_mb`. This is the child's own lifetime high-water mark, so memory is
attributed per `(n, route)` rather than conflated across the sweep. Reading RSS
in-process would require `getrusage` via libc (`unsafe`), which the crate bans,
so the external wrapper is the authoritative source.

## Results

Quiet-machine sweep, `m = 20000`, `k = 0.6 × mean GRM diagonal`, seed 20240617.
Canonical per-stage numbers are in `artifacts/scaling.csv`; the GO/NO-GO
interpretation is in `REPORT.md`. Re-generate everything with
`cargo run --release` (or just `REPORT.md` from existing CSVs with
`cargo run --release -- report`).

### Route A (Cholesky, n+1 cone) — the genomic-scale route

| n | GRM build | Cholesky | Clarabel solve | iters | status | peak RSS |
|---|---|---|---|---|---|---|
| 100 | 0.004 s | 0.000 s | 0.002 s | 9 | Solved | 0.03 GB |
| 500 | 0.029 s | 0.001 s | 0.21 s | 15 | Solved | 0.12 GB |
| 1000 | 0.096 s | 0.004 s | 1.66 s | 16 | Solved | 0.25 GB |
| 2000 | 0.342 s | 0.016 s | 12.84 s | 17 | Solved | 0.61 GB |
| **5000** | **1.91 s** | **0.16 s** | **191.3 s** | 15 | Solved | **2.63 GB** |
| 10000 | 7.38 s | 1.09 s | 1815.5 s | 17 | Solved | 8.80 GB |

Solve grows ≈ `n^2.95` (n=2000→5000) then `n^3.25` (n=5000→10000) — essentially
the cubic cost of the dense KKT factorization Clarabel performs each iteration,
with iteration count flat at 15–17. The GRM build (`ZZᵀ`, `O(n²m)`) and Cholesky
stay sub-10 s even at n=10000; the solve dominates by ~1700× the Cholesky and
~215× the total dense prep (see `REPORT.md` for why that is expected and not a
gate).

### Route B (raw Z, m+1 cone) — does not scale

| n | assemble | solve | iters | peak RSS |
|---|---|---|---|---|
| 100 | 0.004 s | 1.11 s | 10 | 0.24 GB |
| 1000 | 0.069 s | 81.5 s | 10 | 2.26 GB |

At equal n=1000, Route B is ~49× slower (81.5 s vs 1.66 s) and ~9× heavier
(2.26 GB vs 0.25 GB) than Route A, because its conic block is the fully dense
`m×n = 20000×1000` matrix. With `m ≫ n` in the genomic regime, only Route A is
viable; Route B is kept for the head-to-head the brief asked for.

### Interpreting the solve-vs-factorization ratio

The brief lists "solve time comparable to or below the factorization cost" as a
GO signal. Empirically it does **not** hold, and cannot for any interior-point
method here: faer factors the dense Cholesky in milliseconds-to-sub-second,
while Clarabel performs one KKT factorization *per IPM iteration* (~9–20
iterations observed). So the solve necessarily dwarfs the single Cholesky. The
ratio is reported in `REPORT.md` as **context**, not as a verdict gate; the
substantive gates are reliability (`Solved` + feasible), monotone frontier,
bounded conditioning, reaching n=10000 within 36 GB, and sub-`n^3.5` growth.
