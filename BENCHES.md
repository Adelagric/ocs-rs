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

## 2026-07-25 — the matrix-free solver, timed on real panels for the first time

Hardware: Apple M4 Max (14 cores), macOS 25.5.0 arm64. Toolchain: rustc 1.95.0
(release profile), R 4.6.0, optiSel 2.1.0. Genotypes reach the Rust core through
the R binding (`bindings/r`) so support-first and optiSel run in one session on one
instance; the matrix-free times below are the whole `ocs_solve` call including the
R→Rust copy unless stated otherwise.

**Why this section exists.** The manuscript's Table 2 credits its support-first
timings to a NumPy prototype (`research/repro/sf_at_ub.py`) that consumes the dense
kinship matrix `K` and slices `K[S,S]` by indexing. Every export script writes only
`K`, so the *matrix-free* path — the one the crate ships — had never been timed on
real data. It has now.

### optiSel vs the shipped matrix-free solver, same session

Sexed OCS. Synthetic rows use a tight cap (`ub = 1.04 × mean coancestry`), which
forces a large support — the hard case; wheat/mouse use each export script's own
cap (`0.12`/`0.15 × k_greedy`).

| instance | n | m | support | optiSel | matrix-free | speed-up |
|---|---|---|---|---|---|---|
| synthetic (structured) | 1000 | 500 | 55 | 2.10 s | 0.118 s | 18× |
| synthetic (structured) | 2000 | 500 | 105 | 13.46 s | 0.256 s | 53× |
| synthetic (structured) | 5000 | 500 | 152 | 165.8 s | 0.830 s | 200× |
| CIMMYT wheat (real GRM) | 599 | 1279 | 24 | 0.617 s | 0.008 s | 77× |
| HS mouse (real GRM, real sex) | 1814 | 10346 | 19 | 7.04 s | 0.061 s | 125× |

optiSel reproduces its Table 2 numbers (wheat 0.62 vs 0.63, mouse 7.0 vs 6.96), so
the instances match and the whole difference sits in the support-first column. These
matrix-free speed-ups (15×–200×) are lower than Table 2's 90×–2280×, because the
prototype's numbers exclude the cost of *forming* `G` and use dense `K[S,S]` slicing
where the matrix-free path recomputes from `Z`. Exactness is unchanged and
re-confirmed: every row reaches the constraint boundary optiSel's IPM stops just
inside, with a small positive gain gap, budget and sex split exact to 1e-9.

### Charge the dense route for building G

The prototype gets `K` for free. Built end to end, the dense GRM costs (this
hardware): wheat 0.007 s, mouse 0.129 s. So on the mouse the matrix-free solve
(0.032 s, marshalling removed) beats *merely constructing* the dense `G` (0.129 s) —
the 870× in Table 2 came from a 0.008 s that excluded that 0.129 s and is not
reproducible under any honest end-to-end accounting.

### The dense/matrix-free crossover is two-axis, not one

Comparing the whole matrix-free solve against dense `G` construction alone:

- At modest `m` (=1000), matrix-free wins at every `n` measured — it finishes a full
  solve faster than a dense solver can *build* `G` (n=1000: 4.1 ms vs 7 ms build;
  n=30000: 30 ms vs 4.65 s build). `examples/scaling_matrixfree`.
- At marker-rich `m` (=10000, mouse-like), building `G` is cheaper for small `n`
  (n≤~1500), and matrix-free takes over above that. The real panels fall on the
  matrix-free side, but the win is not universal — it is a genuine two-axis
  trade-off (dense build O(n²m) vs a solve whose cost tracks n, m and the support).

Numbers are noisy in the support because support size varies with the cap and the
instance, and it — not n or m alone — sets the per-solve work.

### What made the matrix-free solve fast enough to say this (same hardware)

Two changes, each leaving the optimum bit-identical:

- **Incremental Gram + cached support rows** (`build_gss` → `GramCache`). Rebuilding
  the |S|×|S| Gram every iteration was ~90% of a mouse solve; keeping it and the
  support's rows of `Z` across iterations, and gathering those rows once so later dot
  products scan contiguously, cut the solve roughly in half and more:
  `scaling_matrixfree` (m=1000, sexed) n=1000 0.0083→0.0041 s, n=40000 0.0769→0.0409 s;
  real panels wheat 0.041→0.008 s, mouse 0.624→0.061 s.
- **Bytes marshalling in the Python binding.** Once the solve was milliseconds, pyo3
  extracting a `Vec<f64>` element-by-element dominated a call (n=16000 m=2000: 714 ms).
  Passing `Z.tobytes()` and borrowing `&[u8]`: 714→79 ms; n=1000 44→1.5 ms.

The `G·c` product itself was only ~8% of a mouse solve, not the bottleneck the design
assumed; exploiting the sparsity of `c` there still gives 1.5×–1.8× on the product
(`examples/bench_matvec`), kept because it is free, but it is not where the time was.

### Update — Table 1 re-measured (Clarabel), the language-confound-free comparison

The Clarabel head-to-head had never been re-run after the incremental-Gram work, and
carried no stated provenance, so it sat awkwardly beside a Table 2 that is now
medians with a provenance file. Re-run as one serial sweep on the idle machine
(`compare --n N`, m = 20000, k_frac 0.6, seed 20240617). Both solvers are Rust in the
same process on the same data — this is the comparison that carries **no language
confound**, and it is the honest answer to "you timed Rust against R".

| n | Clarabel solve | Clarabel total¹ | support-first | speed-up (total) | support |
|---|---|---|---|---|---|
| 1000 | 1.602 s | 1.697 s | 0.006 s | 290× | 4 |
| 2000 | 12.308 s | 12.655 s | 0.012 s | 1093× | 5 |
| 5000 | 180.933 s | 182.912 s | 0.018 s | 9917× | 4 |
| 10000 | 1784.170 s | 1792.344 s | 0.024 s | **76179×** | 2 |

¹ forming `G` (7.0 s at n=10000) + Cholesky + cone assembly + solve. Gains agree to
1.2e-9–4.8e-9 on every row.

This **supersedes the earlier Table 1** (126×–37090×): support-first is now about
twice as fast, so the old table understated it. Stability: support-first is
deterministic at this precision (three runs each at n=1000 and n=2000 reproduced
0.006 s and 0.012 s to the millisecond); Clarabel varies a few percent (1.602/1.655/
1.680 s at n=1000). Large-n cells are single runs — Clarabel needs ~30 min at
n=10000.

Caveat carried into the paper: this sweep uses a **loose** cap (k_frac 0.6), so the
support is 2–5 — the easy end. Table 2's synthetic rows are the tight-cap end
(support 59–168) and the factor shrinks there. Both are reported.

### Update — pig measured, and the R binding made zero-copy

The PIC pig panel (n=3534, **52 843 SNP**) was found locally and measured. It exposed
the next bottleneck: the R binding copied the genotype matrix three times (R
`as.numeric` flatten → extendr `Vec<f64>` → faer build), and at 1.5 GB that dominated
the call. Taking the matrix as a borrowed `Robj` (`as_real_slice`, R's own
column-major buffer, one gather into faer) removed two of the three passes. Optimum
bit-identical (28 testthat assertions, all recomputing the answer with R's arithmetic).

Pig, clean, before → after the borrow: end-to-end **1.96 s → 0.61 s**, of which the
copy **1.96 s → 0.276 s**; speed-up vs optiSel **30× → 88×**. The solve itself is
~0.34 s (52k markers streamed per product), not the ~1 ms an earlier noisy
`t_full − t_copy` subtraction suggested — that subtraction was wrong and is not used.

**Clean serial run, all panels, Apple M4 Max, one R session, nothing else running**
(`research/repro/r_binding_bench_all.R`). Shipped matrix-free solver, end to end
(including the R→Rust copy), vs optiSel:

| panel | n | m | optiSel | matrix-free | copy | support | speed-up |
|---|---|---|---|---|---|---|---|
| synthetic (structured) | 1000 | 500 | 2.11 s | 0.120 s | 0.000 s | 59 | 18× |
| synthetic (structured) | 2000 | 500 | 13.56 s | 0.282 s | 0.000 s | 111 | 48× |
| synthetic (structured) | 5000 | 500 | 173.8 s | 0.901 s | 0.004 s | 168 | 193× |
| CIMMYT wheat | 599 | 1279 | 0.634 s | 0.009 s | 0.000 s | 24 | 70× |
| PIC pig | 3534 | 52843 | 53.8 s | 0.611 s | 0.276 s | 28 | 88× |
| HS mouse (real sex) | 1814 | 10346 | 6.94 s | 0.057 s | 0.010 s | 19 | 122× |

optiSel reproduces its manuscript numbers (wheat 0.63, mouse 6.96, pig 54.8), so the
instances match. Range **18×–193×**, every real panel measured. Exactness re-confirmed
on all rows: support-first reaches the constraint boundary optiSel stops inside, with
a small positive gain gap, budget and sex split exact. Synthetic rows are at a tight
cap (large support 59–168, the hard case); the factor grows with `n` as optiSel's
interior-point cost climbs. On the pig the copy (0.276 s) is still the largest single
component of the 0.611 s — the same lesson as the Python `tobytes` fix, one binding
behind: a true zero-copy `MatRef` into the core (no gather) would cut it further.
