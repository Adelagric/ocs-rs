# Table 1 — measured numbers (provenance)

Apple M4 Max (14 cores), macOS 25.5.0 arm64, rustc 1.95.0, release profile. Idle
machine, one serial sweep, nothing else running. Measured 2026-07-28 on commit
`bb8d90e` (the tree before the `Kinship` refactor; see REVISION.md, post-0.3.0 audit).

This comparison carries **no language confound**: both solvers are Rust, called in
the same process on the same generated dataset, by `src/main.rs::cmd_compare`.
Clarabel is the crate's own conic path (build `G` → Cholesky → cone assembly →
interior-point solve, Route A of `src/socp.rs`); support-first works from `Z` and
forms nothing.

## Operating point (the manuscript's Table 1)

`cargo run --release -- compare --n N --m 20000 --k-frac 0.017` — instances
`datagen::generate(n, m=20000, seed=20240617)`, cap `k = 0.017 × mean(diag G)`
(ridge excluded), which puts the support at ≈100 (ΔF ≈ 1 %). This is what
`research/repro/repro.sh` step 3 runs.

| n | m | Clarabel solve | Clarabel total¹ | support-first | speed-up (total) | support | G·c products |
|---|---|---|---|---|---|---|---|
| 1000 | 20000 | 2.109 s | 2.22 s | 1.433 s | 1× | 99 | 151 |
| 2000 | 20000 | 15.871 s | 16.2 s | 1.896 s | 9× | 112 | 165 |
| 5000 | 20000 | 306.4 s | 308 s | 2.358 s | 131× | ≈100² | 158 |
| 10000 | 20000 | 2289.8 s | 2298 s | 3.208 s | **716×** | ≈100² | 173 |

¹ total excluding the shared data generation: forming `G`, its Cholesky, cone
assembly, and the solve. support-first's column is its whole run.
² as tuned (the July log kept only the timing lines for these rows); n=1000 and
n=2000 re-logged on the fixed tree on 2026-09-20 (99 and 112; n=2000: Clarabel
18.3 s, support-first 2.14 s under load, gain Δ 4.6e-9).

At n=1000 the two are at parity: support-first's ~150 products at O(nm) match
Clarabel's O(n³) solve at small n; the gap opens as n³ overtakes n·m.

**Precision of the n=1000 row.** The same Clarabel build runs 2.11–2.14 s from a
binary of `bb8d90e` and 2.41–2.45 s from a binary of the current tree (same
iterations, same gain; a layout effect bisected to the presence of `src/packed.rs`
in the compiled crate — `research/repro/cvxgenrust/RESULTS.md`). The ratio at
n=1000 is therefore known to ±14 % and reads "parity" either way; from n=5000 on
this is immaterial. Table 1 used the faster baseline, the conservative direction.

**Independent baseline.** A Clarabel solver generated from a CVXPY model of the
same cone program by a third-party code generator (`cvxgenrust`) reproduces the
optimum and the per-iteration cost at n=500 and n=1000 —
`research/repro/cvxgenrust/RESULTS.md`.

**Regression note.** Between 2026-07-29 and 2026-09-20 (tag `v0.3.0` included) the
support-first column was not reproducible: a refactor had made the support-Gram
build stride through `Z` (7.5 s instead of 1.6 s at n=1000). Fixed; the fixed tree
reproduces this table (1.65 s at n=1000 under load). REVISION.md has the details.

## Loose cap (superseded headline, kept for the "~76000×" remark)

`compare --n N` at the default `--k-frac 0.6`, same instances. The cap is loose, so
the optimal support is 2–5 — the easy end of the support spectrum, which is why so
few `G·c` products are needed. The manuscript reports this only as the loose-cap
extreme; the operating-point table above is the headline.

| n | m | Clarabel solve | Clarabel total¹ | support-first | speed-up (total) | speed-up (solve) | support | Δgain |
|---|---|---|---|---|---|---|---|---|
| 1000 | 20000 | 1.602 s | 1.697 s | 0.006 s | 290× | 274× | 4 | 1.24e-9 |
| 2000 | 20000 | 12.308 s | 12.655 s | 0.012 s | 1093× | 1063× | 5 | 1.71e-9 |
| 5000 | 20000 | 180.933 s | 182.912 s | 0.018 s | 9917× | 9810× | 4 | 1.93e-9 |
| 10000 | 20000 | 1784.170 s | 1792.344 s | 0.024 s | **76179×** | 75831× | 2 | 4.84e-9 |

¹ total excluding the shared data generation: forming `G` (1.8 s at n=5000, 7.0 s at
n=10000), its Cholesky, cone assembly, and the solve. This is what a dense conic route
actually costs; support-first's column is its whole run, since it builds nothing.

**Stability.** support-first is deterministic at the reported precision: three runs at
n=1000 gave 0.006 s and three at n=2000 gave 0.012 s, identical to the millisecond.
Clarabel varies a few percent (n=1000: 1.602 / 1.655 / 1.680 s; n=2000: 12.308 /
12.717 / 12.746 s). The large-n cells are single runs — Clarabel needs ~30 min at
n=10000 — and their ratio's uncertainty is dominated by that few-percent Clarabel
spread, not by the millisecond support-first side.
