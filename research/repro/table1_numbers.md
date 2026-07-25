# Table 1 — measured numbers (provenance)

Apple M4 Max (14 cores), macOS 25.5.0 arm64, rustc 1.95.0, release profile. Idle
machine, one serial sweep, nothing else running.

Unlike Table 2, this comparison carries **no language confound**: both solvers are
Rust, called in the same process on the same generated dataset, by
`src/main.rs::cmd_compare` (`cargo run --release -- compare --n N`). Clarabel is the
crate's own conic path (build `G` → Cholesky → cone assembly → interior-point solve);
support-first works from `Z` and forms nothing.

Instances: `datagen::generate(n, m, seed=20240617)` with `m = max(2n, 20000)` = 20000
throughout, cap `k = 0.6 × mean(diag G)`.

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

**These supersede the earlier Table 1** (Clarabel 1.40 / 10.6 / 160 / 1579 s,
support-first 0.011 / 0.022 / 0.036 / 0.043 s, 126×–37090×), which predated the
incremental-Gram optimisation and carried no stated provenance. support-first is now
roughly twice as fast, so the earlier table understated it.

**Caveat to carry into the text.** The cap here (`k_frac = 0.6`) is loose, so the
optimal support is 2–5 — the easy end of the support spectrum. That is why only 3–6
`G·c` products are needed. Table 2's synthetic rows sit at the opposite end (tight
cap, support 59–168) and show the factor shrinking accordingly. Both regimes are
reported rather than only the favourable one.
