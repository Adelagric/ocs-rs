# Table 2 — measured numbers (provenance)

Apple M4 Max (14 cores), macOS 25.5.0 arm64; rustc 1.95.0 release; R 4.6.0;
optiSel 2.1.0. Idle machine (no LM Studio / Colima / other apps), serial, nothing
else running. Sexed OCS throughout.

Each timing is a **median over repetitions** taken back to back on the same built
instance: the matrix-free solve (`ocsrs`) and the algorithm prototype (`sf_at_ub.py`,
loaded once) are timed **5×**; optiSel **3×**. Run-to-run spread is sub-% except
where noted. Two prior full runs (one with LM Studio loaded, one idle) agreed with
these medians within noise, so the numbers are robust to machine load.

- **algorithm** column: the active set given a *dense* precomputed `G`, solve only —
  the NumPy prototype `research/repro/sf_at_ub.py`, the regime optiSel runs in (both
  handed `G`). This is what the manuscript's earlier 90×–2280× measured.
- **matrix-free** column: the Rust matrix-free solver via the R binding, end to end
  including the genotype copy, forming no `G` at all
  (`research/repro/r_binding_bench_reps.R`).

Instances identical across columns: prototype and matrix-free return the same gain
and support on every row (verified), so the comparison is like-for-like.

| dataset | n | m | support | algorithm (proto, G given) | matrix-free (shipped, e2e) | optiSel | algo×optiSel | shipped×optiSel |
|---|---|---|---|---|---|---|---|---|
| synthetic (structured) | 1000 | 500 | 59 | 0.023 s | 0.117 s | 2.06 s | 90× | 18× |
| synthetic (structured) | 2000 | 500 | 111 | 0.102 s | 0.280 s | 13.16 s | 129× | 47× |
| synthetic (structured) | 5000 | 500 | 168 | 0.506 s | 0.902 s | 163.8 s | 324× | 182× |
| CIMMYT wheat | 599 | 1279 | 24 | 0.003 s | 0.008 s | 0.595 s | 198× | 74× |
| PIC pig | 3534 | 52843 | 28 | 0.020 s | 0.622 s | 49.9 s | ~2500× | 80× |
| HS mouse (real sex) | 1814 | 10346 | 19 | 0.008 s | 0.057 s | 6.93 s | 866× | 122× |

Spread (5× ocsrs / 5× algo / 3× optiSel), the widest per cell:
- ocsrs: sub-% everywhere except pig [0.603, 0.744] (the 1.5 GB copy jitters).
- algorithm: sub-% except at the ms scale (wheat [0.003, 0.004]).
- optiSel: <1% on all cells including the slow ones (n=5000 [163.6, 164.2], pig
  [49.8, 50.3]).

Headline ranges: shipped × optiSel **18×–182×**; algorithm × optiSel **~90×–2500×**.

Notes:
- optiSel reproduces the manuscript's own numbers (wheat 0.60, mouse 6.93, pig ~50),
  so the instances match.
- The pig algorithm column (~2500×) is the origin of the manuscript's 2280×: the old
  figure was the algorithm comparison (solver vs optiSel, both with `G` formed), a
  legitimate like-for-like — not the shipped solver, which forms no `G` and does 80×
  end to end.
- The algorithm (given `G`) is faster than the matrix-free solve on **every** panel
  (wheat 0.003 vs 0.008, pig 0.020 vs 0.622, mouse 0.008 vs 0.057, …) — expected,
  since given `G` the solve avoids streaming `Z`. The matrix-free solver is not meant
  to beat a solve handed `G`; its point is never building or storing `G` at all. Both
  are one to three orders of magnitude faster than optiSel. (An earlier noisy
  wheat-algorithm timing of 0.011 s had suggested a spurious "matrix-free beats the
  prototype at small m" reversal — the median 0.003 s removes it.)
- An earlier single-sample run reported pig 88× and n=5000 193×; the medians here
  (80×, 182×) supersede them — the difference was optiSel's slow-cell single sample,
  now taken 3×.
- Exactness (all rows): support-first reaches the constraint boundary optiSel's IPM
  stops just inside; gain gaps +0.0001 (mouse) to +0.044 (synthetic n=5000), budget
  and sex split exact to 1e-9.
