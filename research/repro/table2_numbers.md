# Table 2 — measured numbers (provenance)

Apple M4 Max (14 cores), macOS 25.5.0 arm64; rustc 1.95.0 release; R 4.6.0;
optiSel 2.1.0. One serial session, nothing else running (contention inflates
optiSel severalfold). Sexed OCS throughout.

- **algorithm** column: the active set given a *dense* precomputed `G`, timing the
  solve only — the NumPy prototype `research/repro/sf_at_ub.py`, the same regime as
  optiSel (both handed `G`). This is what the manuscript's original 90×–2280×
  measured.
- **shipped** column: the Rust matrix-free solver via the R binding, end to end
  including the genotype copy, forming no `G` at all
  (`research/repro/r_binding_bench_all.R`).

Instances identical across the two columns: prototype and shipped return the same
gain and support on every row (verified), so the comparison is like-for-like.

| dataset | n | m | support | algorithm (proto, G given) | shipped (matrix-free, e2e) | optiSel | algo×optiSel | shipped×optiSel |
|---|---|---|---|---|---|---|---|---|
| synthetic (structured) | 1000 | 500 | 59 | 0.025 s | 0.120 s | 2.11 s | 84× | 18× |
| synthetic (structured) | 2000 | 500 | 111 | 0.105 s | 0.282 s | 13.56 s | 129× | 48× |
| synthetic (structured) | 5000 | 500 | 168 | 0.504 s | 0.901 s | 173.8 s | 345× | 193× |
| CIMMYT wheat | 599 | 1279 | 24 | 0.011 s | 0.009 s | 0.634 s | 58× | 70× |
| PIC pig | 3534 | 52843 | 28 | 0.022 s | 0.611 s | 53.8 s | 2445× | 88× |
| HS mouse (real sex) | 1814 | 10346 | 19 | 0.009 s | 0.057 s | 6.94 s | 771× | 122× |

Notes:
- optiSel reproduces the manuscript's own numbers (wheat 0.63, mouse 6.96, pig 54.8),
  so the instances match and are not contention-inflated.
- The pig algorithm column (2445×) is the origin of the manuscript's 2280×: the old
  figures were the algorithm comparison (solver vs optiSel, both with `G` formed), a
  legitimate like-for-like — not the shipped solver, which forms no `G` and does 88×
  end to end.
- Wheat: the shipped matrix-free solver (0.009 s) beats the dense-G prototype
  (0.011 s) — at small m the matrix-free products are cheaper than indexing a dense
  `G`. At pig's m=52843 the streaming cost reverses that.
- Exactness (all rows): support-first reaches the constraint boundary optiSel's IPM
  stops just inside; gain gaps +0.0001 (mouse) to +0.044 (synthetic n=5000), budget
  and sex split exact to 1e-9.
- Headline speed-up reported in the manuscript = shipped × optiSel (18×–193×). The
  algorithm column (58×–2445×) is reported as the algorithmic result, noting it
  assumes `G` built and stored — which the shipped solver avoids.
