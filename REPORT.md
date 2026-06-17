# OCS × Clarabel — go/no-go verdict

## Verdict: **GO**

| criterion | result | evidence |
|---|---|---|
| Clarabel `Solved` across sweep | ✅ pass | 8/8 points solved |
| Feasible vs original data (±1e-6) | ✅ pass | 8/8 points feasible |
| Frontier gain monotone in k | ✅ pass | max violation 0.00e0 |
| Conditioning (ridge ≤ 1e-3) | ✅ pass | max ridge used 1e-5 |
| Reached n=10000 | ✅ pass | largest solved n=10000 |
| Within 36 GB | ✅ pass | peak RSS 8.80 GB |
| Solve scaling ≤ n^3.5 | ✅ pass | empirical exponent ≈ 3.25 (t ∝ n^p, top two sizes) |
| Solve vs factorization (Route A, context) | — | solve/Cholesky up to 1673×, solve/(GRM+Cholesky) up to 214.5× |

## n=50 correctness

Feasibility margin `k − cᵀGc = 5.998e-6` (≥ 0 required). Solution dumped to `artifacts/correctness/` for independent cross-check (cvxpy / optiSel).

## Scaling sweep

| n | m | route | status | iters | ridge | GRM (s) | factor (s) | assemble (s) | solve (s) | solve/factor | feasible | peak RSS (GB) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 100 | 20000 | A_chol | Solved | 9 | 1e-5 | 0.004 | 0.000 | 0.000 | 0.002 | 35.92× | true | 0.03 |
| 500 | 20000 | A_chol | Solved | 15 | 1e-5 | 0.029 | 0.001 | 0.000 | 0.206 | 161.48× | true | 0.12 |
| 1000 | 20000 | A_chol | Solved | 16 | 1e-5 | 0.096 | 0.004 | 0.001 | 1.656 | 389.83× | true | 0.25 |
| 2000 | 20000 | A_chol | Solved | 17 | 1e-5 | 0.342 | 0.016 | 0.004 | 12.844 | 778.49× | true | 0.61 |
| 5000 | 20000 | A_chol | Solved | 15 | 1e-5 | 1.906 | 0.157 | 0.033 | 191.301 | 1221.40× | true | 2.63 |
| 10000 | 20000 | A_chol | Solved | 17 | 1e-5 | 7.377 | 1.085 | 0.129 | 1815.450 | 1672.64× | true | 8.80 |
| 100 | 20000 | B_raw | Solved | 10 | 0e0 | 0.000 | 0.000 | 0.004 | 1.115 | — | true | 0.24 |
| 1000 | 20000 | B_raw | Solved | 10 | 0e0 | 0.000 | 0.000 | 0.069 | 81.519 | — | true | 2.26 |


## Reading the result

Headline (n=5000, Route A): GRM build 1.91s, Cholesky 0.16s, Clarabel solve 191.30s in 15 iterations (Solved).

All sweep points returned a usable optimum that is feasible against the original G (recomputed from c, not read from solver internals); the gain/diversity frontier is monotone; and the ridge never had to exceed 1e-3. The conic solve does dominate the factorization — up to 1673× the (sub-second) Cholesky — which is expected, since an interior-point method performs one KKT factorization per iteration whereas the Cholesky is a single dense factor; measured against the *total* unavoidable dense prep (GRM build + Cholesky) it is a more modest 214.5×. solve time grows ~n^3.25 (marginally super-cubic) across the top two sizes; n=10000 solved in 1815.5s, 17 iters; peak RSS 8.8 GB < 36 GB. Clarabel's sparse IPM copes with the near-dense conic block at genomic scale: **GO** (offline, once-per-generation use).

## What this spike does NOT prove

- **Synthetic, well-conditioned data.** Independent markers give low relationships, so the GRM stays well-conditioned and ε=1e-5 sufficed at every size (zero ridge escalations). Real, highly-related populations can be far worse-conditioned and may force a larger ridge and/or more iterations — not exercised here.
- **One run per point.** Timings are single-shot wall-clock, not distributions; no warm-up/variance/p99. Order-of-magnitude scaling is the claim, not precise constants.
- **n=10000 took ~31 min.** Comfortable for an offline, once-per-generation decision; not interactive. Throughput is bounded by Clarabel's IPM, which dominates the (sub-second) factorization — not by the linear algebra.
- **Route B does not scale.** The raw-Z (m+1) cone is ~50× slower and far heavier than Route A at equal n (see table); only the Cholesky route (n+1 cone) is viable at genomic scale. This is a property of the formulation, reported rather than a defect.
- **Growth exponent from two points.** The reported exponent uses only the two largest sizes; a richer fit across all sizes would tighten it.
