# Independent conic baseline — cvxgenrust check (provenance)

Question: is the Clarabel baseline of Table 1 hand-assembled in a way that
flatters support-first? Answer it by having a third party build the baseline:
`cvxgenrust` 0.1.0 (Zhu & Boedecker 2026, arXiv:2609.13875) generates a Rust crate
that calls Clarabel natively from a CVXPY model of the same cone program,

```python
c = cp.Variable(n); L = cp.Parameter((n, n)); b = cp.Parameter(n); r = cp.Parameter(nonneg=True)
cp.Problem(cp.Maximize(b @ c), [cp.norm(L.T @ c, 2) <= r, cp.sum(c) == 1, c >= 0])
```

with `L` the Cholesky factor of `G + εI` and `r = √k` — Route A of `src/socp.rs` —
and both solvers are timed on the identical instance (`datagen` seed 20240617,
`m = 20000`, cap `k = 0.017 × mean diag G`, the Table-1 operating point). Same
pinned `clarabel = 0.11.1`, same release profile (`opt-level 3`, thin LTO), default
Clarabel settings on both sides. Script: `run.sh`; the crate's own conic route is
timed by `dump_ocs`, the generated crates by `ocs_bench.rs.tmpl` (dense `L`) and
`ocs_sparse_bench.rs.tmpl` (`L` declared triangular).

Measured 2026-09-20 on the Apple M4 Max (14 cores), rustc 1.98.1, Python 3.13.13,
cvxpy 1.9.3. **Not an idle machine**: unrelated jobs kept the load average at
12–19 throughout. The serial Clarabel timings were only ~3 % above the idle Table-1
values (2.11–2.14 s vs 2.1 s at n=1000); the numbers below are same-session,
interleaved, best of 3 unless stated.

## Same optimum, same per-iteration cost

| n | crate Route A | cvxgenrust (dropzeros) | cvxgenrust (as generated) | gain, crate / generated |
|---|---|---|---|---|
| 500 | 0.251 s, 16 it, 15.7 ms/it | 0.288 s, 19 it, 15.2 ms/it | 1.426 s, 19 it | 1.601139067 / 1.601139056 |
| 1000 | 2.41 s, 21 it, 115 ms/it (current tree)¹ | 2.20 s, 19 it, 116 ms/it | 6.90 s, 19 it | 1.991005810 / 1.991005826 |

Support 100 / 99 on both sides; Σc = 1 to 1e-6; both `Solved`.

- The generated program is the same cone program with one extra scalar variable
  and one extra nonnegative row (CVXPY's canonicalisation of the norm cone), hence
  ±3 iterations; the **per-iteration cost agrees within 3 %** at both sizes
  against a binary of the current tree (15.2 vs 15.7 ms/it; 116 vs 115 ms/it),
  and is 16 % *above* the binary Table 1 was measured with (100 ms/it¹) — the
  hand-assembled baseline is, if anything, the faster one.
- **Explicit zeros — a modelling choice, not a generator defect.** Declared as a
  dense `n×n` parameter (the first spike), `L` reaches Clarabel with the strict
  upper triangle of `Lᵀ` as `n(n−1)/2` explicit zeros (124 750 at n=500, 499 500 at
  n=1000; `A` nnz = n²+2n+2): the generated solver is 4.9× (n=500) / 3.1× (n=1000)
  slower unless Clarabel's `input_sparse_dropzeros = true` is set. **Declared with
  its pattern** — `cp.Parameter((n, n), sparsity=np.tril_indices(n))`
  (`gen_ocs_sparse.py`) — the generator carries only the triangle: `A` nnz
  126 252 / 502 502 = n(n+1)/2+2n+2, zero explicit zeros, and the solve is at parity
  with no setting (n=500: 0.253 s, 19 it; n=1000: 1.98 s, 19 it, 104 ms/it — the
  fastest binary of the set). The setter then expects the entries in the order of
  the sparsity index arrays as given (`np.tril_indices` is row-major over the
  triangle; the column-major packing solves a different problem, |S|=1), and the
  generated `ParameterInfo.shape` reports the pattern length (`[125250]`) rather
  than the matrix shape. The crate's `socp.rs` emits only the triangle, which is
  the same thing done by hand.

¹ The crate's Route A takes 2.11–2.14 s (21 it, 100 ms/it) with a binary built
from commit `bb8d90e`, the tree Table 1 was measured with, and 2.41–2.45 s with a
binary of the current tree — see the next section.

## Build-to-build variance of the baseline (≈14 %)

The same `socp.rs` + `solve.rs` + pinned `clarabel 0.11.1`, same iterations (21),
same gain, but the IPM takes 2.11–2.14 s when `ocs_rs` is built from `bb8d90e`
and 2.41–2.45 s when built from the current tree (fresh target directories;
`lto = thin`, `fat` and `off` all give the latter). Bisected to the mere presence
of `src/packed.rs` in the compiled crate — code that this path never executes and
that is not even linked into the `dump_ocs` binary — so it is a layout/codegen
effect on Clarabel's factorisation, not a change in the work done. Consequences:
the n=1000 ratio of Table 1 is known to ±14 % (it reads "parity" either way);
from n=5000 on (131×, 716×) it is immaterial. Table 1 was measured with the
faster baseline binary, i.e. the conservative direction for the claim.

## Cost of the code-generation route

| n | `L` declared | `generate_code` | peak RSS | generated `src/` | `cargo build --release` | binary |
|---|---|---|---|---|---|---|
| 500 | dense | 0.5 s | 2.3 GB | 10.5 MB | 12 s | — |
| 1000 | dense | 16.7 s | 8.1 GB | 43.5 MB | 19 s (1.9 GB) | 33 MB |
| 500 | triangular | 0.2 s | 1.2 GB | 5.2 MB | — | — |
| 1000 | triangular | 2.3 s | 7.8 GB | 21.7 MB | 94 s (1.0 GB) | — |

Peak memory and source size grow as n² (CVXPY's DPP tensor holds the parameter
entries; the affine maps are emitted as static literals in `data.rs`): ~32 GB at
n=2000, ~200 GB at n=5000 for the dense declaration; the triangular declaration
halves the source but barely the peak memory (7.8 GB at n=1000), which sits in
CVXPY's canonicalisation of `Lᵀc` rather than in the generator. Not attempted above n=1000 on this machine. This is
the embedded design point of the CVXGEN lineage, not a defect; it does mean the
code-generation route is unavailable at the population sizes of Table 1's upper
rows, where the generated solver would in any case run the same O(n³)
factorisation (308 s at n=5000, 2298 s at n=10000).

## What this does and does not establish

- Established: the Table-1 baseline is a correctly assembled Clarabel cone
  program; an independent generator reproduces its optimum and per-iteration
  cost — out of the box once the factor's triangular pattern is declared. The
  active-set speed-ups of Table 1 are not an artefact of the baseline.
- Not established: anything about n ≥ 2000 through the generator (memory), or
  about idle-machine absolute times (see load caveat above).
- Found on the way: the same session exposed a 4.6× regression in the dense
  support-Gram path of `support_first.rs` (commit `e10e2f0`), fixed in the same
  change as this file — see `research/REVISION.md`.
