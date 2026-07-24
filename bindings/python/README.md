# ocs-rs — exact, matrix-free optimum contribution selection

Python bindings for the Rust **support-first** OCS solver.

Optimum contribution selection (OCS) chooses how much each candidate contributes to
the next generation: maximise genetic gain **bᵀc** subject to a cap on the mean
coancestry of the offspring, **cᵀGc ≤ k**, with non-negative contributions summing to
one — or, in the sexed form, to one half per sex.

**Exact, not heuristic.** The solver returns the global optimum of the (convex) OCS
program, certified by the KKT conditions at termination — verified against optiSel, a
conic interior-point solver (agreement 1e-8) and an independent numerical reference
(1.5e-14). Scope: continuous contributions under a single coancestry constraint; it
does not do integer mate allocation.

**Matrix-free.** It never forms the dense n×n relationship matrix: the kinship
products come straight from the genotype matrix as `Gc = εc + Z(Zᵀc)/s`. That is what
lets it solve instances whose dense matrix (≈12 GiB at n = 40 000) will not fit in
memory, and why the cost follows the (tiny) active support and the marker count
rather than n².

On real marker panels it reaches the same optimum as optiSel 90–2280× faster.

## Install

```sh
pip install maturin
maturin develop --release        # from bindings/python/ in the repository
```

(A PyPI release will follow; for now build from source.)

## Quickstart

```python
import numpy as np
import ocs_rs

# Z : (n, m) genotypes centred by twice the allele frequencies
# s : VanRaden scale 2 Σ p(1-p)  —  ocs_rs.vanraden_scale(p) computes it
res = ocs_rs.solve(Z, b, k=0.03, s=s, male=male)

print(res)                # OcsResult(status='Solved', |support|=17, gain=..., ...)
print(res.support)        # the handful of selected candidates
print(res.c[res.support]) # their contributions
```

## API

`solve(Z, b, k, *, s, ridge=1e-5, male=None, caps=None, max_iter=10_000, tol=1e-9)`

| argument | meaning |
|---|---|
| `Z` | (n, m) centred genotype matrix — `G = ZZᵀ/s + ridge·I`, never formed |
| `b` | (n,) breeding values (the selection criterion) |
| `k` | coancestry cap: `cᵀGc ≤ k` |
| `s` | VanRaden scale `2 Σ p(1-p)` |
| `male` | (n,) bool — give it for the sexed form (`Σ_males = Σ_females = ½`) |
| `caps` | (n,) per-candidate upper bounds `c ≤ caps` |

Returns an `OcsResult` with `c`, `support`, `gain`, `quad`, `iterations`, `products`
and `status` (`"Solved"` means KKT-optimal; only then are `c`/`gain` usable).

## Links

- Source, manuscript and reproduction scripts: <https://github.com/Adelagric/ocs-rs>
- Archived release: <https://doi.org/10.5281/zenodo.20746987>

MIT licence.
