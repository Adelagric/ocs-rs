# ocs-rs — exact, matrix-free optimum contribution selection

Python bindings for the Rust **support-first** OCS solver.

Optimum contribution selection (OCS) chooses how much each candidate contributes to
the next generation: maximise genetic gain **bᵀc** subject to a quadratic
relationship bound, **cᵀGc ≤ k**, with non-negative contributions summing to
one — or, in the sexed form, to one half per sex.

The solver targets the global optimum of the convex OCS program and reports
`Solved` when its numerical termination checks pass. Scope: continuous
contributions under a single coancestry constraint; it does not do integer mate
allocation. See the repository's manuscript and validation results for the precise
optimality conditions and numerical tolerances.

**Matrix-free.** It never forms the dense n×n relationship matrix: the kinship
products come straight from the genotype matrix as `Gc = εc + Z(Zᵀc)/s`. A dense
float64 relationship matrix would occupy about 11.9 GiB at n = 40 000; this is
calculated storage, not an observed allocation failure. Working on a small active
support can reduce factorization costs. Small support is an empirical property of
the reported instances; individual contribution caps can force larger supports.

The matrix-free path costs O(nm) per product. The wrappers accept a dense genotype
matrix: Python converts `Z` to contiguous float64, copies it to bytes, and the
native layer copies those entries into its matrix layout. It therefore needs
additional O(nm) memory during a solve. Packed 2-bit PLINK solving is currently
available only through Rust's `read_packed` and `solve_sexed_packed` APIs.

## Install

```sh
pip install ocs-rs
```

Wheels ship for Linux (x86_64, aarch64), macOS (Intel, Apple silicon) and Windows,
`abi3` for Python 3.9+. To build from a **full repository checkout**, install Rust
1.85+ and activate a virtual environment, then run `pip install maturin` and
`maturin develop --release` from `bindings/python/`. The Cargo dependency on the
core crate is relative to that checkout. Current package version: 0.4.0.

## Quickstart

```python
import numpy as np
import ocs_rs

rng = np.random.default_rng(11)
n, m = 12, 20
dosages = rng.binomial(2, rng.uniform(0.1, 0.5, m), size=(n, m))
p = dosages.mean(axis=0) / 2
Z = dosages - 2 * p
s = ocs_rs.vanraden_scale(p)
b = rng.standard_normal(n)
male = np.arange(n) % 2 == 0
caps = np.full(n, 0.15)
ridge = 1e-5  # G = Z @ Z.T / s + ridge * I

# A feasible plan contributes half from each sex, uniformly within that sex.
c0 = np.where(male, 0.5 / male.sum(), 0.5 / (~male).sum())
assert np.all(c0 <= caps)
q0 = float(np.sum((Z.T @ c0) ** 2) / s + ridge * (c0 @ c0))
k = q0 + 0.05  # This cap admits c0; it is not a universal breeding threshold.

res = ocs_rs.solve(Z, b, k=k, s=s, ridge=ridge, male=male, caps=caps)
if res.status != "Solved":
    raise RuntimeError(f"No solved optimum: {res.status}")
print(res)
print(res.support)  # 0-based NumPy indices
print(res.c[res.support])
```

The same example is executable as `python examples/quickstart.py` from this
directory. The positive ridge contributes `ridge * sum(c**2)` to the quadratic relationship
constraint and must also be included when computing a feasible cap.

## API

`solve(Z, b, k, *, s, ridge=1e-5, male=None, caps=None, max_iter=10_000, tol=1e-9)`

| argument | meaning |
|---|---|
| `Z` | (n, m) centred genotype matrix — `G = ZZᵀ/s + ridge·I`, never formed |
| `b` | (n,) breeding values (the selection criterion) |
| `k` | quadratic relationship bound: `cᵀGc ≤ k` |
| `s` | VanRaden scale `2 Σ p(1-p)` |
| `male` | (n,) bool — give it for the sexed form (`Σ_males = Σ_females = ½`) |
| `caps` | (n,) per-candidate upper bounds `c ≤ caps` |

Returns an `OcsResult` with `c`, `support` (0-based indices), `gain`, `quad`,
`iterations`, `products` and `status`. Check `status == "Solved"` before using the
contributions as an optimum; `MaxIter` has no such guarantee.

The API reports `quad = cᵀGc`. Here `G` is an additive relationship matrix;
group coancestry is `quad / 2` using `K = G / 2`. A group-coancestry bound `ub`
therefore corresponds to `k = 2 * ub`. This score does not specify realised
offspring inbreeding without a mating allocation and reference population.

The API accepts **`Z`, not a relationship matrix `G`**. Passing an n×n `G` as `Z`
would solve with `G @ G.T / s + ridge * I`. An arbitrary relationship matrix or
one modified by blending, shrinkage or pedigree adjustment needs a separate
adapter/factorization that reproduces the intended operator; the wrappers do not
provide a general `G` argument.

`b`, `male` and `caps` must have exactly the row order of `Z`. For a PLINK panel,
align them to `panel.ids` before solving. These IDs are the `.fam` **IID column
only**; family IDs (FID) and `.fam` sex are not returned. Preserve that metadata
separately, check that your identifiers are unique, and supply the sex mask
yourself. The reader checks file dimensions and encoding consistency; matching
dimensions do not verify that `.bed`, `.bim` and `.fam` describe the same identities
or ordering. Use a trio produced together, not independently reordered files.

## Links

- Source, manuscript and reproduction scripts: <https://github.com/Adelagric/ocs-rs>
- Archived release: <https://doi.org/10.5281/zenodo.20746987>

MIT licence.
