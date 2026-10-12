# ocs-rs — support-first optimum contribution selection

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20746987.svg)](https://doi.org/10.5281/zenodo.20746987)

An **exact, matrix-free** solver for genomic Optimum Contribution Selection (OCS).
Support-first finds the active contributions by column generation, solves each
fixed support in closed form, and evaluates relationship products without forming
the dense `n×n` matrix. Small supports make this effective on the instances studied;
sparsity is an empirical property, and individual contribution caps can require
many selected candidates. Optimality is checked through KKT conditions within the
solver's numerical tolerances. Comparisons cover public genomic panels and
synthetic populations, with separate results for a given-matrix NumPy prototype
and the distributed Rust matrix-free solver.

> Full write-up — derivation, tables, figures, references — [`research/MANUSCRIPT.tex`](research/MANUSCRIPT.tex)
> (French: [`MANUSCRIPT_fr.tex`](research/MANUSCRIPT_fr.tex)); corrected compiled PDFs:
> [English](research/MANUSCRIPT.pdf) · [French](research/MANUSCRIPT_fr.pdf).
> Reproduction stages, dependencies and historical measurement provenance:
> [`bash research/repro/repro.sh`](research/repro/REPRODUCE.md).

## Results

On real marker panels (CIMMYT wheat n=599, PIC pig n=3534, heterogeneous-stock
mouse n=1814 with recorded sex), including a pig sub-panel of 1194 animals with
EBVs and sex inferred from the pedigree:

- **Numerical optimality checks.** Agrees with the matching Clarabel cone
  program to `1e-8` on tested synthetic instances and with an independent NumPy
  reference to `1.5e-14` on the reported binding instance. Historical optiSel
  comparisons do not always use the requested cap equally, especially on wheat;
  agreement at matched realised coancestry is checked separately. Timing ratios
  alone do not certify identical optima.
- **vs optiSel**, at a contribution-based operating proxy
  `Nₑ = 1/Σcᵢ² = 50`, `ΔF ≈ 1/(2Nₑ) = 1 %`: two distinct results.
  This proxy is not a measured genomic inbreeding rate. **Given `G`**, the NumPy
  active-set prototype is **12–133×** faster than optiSel's interior-point solve — it prices in `O(n²)`,
  independent of the marker count. The **shipped matrix-free** solver (forming no `G`,
  `O(nm)` per iteration) is faster only up to moderate marker density — 3× on wheat,
  ~6× on mouse — and is **slower** than optiSel on the 52k-SNP pig, where forming `G`
  once wins. The matrix-free path is the memory enabler at scale, not a universal
  speed-up (two-column [`Table 2`](research/repro/table2_numbers.md)).
  These historical timing columns were calibrated separately to approximately
  matched contribution proxies. They are not timings at identical exported caps.
  In particular, wheat gain differences reflect different constraint utilization;
  the ratios alone do not establish equal numerical accuracy.
- **vs Clarabel** (a generic conic solver): the tighter-cap synthetic comparison reports
  **716×** at n=10000 (2298 s → 3.21 s). A separate historical loose-cap experiment
  with support 2 reports **76179×** (about 30 minutes → 24 ms). Both compare Rust
  against Rust in one process —
  see [`Table 1`](research/repro/table1_numbers.md).
- **vs AlphaMate** (a heuristic for the *distinct* discrete-mating problem): on the
  continuous relaxation the two share, scored at matched coancestry, the exact
  optimum is no worse on the reported points. This checks the shared relaxation;
  AlphaMate's discrete mate-allocation task is distinct. Its emulated run time is
  not used for a timing comparison.
- **Scaling.** The reported structured synthetic sweep fixes 1000 markers and
  the contribution proxy `Nₑ = 50`. At 40000 candidates, the support remains near
  90 and the measured matrix-free solve takes approximately **0.65 s**. A dense
  double-precision `G` alone would occupy **11.9 GiB**, compared with **0.30 GiB**
  for dense `Z`; these are calculated array footprints. Dense construction was
  deliberately omitted above 30000 candidates, so the experiment establishes
  reduced storage requirements rather than a demonstrated out-of-memory failure
  (Figure 1: [`research/fig_scaling.pdf`](research/fig_scaling.pdf)).

## The model

OCS picks how much each of `n` candidates contributes to the next generation:
maximise genetic gain subject to a cap on average coancestry (to preserve
diversity).

```
maximize    bᵀc                  b = genomic estimated breeding values (GEBV)
subject to  A c = d              budget: Σcᵢ = 1, or sexed Σ_males = Σ_females = ½
            0 ≤ c ≤ u            optional per-candidate contribution cap (u = 1 ⇒ off)
            cᵀ G c ≤ k           G = VanRaden additive relationship matrix, k = quadratic score bound
```

`G = ZZᵀ/s + εI` is symmetric positive definite for `ε > 0`; `cᵀGc ≤ k` is a
convex quadratic constraint (a second-order cone). With this additive relationship
convention, group coancestry is `cᵀ(G/2)c`: the API returns `quad = cᵀGc` and
uses twice the group-coancestry bound as `k`. Offspring inbreeding also depends
on mating allocation and the reference population. When parents are separated
by sex, the *sexed* form replaces the simplex with two equality rows, each sex
contributing half. The simplex form remains available for other selection models.

## How it works

- **Active-set column generation.** Seed the support with the best male and
  female; price candidates by reduced cost and add the best, toward a *single*
  fixed quadratic relationship cap `k`; drop negatives. A feasible point with no positive
  reduced cost is KKT-optimal — exact, not heuristic.
- **Closed form per support.** Eliminate the equality constraints through a
  `q×q` reduction `P = A_S G_S⁻¹ A_Sᵀ` (q ≤ 2) and read the multiplier off a
  scalar quadratic — one Cholesky, no inner iteration.
- **Matrix-free.** `Gc = εc + Z(Zᵀc)/s`; the dense `n×n` `G` is never formed.
- **Per-candidate caps `c ≤ u`.** The same active set absorbs upper bounds with no
  loss of the closed form — a candidate at its cap moves to an *upper* set and
  enters the restricted solve as a constant offset; verified against Clarabel on
  the matching `c ≤ u` cone program.

Derivation in [`research/MANUSCRIPT.tex`](research/MANUSCRIPT.tex) (Methods);
implementation in [`src/support_first.rs`](src/support_first.rs) (`solve` /
`solve_sexed`, and `solve_capped` / `solve_sexed_capped` for `c ≤ u`).

## Running

```sh
# support-first head-to-head and benchmarks
cargo run --release -- compare --n 5000          # Clarabel vs support-first, same optimum
cargo run --release --example bench_sexed         # sexed solver: release timing + cross-language export
cargo run --release --example scaling_matrixfree  # dense-build timing and storage data
# Figure 1's operational solve/support data require AlphaSimR and the R binding:
Rscript research/repro/scaling_structured.R
python3 research/fig_scaling.py

# the original Clarabel evaluation harness
cargo run --release                               # = `all`: correctness + frontier + scaling -> REPORT.md
cargo run --release -- scale --n 5000 --route a   # one timing point (CSV row on stdout)
```

`stdout` is data; progress goes to `stderr`. Flags:
`--n --m --seed --k-frac --route a|b --max-iter --time-limit --max-n`.

## Genotypes in

Real panels are stored as PLINK binaries, so the solver reads them directly — no
conversion step, and the dense `G` still never exists:

```rust
let panel = ocs_rs::plink::read_panel(Path::new("cattle"))?;   // cattle.bed/.bim/.fam
```

```python
panel = ocs_rs.read_plink("cattle")  # requires your cattle.bed/.bim/.fam files
# b must already be aligned to panel.ids; the standalone example below needs no files.
res = ocs_rs.solve(panel.Z, b, k=0.03, s=panel.s)
if res.status != "Solved":
    raise RuntimeError(res.status)
```

Genotypes are centred by twice the allele frequency estimated from the non-missing
calls (missing ones imputed to the marker mean), so `panel.Z` and `panel.s` feed
`solve` as they are. SNP-major `.bed` only — PLINK ≥ 1.9's default; the trio is
checked for compatible dimensions. This does not detect an unrelated trio with
the same dimensions. `panel.ids` contains the `.fam` IID field in file order;
family IDs, sex and phenotypes are not returned. Align breeding values, sex flags
and caps explicitly to that order, and check IID uniqueness before joining.
VCF is not read: convert with `plink --vcf`.

### Packed storage through the Rust API

`.bed` stores genotypes two bits each, four per byte, column-major — the layout the
solver's packed store uses. `read_packed` therefore remaps the file's bytes and
solves from those, expanding nothing:

```rust
let panel = ocs_rs::plink::read_packed(Path::new("cattle"))?;  // 2 bits per genotype
let out = solve_sexed_packed(&panel.geno, 1e-5, &b, &male, k, 10_000, 1e-9);
```

On a 5000 × 50000 panel (`cargo run --release --example plink_packed`), peak resident
memory is **69 MB against 2.02 GB** for the dense route, on a 62.5 MB file — the same
optimum, to 7e-14 in gain, missing calls included (a missing genotype takes the fourth
2-bit code and decodes to the marker mean exactly). The dense `f64` route stays the
faster one per solve; packing reduces the genotype storage requirement. These
packed entry points are currently Rust-only. The R and Python interfaces read
and solve with a dense centered `Z`, whose storage remains `O(nm)`.

## Python

The people who run OCS work in R and Python, not Rust — so the solver is callable
from Python, numpy in and numpy out:

```sh
pip install ocs-rs
```

```python
import numpy as np
import ocs_rs

rng = np.random.default_rng(7)
M = rng.binomial(2, 0.5, size=(16, 12))
p = M.mean(axis=0) / 2
Z = M - 2 * p
s = ocs_rs.vanraden_scale(p)
b = rng.normal(size=16)           # illustrative selection values
male = np.arange(16) % 2 == 0
ridge = 1e-5
c0 = np.full(16, 1 / 16)         # feasible half-per-sex plan
k = np.sum((Z.T @ c0) ** 2) / s + ridge * (c0 @ c0) + 0.1
res = ocs_rs.solve(Z, b, k=k, s=s, ridge=ridge, male=male)
if res.status != "Solved":
    raise RuntimeError(res.status)
print(res.support, res.gain)      # Python support indices are 0-based
```

Wheels are published for Linux (x86_64, aarch64), macOS (Intel, Apple silicon) and
Windows; they are `abi3`, so one wheel per platform serves Python 3.9 onward. To build
from a full checkout instead, use
`pip install maturin` then `cd bindings/python && maturin develop --release`.

Details in [`bindings/python/README.md`](bindings/python/README.md). The bindings sit in
their own crate, so the core keeps the pure-Rust dependency surface above. Correctness
is checked against an independent optimiser (SciPy SLSQP), not just against the Rust
solver. [`.github/workflows/wheels.yml`](.github/workflows/wheels.yml) builds the
wheel targets and tests installed native wheels on Linux x86_64, Apple silicon and
Windows, with Python 3.9 and the latest supported Python. Linux aarch64 and Intel
macOS have build coverage only. Publication on a `v*` tag is gated on those tests
and the core MSRV/stable checks, using PyPI Trusted Publishing.

## R

The domain's own tool is an R package, so the solver is an R package too:

```sh
R CMD INSTALL bindings/r   # run from the full repository checkout
```

```r
library(ocsrs)
set.seed(7)
M <- matrix(rbinom(16 * 12, 2, 0.5), nrow = 16)
p <- colMeans(M) / 2
Z <- sweep(M, 2, 2 * p)
s <- vanraden_scale(p)
b <- rnorm(16)                  # illustrative selection values
male <- seq_len(16) %% 2 == 1
ridge <- 1e-5
c0 <- rep(1 / 16, 16)           # feasible half-per-sex plan
k <- sum(crossprod(Z, c0)^2) / s + ridge * sum(c0^2) + 0.1
res <- ocs_solve(Z, b, k = k, s = s, ridge = ridge, male = male)
stopifnot(res$status == "Solved")
print(res$support)              # R support indices are 1-based
```

Source builds require Rust 1.85 or later; the R package requires R 4.2 or later.
The R Rust crate references the repository core by a relative path, so keep the
full checkout when installing. Both interfaces accept centered `Z` and its scale,
not an arbitrary supplied `G`; a modified relationship matrix needs a compatible
factorization and a separately validated adapter.

R stores matrices column-major and so does the solver's linear algebra, so the
double-precision genotype matrix is borrowed: the Rust side views R's buffer in place
(1.5 GB left where it is, on a 52k-marker panel). Supports come back 1-based — that
and the layout are the two conventions that would otherwise fail silently, and both
are covered by tests that recompute the answer with R's own arithmetic. Integer
matrices are converted to doubles. Python's wrapper copies `Z` into native dense
storage; its memory use differs from both this borrowed R path and Rust packing.
[`research/repro/r_binding_optisel.R`](research/repro/r_binding_optisel.R) runs
support-first and optiSel in one session on one instance with the conventions
aligned exactly (optiSel constrains `c' sKin c ≤ ub` on `sKin = G/2 + εI`, hence
`k = 2·ub` and `ridge = 2ε`), and prints gain, realised coancestry, support and time
for both — reproducible rather than asserted. The FFI lives in
[`bindings/r`](bindings/r), never in the core.

## Reproduction

```sh
bash research/repro/repro.sh
```

Runs the documented stages available in the local environment and reports missing
dependencies or datasets as omissions. Executed-stage failures produce a nonzero
exit status; a partial run does not reproduce all manuscript results. Operational
comparisons require `BGLR`, `optiSel` and the freshly built R binding; structured
scaling and multigeneration simulations also require `AlphaSimR`.
The PIC pig panel needs a
one-time manual download (URL + layout in
[`research/repro/REPRODUCE.md`](research/repro/REPRODUCE.md)); AlphaMate runs from
a Linux binary under Colima + Rosetta (recipe in the same file).

## Verification

```sh
cargo test                                 # KKT certificates, feasibility, monotonicity, CSC properties
cargo clippy --all-targets -- -D warnings
cargo fmt --check
```

A `Solved` result reports convergence within the configured tolerances. Check the
status before using contributions; `MaxIter` does not certify an optimum or
feasibility. Tests verify budget, sex totals, individual caps and kinship across a
range of instances. The optimum is cross-checked
independently of this crate's arithmetic: against Clarabel (conic IPM), a NumPy
reference (`1.5e-14`), SciPy SLSQP, and optiSel on real data. The Clarabel
baseline itself is reconstructed by a third party's code generator, cvxgenrust,
from a CVXPY model of the cone program — same optimum, same per-iteration cost
([`research/repro/cvxgenrust/RESULTS.md`](research/repro/cvxgenrust/RESULTS.md)).

## Origin — the Clarabel spike

This began as a go/no-go on whether the [Clarabel](https://clarabel.org) conic
interior-point solver handles genomic-scale OCS (verdict: **GO**, see
[`REPORT.md`](REPORT.md), [`BENCHES.md`](BENCHES.md), [`DECISIONS.md`](DECISIONS.md)).
That spike exposed the opening: the conic IPM pays `O(n³)` per iteration to
describe a solution that activates a handful of candidates. support-first
exploits exactly that. Clarabel is kept as an independent cross-check oracle.

## Honest caveats

Wheat and mouse use recorded phenotypes as selection criteria. The pig panels use
real EBVs; the 1194-animal sub-panel combines EBVs with sex inferred from the
pedigree, but includes animals that became parents and is a post-selection subset.
Mouse sex is recorded, while wheat and the full pig panel use arbitrary balanced
splits. A GBLUP-derived criterion on the pig sub-panel is an additional within-panel
check, not independent predictive validation. The given-`G` optiSel comparison
times a NumPy prototype against R; the Rust matrix-free comparison is separate.
Small support is empirical: the no-ridge theorem guarantees existence of an
optimum with at most `q+r+1` nonzeros only without restrictive individual caps.
With caps it bounds free coordinates, and many positive capped coordinates may
be necessary. No universal bound depending only on `q` and the un-ridged rank `r`
survives a positive ridge
([`research/support_bound_sketch.md`](research/support_bound_sketch.md)).
The solver handles a single quadratic constraint and continuous contributions;
per-candidate caps `c ≤ u` are supported, but multiple quadratic constraints and
integer mate allocation are not. Large-candidate benchmarks are synthetic;
these results do not establish performance inside a production breeding evaluation.
These limits are stated in the manuscript's Discussion.

## Archive versions

The DOI badge identifies all versions and resolves to the latest published Zenodo
record. A repository edit does not replace archived software or PDFs. See
[`docs/ARCHIVING.md`](docs/ARCHIVING.md) for version-specific citations and the
process for archiving substantive corrections.

## Credits

The results stand on other people's work, cited in the manuscript and named here
so the credit is visible from the README alone.

- **cvxgenrust** — Hao Zhu & Joschka Boedecker (University of Freiburg),
  [github.com/dxogrp/cvxgenrust](https://github.com/dxogrp/cvxgenrust),
  [arXiv:2609.13875](https://doi.org/10.48550/arXiv.2609.13875). Used to
  independently reconstruct the Clarabel baseline from the CVXPY formulation;
  what we learned doing so is reported upstream as
  [dxogrp/cvxgenrust#1](https://github.com/dxogrp/cvxgenrust/issues/1), addressed
  in cvxgenrust 0.1.1.
- **Clarabel** — Paul Goulart & Yuwen Chen, [clarabel.org](https://clarabel.org).
  The conic interior-point oracle used in independent numerical validation.
- **faer** — Sarah Quiñones, [github.com/sarah-quinones/faer-rs](https://github.com/sarah-quinones/faer-rs).
  The pure-Rust dense linear algebra underneath.
- **optiSel** — Robin Wellmann; **AlphaMate** — Gregor Gorjanc & John Hickey;
  **AlphaSimR** — Chris Gaynor et al. The domain baselines and the simulator of
  the multi-generation validation.
- **Data** — the wheat and mouse panels ship with **BGLR** (Paulino Pérez &
  Gustavo de los Campos); the pig panel is the PIC common dataset (Cleveland,
  Hickey & Forni 2012).

## Constraints honoured

Pure-Rust dependencies only (`clarabel`, `faer`, `rand`, `rand_distr`) — no
system BLAS/LAPACK, no FFI, no `unsafe`. Versions pinned exactly. Single-author.
