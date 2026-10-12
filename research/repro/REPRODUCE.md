# Reproducing the support-first results

The current automated pipeline is `repro.sh`. It produces fresh run evidence in
`artifacts/repro_runs/<UTC timestamp>-<pid>/`, including stage logs, a status table,
named panel exports, measured CSVs and figures. A completed stage copies its CSV
to `artifacts/`; stale CSVs for attempted stages are removed first. Dependent
stages read only the new run directory. An omitted or failed stage cannot consume
an old export from `/tmp` or an earlier run.

The manuscript's recorded measurements are separately archived under
`research/repro/reported/v0.4.0/`. Their SHA256 manifest and README describe the
source and its limitations. These are original summary CSVs, not a new benchmark
run and not proof that every original solver result passed today's checks.

Fresh corrected-source simulation evidence and reusable founder checkpoints from
11 October 2026 are provided in [`verified/2026-10-11/`](verified/2026-10-11/README.md).
That run reproduces diversity retention, while its final gain slightly favours
truncation; the historical late OCS gain advantage is not universal.

## Tables and figures

| Stage name | Script | Output |
|---|---|---|
| `table1` | Rust `compare --m 20000 --k-frac 0.017` | tighter-cap synthetic Table 1 logs (target support, no verified Ne_c); n=1000,2000,5000 |
| `cvxgenrust` | `cvxgenrust/run.sh 500 1000` | independent generated Clarabel baseline logs |
| `scaling_matrixfree` | `examples/scaling_matrixfree.rs` | measured dense-build curve on independent-locus panels at m=1000 |
| `operational_realpanels` | `r_operational_ne.R both` | wheat/mouse cap calibration and matrix-free/optiSel results |
| `operational_pig` | `r_operational_ne_pig.R` | PIC sexed-panel caps and results; exact dense panel exported |
| `exports` | `export_K_panels.R wheat mouse` | validation of the named benchmark exports |
| `algorithm_realpanels` | `sf_operational_ne.py wheat mouse` | verified dense prototype at those exact R caps |
| `algorithm_pig` | `sf_operational_ne_pig.py` | verified dense pig prototype at those exact R caps |
| `table2` | `combine_operational.py` | complete matched-cap `table2_operational.csv` |
| `structured` | `scaling_structured.R` | AlphaSimR m=1000, n=1000 through 40000, contribution Ne_c near 50 |
| `multigen` | `multigen_validation.R` | three replicates, two methods, twelve generations |
| `fig_scaling` | `research/fig_scaling.py` | measured build/structured solve and theoretical footprints |
| `fig_frontier` | `research/fig_frontier.py` | all three panels and three contribution-proxy targets |
| `fig_multigen` | `research/fig_multigen.py` | recorded gain and observed heterozygosity trajectories |

The independent-locus scaling example uses a historical loose cap for its own
solve. That solve is **not** substituted for the structured solve curve. The
scaling figure requires equal positive marker counts, the recorded target Ne_c=50
and achieved Ne_c near 50, and `Solved` status in fresh CSVs. The n=40000 dense build
is `not_measured`: it exceeds the example's n=30000 allocation policy. A legacy
`infeasible` cell is interpreted explicitly as the same policy omission. The
open marker/dashed segment is an n-squared projection from the last measured
build. The 8 n² and 8 nm memory curves are theoretical float64 array footprints;
they exclude solver workspace and say nothing universal about available RAM.

Contribution Ne_c=1/sum(c_i²) measures the concentration of a continuous allocation.
The displayed proxy 50/Ne_c percent is derived from 1/(2Ne_c). It is not an observed
realized genomic inbreeding rate. Observed heterozygosity in the multigeneration
experiment is a separate measurement. AlphaSimR panels retain family structure,
but re-centering markers per panel still drives the uniform-allocation kinship
floor toward zero; the scaling experiment does not establish a positive floor.

MaCS founder generation in the installed AlphaSimR 2.1.0 did not repeat the same
haplotypes after `set.seed`, even at one thread. The upstream maintainers document
this seed limitation and a later fix in
[their reproducibility discussion](https://github.com/gaynorr/AlphaSimR/discussions/219).
The installed `runMacs` has no explicit `seed` argument; it samples an internal
seed vector from R, which does not eliminate the observed limitation.

The simulations now default to `REPRO_SIM_THREADS=1`, passed to `runMacs` and
`SimParam$nThreads`. They persist founder populations and the post-generation R
RNG state in `artifacts/simulation_state/`; `REPRO_SIM_STATE_DIR` selects a separate
checkpoint directory. Replay validates generator parameters, seed, thread count
and AlphaSimR version, restores the stored RNG state, and records checkpoint
fingerprints in each run's `simulation_state_manifest.csv`. Parameter/version
mismatch fails rather than silently substituting a cached population. These are
input checkpoints, not reused benchmark results. Share the checkpoints alongside
run evidence to replay the same panels. De novo founder generation may still
differ across runs/installations; these checkpoints do not make MaCS generation
from a seed deterministic. A new checkpoint directory generates a new population.
The small repeated smoke test verifies identical founder/candidate genotypes,
criterion, sex, native contributions, support and gain on checkpoint replay.
No independent-locus generator is substituted for the coalescent family model.

The additive relationship quantity is G_eps=ZZ'/s+2e-5I. Native `quad`, its cap k,
and the structured CSV's legacy `coancestry` column refer to c'G_eps c. Biological
group coancestry uses K=G_eps/2 and is c'Kc=quad/2, with cap ub=k/2. The real-panel
operational CSV's `coancestry` column uses K. The historical structured column
excluded the tiny ridge term; its source convention is preserved in the archive.
The legacy `k_floor_uniform` is a feasible uniform-reference value, not proof of
an exact minimum or a positive persistent kinship floor.

Both operational R scripts use K=ZZ'/(2s)+1e-5I and native ridge 2e-5 with k=2ub.
Structured scaling now uses that same ridge convention, including it in both
calibration brackets. The historical structured cross-check used different
ridges and must not establish same-optimum exactness. Every fresh R solve in the
operational, structured and multigeneration stages requires `Solved` and primal
feasibility. The independent-locus Rust build-timing example checks `Solved`.
The dense prototype requires convergence,
primal/dual feasibility, stationarity, reduced costs and complementarity before
accepting timing output. Positive ridge makes the exported genomic kinship matrix
positive definite, the condition under which these KKT checks certify a global
optimum. The Python loader verifies this once by Cholesky outside the timed
solve. Symmetric indefinite inputs are rejected. Exhausting iterations is a failure.

R calibrates each cap once and exports `<panel>_K.csv`, `<panel>_bc.csv` and panel
metadata under the run's `inputs/`. MD5 fingerprints bind the criterion and sex
ordering plus the kinship matrix to each cap row. The Python prototype requires native Solved status and finite target/achieved Ne,
cap and gain, then reuses those exact caps and checks both fingerprints and the native gain. It does not
independently bisect a second cap. The combined fresh table requires all nine
operating points; missing pig data omits that table and frontier figure.

## Dependencies and data

Use a fresh source-built `ocsrs` R package (`R CMD INSTALL bindings/r`), rather
than an older package already on the library path. Rust must meet the crate's
`rust-version`; Python needs numpy, scipy and matplotlib. R needs AlphaSimR for
the simulations, BGLR for wheat/mouse, optiSel for dense baselines, and data.table
for PIC pig. The archived experiment used optiSel 2.1.0; that source is available
from the CRAN archive:

```r
install.packages(c("BGLR", "data.table", "AlphaSimR", "optiSolve", "ECOSolveR", "quadprog"))
install.packages("https://cran.r-project.org/src/contrib/Archive/optiSel/optiSel_2.1.0.tar.gz",
                 repos = NULL, type = "source")
```

Wheat and mouse data ship in BGLR. Mouse has recorded sex; wheat uses an arbitrary
balanced split. PIC data are a manual input: Cleveland, Hickey and Forni (2012),
*A common dataset for genomic analysis of livestock populations*, G3 2:429–435,
[DOI 10.1534/g3.111.001453](https://doi.org/10.1534/g3.111.001453), Supporting
Information `FileS1.zip`. Put `genotypes.txt`, `ebvs.txt` and `pedigree.txt` in
`/tmp/pig/FileS1/`, or set `PIG_DATA_DIR` to that extracted directory. Pedigree
membership resolves sex for 1194 of 3534 genotyped animals; these are animals that
became parents, a selected subset. No download or redistribution is automated.

## Run modes

```bash
bash research/repro/repro.sh
REPRO_STRICT=1 bash research/repro/repro.sh
REPRO_FULL=1 bash research/repro/repro.sh
REPRO_STEPS=scaling_matrixfree,structured,fig_scaling bash research/repro/repro.sh
REPRO_STEPS=operational_realpanels,exports,algorithm_realpanels bash research/repro/repro.sh
python3 research/repro/test_reproduction.py
Rscript research/repro/test_simulation_reproducibility.R
```

`REPRO_STEPS` selects a subset in fixed pipeline order; it does not reorder execution.
Empty or unknown stage names fail explicitly. For dependencies, include the producing
stages in the same run. Set `REPRO_PYTHON` to an isolated Python executable and
`R_LIBS_USER` to the freshly built R library when needed. Actual execution errors
return nonzero even in partial mode. Missing dependencies/data are separately
reported omissions; `REPRO_STRICT=1` returns2 for them. The intentional expensive
Clarabel n=10000 omission is included only with `REPRO_FULL=1`; its historical
tighter-cap run took about 38 minutes. AlphaMate is manual in every mode.

New figures remain in the run's `figures/` by default. Existing reported manuscript
figures are preserved on omissions and failures. Set `REPRO_UPDATE_FIGURES=1` to
copy a successfully produced figure into `research/`. Promotion checks every
copy, prepares both formats before replacement and restores prior files if a
replacement fails; copy/replacement errors are actual failures. Benchmark on an idle
machine, serially, before using new timings to update manuscript cells.

To redraw **historical** figures from the original summaries, without executing
or claiming new benchmarks:

```bash
python3 research/fig_scaling.py --reported
python3 research/fig_frontier.py --reported
python3 research/fig_multigen.py --reported
```

These commands validate the archived hashes and label the plots historical.
They do not fabricate missing original solver statuses or prove matched caps.
The previous seven-panel loose-cap timing notes are retained only in
`table2_loosecap_legacy.md`; they are not the current operational Table 2.

## AlphaMate manual provenance

AlphaMate solves discrete mate allocation. Its continuous contribution vector
can be rescored at matched coancestry as a relaxation check; a wall-time comparison
with the native continuous solver is not supported by the x86 emulated run.
The recorded execution used the Linux x86-64 binary under Rosetta-backed Colima.
The binary/source and generated contributor files are not included in the
recovered CSV bundle. Therefore the historical AlphaMate table cannot be rebuilt
from that bundle alone.

The documented recipe uses the full parent set, fewer matings than candidates,
a positively shifted criterion, seed 15 (the input generator default), the matrix
scaling convention below, and one output directory per target degree. The exact
old contributor file hashes, binary hash and complete specification are absent;
this is a rerun recipe, not certified recovery of the historical experiment.

```bash
colima start --vm-type vz --vz-rosetta --cpu 4 --memory 6
python3 alphamate/bench/make_alphamate_inputs.py --matrix G.txt --scale 2 \
    --ebv b_pos.txt --sex sex.txt --outdir run --matings 300 \
    --male-parents <nM> --female-parents <nF> --seed 15 --target-degree 45
docker run --rm --platform linux/amd64 --ulimit stack=-1 \
    -e OMP_NUM_THREADS=1 -e MKL_NUM_THREADS=1 \
    -v "$PWD/alphamate":/work -w /work/run debian:bookworm-slim \
    sh -c 'ulimit -s unlimited; exec /work/AlphaMate-src/binaries/AlphaMate_Unix AlphaMateSpec.txt'
python3 alphamate/bench/eval_contrib.py <n> run/Contributors*.txt
python3 research/repro/sf_at_ub.py <n> <rescored-cap> --legacy
```

`eval_contrib.py` uses the old `/tmp/bench_*_<n>` convention. For a manual rerun,
retain and hash the exact instance files, all input specifications, generated
contributor files, the original unshifted criterion for scoring, the binary and
environment versions. Those n-keyed scratch files are deliberately excluded
from the current automated operational pipeline.
