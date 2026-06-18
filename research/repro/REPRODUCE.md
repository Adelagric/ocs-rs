# Reproducing the support-first results

Every number, table, and figure in the manuscript is regenerated from this
directory plus the Rust crate. The automatable pipeline is driven by `repro.sh`;
two inputs need a one-time manual step (the PIC pig download, and — optionally —
AlphaMate), documented below.

## What lives where

| Path | Produces |
|---|---|
| `examples/bench_sexed.rs` (crate) | Rust release timing of `solve_sexed` + cross-language export |
| `examples/scaling_matrixfree.rs` (crate) | the scaling sweep (Figure 1 data: support, solve, G build, G/Z memory) |
| `research/fig_scaling.py` | **Figure 1** (`research/fig_scaling.{pdf,png}`) |
| `research/support_first.py` | pool-unique (simplex) support-first prototype |
| `research/support_first_sex.py` | sexed support-first prototype (the NumPy oracle) |
| `research/repro/sf_at_ub.py` | sexed support-first at an arbitrary kinship cap (frontier points) |
| `research/repro/{wheat,mouse,pig}_export.R` | VanRaden GRM + EBV + sex → `/tmp/bench_{K,bc,ub}_{n}.{csv,txt}` |
| `research/repro/synthetic_optisel.R` | structured-population synthetic optiSel benchmark |
| `research/optisel_benchmark.R` | optiSel timing harness |
| `alphamate/bench/make_alphamate_inputs.py` | OCS instance → AlphaMate input files |
| `alphamate/bench/eval_contrib.py` | score an AlphaMate contribution vector in our metric |

Export scripts use `/tmp` as scratch: each writes `bench_K_{n}.csv` (the kinship
matrix), `bench_bc_{n}.csv` (`bv,oc,sex`), and `bench_ub_{n}.txt`; the Python
benchmarks read those by `n`.

## Dependencies

- **Rust** (stable) — the crate; `cargo run --release` builds the examples.
- **Python 3** with `numpy`, `scipy`, `matplotlib`.
- **R** with `BGLR` (wheat + mouse data) and `optiSel` (the exact baseline).
  `install.packages(c("BGLR","optiSel","data.table"))`.
- *Optional:* Docker + Colima (Apple Silicon: `colima start --vm-type vz
  --vz-rosetta`) for AlphaMate.

## Datasets

- **Wheat** and **mouse** ship inside `BGLR` (`data(wheat)`, `data(mice)`) — fetched
  automatically by the export scripts. The mouse panel carries a real recorded sex
  (934 males, 880 females); wheat is given an arbitrary balanced split.
- **PIC pig** — a one-time manual download (not redistributable here):
  - Cleveland, M.A., Hickey, J.M. & Forni, S. (2012). *A common dataset for
    genomic analysis of livestock populations.* G3 2(4):429–435.
    **DOI 10.1534/g3.111.001453.**
  - Obtain `FileS1.zip` from the article's **Supporting Information** (G3 / GSA
    figshare): <https://doi.org/10.1534/g3.111.001453> → *Supporting Information*.
  - `FileS1.zip` contains `genotypes.txt` (3534 individuals × 52,843 SNP),
    `ebvs.txt` (EBVs; trait 3 is used as **b**), and `phenotypes.txt`.
  - Unzip so the files sit at **`/tmp/pig/FileS1/`** (i.e. `/tmp/pig/FileS1/genotypes.txt`),
    matching `pig_export.R`. The PIC panel ships no usable sex (sex chromosomes
    removed), so an arbitrary balanced split is used.

## Run it

```bash
# from the crate root
bash research/repro/repro.sh
```

`repro.sh` checks each toolchain, runs everything it can, and skips (with a clear
message) any step whose dependency or dataset is absent — so a partial environment
still reproduces the parts it can. The PIC pig and AlphaMate steps run only if
their inputs are present.

## AlphaMate (optional)

No macOS build exists and the source is locked to an Intel toolchain; the Linux
x86-64 binary runs under Rosetta-backed Colima:

```bash
colima start --vm-type vz --vz-rosetta --cpu 4 --memory 6
# convert an instance, then (mount under $HOME, not /tmp):
python3 alphamate/bench/make_alphamate_inputs.py --matrix G.txt --scale 2 \
    --ebv b_pos.txt --sex sex.txt --outdir run --matings 300 \
    --male-parents <nM> --female-parents <nF> --target-degree 45
docker run --rm --platform linux/amd64 --ulimit stack=-1 \
    -e OMP_NUM_THREADS=1 -e MKL_NUM_THREADS=1 \
    -v "$PWD/alphamate":/work -w /work/run debian:bookworm-slim \
    sh -c 'ulimit -s unlimited; exec /work/AlphaMate-src/binaries/AlphaMate_Unix AlphaMateSpec.txt'
python3 alphamate/bench/eval_contrib.py <n> run/Contributors*.txt
```

Note the three input requirements found the hard way: matings < n; the **full**
parent set (a reduced count triggers a setup segfault); and a **positively shifted**
selection criterion (the value/max objective inverts on negative EBVs). The
criterion shift is argmax-invariant, so it does not change the optimum.
