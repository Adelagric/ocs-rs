# The real pig panel through the packed route

Reproduces the two numbers the manuscript's Methods quotes for the 52k-SNP PIC pig
panel: what the packed store costs against the dense `f64` matrix, and what rounding
this panel's imputed dosages to hard calls does to the optimum.

The panel is the Cleveland, Hickey & Forni (2012) common dataset
(DOI 10.1534/g3.111.001453); `FileS1` holds `genotypes.txt` (934 MB, imputed dosages).

```sh
cd research/repro/pig_to_bed
cargo run --release -- /path/to/FileS1 0.05
```

One pass over the file builds the raw-dosage matrix and writes a PLINK trio with the
rounded hard calls (what `plink --make-bed` would produce from them), then solves three
ways at one cap: raw dosages dense, rounded dense, rounded packed.

Measured 2026-09-30, Apple M4 Max, release profile (load average ~7, so timings are
indicative; the sizes and optima are exact):

```
PIC pig panel: n=3534 m=52843  read in 2.2s
  imputed dosages: 237563 of 186747162 values are non-integer (0.127%)
  cap k = 0.05 x mean diag G = 0.051679
  raw dosages, dense f64 :   1494.0 MB  |S|=23  gain 0.997436648
  rounded, dense f64     :   1494.0 MB  |S|=23  gain 0.997436469
  rounded, 2-bit packed  :     46.7 MB  |S|=23  gain 0.997436469
  packed vs dense-from-bed: support match true  Δgain 4.07e-14
  rounding cost (raw vs rounded): support 23 vs 23, shared 23, Δgain 1.79e-7
  (a dense n x n G for this panel: 99.9 MB)
```

So: 32× less memory for the same optimum, and the representation's own cost — rounding
0.127 % of the values — is 1.8e-7 in gain with the support unchanged. The selection
criterion here is a deterministic stand-in, not the panel's EBV: this compares two
representations of one panel, not two breeding recommendations.
