# Table 2 — reported operational measurements

Current manuscript operating point: contribution Ne_c=1/sum(c_i^2) near 50,
proxy ΔF=1/(2Ne_c) near 1%. This proxy is not observed genomic inbreeding.
Values are from the recovered historical `reported/v0.4.0/table2_operational.csv`,
protected by `manifest.json`; they were not rerun during the correction.

| Panel | n | m | support | Given-G prototype | Matrix-free R binding | optiSel | Given-G/optiSel speed-up | Matrix-free/optiSel speed-up |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| CIMMYT wheat |599|1279|109|0.0553 s|0.214 s|0.667 s|12×|3×|
| HS mouse, recorded sex |1814|10346|89|0.0605 s|1.455 s|8.056 s|133×|~6×|
| PIC pig, EBV and pedigree-derived sex |1194|52843|110|0.078 s|10.063 s|2.768 s|~35×|0.3×|

Ratios above are recomputed from these stored values before rounding; the original
combined CSV stores some rounded ratios as integers (mouse 133, pig 35). The
manuscript displays rounded timing cells. The matrix-free wheat/mouse numbers
are single timed runs; pig and given-G provenance describes medians over three
runs. The archive lacks all raw repetition logs and does not establish a median
for every displayed cell. The recorded machine description is Apple M4 Max,
macOS 25.5.0, R 4.6.0 with OpenBLAS 0.3.33, optiSel 2.1.0 and release Rust 1.95.0.
These are original experiment descriptions, not a reconstructed environment.

The historical given-G scripts bisected their own caps to target Ne; no matched
cap/fingerprint records establish identical instances with the matrix-free
column. The recovered wheat/mouse matrix-free caps and gain gaps are in
`reported/v0.4.0/operational_ne_realpanels.csv`: at target 50 the wheat gain gap
against optiSel is 0.046743 and mouse 0.0008653. Especially for wheat, claiming all
methods attained an identical optimum is unsupported.

Fresh reproduction now calibrates once in R, exports the exact named panel and
records K/criterion/sex fingerprints. The given-G scripts reuse those same caps,
validate primal/KKT convergence and compare gain before timing output is accepted.
`combine_operational.py` requires all three panels, all three target Ne points,
identical caps and fingerprints; it creates a separate fresh table in the run
artifact directory. Missing pig data produces a recorded omission rather than
a partial table or an empty legend. Timings from a loaded machine must not be
substituted for these historical reported cells without a new provenance note.

The earlier seven-panel loose-cap notes are archived in
`table2_loosecap_legacy.md`; their 12–180× matrix-free claims do not describe this
operating point. Figure 1's independent-locus build curve is also separate from
the structured contribution-Ne solve curve.
