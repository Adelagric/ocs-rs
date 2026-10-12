# Reported v0.4.0 measurements

These five summary CSVs were recovered byte-for-byte from the local experiment
checkout at commit e6790c99. The SHA256 manifest detects accidental edits. They
are historical reported measurements, not measurements obtained during the
correction. Original raw solver logs, environment snapshots and AlphaMate
contributor files were not present in this summary bundle.

`table2_operational.csv` combines the three panels at contribution targets
Ne_c=25,50,100. Its given-G prototype calibrated its own cap independently; the
bundle does not establish equal caps across the algorithm and matrix-free
columns. `operational_ne_realpanels.csv` records the matrix-free caps and optiSel
gain gaps, including wheat gaps of 0.0135, 0.0467 and 0.2482. These are not proof
that optiSel reached the same optimum. The pig's achieved contribution Ne_c at
target 50 is rounded to 49 in the combined summary.

`scaling_structured.csv` records m=1000 and achieved contribution Ne_c near 50,
including 0.646 s at n=40000. Its producing script used solver ridge 1e-5 whereas
the optiSel cross-check used G/2+1e-5I; those two formulations differed by a factor
two in the ridge. This archive does not repair or revalidate that comparison.
The historical CSV has no solver-status field. Explicit `--reported` plotting
checks the archived schema and hash while labelling the source historical;
fresh plotting requires recorded target and Solved status.

`scaling_matrixfree.csv` used independent-locus panels and a loose fixed cap.
Only its dense-build timing curve is used in the scaling plot. The legacy
`infeasible` at n=40000 was an allocation-policy omission above n=30000. The
11.92 GiB G and 0.298 GiB Z values are theoretical float64 array footprints,
excluding workspace, not measured peak memory or a universal machine limit.

`multigen_validation.csv` contains 78 population summaries (three replicates,
two methods, thirteen generations). Observed heterozygosity is measured; the
selection target50/Ne_c is a contribution-based proxy, not observed genomic
inbreeding. Regeneration from this archive reproduces the recorded trajectories
without asserting that a new simulation was run.
