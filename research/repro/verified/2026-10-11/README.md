# Verification run — 11 October 2026

This is fresh verification evidence from the corrected scripts, separate from the
historical v0.4.0 summaries used in the figures. Both selected stages completed:
structured scaling (1,000 to 40,000 candidates, 1,000 markers) and three replicates
of twelve-generation OCS/truncation selection. Every R solve passed the status
and primal acceptance checks. Timings were obtained during other verification
work and are not replacement headline benchmarks.

At n=40,000, the recorded support is 83 and the contribution participation ratio
is approximately 49.7. In the multigeneration run, generation-12 mean gain is
11.46 for OCS and 11.51 for truncation; observed heterozygosity is 0.234 and 0.213,
respectively. Diversity retention is reproduced; a late gain advantage is not.

The input checkpoints contain founder populations and post-generation R RNG
states, not solver results. The scripts verify AlphaSimR version and generation
parameters before reuse. They were generated using AlphaSimR 2.1.0, R 4.6.0 and
one simulation thread, with the corrected ocsrs 0.4.0 local source build.
Use the same dependency versions and run from the repository root:

```sh
REPRO_SIM_STATE_DIR=research/repro/verified/2026-10-11/simulation_state \
REPRO_SIM_THREADS=1 REPRO_STEPS=structured,multigen \
bash research/repro/repro.sh
```

The native and BLAS runtime can affect timing and floating-point details. Stored
founders address the observed MaCS seed limitation; a seed alone does not promise
identical de novo founders. `manifest.json` records SHA256 hashes for this evidence.
The complete local logs remain in the corresponding run directory; only portable
summary evidence and input checkpoints are included here.
