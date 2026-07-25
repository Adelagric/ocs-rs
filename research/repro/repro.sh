#!/usr/bin/env bash
# Reproduce the support-first results. Run from anywhere; it cds to the crate
# root. Each step is guarded: a missing toolchain or dataset is reported and
# skipped, so a partial environment still reproduces the parts it can.
#   bash research/repro/repro.sh
#
# Steps are labelled with the table or figure they produce, so a reader can check
# any single number. The one long step (Clarabel at n=10000, ~26 min) is behind
# REPRO_FULL=1; everything else runs in minutes:
#   REPRO_FULL=1 bash research/repro/repro.sh
#
# Timings in the paper were taken on an idle machine, serially. Running anything
# else concurrently inflates the interior-point baselines severalfold.
set -u
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
cd "$ROOT"

have() { command -v "$1" >/dev/null 2>&1; }
pyok() { python3 -c "import $1" >/dev/null 2>&1; }
rok() { Rscript -e "suppressMessages(library($1))" >/dev/null 2>&1; }
banner() { printf '\n=== %s ===\n' "$1"; }

banner "1/7  Rust release timings + scaling sweep (Figure 1 data)"
if have cargo; then
  cargo run --release --example bench_sexed
  cargo run --release --example scaling_matrixfree
else
  echo "SKIP: cargo not found"
fi

banner "2/7  Figure 1"
if have python3 && pyok matplotlib; then
  python3 research/fig_scaling.py
else
  echo "SKIP: python3 + matplotlib required"
fi

banner "3/7  Table 1 — support-first vs Clarabel (Rust vs Rust, one process)"
if have cargo; then
  ns="1000 2000 5000"
  if [ "${REPRO_FULL:-0}" = "1" ]; then
    ns="$ns 10000"
  else
    echo "  (n=10000 skipped: Clarabel takes ~26 min there; set REPRO_FULL=1 to include it)"
  fi
  for n in $ns; do
    cargo run --release --quiet -- compare --n "$n" || echo "  compare n=$n failed"
  done
else
  echo "SKIP: cargo not found"
fi

banner "4/7  Real-data GRM export (wheat, mouse via BGLR; pig if downloaded)"
if have Rscript && rok BGLR && rok optiSel; then
  Rscript research/repro/wheat_export.R || echo "  wheat export failed"
  Rscript research/repro/mouse_export.R || echo "  mouse export failed"
  if [ -f /tmp/pig/FileS1/genotypes.txt ]; then
    Rscript research/repro/pig_export.R || echo "  pig export failed"
  else
    echo "  SKIP pig: /tmp/pig/FileS1/genotypes.txt absent — see REPRODUCE.md for the download"
  fi
else
  echo "SKIP: Rscript + BGLR + optiSel required"
fi

banner "5/7  Table 2, 'matrix-free' column — shipped solver vs optiSel, same R session"
if have Rscript && rok BGLR && rok optiSel && rok ocsrs; then
  # Medians over repetitions, all six panels; the pig row runs only if its
  # genotypes are present (the script skips it otherwise).
  Rscript research/repro/r_binding_bench_reps.R 5 || echo "  bench reps failed"
  # optiSel repeated on the two slow cells (n=5000, pig), which the reps script
  # samples once.
  if [ -f /tmp/pig/FileS1/genotypes.txt ]; then
    Rscript research/repro/r_optisel_slow_reps.R || echo "  slow-cell reps failed"
  fi
else
  echo "SKIP: needs Rscript + BGLR + optiSel + the ocsrs package"
  echo "      build it with:  R CMD INSTALL bindings/r"
fi

banner "6/7  Table 2, 'algorithm' column — active set given a dense G"
if have Rscript && rok BGLR; then
  # Export K at the same caps the matrix-free run used, so both columns solve the
  # identical instances (verified by matching gain and support).
  Rscript research/repro/r_export_K.R || echo "  K export failed"
fi
if have python3 && pyok numpy && pyok scipy; then
  for n in 1000 2000 5000 599 1814 3534; do
    ub_file="/tmp/bench_ub_${n}.txt"
    if [ -f "$ub_file" ] && [ -f "/tmp/bench_K_${n}.csv" ]; then
      ub="$(cat "$ub_file")"
      printf 'n=%-6s ' "$n"
      # five repetitions, K loaded once; the paper reports the median
      python3 research/repro/sf_at_ub.py "$n" "$ub" "$ub" "$ub" "$ub" "$ub"
    else
      echo "  SKIP n=${n}: exports absent (run step 6's R part first)"
    fi
  done
else
  echo "SKIP: python3 + numpy + scipy required"
fi

banner "7/7  AlphaMate (optional, not run automatically)"
echo "  Linux binary under Colima + Rosetta — see REPRODUCE.md (AlphaMate section)."

banner "done"
