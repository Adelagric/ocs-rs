#!/usr/bin/env bash
# Reproduce the support-first results. Run from anywhere; it cds to the crate
# root. Each step is guarded: a missing toolchain or dataset is reported and
# skipped, so a partial environment still reproduces the parts it can.
#   bash research/repro/repro.sh
set -u
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
cd "$ROOT"

have() { command -v "$1" >/dev/null 2>&1; }
pyok() { python3 -c "import $1" >/dev/null 2>&1; }
rok() { Rscript -e "suppressMessages(library($1))" >/dev/null 2>&1; }
banner() { printf '\n=== %s ===\n' "$1"; }

banner "1/5  Rust release timings + scaling sweep (Figure 1 data)"
if have cargo; then
  cargo run --release --example bench_sexed
  cargo run --release --example scaling_matrixfree
else
  echo "SKIP: cargo not found"
fi

banner "2/5  Figure 1"
if have python3 && pyok matplotlib; then
  python3 research/fig_scaling.py
else
  echo "SKIP: python3 + matplotlib required"
fi

banner "3/5  Real-data GRM export (wheat, mouse via BGLR; pig if downloaded)"
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

banner "4/5  support-first vs optiSel on the exported instances"
if have python3 && pyok numpy && pyok scipy; then
  for n in 599 1814 3534; do
    if [ -f "/tmp/bench_K_${n}.csv" ]; then
      python3 research/support_first_sex.py "$n"
    else
      echo "  SKIP n=${n}: /tmp/bench_K_${n}.csv absent (run step 3 first)"
    fi
  done
else
  echo "SKIP: python3 + numpy + scipy required"
fi

banner "5/5  AlphaMate (optional, not run automatically)"
echo "  Linux binary under Colima + Rosetta — see REPRODUCE.md (AlphaMate section)."

banner "done"
