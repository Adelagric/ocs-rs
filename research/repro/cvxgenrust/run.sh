#!/usr/bin/env bash
# Independent check of the Table-1 conic baseline: generate a Clarabel solver for
# the same cone program with cvxgenrust (Zhu & Boedecker 2026, a third-party code
# generator) from a six-line CVXPY model, and time it on the identical instance
# next to the crate's own hand-assembled conic route.
#
#   bash research/repro/cvxgenrust/run.sh [n ...]      (default: 500 1000)
#
# Needs: uv (Python >= 3.12; installs cvxgenrust 0.1.0 into a local venv), cargo.
# Generation is O(n^2) in memory (8 GB at n=1000, ~32 GB at n=2000): stay at
# n <= 1000 on a laptop. Everything lands in ./artifacts (git-ignored).
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
ART="$HERE/artifacts"
mkdir -p "$ART"
cd "$HERE"

if ! command -v uv >/dev/null 2>&1; then
  echo "SKIP: uv not found (needed to install cvxgenrust)"; exit 0
fi
[ -d "$ART/.venv" ] || uv venv -q --python 3.13 "$ART/.venv"
# shellcheck disable=SC1091
source "$ART/.venv/bin/activate"
uv pip install -q "cvxgenrust==0.1.0"

( cd dump_ocs && cargo build --release -q )

for n in "${@:-500 1000}"; do
  echo "=== n=$n ==="
  ( cd "$ART" && "$HERE/dump_ocs/target/release/dump_ocs" "$n" )
  if [ "$(uname)" = Darwin ]; then TIME="/usr/bin/time -l"; else TIME="/usr/bin/time -v"; fi
  ( cd "$ART" && $TIME python "$HERE/gen_ocs.py" "$n" 2>&1 | grep -E "^\[gen\]|maximum resident|Maximum resident" )
  crate="$ART/ocs_n${n}_cgr"
  sed "s/MODULE/ocs_n${n}/g" ocs_bench.rs.tmpl > "$crate/examples/ocs_bench.rs"
  # Apples to apples with the crate: the same pinned clarabel and release profile.
  sed -i '' 's/^clarabel = "0.11.1"/clarabel = "=0.11.1"/' "$crate/Cargo.toml"
  grep -q '^\[profile.release\]' "$crate/Cargo.toml" || printf '\n[profile.release]\nopt-level = 3\nlto = "thin"\n' >> "$crate/Cargo.toml"
  # The generated crate sits inside this repository: keep it out of the crate's workspace.
  grep -q '^\[workspace\]' "$crate/Cargo.toml" || printf '\n[workspace]\n' >> "$crate/Cargo.toml"
  ( cd "$crate" && cargo build --release -q --example ocs_bench )
  "$crate/target/release/examples/ocs_bench" "$ART/ocs_n${n}.bin" 3

  # Same family with L declared lower-triangular (sparsity=): the generator then
  # carries only the triangle and needs no Clarabel setting to be at parity.
  ( cd "$ART" && $TIME python "$HERE/gen_ocs_sparse.py" "$n" 2>&1 | grep -E "^\[gen-sparse\]|maximum resident|Maximum resident" )
  scrate="$ART/ocs_sparse_n${n}_cgr"
  sed "s/MODULE/ocs_sparse_n${n}/g" ocs_sparse_bench.rs.tmpl > "$scrate/examples/ocs_sparse_bench.rs"
  sed -i '' 's/^clarabel = "0.11.1"/clarabel = "=0.11.1"/' "$scrate/Cargo.toml"
  grep -q '^\[profile.release\]' "$scrate/Cargo.toml" || printf '\n[profile.release]\nopt-level = 3\nlto = "thin"\n' >> "$scrate/Cargo.toml"
  grep -q '^\[workspace\]' "$scrate/Cargo.toml" || printf '\n[workspace]\n' >> "$scrate/Cargo.toml"
  ( cd "$scrate" && cargo build --release -q --example ocs_sparse_bench )
  "$scrate/target/release/examples/ocs_sparse_bench" "$ART/ocs_n${n}.bin"
done
