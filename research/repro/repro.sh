#!/usr/bin/env bash
# Fresh automated reproduction; partial environments report omissions separately
# from actual failures. REPRO_STRICT=1 fails on missing automated dependencies.
# REPRO_FULL=1 includes expensive Clarabel n=10000. REPRO_STEPS is a comma-separated
# subset in fixed pipeline order, e.g. scaling_matrixfree,structured,fig_scaling.
set -u
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
cd "$ROOT"
PYTHON="${REPRO_PYTHON:-python3}"
RUN_DIR="$ROOT/artifacts/repro_runs/$(date -u +%Y%m%dT%H%M%SZ)-$$"
export REPRO_OUTPUT_DIR="$RUN_DIR"
export REPRO_INPUT_DIR="$RUN_DIR/inputs"
export REPRO_FIGURE_DIR="$RUN_DIR/figures"
export MPLCONFIGDIR="$RUN_DIR/plot-cache"
export XDG_CACHE_HOME="$RUN_DIR/cache"
mkdir -p "$RUN_DIR/logs" "$REPRO_INPUT_DIR" "$REPRO_FIGURE_DIR"
{
  printf 'started_utc=%s\nsteps=%s\nfull=%s\nstrict=%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "${REPRO_STEPS:-all}" "${REPRO_FULL:-0}" "${REPRO_STRICT:-0}"
  if command -v git >/dev/null 2>&1; then git rev-parse HEAD; git diff --stat; fi
  if command -v rustc >/dev/null 2>&1; then rustc --version; fi
  if command -v "$PYTHON" >/dev/null 2>&1; then "$PYTHON" --version; fi
  if command -v Rscript >/dev/null 2>&1; then Rscript -e 'print(sessionInfo()); for(p in c("ocsrs","AlphaSimR","optiSel","BGLR")) if(requireNamespace(p, quietly=TRUE)) cat(p,as.character(packageVersion(p)),"\n")'; fi
} > "$RUN_DIR/provenance.txt" 2>&1
failures=0
omissions=0
reason=""
selected() { [ -z "${REPRO_STEPS:-}" ] || case ",${REPRO_STEPS}," in *",$1,"*) true;; *) false;; esac; }
have() { command -v "$1" >/dev/null 2>&1; }
pydeps() { have "$PYTHON" && "$PYTHON" -c "import $1" >/dev/null 2>&1; }
rdeps() { have Rscript && Rscript -e "for (p in c($1)) suppressMessages(library(p, character.only=TRUE))" >/dev/null 2>&1; }
need_cargo() { reason="cargo unavailable"; have cargo; }
need_cvgen() { reason="uv and cargo required"; have uv && have cargo; }
need_real() { reason="R packages BGLR, optiSel and fresh ocsrs required"; rdeps '"BGLR","optiSel","ocsrs"'; }
need_pig() {
  reason="R packages data.table, optiSel, ocsrs and PIC genotypes/EBVs/pedigree required"
  rdeps '"data.table","optiSel","ocsrs"' || return 1
  local p="${PIG_DATA_DIR:-/tmp/pig/FileS1}"
  [ -f "$p/genotypes.txt" ] && [ -f "$p/ebvs.txt" ] && [ -f "$p/pedigree.txt" ]
}
need_exports() { reason="fresh operational wheat/mouse exports required"; have Rscript && [ -f "$RUN_DIR/operational_ne_realpanels.csv" ]; }
need_algorithm_real() { reason="numpy/scipy and successful same-run wheat/mouse caps/exports required"; pydeps 'numpy, scipy' && [ -f "$RUN_DIR/operational_ne_realpanels.csv" ] && [ -f "$REPRO_INPUT_DIR/wheat_K.csv" ] && [ -f "$REPRO_INPUT_DIR/mouse_K.csv" ]; }
need_algorithm_pig() { reason="numpy/scipy and successful same-run pig caps/exports required"; pydeps 'numpy, scipy' && [ -f "$RUN_DIR/operational_ne_pig.csv" ] && [ -f "$REPRO_INPUT_DIR/pig_K.csv" ]; }
need_table2() { reason="all three current operational panels required"; pydeps csv && [ -f "$RUN_DIR/algorithm_operational_wheat.csv" ] && [ -f "$RUN_DIR/algorithm_operational_mouse.csv" ] && [ -f "$RUN_DIR/algorithm_operational_pig.csv" ]; }
need_structured() { reason="R packages AlphaSimR, ocsrs and optiSel required"; rdeps '"AlphaSimR","ocsrs","optiSel"'; }
need_multigen() { reason="R packages AlphaSimR and ocsrs required"; rdeps '"AlphaSimR","ocsrs"'; }
need_scaling_fig() { reason="matplotlib/numpy and both successful same-run scaling CSVs required"; pydeps 'matplotlib, numpy' && [ -f "$RUN_DIR/scaling_matrixfree.csv" ] && [ -f "$RUN_DIR/scaling_structured.csv" ]; }
need_frontier_fig() { reason="matplotlib/numpy and complete same-run operational table required"; pydeps 'matplotlib, numpy' && [ -f "$RUN_DIR/table2_operational.csv" ]; }
need_multigen_fig() { reason="matplotlib/numpy and successful same-run multigeneration CSV required"; pydeps 'matplotlib, numpy' && [ -f "$RUN_DIR/multigen_validation.csv" ]; }
promote_outputs() {
  local outputs="$1" output source destination i count=0 ready=1
  local targets=() temporary=() previous=() had_previous=() keep_backup=()
  for output in $outputs; do
    case "$output" in
      fig_*)
        [ "${REPRO_UPDATE_FIGURES:-0}" = 1 ] || continue
        source="$REPRO_FIGURE_DIR/$output"; destination="$ROOT/research/$output";;
      *) source="$RUN_DIR/$output"; destination="$ROOT/artifacts/$output";;
    esac
    targets[$count]="$destination"
    temporary[$count]="$destination.repro-$$.tmp"
    previous[$count]="$destination.repro-$$.previous.tmp"
    had_previous[$count]=0
    keep_backup[$count]=0
    if [ -e "$destination" ]; then had_previous[$count]=1; fi
    count=$((count + 1))
    if ! cp "$source" "${temporary[$((count - 1))]}"; then ready=0; break; fi
    if [ "${had_previous[$((count - 1))]}" = 1 ] &&
       ! cp "$destination" "${previous[$((count - 1))]}"; then ready=0; break; fi
  done
  if [ "$ready" = 1 ]; then
    for ((i=0; i<count; i++)); do
      if ! mv "${temporary[$i]}" "${targets[$i]}"; then ready=0; break; fi
    done
    if [ "$ready" = 0 ]; then
      for ((i=0; i<count; i++)); do
        if [ "${had_previous[$i]}" = 1 ]; then
          if ! mv "${previous[$i]}" "${targets[$i]}"; then
            keep_backup[$i]=1
            printf 'Could not restore %s; recovery copy retained at %s\n' "${targets[$i]}" "${previous[$i]}" >&2
          fi
        else
          rm -f "${targets[$i]}"
        fi
      done
    fi
  fi
  for ((i=0; i<count; i++)); do
    rm -f "${temporary[$i]}"
    if [ "${keep_backup[$i]}" = 0 ]; then rm -f "${previous[$i]}"; fi
  done
  [ "$ready" = 1 ]
}
run_stage() {
  local id="$1" outputs="$2" needs="$3"
  shift 3
  selected "$id" || return 0
  printf '\n=== %s ===\n' "$id"
  local output
  for output in $outputs; do
    case "$output" in
      fig_*) :;;
      *) rm -f "$ROOT/artifacts/$output";;
    esac
  done
  if ! "$needs"; then
    printf 'OMITTED %s: %s\n' "$id" "$reason"
    printf '%s\tomitted\t%s\n' "$id" "$reason" >> "$RUN_DIR/status.tsv"
    omissions=$((omissions + 1))
    return 0
  fi
  if "$@" > "$RUN_DIR/logs/$id.log" 2>&1; then
    local valid=1
    for output in $outputs; do
      case "$output" in fig_*) [ -f "$REPRO_FIGURE_DIR/$output" ] || valid=0;; *) [ -f "$RUN_DIR/$output" ] || valid=0;; esac
    done
    if [ "$valid" = 1 ]; then
      cat "$RUN_DIR/logs/$id.log"
      if promote_outputs "$outputs"; then
        printf '%s\tsuccess\tlogs/%s.log\n' "$id" "$id" >> "$RUN_DIR/status.tsv"
        return 0
      fi
      printf 'could not promote stage outputs\n' >> "$RUN_DIR/logs/$id.log"
    else
      printf 'stage exited zero but a declared output is absent\n' >> "$RUN_DIR/logs/$id.log"
    fi
  fi
  cat "$RUN_DIR/logs/$id.log"
  for output in $outputs; do
    case "$output" in fig_*) rm -f "$REPRO_FIGURE_DIR/$output";; *) rm -f "$RUN_DIR/$output";; esac
  done
  printf 'FAILED %s; see %s\n' "$id" "$RUN_DIR/logs/$id.log"
  printf '%s\tfailed\tlogs/%s.log\n' "$id" "$id" >> "$RUN_DIR/status.tsv"
  failures=$((failures + 1))
}
table1() {
  local n bad=0
  for n in 1000 2000 5000; do
    cargo run --release --quiet -- compare --n "$n" --m 20000 --k-frac 0.017 || bad=1
  done
  if [ "${REPRO_FULL:-0}" = 1 ]; then
    cargo run --release --quiet -- compare --n 10000 --m 20000 --k-frac 0.017 || bad=1
  else
    printf 'POLICY OMISSION: Clarabel n=10000; enable REPRO_FULL=1\n'
  fi
  return "$bad"
}
if [ -n "${REPRO_STEPS:-}" ]; then
  selection="${REPRO_STEPS},"
  while [ -n "$selection" ]; do
    stage="${selection%%,*}"
    selection="${selection#*,}"
    case "$stage" in table1|cvxgenrust|scaling_matrixfree|operational_realpanels|operational_pig|exports|algorithm_realpanels|algorithm_pig|table2|structured|multigen|fig_scaling|fig_frontier|fig_multigen) :;; *) printf 'Unknown stage: %s\n' "$stage" >&2; exit 1;; esac
  done
fi
run_stage table1 "" need_cargo table1
run_stage cvxgenrust "" need_cvgen bash research/repro/cvxgenrust/run.sh 500 1000
run_stage scaling_matrixfree "scaling_matrixfree.csv" need_cargo cargo run --release --example scaling_matrixfree
run_stage operational_realpanels "operational_ne_realpanels.csv" need_real Rscript research/repro/r_operational_ne.R both
run_stage operational_pig "operational_ne_pig.csv" need_pig Rscript research/repro/r_operational_ne_pig.R
run_stage exports "" need_exports Rscript research/repro/export_K_panels.R wheat mouse
run_stage algorithm_realpanels "algorithm_operational_wheat.csv algorithm_operational_mouse.csv" need_algorithm_real "$PYTHON" research/repro/sf_operational_ne.py wheat mouse
run_stage algorithm_pig "algorithm_operational_pig.csv" need_algorithm_pig "$PYTHON" research/repro/sf_operational_ne_pig.py
run_stage table2 "table2_operational.csv" need_table2 "$PYTHON" research/repro/combine_operational.py
run_stage structured "scaling_structured.csv" need_structured Rscript research/repro/scaling_structured.R
run_stage multigen "multigen_validation.csv" need_multigen Rscript research/repro/multigen_validation.R
run_stage fig_scaling "fig_scaling.pdf fig_scaling.png" need_scaling_fig "$PYTHON" research/fig_scaling.py
run_stage fig_frontier "fig_frontier.pdf fig_frontier.png" need_frontier_fig "$PYTHON" research/fig_frontier.py
run_stage fig_multigen "fig_multigen.pdf fig_multigen.png" need_multigen_fig "$PYTHON" research/fig_multigen.py
printf '\nActual failures: %d; dependency/data omissions: %d\nRun evidence: %s\n' "$failures" "$omissions" "$RUN_DIR"
printf 'AlphaMate is manual; its historical rows require archived contributor files (see REPRODUCE.md).\n'
[ "$failures" -eq 0 ] || exit 1
if [ "${REPRO_STRICT:-0}" = 1 ] && [ "$omissions" -gt 0 ]; then exit 2; fi
