#!/usr/bin/env bash
# Run ONE CSP2 pipeline invocation, capture timing, and normalize its outputs into an
# equivalence manifest. Used to (a) freeze the baseline on the current MUMmer pipeline and
# (b) validate that a refactor PR produces the SAME manifest.
#
# Usage:
#   benchmark/run.sh <name> <snp|screen|align> <dataset_dir> [extra nextflow/CSP2 args...]
#
# Examples:
#   benchmark/run.sh baseline_small snp benchmark/data/refsel_small --n_ref 1
#   benchmark/run.sh baseline_full  snp benchmark/data/refsel_full  --n_ref 3
#
# Output: benchmark/runs/<name>/{out/,manifest.json,trace.txt,wall_time.txt,nextflow.log}
set -euo pipefail

if [[ $# -lt 3 ]]; then
    grep '^#' "$0" | sed 's/^# \{0,1\}//'; exit 2
fi
NAME="$1"; MODE="$2"; DATASET="$3"; shift 3
EXTRA=("$@")

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

ENV="${CSP2_BENCH_ENV:-csp2-bench}"
RUNDIR="benchmark/runs/${NAME}"
DATA_ABS="$(cd "$DATASET" && pwd)"

# Activate the env by PATH (mamba run's activation wrapper mangles `-c`/`--` flags).
ENV_PREFIX="$(mamba env list | awk -v e="$ENV" '$1==e{print $NF}')"
[[ -n "$ENV_PREFIX" && -x "$ENV_PREFIX/bin/nextflow" ]] || {
    echo "!! env '$ENV' not found or nextflow missing; run benchmark/setup.sh" >&2; exit 1; }
export PATH="$ENV_PREFIX/bin:$PATH"

rm -rf "$RUNDIR"
mkdir -p "$RUNDIR"

# SSL: conda proxy is a Zscaler MITM; disable verification for any conda ops the run
# triggers. NXF_OFFLINE stops nextflow from fetching its plugin index (whose TLS would
# fail against the JDK truststore) — CSP2 uses no plugins.
export CONDA_SSL_VERIFY=false MAMBA_SSL_VERIFY=false NXF_ANSI_LOG=false \
       NXF_OFFLINE=true NXF_PLUGINS_DEFAULT=false NXF_OPTS='-Xmx700m'

echo ">> run '$NAME'  mode=$MODE  dataset=$DATASET  extra=${EXTRA[*]:-(none)}"
START=$(date +%s)
set +e
nextflow \
    -c benchmark/nextflow_bench.config \
    run CSP2.nf -profile standard \
    --runmode "$MODE" \
    --fasta "$DATA_ABS" \
    --outroot "$RUNDIR" --out out \
    --notree skip \
    -work-dir "$RUNDIR/work" \
    -with-trace "$RUNDIR/trace.txt" \
    -with-report "$RUNDIR/report.html" \
    -with-timeline "$RUNDIR/timeline.html" \
    "${EXTRA[@]}" 2>&1 | tee "$RUNDIR/nextflow.log"
STATUS=${PIPESTATUS[0]}
set -e
END=$(date +%s)
echo "wall_seconds=$((END - START))" | tee "$RUNDIR/wall_time.txt"

if [[ $STATUS -ne 0 ]]; then
    echo "!! nextflow exited $STATUS — see $RUNDIR/nextflow.log" >&2
    exit $STATUS
fi

# CSP2 appends a timestamp to the output dir if 'out' already existed; resolve the real one.
REAL_OUT="$RUNDIR/out"
if [[ ! -d "$REAL_OUT" ]]; then
    REAL_OUT="$(find "$RUNDIR" -maxdepth 1 -type d -name 'out*' | sort | head -1)"
fi
echo ">> normalizing $REAL_OUT"
python benchmark/normalize.py "$REAL_OUT" -o "$RUNDIR/manifest.json"

echo ">> done: $RUNDIR/manifest.json  (wall $((END - START))s)"
