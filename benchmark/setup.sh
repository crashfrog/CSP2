#!/usr/bin/env bash
# One-time provisioning for the CSP2 equivalence benchmark on this machine.
#   1. create the pinned conda env (bio tools + python libs + java)
#   2. install nextflow INTO that env
#   3. verify every tool the pipeline shells out to is present
#
# Idempotent: safe to re-run. Network here sits behind a Zscaler MITM proxy, so conda SSL
# verification is disabled (CA bundle at ~/.certs is used elsewhere but conda downloads
# flake against it); this is acceptable on a controlled dev box.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"
ENV="${CSP2_BENCH_ENV:-csp2-bench}"
export CONDA_SSL_VERIFY=false

echo "== [1/3] conda env '$ENV' =="
if mamba env list | grep -qE "^\s*${ENV}\s"; then
    echo "   exists; updating from environment.yml"
    mamba env update -n "$ENV" -f benchmark/environment.yml
else
    mamba env create -y -n "$ENV" -f benchmark/environment.yml
fi

echo "== [2/3] nextflow =="
if mamba run -n "$ENV" bash -lc 'command -v nextflow' >/dev/null 2>&1; then
    echo "   already installed: $(mamba run -n "$ENV" nextflow -v 2>&1 | head -1)"
else
    # Prefer a bioconda LTS close to CSP2's dev era (22.10.7). Fall back progressively.
    installed=""
    for ver in "nextflow=24.10.*" "nextflow=24.*" "nextflow"; do
        if mamba install -y -n "$ENV" -c bioconda -c conda-forge $ver 2>/dev/null; then
            installed="$ver"; break
        fi
        echo "   '$ver' unavailable, trying next..."
    done
    [[ -n "$installed" ]] || { echo "!! could not install nextflow from bioconda" >&2; exit 1; }
    echo "   installed via '$installed'"
fi

echo "== [2b/3] deterministic nucmer shadow =="
# nucmer4 defaults to 2 threads and its multithreaded anchoring is NONDETERMINISTIC:
# back-to-back dnadiff runs on the same inputs yield different SNP/indel/coord counts,
# which would make baseline-vs-candidate equivalence impossible. dnadiff invokes nucmer
# from its OWN bin dir (ignores PATH), so we shadow the env's nucmer with a -t 1 wrapper.
# This affects ONLY the benchmark env; baseline and candidate both use it, so the
# equivalence comparison stays valid. (Verified: -t 1 is bit-identical across runs.)
ENV_PREFIX="$(mamba env list | awk -v e="$ENV" '$1==e{print $NF}')"
NUCMER="$ENV_PREFIX/bin/nucmer"
if [[ -f "$NUCMER" && ! -f "$ENV_PREFIX/bin/nucmer.real" ]]; then
    mv "$NUCMER" "$ENV_PREFIX/bin/nucmer.real"
    cat > "$NUCMER" <<'EOF'
#!/usr/bin/env bash
# Benchmark shadow: force nucmer single-threaded for reproducible alignments.
exec "$(dirname "$0")/nucmer.real" -t 1 "$@"
EOF
    chmod +x "$NUCMER"
    echo "   installed single-thread nucmer shadow"
else
    echo "   shadow already present (or nucmer missing)"
fi

echo "== [3/3] tool check =="
mamba run -n "$ENV" bash -lc '
    fail=0
    for t in nextflow java python dnadiff nucmer show-coords show-snps bedtools mash \
             kmercountexact.sh; do
        if command -v "$t" >/dev/null 2>&1; then printf "  ok   %s\n" "$t";
        else printf "  MISS %s\n" "$t"; fail=1; fi
    done
    python - <<PY
import importlib, sys
miss=[m for m in ("pandas","numpy","scipy","sklearn","pybedtools","Bio") if not importlib.util.find_spec(m)]
print("  py-libs MISSING:", miss) if miss else print("  ok   python libs (pandas/numpy/scipy/sklearn/pybedtools/Bio)")
sys.exit(1 if miss else 0)
PY
    exit $fail
'
echo "== setup complete =="
