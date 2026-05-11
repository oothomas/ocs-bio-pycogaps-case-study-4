#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"

THREADS="${1:-}"
MODE="${2:-pilot}"
PREPROCESSED_H5AD="${PREPROCESSED_H5AD:-${ROOT_DIR}/data/processed/preprocessed_cells_hvg3000.h5ad}"
RUNNER="${ROOT_DIR}/scripts/cogaps_run_one_singleprocess_r.R"

if [[ -z "${THREADS}" ]]; then
  cat >&2 <<'EOF'
Usage:
  bash scripts/run_r_nosparse_multithread_experiment.sh <threads> [pilot|full]

Examples:
  bash scripts/run_r_nosparse_multithread_experiment.sh 2 pilot
  bash scripts/run_r_nosparse_multithread_experiment.sh 4 full
EOF
  exit 1
fi

if [[ "${MODE}" != "pilot" && "${MODE}" != "full" ]]; then
  echo "ERROR: mode must be 'pilot' or 'full'" >&2
  exit 1
fi

case "${THREADS}" in
  ''|*[!0-9]*)
    echo "ERROR: threads must be a positive integer" >&2
    exit 1
    ;;
esac

if (( THREADS < 1 )); then
  echo "ERROR: threads must be >= 1" >&2
  exit 1
fi

if [[ ! -f "${PREPROCESSED_H5AD}" ]]; then
  echo "ERROR: preprocessed H5AD not found at ${PREPROCESSED_H5AD}" >&2
  exit 1
fi

export OMP_NUM_THREADS="${THREADS}"
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

if [[ "${MODE}" == "pilot" ]]; then
  OUTDIR="${ROOT_DIR}/data/results_r_nosparse_mt_t${THREADS}_pilot"
  NITER=1000
  CHECKPOINT_INTERVAL=250
  SNAPSHOTS=0
  TAKE_PUMP_ARGS=()
else
  OUTDIR="${ROOT_DIR}/data/results_r_nosparse_mt_t${THREADS}_full"
  NITER=2000
  CHECKPOINT_INTERVAL=500
  SNAPSHOTS=10
  TAKE_PUMP_ARGS=(--take-pump-samples)
fi

mkdir -p "${OUTDIR}"

echo "Launching no-sparse R CoGAPS experiment"
echo "  mode: ${MODE}"
echo "  threads: ${THREADS}"
echo "  output: ${OUTDIR}"
echo "  OMP_NUM_THREADS=${OMP_NUM_THREADS}"
echo "  OPENBLAS_NUM_THREADS=${OPENBLAS_NUM_THREADS}"
echo "  MKL_NUM_THREADS=${MKL_NUM_THREADS}"
echo "  NUMEXPR_NUM_THREADS=${NUMEXPR_NUM_THREADS}"

exec Rscript "${RUNNER}" \
  --preprocessed-h5ad "${PREPROCESSED_H5AD}" \
  --outdir "${OUTDIR}" \
  --k 7 \
  --seed 2 \
  --n-iter "${NITER}" \
  --cogaps-threads "${THREADS}" \
  --async-updates \
  --no-sparse-opt \
  --output-frequency 100 \
  --checkpoint-interval "${CHECKPOINT_INTERVAL}" \
  --n-snapshots "${SNAPSHOTS}" \
  --snapshot-phase all \
  "${TAKE_PUMP_ARGS[@]}" \
  --force-rerun
