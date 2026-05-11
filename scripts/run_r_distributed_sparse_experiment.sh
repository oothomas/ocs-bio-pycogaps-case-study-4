#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

PROFILE="${1:-full}"
SPARSITY_MODE="${2:-on}"

PREPROCESSED_H5AD="${ROOT}/data/processed/preprocessed_cells_hvg3000.h5ad"
RUNNER="${ROOT}/scripts/cogaps_run_distributed_r.R"

K=7
SEED=2
DISTRIBUTED_MODE="single-cell"
NSETS=4
CUT=7
MIN_NS=2
MAX_NS=6
BP_WORKERS=4
STRATA_COLS="cell_type,condition,replicate"

case "${PROFILE}" in
  smoke)
    NITER=50
    OUTPUT_FREQUENCY=25
    SUFFIX="smoke"
    ;;
  full)
    NITER=2000
    OUTPUT_FREQUENCY=100
    SUFFIX="full"
    ;;
  *)
    echo "Usage: bash scripts/run_r_distributed_sparse_experiment.sh [smoke|full] [on|off]" >&2
    exit 1
    ;;
esac

case "${SPARSITY_MODE}" in
  on)
    SPARSE_FLAG="--use-sparse-opt"
    OUTDIR="${ROOT}/data/results_r_distributed_sparse_on_${SUFFIX}"
    ;;
  off)
    SPARSE_FLAG="--no-sparse-opt"
    OUTDIR="${ROOT}/data/results_r_distributed_sparse_off_${SUFFIX}"
    ;;
  *)
    echo "Usage: bash scripts/run_r_distributed_sparse_experiment.sh [smoke|full] [on|off]" >&2
    exit 1
    ;;
esac

SHARED_DIR="${ROOT}/data/results_r_distributed_shared"
mkdir -p "${SHARED_DIR}"
SHARED_EXPLICIT_SETS="${SHARED_DIR}/cogaps_K${K}_seed${SEED}_iter${NITER}_single_cell.explicit_sets.rds"

if [[ ! -f "${SHARED_EXPLICIT_SETS}" ]]; then
  Rscript "${RUNNER}" \
    --preprocessed-h5ad "${PREPROCESSED_H5AD}" \
    --outdir "${SHARED_DIR}" \
    --k "${K}" \
    --seed "${SEED}" \
    --n-iter "${NITER}" \
    --distributed-mode "${DISTRIBUTED_MODE}" \
    --n-sets "${NSETS}" \
    --cut "${CUT}" \
    --min-ns "${MIN_NS}" \
    --max-ns "${MAX_NS}" \
    --bp-workers "${BP_WORKERS}" \
    --strata-cols "${STRATA_COLS}" \
    --output-frequency "${OUTPUT_FREQUENCY}" \
    --write-sets-only
fi

Rscript "${RUNNER}" \
  --preprocessed-h5ad "${PREPROCESSED_H5AD}" \
  --outdir "${OUTDIR}" \
  --k "${K}" \
  --seed "${SEED}" \
  --n-iter "${NITER}" \
  --distributed-mode "${DISTRIBUTED_MODE}" \
  --n-sets "${NSETS}" \
  --cut "${CUT}" \
  --min-ns "${MIN_NS}" \
  --max-ns "${MAX_NS}" \
  --bp-workers "${BP_WORKERS}" \
  --strata-cols "${STRATA_COLS}" \
  --output-frequency "${OUTPUT_FREQUENCY}" \
  --explicit-sets-rds "${SHARED_EXPLICIT_SETS}" \
  ${SPARSE_FLAG}
