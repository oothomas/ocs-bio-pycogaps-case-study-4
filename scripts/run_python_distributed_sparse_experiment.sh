#!/usr/bin/env bash
set -euo pipefail

PROFILE="${1:-full}"
SPARSE_MODE="${2:-on}"

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
IMAGE="${COGAPS_RUNTIME_IMAGE:-othomas2/pycogaps-runtime-guide:0.3.0}"
DOCKER_TTY=()
if [[ -t 1 ]]; then
  DOCKER_TTY=(-t)
fi

PREPROCESSED_H5AD="${ROOT}/data/processed/preprocessed_cells_hvg3000.h5ad"
COGAPS_INPUT_H5AD="${ROOT}/data/processed/cogaps_input_genesxcells_hvg3000_float64.h5ad"
SHARED_DIR="${ROOT}/data/results_r_distributed_shared"
EXPLICIT_SETS_RDS="${SHARED_DIR}/cogaps_K7_seed2_iter2000_single_cell.explicit_sets.rds"
EXPLICIT_SETS_JSON="${SHARED_DIR}/cogaps_K7_seed2_iter2000_single_cell.explicit_sets.python.json"
EXPLICIT_SETS_SUMMARY="${SHARED_DIR}/cogaps_K7_seed2_iter2000_single_cell.explicit_sets.python.summary.csv"

if [[ ! -f "${PREPROCESSED_H5AD}" ]]; then
  echo "Missing preprocessed H5AD: ${PREPROCESSED_H5AD}" >&2
  exit 1
fi

if [[ ! -f "${COGAPS_INPUT_H5AD}" ]]; then
  echo "Missing CoGAPS input H5AD: ${COGAPS_INPUT_H5AD}" >&2
  exit 1
fi

if [[ ! -f "${EXPLICIT_SETS_RDS}" ]]; then
  echo "Missing explicit set RDS: ${EXPLICIT_SETS_RDS}" >&2
  exit 1
fi

docker run --rm "${DOCKER_TTY[@]}" --platform linux/amd64 \
  -v "${ROOT}:/workspace/case-study" \
  -w /workspace/case-study \
  "${IMAGE}" \
  bash -lc "Rscript scripts/export_distributed_explicit_sets_for_python.R \
    --explicit-sets-rds data/results_r_distributed_shared/cogaps_K7_seed2_iter2000_single_cell.explicit_sets.rds \
    --preprocessed-h5ad data/processed/preprocessed_cells_hvg3000.h5ad \
    --out-json data/results_r_distributed_shared/cogaps_K7_seed2_iter2000_single_cell.explicit_sets.python.json \
    --subset-summary-csv data/results_r_distributed_shared/cogaps_K7_seed2_iter2000_single_cell.explicit_sets.python.summary.csv \
    --label distributed_single_cell_k7_seed2"

case "${PROFILE}" in
  smoke)
    N_ITER=50
    N_SNAPSHOTS=0
    TAKE_PUMP=""
    OUTPUT_FREQ=100
    SNAPSHOT_PHASE="sampling"
    SUFFIX="smoke"
    ;;
  full)
    N_ITER=2000
    N_SNAPSHOTS=0
    TAKE_PUMP=""
    OUTPUT_FREQ=100
    SNAPSHOT_PHASE="sampling"
    SUFFIX="light_full"
    ;;
  *)
    echo "Usage: bash scripts/run_python_distributed_sparse_experiment.sh [full|smoke] [on|off]" >&2
    exit 1
    ;;
esac

case "${SPARSE_MODE}" in
  on)
    OUTDIR="data/results_python_distributed_sparse_on_${SUFFIX}"
    SPARSE_ARG="--use-sparse-opt"
    ;;
  off)
    OUTDIR="data/results_python_distributed_sparse_off_${SUFFIX}"
    SPARSE_ARG="--no-sparse-opt"
    ;;
  *)
    echo "Sparse mode must be 'on' or 'off'." >&2
    exit 1
    ;;
esac

docker run --rm "${DOCKER_TTY[@]}" --platform linux/amd64 \
  -v "${ROOT}:/workspace/case-study" \
  -w /workspace/case-study \
  "${IMAGE}" \
  bash -lc "PYTHONWARNINGS=ignore::FutureWarning python scripts/cogaps_run_distributed_python.py \
    --cogaps-input-h5ad data/processed/cogaps_input_genesxcells_hvg3000_float64.h5ad \
    --preprocessed-h5ad data/processed/preprocessed_cells_hvg3000.h5ad \
    --worker-source-h5ad data/processed/preprocessed_cells_hvg3000.h5ad \
    --worker-source-layout cellsxgenes \
    --explicit-sets-json data/results_r_distributed_shared/cogaps_K7_seed2_iter2000_single_cell.explicit_sets.python.json \
    --outdir ${OUTDIR} \
    --k 7 \
    --seed 2 \
    --n-iter ${N_ITER} \
    --top-genes 50 \
    --stim-label stim \
    --distributed-mode single-cell \
    --n-sets 4 \
    --cut 7 \
    --min-ns 2 \
    --max-ns 6 \
    --pool-workers 4 \
    --messages \
    --output-frequency ${OUTPUT_FREQ} \
    --n-snapshots ${N_SNAPSHOTS} \
    --snapshot-phase ${SNAPSHOT_PHASE} \
    ${SPARSE_ARG} \
    ${TAKE_PUMP}"
