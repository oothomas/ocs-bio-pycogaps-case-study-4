#!/usr/bin/env bash
set -euo pipefail

MODE="${1:-sparse}"

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
IMAGE="${COGAPS_RUNTIME_IMAGE:-othomas2/pycogaps-runtime-guide:0.3.0}"
DOCKER_TTY=()
if [[ -t 1 ]]; then
  DOCKER_TTY=(-t)
fi

INPUT_H5AD="${ROOT}/data/processed/cogaps_input_genesxcells_hvg3000_float64.h5ad"
PREPROCESSED_H5AD="${ROOT}/data/processed/preprocessed_cells_hvg3000.h5ad"

COMMON_ARGS=(
  --cogaps-input-h5ad "data/processed/cogaps_input_genesxcells_hvg3000_float64.h5ad"
  --preprocessed-h5ad "data/processed/preprocessed_cells_hvg3000.h5ad"
  --k 7
  --seed 2
  --n-iter 2000
  --top-genes 50
  --stim-label stim
  --blas-threads 1
  --cogaps-threads 1
  --output-frequency 100
  --checkpoint-interval 0
  --n-snapshots 0
  --snapshot-phase sampling
)

case "${MODE}" in
  sparse)
    OUTDIR="data/results_python_sparse_light_local"
    SPARSE_ARG="--use-sparse-opt"
    ;;
  nosparse)
    OUTDIR="data/results_python_nosparse_light_local"
    SPARSE_ARG="--no-sparse-opt"
    ;;
  *)
    echo "Usage: bash scripts/run_python_singleprocess_reproduction.sh [sparse|nosparse]" >&2
    exit 1
    ;;
esac

if [[ ! -f "${INPUT_H5AD}" ]]; then
  echo "Missing input H5AD: ${INPUT_H5AD}" >&2
  exit 1
fi

docker run --rm "${DOCKER_TTY[@]}" --platform linux/amd64 \
  -v "${ROOT}:/workspace/case-study" \
  -w /workspace/case-study \
  "${IMAGE}" \
  bash -lc "PYTHONWARNINGS=ignore::FutureWarning python scripts/cogaps_run_one_singleprocess.py ${SPARSE_ARG} --outdir ${OUTDIR} ${COMMON_ARGS[*]}"
