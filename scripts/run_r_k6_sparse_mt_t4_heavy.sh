#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
IMAGE="${COGAPS_RUNTIME_IMAGE:-othomas2/pycogaps-runtime-guide:0.3.0}"
OUTDIR="data/results_r_k6_sparse_mt_t4_heavy"
TAG="K6_seed2_iter2000"

cd "${ROOT}"
mkdir -p "${OUTDIR}"

# Clear only files from this interrupted benchmark before relaunching.
rm -f \
  "${OUTDIR}/cogaps_${TAG}.rds" \
  "${OUTDIR}/cogaps_${TAG}.metrics.json" \
  "${OUTDIR}/cogaps_${TAG}.checkpoint.out"

docker run --rm --platform linux/amd64 \
  -v "${ROOT}:/workspace/case-study" \
  -w /workspace/case-study \
  "${IMAGE}" \
  bash -lc 'set -euo pipefail
    export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
    Rscript -e "cat(\"R CoGAPS OpenMP support:\", CoGAPS::compiledWithOpenMPSupport(), \"\n\")"
    Rscript scripts/cogaps_run_one_singleprocess_r.R \
      --preprocessed-h5ad data/processed/preprocessed_cells_hvg3000.h5ad \
      --outdir data/results_r_k6_sparse_mt_t4_heavy \
      --k 6 \
      --seed 2 \
      --n-iter 2000 \
      --top-genes 50 \
      --stim-label stim \
      --cogaps-threads 4 \
      --async-updates \
      --use-sparse-opt \
      --output-frequency 100 \
      --checkpoint-interval 500 \
      --n-snapshots 10 \
      --snapshot-phase all \
      --take-pump-samples \
      --force-rerun'
