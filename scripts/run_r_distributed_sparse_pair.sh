#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

PAIR_LOG_DIR="${ROOT}/data/results_r_distributed_pair_logs"
mkdir -p "${PAIR_LOG_DIR}"

TS="${RUN_ID:-$(date +%Y%m%d_%H%M%S)}"
PAIR_LOG="${PAIR_LOG_DIR}/distributed_sparse_pair_${TS}.log"
PAIR_STATUS="${PAIR_LOG_DIR}/distributed_sparse_pair_${TS}.status"
LATEST_LOG_LINK="${PAIR_LOG_DIR}/distributed_sparse_pair_latest.log"
LATEST_STATUS_LINK="${PAIR_LOG_DIR}/distributed_sparse_pair_latest.status"

ln -sfn "${PAIR_LOG}" "${LATEST_LOG_LINK}"
ln -sfn "${PAIR_STATUS}" "${LATEST_STATUS_LINK}"

{
  echo "[START] $(date '+%Y-%m-%d %H:%M:%S %Z')"
  echo "[ROOT] ${ROOT}"
  echo "[STEP] full sparse-on"
  bash "${ROOT}/scripts/run_r_distributed_sparse_experiment.sh" full on
  echo "[DONE] full sparse-on $(date '+%Y-%m-%d %H:%M:%S %Z')"

  echo "[STEP] full sparse-off"
  bash "${ROOT}/scripts/run_r_distributed_sparse_experiment.sh" full off
  echo "[DONE] full sparse-off $(date '+%Y-%m-%d %H:%M:%S %Z')"

  echo "[FINISH] $(date '+%Y-%m-%d %H:%M:%S %Z')"
  echo "ok" > "${PAIR_STATUS}"
} 2>&1 | tee "${PAIR_LOG}"
