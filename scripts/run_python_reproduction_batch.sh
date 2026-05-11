#!/usr/bin/env bash
set -uo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
LOG_DIR="${ROOT}/data/results_python_batch_logs"
mkdir -p "${LOG_DIR}"

RUN_ID="${RUN_ID:-$(date +%Y%m%d_%H%M%S)}"
MASTER_LOG="${LOG_DIR}/python_light_local_batch_${RUN_ID}.log"
STATUS_TSV="${LOG_DIR}/python_light_local_batch_${RUN_ID}.status.tsv"

{
  echo -e "step\tstatus\tstarted_at\tfinished_at"
} > "${STATUS_TSV}"

timestamp() {
  date "+%Y-%m-%d %H:%M:%S %Z"
}

log() {
  printf "[%s] %s\n" "$(timestamp)" "$*" | tee -a "${MASTER_LOG}"
}

run_step() {
  local step_name="$1"
  shift
  local started finished rc

  started="$(timestamp)"
  log "START ${step_name}"
  "$@"
  rc=$?
  finished="$(timestamp)"

  if [[ ${rc} -eq 0 ]]; then
    log "DONE  ${step_name}"
    printf "%s\tok\t%s\t%s\n" "${step_name}" "${started}" "${finished}" >> "${STATUS_TSV}"
  else
    log "FAIL  ${step_name} (exit ${rc})"
    printf "%s\tfailed(%s)\t%s\t%s\n" "${step_name}" "${rc}" "${started}" "${finished}" >> "${STATUS_TSV}"
  fi

  return 0
}

log "Python CoGAPS light-local batch starting"
log "ROOT=${ROOT}"
log "MASTER_LOG=${MASTER_LOG}"
log "STATUS_TSV=${STATUS_TSV}"

cd "${ROOT}" || exit 1

run_step "python_single_sparse_light" bash scripts/run_python_singleprocess_reproduction.sh sparse
run_step "python_distributed_sparse_on_light" bash scripts/run_python_distributed_sparse_experiment.sh full on
run_step "python_distributed_sparse_off_light" bash scripts/run_python_distributed_sparse_experiment.sh full off

log "Python CoGAPS light-local batch finished"
