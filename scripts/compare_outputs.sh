#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
JOBS_DIR="${ROOT_DIR}/outputs"
JOB_CPP="${1:-job_run_cpp}"
JOB_RS="${2:-job_run_rs}"

LOG_CPP="${JOBS_DIR}/${JOB_CPP}/errat.logf"
LOG_RS="${JOBS_DIR}/${JOB_RS}/errat.logf"
PS_CPP="${JOBS_DIR}/${JOB_CPP}/errat.ps"
PS_RS="${JOBS_DIR}/${JOB_RS}/errat.ps"

if [[ ! -f "${LOG_CPP}" || ! -f "${LOG_RS}" || ! -f "${PS_CPP}" || ! -f "${PS_RS}" ]]; then
  echo "[compare] missing outputs. expected:" >&2
  echo "  ${LOG_CPP}" >&2
  echo "  ${LOG_RS}" >&2
  echo "  ${PS_CPP}" >&2
  echo "  ${PS_RS}" >&2
  exit 1
fi

diff -u "${LOG_CPP}" "${LOG_RS}"
diff -u "${PS_CPP}" "${PS_RS}"

echo "[compare] logf and ps are identical"
