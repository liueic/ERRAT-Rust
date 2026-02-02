#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
JOBS_DIR="${ROOT_DIR}/outputs"
RUST_BIN="${ROOT_DIR}/target/release/errat"
CPP_BIN="${ROOT_DIR}/errat_cpp"
PDB_SAMPLE="${ROOT_DIR}/data/Hpyr004913.1-R65830.mRNA_relaxed_rank_001_alphafold2_ptm_model_5_seed_000.pdb"

RUNS="${1:-10}"
PROTEIN_ID="${2:-Hpyr004913}"
JOB_CPP="${3:-job_run_cpp}"
JOB_RS="${4:-job_run_rs}"

mkdir -p "${JOBS_DIR}/${JOB_CPP}" "${JOBS_DIR}/${JOB_RS}"
cp -f "${PDB_SAMPLE}" "${JOBS_DIR}/${JOB_CPP}/errat.pdb"
cp -f "${PDB_SAMPLE}" "${JOBS_DIR}/${JOB_RS}/errat.pdb"

CPP_TIMES="${ROOT_DIR}/cpp_times_${RUNS}.txt"
RS_TIMES="${ROOT_DIR}/rs_times_${RUNS}.txt"
rm -f "${CPP_TIMES}" "${RS_TIMES}"

if [[ -x "${CPP_BIN}" ]]; then
  for _ in $(seq 1 "${RUNS}"); do
    ulimit -s unlimited
    /usr/bin/time -f "%e" "${CPP_BIN}" "${PROTEIN_ID}" "${JOB_CPP}" 1>/dev/null 2>>"${CPP_TIMES}"
  done
else
  echo "[bench] C++ binary not found: ${CPP_BIN} (skip C++ benchmark)" >&2
fi

if [[ -x "${RUST_BIN}" ]]; then
  for _ in $(seq 1 "${RUNS}"); do
    ERRAT_JOBS_PATH="${JOBS_DIR}" /usr/bin/time -f "%e" "${RUST_BIN}" "${PROTEIN_ID}" "${JOB_RS}" 1>/dev/null 2>>"${RS_TIMES}"
  done
else
  echo "[bench] Rust binary not found: ${RUST_BIN}" >&2
  exit 1
fi

echo "C++ mean:"; if [[ -f "${CPP_TIMES}" ]]; then awk '{s+=$1} END{printf "%.6f\n", s/NR}' "${CPP_TIMES}"; fi

echo "Rust mean:"; awk '{s+=$1} END{printf "%.6f\n", s/NR}' "${RS_TIMES}"

if [[ -f "${CPP_TIMES}" ]]; then
  CPP_MEAN=$(awk '{s+=$1} END{printf "%.6f", s/NR}' "${CPP_TIMES}")
  RS_MEAN=$(awk '{s+=$1} END{printf "%.6f", s/NR}' "${RS_TIMES}")
  awk -v c="${CPP_MEAN}" -v r="${RS_MEAN}" 'BEGIN{printf "speedup=%.2fx\n", c/r; printf "time_reduction=%.1f%%\n", (1 - r/c)*100}'
fi
