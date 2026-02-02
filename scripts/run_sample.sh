#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
JOBS_DIR="${ROOT_DIR}/outputs"
RUST_BIN="${ROOT_DIR}/target/release/errat"
PDB_SAMPLE="${ROOT_DIR}/data/Hpyr004913.1-R65830.mRNA_relaxed_rank_001_alphafold2_ptm_model_5_seed_000.pdb"

JOB_ID="${1:-job_run}"
PROTEIN_ID="${2:-Hpyr004913}"

mkdir -p "${JOBS_DIR}/${JOB_ID}"
cp -f "${PDB_SAMPLE}" "${JOBS_DIR}/${JOB_ID}/errat.pdb"

ERRAT_JOBS_PATH="${JOBS_DIR}" "${RUST_BIN}" "${PROTEIN_ID}" "${JOB_ID}"

echo "[run] output: ${JOBS_DIR}/${JOB_ID}/errat.logf"
echo "[run] output: ${JOBS_DIR}/${JOB_ID}/errat.ps"
