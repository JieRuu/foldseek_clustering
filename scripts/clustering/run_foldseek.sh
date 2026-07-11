#!/bin/bash
set -euo pipefail

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"

INPUT_DIR="${PROJECT_DIR}/results/antigen_only_pdb"
OUTPUT_DIR="${PROJECT_DIR}/results/foldseek"

DB_DIR="${OUTPUT_DIR}/db"
TMP_ROOT="${OUTPUT_DIR}/tmp"

mkdir -p "${OUTPUT_DIR}"

echo "========================================"
echo "Cleaning old Foldseek results"
echo "========================================"

rm -rf "${OUTPUT_DIR}"
mkdir -p "${OUTPUT_DIR}"
mkdir -p "${DB_DIR}"
mkdir -p "${TMP_ROOT}"

COVERAGE=0.90
ALIGNMENT_TYPE=1

THRESHOLDS=(
    0.70
    0.75
    0.80
    0.85
    0.90
)

DB="${DB_DIR}/antigen_only_db"

echo "========================================"
echo "Checking input"
echo "========================================"

if [ ! -d "${INPUT_DIR}" ]; then
    echo "ERROR: input directory not found:"
    echo "${INPUT_DIR}"
    exit 1
fi

N_INPUT=$(find "${INPUT_DIR}" -name "*.pdb" | wc -l)

echo "Input dir: ${INPUT_DIR}"
echo "Number of PDB files: ${N_INPUT}"

if [ "${N_INPUT}" -eq 0 ]; then
    echo "ERROR: no PDB files found in input directory."
    exit 1
fi

echo ""
echo "========================================"
echo "Creating Foldseek database"
echo "========================================"

foldseek createdb "${INPUT_DIR}" "${DB}"

for THRESHOLD in "${THRESHOLDS[@]}"
do
    echo ""
    echo "========================================"
    echo "Running threshold ${THRESHOLD}"
    echo "========================================"

    TMP_DIR="${TMP_ROOT}/t${THRESHOLD}"
    mkdir -p "${TMP_DIR}"

    CLUSTER_DB="${DB_DIR}/cluster_t${THRESHOLD}"
    CLUSTER_TSV="${OUTPUT_DIR}/cluster_t${THRESHOLD}.tsv"

    foldseek cluster \
        "${DB}" \
        "${CLUSTER_DB}" \
        "${TMP_DIR}" \
        --min-seq-id "${THRESHOLD}" \
        -c "${COVERAGE}" \
        --alignment-type "${ALIGNMENT_TYPE}"

    foldseek createtsv \
        "${DB}" \
        "${DB}" \
        "${CLUSTER_DB}" \
        "${CLUSTER_TSV}"

    echo "Saved cluster TSV:"
    echo "${CLUSTER_TSV}"
done

echo ""
echo "========================================"
echo "All Foldseek clustering runs finished"
echo "========================================"

echo "Cluster TSV files:"
ls -lh "${OUTPUT_DIR}"/*.tsv

echo ""
echo "Foldseek DB files are stored in:"
echo "${DB_DIR}"

echo ""
echo "Temporary files are stored in:"
echo "${TMP_ROOT}"