#!/bin/bash
#SBATCH -t 12:00:00
#SBATCH --mem=64G
#SBATCH --account=hansengroup
#SBATCH -J transdecoder_pmlb
#SBATCH -p all
#SBATCH -c 4
#SBATCH -N 1
#SBATCH --mail-user=vivian.chu@mail.utoronto.ca
#SBATCH --mail-type=ALL
#SBATCH -o /cluster/projects/livingbank/workspace/vivian/neo/novel/logs/transdecoder_%A_%a.out

# Run TransDecoder ORF prediction as a SLURM array job.
# Submit with: sbatch --array=0-<N-1> run_transdecoder.sh
# where N = number of *_novel_sequences.fasta files in novel-dir.
#
# Usage:
#   sbatch --array=0-<N-1> run_transdecoder.sh [--novel-dir <path>] [--out-dir <path>] [--min-aa <int>] [--dry-run]
#
# Defaults:
#   --novel-dir  /cluster/projects/livingbank/workspace/vivian/neo/novel/novel_transcripts
#   --out-dir    /cluster/projects/livingbank/workspace/vivian/neo/novel/novel_translated
#   --min-aa     25

set -euo pipefail

echo "======================================"
echo "Job started: $(date)"
echo "Host: $(hostname)"
echo "Array Task ID: ${SLURM_ARRAY_TASK_ID:-0}"
echo "======================================"

#########################
# Defaults
#########################

NOVEL_DIR="/cluster/projects/livingbank/workspace/vivian/neo/novel/novel_transcripts"
OUT_DIR="/cluster/projects/livingbank/workspace/vivian/neo/novel/novel_translated"
MIN_AA=25
DRY_RUN=false

#########################
# Argument parsing
#########################

while [[ $# -gt 0 ]]; do
    case "$1" in
        --novel-dir) NOVEL_DIR="$2"; shift 2 ;;
        --out-dir)   OUT_DIR="$2";   shift 2 ;;
        --min-aa)    MIN_AA="$2";    shift 2 ;;
        --dry-run)   DRY_RUN=true;   shift ;;
        *) echo "Unknown argument: $1" >&2; exit 1 ;;
    esac
done

if [[ ! -d "${NOVEL_DIR}" ]]; then
    echo "ERROR: novel-dir not found: ${NOVEL_DIR}" >&2
    exit 1
fi

mkdir -p "${OUT_DIR}"
mkdir -p "/cluster/projects/livingbank/workspace/vivian/neo/novel/logs"

#########################
# Select this task's FASTA
#
# Supports both:
#   NOVEL_DIR/<SAMPLE>_novel_sequences.fasta
#   NOVEL_DIR/<SAMPLE>/<SAMPLE>_novel_sequences.fasta
#########################

FASTAS=()
while IFS= read -r fasta; do
    FASTAS+=("${fasta}")
done < <(
    find "${NOVEL_DIR}" -type f -name "*_novel_sequences.fasta" \
        ! -path "${OUT_DIR}/*" \
        ! -path "${NOVEL_DIR}/logs/*" \
        | sort
)

if [[ ${#FASTAS[@]} -eq 0 ]]; then
    echo "ERROR: No *_novel_sequences.fasta files found under ${NOVEL_DIR}" >&2
    exit 1
fi

TASK_ID="${SLURM_ARRAY_TASK_ID:-0}"
if [[ ${TASK_ID} -ge ${#FASTAS[@]} ]]; then
    echo "ERROR: TASK_ID ${TASK_ID} out of range (${#FASTAS[@]} samples)" >&2
    exit 1
fi

FASTA="${FASTAS[${TASK_ID}]}"
SAMPLE=$(basename "${FASTA}" _novel_sequences.fasta)

echo "novel-dir:  ${NOVEL_DIR}"
echo "out-dir:    ${OUT_DIR}"
echo "min-aa:     ${MIN_AA}"
echo "sample:     ${SAMPLE} (task ${TASK_ID} of ${#FASTAS[@]})"
echo "dry-run:    ${DRY_RUN}"
echo ""

#########################
# Run TransDecoder
#########################

SEQ_COUNT=$(grep -c "^>" "${FASTA}" 2>/dev/null) || SEQ_COUNT=0

if [[ ${SEQ_COUNT} -eq 0 ]]; then
    echo "[SKIP] ${SAMPLE}: empty FASTA"
    exit 0
fi

echo "--- ${SAMPLE} (${SEQ_COUNT} sequences) ---"

WORKDIR="${OUT_DIR}/${SAMPLE}_transdecoder_wd"
mkdir -p "${WORKDIR}"

ABS_FASTA=$(realpath "${FASTA}")
BASENAME=$(basename "${FASTA}")
ABS_PEP="${OUT_DIR}/${BASENAME}.transdecoder.pep"

if [[ "${DRY_RUN}" == "true" ]]; then
    echo "[dry-run] Would run TransDecoder on ${ABS_FASTA}"
    echo "[dry-run] Would write ${ABS_PEP}"
    exit 0
fi

(
    cd "${WORKDIR}"

    TransDecoder.LongOrfs \
        -t "${ABS_FASTA}" \
        -m "${MIN_AA}" \
        --output_dir . \
        2>&1

    TransDecoder.Predict \
        -t "${ABS_FASTA}" \
        --single_best_only \
        --output_dir . \
        2>&1

    LOCAL_PEP="${BASENAME}.transdecoder.pep"

    if [[ -f "${LOCAL_PEP}" ]]; then
        cp "${LOCAL_PEP}" "${ABS_PEP}"
        ORF_COUNT=$(grep -c "^>" "${ABS_PEP}" 2>/dev/null) || ORF_COUNT=0
        echo "[OK]   ${SAMPLE}: ${ORF_COUNT} ORFs → $(basename "${ABS_PEP}")"
    else
        echo "[WARN] ${SAMPLE}: TransDecoder produced no .pep output"
    fi
)

echo ""
echo "======================================"
echo "Job finished: $(date)"
echo "======================================"
