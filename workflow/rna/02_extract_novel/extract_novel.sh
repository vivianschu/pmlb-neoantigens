#!/bin/bash
#SBATCH -t 12:00:00
#SBATCH --mem=16G
#SBATCH --account=hansengroup
#SBATCH -J extract_novel_pmlb
#SBATCH -p all
#SBATCH -c 4
#SBATCH -N 1
#SBATCH --mail-user=vivian.chu@mail.utoronto.ca
#SBATCH --mail-type=ALL
#SBATCH -o /cluster/projects/livingbank/workspace/vivian/neo/novel/logs/novel_%A_%a.out

set -euo pipefail

echo "======================================"
echo "Job started: $(date)"
echo "Host: $(hostname)"
echo "Array Task ID: ${SLURM_ARRAY_TASK_ID:-0}"
echo "======================================"

#########################
# Configuration
#########################

module load gffread
module load gffcompare

PATH="${PATH}:/cluster/home/t117036uhn/local_soft/stringtie2:/cluster/home/t117036uhn/local_soft/gffcompare"

# Reference
REF_DIR="/cluster/projects/livingbank/workspace/references/hg38_ek12"
REF_FASTA="${REF_DIR}/hg38_ek12.fa"
REF_GTF="/cluster/projects/livingbank/workspace/references/hg_38/gencode.v25.annotation.gtf"

# Directories
INPUT_BASE="/cluster/projects/livingbank/workspace/vivian/RNA/process"
OUTPUT_BASE="/cluster/projects/livingbank/workspace/vivian/neo/novel"

# Class codes to extract:
#   u = intergenic (novel locus, no overlap with reference)
#   n = intronic (novel transcript within known gene intron)
#   j = novel junction (new splice variant of known gene)
CLASS_CODES="u|n|j"

# Load samples dynamically from manifest
MANIFEST="/cluster/projects/livingbank/workspace/vivian/neo/260313_manifest_final.tsv"
if [[ ! -f "${MANIFEST}" ]]; then
    echo "ERROR: Manifest not found: ${MANIFEST}" >&2
    exit 1
fi
mapfile -t SAMPLES < <(tail -n +2 "${MANIFEST}" | cut -f1)

TASK_ID=${SLURM_ARRAY_TASK_ID:-0}
SAMPLE="${SAMPLES[${TASK_ID}]}"

echo "Processing sample ${TASK_ID}/${#SAMPLES[@]}: ${SAMPLE}"

#########################
# Setup paths
#########################

SAMPLE_INPUT="${INPUT_BASE}/${SAMPLE}"
SAMPLE_OUTPUT="${OUTPUT_BASE}/${SAMPLE}"
STRINGTIE_GTF="${SAMPLE_INPUT}/${SAMPLE}_annotation.gtf"

# Validate input exists
if [[ ! -f "${STRINGTIE_GTF}" ]]; then
    echo "ERROR: StringTie output not found: ${STRINGTIE_GTF}" >&2
    echo "Run process_rnaseq.sh first to generate this file." >&2
    exit 1
fi

mkdir -p "${SAMPLE_OUTPUT}"

echo "Input GTF: ${STRINGTIE_GTF}"
echo "Reference GTF: ${REF_GTF}"
echo "Reference FASTA: ${REF_FASTA}"
echo "Extracting class codes: ${CLASS_CODES}"

#########################
# Run gffcompare
#########################

GFFCOMPARE_PREFIX="${SAMPLE_OUTPUT}/${SAMPLE}_gffcmp"
GFFCOMPARE_GTF="${GFFCOMPARE_PREFIX}.annotated.gtf"

echo ""
echo "Running gffcompare to annotate transcripts..."

gffcompare -r "${REF_GTF}" -o "${GFFCOMPARE_PREFIX}" "${STRINGTIE_GTF}"

echo "gffcompare output: ${GFFCOMPARE_GTF}"

#########################
# Filter novel transcripts
#########################

FILTERED_GTF="${SAMPLE_OUTPUT}/${SAMPLE}_novel_transcripts.gtf"

echo ""
echo "Filtering transcripts with class_code (${CLASS_CODES})..."

grep -P "class_code \"(${CLASS_CODES})\"" "${GFFCOMPARE_GTF}" > "${FILTERED_GTF}" || true

# Count transcripts per class code
# Note: grep -c returns exit code 1 when no matches, so use || assignment not || echo
NOVEL_COUNT=$(grep -c 'class_code "u"' "${FILTERED_GTF}" 2>/dev/null) || NOVEL_COUNT=0
INTRONIC_COUNT=$(grep -c 'class_code "n"' "${FILTERED_GTF}" 2>/dev/null) || INTRONIC_COUNT=0
JUNCTION_COUNT=$(grep -c 'class_code "j"' "${FILTERED_GTF}" 2>/dev/null) || JUNCTION_COUNT=0
TOTAL_COUNT=$((NOVEL_COUNT + INTRONIC_COUNT + JUNCTION_COUNT))

echo "Found transcripts:"
echo "  - Class 'u' (intergenic):     ${NOVEL_COUNT}"
echo "  - Class 'n' (intronic):       ${INTRONIC_COUNT}"
echo "  - Class 'j' (novel junction): ${JUNCTION_COUNT}"
echo "  - Total:                      ${TOTAL_COUNT}"

if [[ ${TOTAL_COUNT} -eq 0 ]]; then
    echo "WARNING: No novel transcripts found for ${SAMPLE}"
    touch "${SAMPLE_OUTPUT}/${SAMPLE}_novel_sequences.fasta"
    exit 0
fi

#########################
# Extract sequences
#########################

NOVEL_FASTA="${SAMPLE_OUTPUT}/${SAMPLE}_novel_sequences.fasta"

echo "Extracting transcript sequences..."

gffread -w "${NOVEL_FASTA}" -g "${REF_FASTA}" "${FILTERED_GTF}"

# Count extracted sequences
SEQ_COUNT=$(grep -c "^>" "${NOVEL_FASTA}" 2>/dev/null) || SEQ_COUNT=0

echo "Extracted ${SEQ_COUNT} sequences to ${NOVEL_FASTA}"

if [[ ${SEQ_COUNT} -eq 0 && ${TOTAL_COUNT} -gt 0 ]]; then
    echo "ERROR: gffread extracted 0 sequences from ${TOTAL_COUNT} transcripts." >&2
    echo "  Check chromosome name compatibility between:" >&2
    echo "    GTF:   $(cut -f1 "${FILTERED_GTF}" | sort -u | tr '\n' ' ')" >&2
    echo "    FASTA: $(grep '^>' "${REF_FASTA}" | head -3 | sed 's/^>//' | cut -d' ' -f1 | tr '\n' ' ')" >&2
    exit 1
fi

#########################
# Summary
#########################

echo "======================================"
echo "Job completed: $(date)"
echo "Output files:"
ls -lh "${FILTERED_GTF}" "${NOVEL_FASTA}"
echo "Sequences: ${SEQ_COUNT}"
echo "======================================"
echo "NOTE: Run workflow/novel/run_transdecoder.sh locally after transferring novel_transcripts/"
