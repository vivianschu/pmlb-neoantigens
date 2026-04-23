#!/bin/bash
#SBATCH -t 06:00:00
#SBATCH --mem=15G
#SBATCH --account=hansengroup
#SBATCH -J hla_pmlb
#SBATCH -p all
#SBATCH -c 12
#SBATCH -N 1
#SBATCH --mail-user=vivian.chu@mail.utoronto.ca
#SBATCH --mail-type=ALL
#SBATCH -o /cluster/projects/livingbank/workspace/vivian/neo/hla/logs/hla_%A_%a.out

set -euo pipefail

echo "======================================"
echo "Job started: $(date)"
echo "Host: $(hostname)"
echo "Array Task ID: ${SLURM_ARRAY_TASK_ID:-0}"
echo "======================================"

#########################
# Configuration
#########################

module load kallisto
module load samtools

# Input: processed BAM files from process.sh
SAMPLES_DIR="/cluster/projects/livingbank/workspace/vivian/RNA/process"
OUTPUT_DIR="/cluster/projects/livingbank/workspace/vivian/neo/hla"
THREADS=12

# Load samples dynamically from manifest
MANIFEST="/cluster/projects/livingbank/workspace/vivian/neo/260313_manifest_final.tsv"
if [[ ! -f "${MANIFEST}" ]]; then
    echo "ERROR: Manifest not found: ${MANIFEST}" >&2
    exit 1
fi
mapfile -t SAMPLES < <(tail -n +2 "${MANIFEST}" | cut -f1)

mkdir -p "${OUTPUT_DIR}"
mkdir -p "${OUTPUT_DIR}/logs"

#########################
# Get sample for this task
#########################

TASK_ID=${SLURM_ARRAY_TASK_ID:-0}

if [[ ${TASK_ID} -ge ${#SAMPLES[@]} ]]; then
    echo "ERROR: Task ID ${TASK_ID} exceeds number of samples (${#SAMPLES[@]})" >&2
    exit 1
fi

SAMPLE="${SAMPLES[${TASK_ID}]}"

echo "Processing sample ${TASK_ID}/${#SAMPLES[@]}: ${SAMPLE}"

#########################
# Setup paths
#########################

SAMPLE_INPUT="${SAMPLES_DIR}/${SAMPLE}"
SAMPLE_OUTPUT="${OUTPUT_DIR}/${SAMPLE}"
INPUT_BAM="${SAMPLE_INPUT}/${SAMPLE}_Aligned.sortedByCoord.out.bam"
GENOTYPE_JSON="${SAMPLE_OUTPUT}/${SAMPLE}.genotype.json"

# Check if already completed
if [[ -f "${GENOTYPE_JSON}" ]]; then
    echo "SKIP: Genotype file already exists"
    exit 0
fi

# Validate input
if [[ ! -f "${INPUT_BAM}" ]]; then
    echo "ERROR: BAM file not found: ${INPUT_BAM}" >&2
    exit 1
fi

mkdir -p "${SAMPLE_OUTPUT}"

#########################
# Step 1: Extract HLA reads
#########################

HLA_CLASS_I="A,B,C,E,F,G"
HLA_CLASS_II="DRA,DPA1,DPB1,DQA1,DQB1,DRB1,DRB3,DRB4,DRB5,DMA,DMB,DOA,DOB"
HLA_GENES="${HLA_CLASS_I},${HLA_CLASS_II}"

echo "Extracting HLA reads..."
arcasHLA extract "${INPUT_BAM}" -o "${SAMPLE_OUTPUT}" -t "${THREADS}" -v

#########################
# Step 2: Find extracted files
#########################

EXTRACTED_1=$(find "${SAMPLE_OUTPUT}" -name "*.extracted.1.fq.gz" -type f | head -1)
EXTRACTED_2=$(find "${SAMPLE_OUTPUT}" -name "*.extracted.2.fq.gz" -type f | head -1)

#########################
# Step 3: Genotyping
#########################

if [[ -n "${EXTRACTED_1}" && -n "${EXTRACTED_2}" ]]; then
    echo "Running paired-end genotyping..."
    arcasHLA genotype "${EXTRACTED_1}" "${EXTRACTED_2}" \
        --threads "${THREADS}" \
        --genes "${HLA_GENES}" \
        -o "${SAMPLE_OUTPUT}" \
        -v
else
    # Try single-end
    EXTRACTED_SINGLE=$(find "${SAMPLE_OUTPUT}" -name "*.extracted.fq.gz" -type f | head -1)
    if [[ -n "${EXTRACTED_SINGLE}" ]]; then
        echo "Running single-end genotyping..."
        arcasHLA genotype "${EXTRACTED_SINGLE}" \
            --threads "${THREADS}" \
            --genes "${HLA_GENES}" \
            -o "${SAMPLE_OUTPUT}" \
            -v
    else
        echo "ERROR: No extracted FASTQ files found" >&2
        exit 1
    fi
fi

#########################
# Step 4: Verify and rename output
#########################

GENOTYPE_FILE=$(find "${SAMPLE_OUTPUT}" -name "*.genotype.json" -type f | head -1)

if [[ -z "${GENOTYPE_FILE}" || ! -f "${GENOTYPE_FILE}" ]]; then
    echo "ERROR: Genotype file not generated" >&2
    exit 1
fi

# Rename to consistent naming if needed
if [[ "${GENOTYPE_FILE}" != "${GENOTYPE_JSON}" ]]; then
    mv "${GENOTYPE_FILE}" "${GENOTYPE_JSON}"
fi

echo "SUCCESS: Genotype complete"
cat "${GENOTYPE_JSON}"

# Convert to TSV
arcasHLA convert --infile "${GENOTYPE_JSON}" \
    --outfile "${SAMPLE_OUTPUT}/${SAMPLE}.genotype.tsv" \
    --format tsv 2>/dev/null || true

# Clean up intermediate files
rm -f "${SAMPLE_OUTPUT}"/*.extracted.*.fq.gz

echo "Job finished: $(date)"
