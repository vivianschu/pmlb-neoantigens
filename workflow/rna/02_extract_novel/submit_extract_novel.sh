#!/bin/bash
# Submit novel transcript extraction jobs for neo-SCLC samples
# Usage: ./submit_extract_novel_neosclc.sh [--dry-run] [--merge-only]

set -euo pipefail

SAMPLES_DIR="/cluster/projects/livingbank/workspace/vivian/RNA/process"
OUTPUT_DIR="/cluster/projects/livingbank/workspace/vivian/neo/novel"
SCRIPT_PATH="$(dirname "$(realpath "$0")")/extract_novel.sh"
LOG_DIR="${OUTPUT_DIR}/logs"

# Load samples dynamically from manifest
MANIFEST="/cluster/projects/livingbank/workspace/vivian/neo/260313_manifest_final.tsv"
if [[ ! -f "${MANIFEST}" ]]; then
    echo "ERROR: Manifest not found: ${MANIFEST}" >&2
    exit 1
fi
mapfile -t SAMPLES < <(tail -n +2 "${MANIFEST}" | cut -f1)

mkdir -p "${OUTPUT_DIR}" "${LOG_DIR}"

#########################
# Merge function
#########################

merge_results() {
    echo "Collecting novel transcript results..."
    echo ""

    local SUMMARY_FILE="${OUTPUT_DIR}/novel_transcripts_summary.csv"
    local NOVEL_TRANSCRIPTS_DIR="${OUTPUT_DIR}/novel_transcripts"

    mkdir -p "${NOVEL_TRANSCRIPTS_DIR}"

    # Generate summary and copy per-sample files
    {
        echo "sample,total_novel,class_u,class_n,class_j,sequences_extracted,orfs_predicted"

        for sample in "${SAMPLES[@]}"; do
            sample_dir="${OUTPUT_DIR}/${sample}"

            [[ ! -d "${sample_dir}" ]] && continue

            novel_gtf="${sample_dir}/${sample}_novel_transcripts.gtf"
            novel_fasta="${sample_dir}/${sample}_novel_sequences.fasta"
            td_pep="${novel_fasta}.transdecoder.pep"

            if [[ -f "${novel_gtf}" ]]; then
                class_u=$(grep -c 'class_code "u"' "${novel_gtf}" 2>/dev/null) || class_u=0
                class_n=$(grep -c 'class_code "n"' "${novel_gtf}" 2>/dev/null) || class_n=0
                class_j=$(grep -c 'class_code "j"' "${novel_gtf}" 2>/dev/null) || class_j=0
                total=$(grep -c $'\ttranscript\t' "${novel_gtf}" 2>/dev/null) || total=0

                seqs=0
                [[ -f "${novel_fasta}" ]] && { seqs=$(grep -c "^>" "${novel_fasta}" 2>/dev/null) || seqs=0; }

                orfs=0
                [[ -f "${td_pep}" ]] && { orfs=$(grep -c "^>" "${td_pep}" 2>/dev/null) || orfs=0; }

                # Copy per-sample transcript files to novel_transcripts/
                cp "${novel_gtf}" "${NOVEL_TRANSCRIPTS_DIR}/${sample}_novel_transcripts.gtf"
                [[ -f "${novel_fasta}" ]] && cp "${novel_fasta}" "${NOVEL_TRANSCRIPTS_DIR}/${sample}_novel_sequences.fasta"
                [[ -f "${td_pep}" ]] && cp "${td_pep}" "${NOVEL_TRANSCRIPTS_DIR}/${sample}_novel_sequences.fasta.transdecoder.pep"

                echo "${sample},${total},${class_u},${class_n},${class_j},${seqs},${orfs}"
            fi
        done
    } > "${SUMMARY_FILE}"

    # Calculate totals
    local TOTAL_SAMPLES=$(tail -n +2 "${SUMMARY_FILE}" | wc -l)
    local TOTAL_TRANSCRIPTS=$(tail -n +2 "${SUMMARY_FILE}" | awk -F',' '{sum+=$2} END {print sum}')
    local TOTAL_SEQS=$(tail -n +2 "${SUMMARY_FILE}" | awk -F',' '{sum+=$6} END {print sum}')
    local TOTAL_ORFS=$(tail -n +2 "${SUMMARY_FILE}" | awk -F',' '{sum+=$7} END {print sum}')
    local TD_SAMPLES=$(tail -n +2 "${SUMMARY_FILE}" | awk -F',' '$7>0' | wc -l)

    echo "======================================"
    echo "Novel Transcript Summary (PMLB)"
    echo "======================================"
    echo "Samples processed:    ${TOTAL_SAMPLES}"
    echo "Total transcripts:    ${TOTAL_TRANSCRIPTS}"
    echo "Total sequences:      ${TOTAL_SEQS}"
    echo "Total ORFs predicted: ${TOTAL_ORFS} (${TD_SAMPLES} samples with TransDecoder output)"
    echo ""
    echo "Class code breakdown (across all samples):"
    echo "  u (intergenic):     $(tail -n +2 "${SUMMARY_FILE}" | awk -F',' '{sum+=$3} END {print sum}')"
    echo "  n (intronic):       $(tail -n +2 "${SUMMARY_FILE}" | awk -F',' '{sum+=$4} END {print sum}')"
    echo "  j (novel junction): $(tail -n +2 "${SUMMARY_FILE}" | awk -F',' '{sum+=$5} END {print sum}')"
    echo ""
    echo "Output files:"
    echo "  Summary:              ${SUMMARY_FILE}"
    echo "  Novel transcripts:    ${NOVEL_TRANSCRIPTS_DIR}/ (${TOTAL_SAMPLES} samples)"
    echo "======================================"
}

#########################
# Main
#########################

if [[ "${1:-}" == "--merge-only" ]]; then
    merge_results
    exit 0
fi

SAMPLE_COUNT=${#SAMPLES[@]}
MAX_IDX=$((SAMPLE_COUNT - 1))

echo "PMLB Novel Transcript Extraction"
echo "====================================="
echo "Samples: ${SAMPLES[*]}"
echo "Sample count: ${SAMPLE_COUNT}"
echo ""

# Check for StringTie outputs (prerequisite for gffcompare)
READY_COUNT=0
for sample in "${SAMPLES[@]}"; do
    if [[ -f "${SAMPLES_DIR}/${sample}/${sample}_annotation.gtf" ]]; then
        ((READY_COUNT++)) || true
    fi
done

echo "Samples with StringTie output: ${READY_COUNT}/${SAMPLE_COUNT}"

if [[ ${READY_COUNT} -eq 0 ]]; then
    echo "ERROR: No StringTie outputs found. Run process_rnaseq.sh first." >&2
    exit 1
fi

if [[ "${1:-}" == "--dry-run" ]]; then
    echo ""
    echo "Dry run - would execute:"
    echo "sbatch --array=0-${MAX_IDX} ${SCRIPT_PATH}"
else
    echo ""
    echo "Submitting job array..."
    JOB_ID=$(sbatch --array=0-${MAX_IDX} "${SCRIPT_PATH}" | awk '{print $4}')
    echo "Submitted job array: ${JOB_ID}"
    echo ""
    echo "Monitor with: squeue -j ${JOB_ID}"
    echo "After completion, run: $0 --merge-only"
fi
