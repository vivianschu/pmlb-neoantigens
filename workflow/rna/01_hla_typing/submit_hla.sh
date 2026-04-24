#!/bin/bash
# Submit HLA typing array jobs and optionally merge results
# Usage: ./submit_hla.sh [--dry-run] [--merge-only]

set -euo pipefail

# Input: processed BAM files from process.sh
SAMPLES_DIR="/cluster/projects/livingbank/workspace/vivian/RNA/process"
OUTPUT_DIR="/cluster/projects/livingbank/workspace/vivian/neo/hla"
SCRIPT_PATH="$(dirname "$(realpath "$0")")/hla_typing_array.sh"
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
    echo "Merging HLA typing results..."
    echo ""

    local SUMMARY_FILE="${OUTPUT_DIR}/hla_summary.csv"
    local FAILED_FILE="${OUTPUT_DIR}/failed_samples.txt"
    local MERGED_FILE="${OUTPUT_DIR}/all_genotypes.tsv"

    local -a HLA_CLASS_I=(A B C E F G)
    local -a HLA_CLASS_II=(DRA DPA1 DPB1 DQA1 DQB1 DRB1 DRB3 DRB4 DRB5 DMA DMB DOA DOB)
    local -a ALL_GENES=("${HLA_CLASS_I[@]}" "${HLA_CLASS_II[@]}")

    # Arrays to track status
    local -a SUCCESSFUL_SAMPLES=()
    local -a FAILED_SAMPLES=()
    local -a FAILED_REASONS=()

    # Build CSV header and NA placeholder from gene list
    local HEADER="sample,status,reason"
    local NA_COLS=""
    for gene in "${ALL_GENES[@]}"; do
        gene_lower=$(echo "${gene}" | tr '[:upper:]' '[:lower:]')
        HEADER="${HEADER},${gene_lower}_1,${gene_lower}_2"
        NA_COLS="${NA_COLS},NA,NA"
    done

    # Generate summary CSV and collect status
    {
        echo "${HEADER}"

        for sample in "${SAMPLES[@]}"; do
            sample_output="${OUTPUT_DIR}/${sample}"
            genotype_json="${sample_output}/${sample}.genotype.json"
            input_bam="${SAMPLES_DIR}/${sample}/${sample}_Aligned.sortedByCoord.out.bam"

            if [[ -f "${genotype_json}" ]]; then
                allele_vals=""
                for gene in "${ALL_GENES[@]}"; do
                    vals=$(jq -r --arg g "${gene}" \
                        '[(.[$g] // [])[0] // "NA", (.[$g] // [])[1] // "NA"] | join(",")' \
                        "${genotype_json}" 2>/dev/null || echo "NA,NA")
                    allele_vals="${allele_vals},${vals}"
                done
                echo "${sample},success,completed${allele_vals}"
                SUCCESSFUL_SAMPLES+=("${sample}")
            else
                # Diagnose failure reason
                local reason="unknown"

                if [[ ! -f "${input_bam}" ]]; then
                    reason="missing_bam"
                elif [[ ! -d "${sample_output}" ]]; then
                    reason="not_started"
                elif ls "${sample_output}"/*.extracted.*.fq.gz &>/dev/null; then
                    reason="genotyping_failed"
                elif [[ -f "${sample_output}/${sample}.alignment.p" ]]; then
                    reason="genotyping_failed"
                else
                    # Check slurm logs for errors
                    local log_pattern="${LOG_DIR}/hla_*_*.out"
                    if compgen -G "${log_pattern}" > /dev/null; then
                        if grep -l "ERROR.*${sample}" ${log_pattern} &>/dev/null; then
                            reason="error_in_log"
                        elif grep -l "OOM\|memory" ${log_pattern} &>/dev/null; then
                            reason="out_of_memory"
                        fi
                    fi
                    [[ "${reason}" == "unknown" ]] && reason="extraction_failed"
                fi

                echo "${sample},failed,${reason}${NA_COLS}"
                FAILED_SAMPLES+=("${sample}")
                FAILED_REASONS+=("${reason}")
            fi
        done
    } > "${SUMMARY_FILE}"

    # Print summary statistics
    local TOTAL=$((${#SUCCESSFUL_SAMPLES[@]} + ${#FAILED_SAMPLES[@]}))
    echo "======================================"
    echo "HLA Typing Summary"
    echo "======================================"
    echo "Total samples:  ${TOTAL}"
    echo "Successful:     ${#SUCCESSFUL_SAMPLES[@]}"
    echo "Failed/Pending: ${#FAILED_SAMPLES[@]}"
    echo ""

    # List failed samples with reasons
    if [[ ${#FAILED_SAMPLES[@]} -gt 0 ]]; then
        echo "Failed/Pending Samples:"
        echo "--------------------------------------"
        printf "%-30s %s\n" "SAMPLE" "REASON"
        echo "--------------------------------------"

        # Also write to file
        > "${FAILED_FILE}"

        for i in "${!FAILED_SAMPLES[@]}"; do
            printf "%-30s %s\n" "${FAILED_SAMPLES[$i]}" "${FAILED_REASONS[$i]}"
            echo "${FAILED_SAMPLES[$i]},${FAILED_REASONS[$i]}" >> "${FAILED_FILE}"
        done

        echo "--------------------------------------"
        echo ""
        echo "Failed samples list written to: ${FAILED_FILE}"

        # Group by reason
        echo ""
        echo "Failure breakdown:"
        sort "${FAILED_FILE}" | cut -d',' -f2 | uniq -c | sort -rn
    fi

    echo ""
    echo "Summary written to: ${SUMMARY_FILE}"

    # Try arcasHLA merge
    if command -v arcasHLA &>/dev/null && [[ ${#SUCCESSFUL_SAMPLES[@]} -gt 0 ]]; then
        arcasHLA merge --indir "${OUTPUT_DIR}" --outfile "${MERGED_FILE}" 2>/dev/null && \
            echo "Merged file written to: ${MERGED_FILE}" || \
            echo "Note: arcasHLA merge failed"
    fi
}

#########################
# Main
#########################

if [[ "${1:-}" == "--merge-only" ]]; then
    merge_results
    exit 0
fi

# Sample count from fixed list
SAMPLE_COUNT=${#SAMPLES[@]}
MAX_IDX=$((SAMPLE_COUNT - 1))

echo "Found ${SAMPLE_COUNT} PMLB samples (array indices 0-${MAX_IDX})"
echo "Samples: ${SAMPLES[*]}"

# Check for already completed samples
COMPLETED=0
for sample in "${SAMPLES[@]}"; do
    if [[ -f "${OUTPUT_DIR}/${sample}/${sample}.genotype.json" ]]; then
        COMPLETED=$((COMPLETED + 1))
    fi
done
echo "Already completed: ${COMPLETED}/${SAMPLE_COUNT}"

# Check for input BAM files
echo ""
echo "Checking input BAM files:"
for sample in "${SAMPLES[@]}"; do
    bam="${SAMPLES_DIR}/${sample}/${sample}_Aligned.sortedByCoord.out.bam"
    if [[ -f "${bam}" ]]; then
        echo "  [OK] ${sample}"
    else
        echo "  [MISSING] ${sample}: ${bam}"
    fi
done

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
