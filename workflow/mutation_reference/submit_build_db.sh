#!/bin/bash
# =============================================================================
# submit_build_db.sh
#
# Top-level submission wrapper for the FragPipe reference database pipeline.
#
# Stages and their SLURM dependencies:
#
#   [01] build_mutated_db    single job  → neo/mutation_reference/mutated_proteins/
#   [02] novel_orfs          array job   → neo/mutation_reference/novel_proteins/
#                            (runs in parallel with 01; both depend on nothing)
#   [03] merge_fragpipe_db   single job  → neo/mutation_reference/final/
#                            (depends on 01 AND 02)
#
# Modes:
#   default          submit 01 only (mutated proteins)
#   --with-novel     submit 01 + 02 + 03 as a chained pipeline
#   --novel-only     submit 02 + 03 (skip mutation calling)
#   --merge-only     submit 03 only (both 01 and 02 already done)
#   --dry-run        print sbatch commands without submitting
#
# Usage:
#   bash submit_build_db.sh --with-novel
#   bash submit_build_db.sh --with-novel \
#       --sample-manifest /cluster/projects/livingbank/workspace/vivian/neo/260313_manifest_final.tsv
#   bash submit_build_db.sh --merge-only
#   bash submit_build_db.sh --dry-run --with-novel
#
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# ---------------------------------------------------------------------------
# Cluster paths  (edit ONLY these lines)
# ---------------------------------------------------------------------------
NEO_BASE="/cluster/projects/livingbank/workspace/vivian/neo"
NOVEL_DIR="${NEO_BASE}/novel/novel_transcripts"
NOVEL_PROTEINS_DIR="${NEO_BASE}/mutation_reference/novel_proteins"
NOVEL_PARTITION="himem"
NOVEL_MEM="128G"
NOVEL_TIME="24:00:00"
FINAL_DIR="${NEO_BASE}/mutation_reference/final"

ORF_FASTA="${NEO_BASE}/mutation_reference/orf/cdna_orfs.fasta"
MUTATIONS="/cluster/projects/livingbank/Project/Pan-organoid/Mutation/Org_exome.data_mutations_extended.gt4.202507.txt"
MUTATED_DIR="${NEO_BASE}/mutation_reference/mutated_proteins"

# ---------------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------------
SAMPLE_MANIFEST=""
EXPRESSION_TSV=""
RNA_VARIANTS_TSV=""
HLA_TSV=""
BINDING_TSV=""
LABELS_TSV=""
WITH_NOVEL=false
NOVEL_ONLY=false
MERGE_ONLY=false
DRY_RUN=false
EXTRA_ARGS=""

# ---------------------------------------------------------------------------
# Parse arguments
# ---------------------------------------------------------------------------
while [[ $# -gt 0 ]]; do
    case "$1" in
        --orf-fasta)          ORF_FASTA="$2";        shift 2 ;;
        --mutations)          MUTATIONS="$2";         shift 2 ;;
        --sample-manifest)    SAMPLE_MANIFEST="$2";   shift 2 ;;
        --expression-tsv)     EXPRESSION_TSV="$2"; EXTRA_ARGS="${EXTRA_ARGS} --expression-tsv $2"; shift 2 ;;
        --expression-sample-column) EXTRA_ARGS="${EXTRA_ARGS} --expression-sample-column $2"; shift 2 ;;
        --expression-gene-column) EXTRA_ARGS="${EXTRA_ARGS} --expression-gene-column $2"; shift 2 ;;
        --expression-transcript-column) EXTRA_ARGS="${EXTRA_ARGS} --expression-transcript-column $2"; shift 2 ;;
        --expression-tpm-column) EXTRA_ARGS="${EXTRA_ARGS} --expression-tpm-column $2"; shift 2 ;;
        --rna-variants-tsv)   RNA_VARIANTS_TSV="$2"; EXTRA_ARGS="${EXTRA_ARGS} --rna-variants-tsv $2"; shift 2 ;;
        --hla-tsv)            HLA_TSV="$2"; EXTRA_ARGS="${EXTRA_ARGS} --hla-tsv $2"; shift 2 ;;
        --hla-sample-column)  EXTRA_ARGS="${EXTRA_ARGS} --hla-sample-column $2"; shift 2 ;;
        --hla-allele-column)  EXTRA_ARGS="${EXTRA_ARGS} --hla-allele-column $2"; shift 2 ;;
        --binding-tsv)        BINDING_TSV="$2"; EXTRA_ARGS="${EXTRA_ARGS} --binding-tsv $2"; shift 2 ;;
        --labels-tsv)         LABELS_TSV="$2"; EXTRA_ARGS="${EXTRA_ARGS} --labels-tsv $2"; shift 2 ;;
        --min-tpm)            EXTRA_ARGS="${EXTRA_ARGS} --min-tpm $2"; shift 2 ;;
        --min-rna-vaf)        EXTRA_ARGS="${EXTRA_ARGS} --min-rna-vaf $2"; shift 2 ;;
        --min-rna-depth)      EXTRA_ARGS="${EXTRA_ARGS} --min-rna-depth $2"; shift 2 ;;
        --require-rna-filters) EXTRA_ARGS="${EXTRA_ARGS} --require-rna-filters"; shift ;;
        --binder-ic50-threshold) EXTRA_ARGS="${EXTRA_ARGS} --binder-ic50-threshold $2"; shift 2 ;;
        --binder-rank-threshold) EXTRA_ARGS="${EXTRA_ARGS} --binder-rank-threshold $2"; shift 2 ;;
        --compress-tables)    EXTRA_ARGS="${EXTRA_ARGS} --compress-tables"; shift ;;
        --with-novel)         WITH_NOVEL=true;        shift   ;;
        --novel-only)         NOVEL_ONLY=true;        shift   ;;
        --merge-only)         MERGE_ONLY=true;        shift   ;;
        --dry-run)            DRY_RUN=true;           shift   ;;
        --include-reference)  EXTRA_ARGS="${EXTRA_ARGS} --include-reference"; shift ;;
        --allow-ref-mismatch) EXTRA_ARGS="${EXTRA_ARGS} --allow-ref-mismatch"; shift ;;
        *) echo "Unknown argument: $1" >&2; exit 1 ;;
    esac
done

mkdir -p logs
DEFAULT_HLA_SUMMARY="$(cd "${SCRIPT_DIR}/../.." && pwd)/data/rna/hla_summary.csv"
if [[ -z "${HLA_TSV}" && -f "${DEFAULT_HLA_SUMMARY}" ]]; then
    HLA_TSV="${DEFAULT_HLA_SUMMARY}"
    EXTRA_ARGS="${EXTRA_ARGS} --hla-tsv ${HLA_TSV} --hla-sample-column sample"
fi

# ---------------------------------------------------------------------------
# Helper: run or echo an sbatch command; returns the job ID
# ---------------------------------------------------------------------------
_sbatch() {
    if [[ "${DRY_RUN}" == "true" ]]; then
        echo "[dry-run] sbatch $*" >&2
        echo "DRY_RUN_JOB_ID"
        return 0
    fi
    local job_name
    job_name=$(printf '%s\n' "$@" | grep -m1 '^--job-name=' | cut -d= -f2)
    job_name="${job_name:-unknown}"
    if ! sbatch --parsable "$@"; then
        echo "[submit_build_db.sh] ERROR: sbatch submission failed for job '${job_name}'" >&2
        exit 1
    fi
}

# ---------------------------------------------------------------------------
# Validate common inputs (skip for merge-only / novel-only)
# ---------------------------------------------------------------------------
if [[ "${DRY_RUN}" == "false" && "${MERGE_ONLY}" == "false" && "${NOVEL_ONLY}" == "false" ]]; then
    [[ ! -f "${ORF_FASTA}"  ]] \
        && echo "[submit_build_db.sh] ERROR: ORF FASTA not found: ${ORF_FASTA}" >&2 \
        && echo "  → Run 00_fetch_cdna_orfs.sh first." >&2 && exit 1
    [[ ! -f "${MUTATIONS}"  ]] \
        && echo "[submit_build_db.sh] ERROR: Mutations file not found: ${MUTATIONS}" >&2 && exit 1
fi
if [[ "${DRY_RUN}" == "false" && -n "${SAMPLE_MANIFEST}" && ! -f "${SAMPLE_MANIFEST}" ]]; then
    echo "[submit_build_db.sh] ERROR: Manifest not found: ${SAMPLE_MANIFEST}" >&2; exit 1
fi
for optional_file in "${EXPRESSION_TSV}" "${RNA_VARIANTS_TSV}" "${HLA_TSV}" "${BINDING_TSV}" "${LABELS_TSV}"; do
    if [[ "${DRY_RUN}" == "false" && -n "${optional_file}" && ! -f "${optional_file}" ]]; then
        echo "[submit_build_db.sh] ERROR: optional input not found: ${optional_file}" >&2; exit 1
    fi
done

echo "======================================"
echo " FragPipe reference DB pipeline"
[[ "${WITH_NOVEL}"  == "true" ]] && echo " Mode: mutated + novel (full)"
[[ "${NOVEL_ONLY}"  == "true" ]] && echo " Mode: novel only"
[[ "${MERGE_ONLY}"  == "true" ]] && echo " Mode: merge only"
[[ "${WITH_NOVEL}"  == "false" && "${NOVEL_ONLY}" == "false" && "${MERGE_ONLY}" == "false" ]] \
    && echo " Mode: mutated only"
[[ "${DRY_RUN}"     == "true" ]] && echo " *** DRY RUN ***"
echo "======================================"

# ---------------------------------------------------------------------------
# Count novel samples for array range
# ---------------------------------------------------------------------------
_count_novel_samples() {
    find "${NOVEL_DIR}" -type f -name "*_novel_sequences.fasta" \
        ! -path "${NOVEL_DIR}/logs/*" \
        | wc -l
}

MANIFEST_ARG=""
[[ -n "${SAMPLE_MANIFEST}" ]] && MANIFEST_ARG="--sample-manifest ${SAMPLE_MANIFEST}"

# ===========================================================================
# Merge-only mode
# ===========================================================================
if [[ "${MERGE_ONLY}" == "true" ]]; then
    JOB_ID=$(_sbatch \
        --job-name=merge_fragpipe_db \
        --output="${NEO_BASE}/mutation_reference/logs/merge_fragpipe_db_%A_%a.out" \
        --account=hansengroup -p all -N 1 \
        "${SCRIPT_DIR}/03_merge_fragpipe_db.sh" \
            --mutated-dir "${MUTATED_DIR}/per_sample" \
            --novel-dir   "${NOVEL_PROTEINS_DIR}" \
            --outdir      "${FINAL_DIR}" \
            --mode        combined)
    echo "Merge job: ${JOB_ID}"
    echo "  Output → ${FINAL_DIR}/all_samples_combined.fasta"
    exit 0
fi

# ===========================================================================
# Novel-only mode
# ===========================================================================
if [[ "${NOVEL_ONLY}" == "true" ]]; then
    N_NOVEL=$(_count_novel_samples)
    if [[ "${N_NOVEL}" -eq 0 ]]; then
        echo "ERROR: No sample subdirectories found in ${NOVEL_DIR}" >&2; exit 1
    fi

    NOVEL_JOB=$(_sbatch \
        --job-name=novel_orfs \
        --array="0-$((N_NOVEL - 1))" \
        --output="${NEO_BASE}/mutation_reference/logs/novel_orfs_%A_%a.out" \
        --account=hansengroup -p "${NOVEL_PARTITION}" -N 1 \
        --time="${NOVEL_TIME}" \
        --mem="${NOVEL_MEM}" \
        "${SCRIPT_DIR}/02_novel_orfs.sh" \
            --novel-dir      "${NOVEL_DIR}" \
            --out-dir        "${NOVEL_PROTEINS_DIR}")
    echo "Novel-orfs array job: ${NOVEL_JOB} (${N_NOVEL} tasks)"

    MERGE_JOB=$(_sbatch \
        --job-name=merge_fragpipe_db \
        --dependency="afterok:${NOVEL_JOB}" \
        --output="${NEO_BASE}/mutation_reference/logs/merge_fragpipe_db_%A_%a.out" \
        --account=hansengroup -p all -N 1 \
        "${SCRIPT_DIR}/03_merge_fragpipe_db.sh" \
            --mutated-dir "${MUTATED_DIR}/per_sample" \
            --novel-dir   "${NOVEL_PROTEINS_DIR}" \
            --outdir      "${FINAL_DIR}" \
            --mode        novel-only)
    echo "Merge job: ${MERGE_JOB} (depends on ${NOVEL_JOB})"
    exit 0
fi

# ===========================================================================
# Default: submit 01 (mutated proteins)
# ===========================================================================
MUTATED_JOB=$(_sbatch \
    --job-name=build_mutated_db \
    --output="${NEO_BASE}/mutation_reference/logs/build_mutated_db_%A_%a.out" \
    --account=hansengroup -p all -N 1 \
    "${SCRIPT_DIR}/01_build_mutated_db.sh" \
        --script-dir "${SCRIPT_DIR}" \
        --orf-fasta  "${ORF_FASTA}"  \
        --mutations  "${MUTATIONS}"  \
        --outdir     "${MUTATED_DIR}" \
        ${MANIFEST_ARG} \
        ${EXTRA_ARGS})
echo "Mutated-DB job: ${MUTATED_JOB}"
echo "  Output → ${MUTATED_DIR}/"

if [[ "${WITH_NOVEL}" == "false" ]]; then
    exit 0
fi

# ===========================================================================
# --with-novel: also submit 02 (runs in parallel with 01) then 03
# ===========================================================================
N_NOVEL=$(_count_novel_samples)
if [[ "${N_NOVEL}" -eq 0 ]]; then
    echo "WARNING: --with-novel set but no sample directories found in ${NOVEL_DIR}" >&2
    exit 1
fi

NOVEL_JOB=$(_sbatch \
    --job-name=novel_orfs \
    --array="0-$((N_NOVEL - 1))" \
    --output="${NEO_BASE}/mutation_reference/logs/novel_orfs_%A_%a.out" \
    --account=hansengroup -p "${NOVEL_PARTITION}" -N 1 \
    --time="${NOVEL_TIME}" \
    --mem="${NOVEL_MEM}" \
    "${SCRIPT_DIR}/02_novel_orfs.sh" \
        --novel-dir      "${NOVEL_DIR}" \
        --out-dir        "${NOVEL_PROTEINS_DIR}")
echo "Novel-orfs array job: ${NOVEL_JOB} (${N_NOVEL} tasks)"
echo "  Output → ${NOVEL_PROTEINS_DIR}/"

# 03 waits for BOTH 01 and 02
MERGE_JOB=$(_sbatch \
    --job-name=merge_fragpipe_db \
    --dependency="afterok:${MUTATED_JOB}:${NOVEL_JOB}" \
    --output="${NEO_BASE}/mutation_reference/logs/merge_fragpipe_db_%A_%a.out" \
    --account=hansengroup -p all -N 1 \
    "${SCRIPT_DIR}/03_merge_fragpipe_db.sh" \
        --mutated-dir "${MUTATED_DIR}/per_sample" \
        --novel-dir   "${NOVEL_PROTEINS_DIR}" \
        --outdir      "${FINAL_DIR}" \
        --mode        combined)
echo "Merge job: ${MERGE_JOB} (depends on ${MUTATED_JOB} and ${NOVEL_JOB})"
echo "  Final DB → ${FINAL_DIR}/all_samples_combined.fasta"
