#!/bin/bash
# =============================================================================
# 01_build_mutated_db.sh
#
# SLURM job: map sample-specific somatic mutations to cDNA ORFs, translate
# compact variant-centered protein contexts, and write per-sample + combined
# mutant FASTAs for FragPipe.
#
# Inputs (all local to the repo):
#   /cluster/projects/livingbank/workspace/vivian/neo/mutation_reference/orf/cdna_orfs.fasta
#   /cluster/projects/livingbank/Project/Pan-organoid/Mutation/Org_exome.data_mutations_extended.gt4.202507.txt
#
# Outputs:
#   neo/mutation_reference/mutated_proteins/
#   ├── all_samples_mutated_proteins.fasta
#   ├── per_sample/<SAMPLE>.mutated_proteins.fasta
#   ├── mutation_application_report.tsv
#   ├── translation_report.tsv
#   └── validation_report.tsv
#
# Usage:
#   sbatch 01_build_mutated_db.sh
#   sbatch 01_build_mutated_db.sh --sample-manifest \
#       /cluster/projects/livingbank/workspace/vivian/neo/260313_manifest_final.tsv
#   bash   01_build_mutated_db.sh --dry-run
#
# =============================================================================

#SBATCH --job-name=build_mutated_db
#SBATCH --account=hansengroup
#SBATCH -p all
#SBATCH -N 1
#SBATCH --mail-user=vivian.chu@mail.utoronto.ca
#SBATCH --mail-type=ALL
#SBATCH --time=06:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH -o /cluster/projects/livingbank/workspace/vivian/neo/mutation_reference/logs/build_mutated_db_%A_%a.out


set -euo pipefail

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
# SCRIPT_DIR is passed explicitly from submit_build_db.sh via --script-dir.
# Direct sbatch submissions are resolved below from SLURM_SUBMIT_DIR because
# Slurm may execute a copy of this script from /var/spool/slurmd/.
SCRIPT_DIR=""
NEO_BASE="/cluster/projects/livingbank/workspace/vivian/neo"

ORF_FASTA="${NEO_BASE}/mutation_reference/orf/cdna_orfs.fasta"
MUTATIONS="/cluster/projects/livingbank/Project/Pan-organoid/Mutation/Org_exome.data_mutations_extended.gt4.202507.txt"
OUTDIR="${NEO_BASE}/mutation_reference/mutated_proteins"
SAMPLE_MANIFEST=""
EXPRESSION_TSV=""
RNA_VARIANTS_TSV=""
HLA_TSV=""
BINDING_TSV=""
LABELS_TSV=""
DRY_RUN=false
EXTRA_ARGS=""

# ---------------------------------------------------------------------------
# Parse arguments
# ---------------------------------------------------------------------------
while [[ $# -gt 0 ]]; do
    case "$1" in
        --script-dir)         SCRIPT_DIR="$2";  shift 2 ;;
        --orf-fasta)          ORF_FASTA="$2";   shift 2 ;;
        --mutations)          MUTATIONS="$2";    shift 2 ;;
        --outdir)             OUTDIR="$2";       shift 2 ;;
        --sample-manifest)    SAMPLE_MANIFEST="$2"; shift 2 ;;
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
        --variant-context-flank-aa) EXTRA_ARGS="${EXTRA_ARGS} --variant-context-flank-aa $2"; shift 2 ;;
        --frameshift-tail-aa) EXTRA_ARGS="${EXTRA_ARGS} --frameshift-tail-aa $2"; shift 2 ;;
        --compress-tables)    EXTRA_ARGS="${EXTRA_ARGS} --compress-tables"; shift ;;
        --include-reference)  EXTRA_ARGS="${EXTRA_ARGS} --include-reference"; shift ;;
        --allow-ref-mismatch) EXTRA_ARGS="${EXTRA_ARGS} --allow-ref-mismatch"; shift ;;
        --dry-run)            DRY_RUN=true;      shift ;;
        *) echo "Unknown argument: $1" >&2; exit 1 ;;
    esac
done

_resolve_script_dir() {
    local candidates=()
    local bash_dir=""

    [[ -n "${SCRIPT_DIR}" ]] && candidates+=("${SCRIPT_DIR}")
    [[ -n "${MUTATION_REFERENCE_SCRIPT_DIR:-}" ]] && candidates+=("${MUTATION_REFERENCE_SCRIPT_DIR}")

    bash_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    candidates+=("${bash_dir}")

    if [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then
        candidates+=("${SLURM_SUBMIT_DIR}")
        candidates+=("${SLURM_SUBMIT_DIR}/workflow/mutation_reference")
    fi

    candidates+=("${PWD}")
    candidates+=("${PWD}/workflow/mutation_reference")

    for candidate in "${candidates[@]}"; do
        [[ -n "${candidate}" ]] || continue
        if [[ -f "${candidate}/build_mutated_protein_db.py" ]]; then
            SCRIPT_DIR="$(cd "${candidate}" && pwd)"
            return 0
        fi
    done

    echo "[01_build_mutated_db.sh] ERROR: Could not locate build_mutated_protein_db.py." >&2
    echo "  Submit from the repository root or workflow/mutation_reference directory," >&2
    echo "  or pass --script-dir /path/to/workflow/mutation_reference." >&2
    echo "  SLURM_SUBMIT_DIR=${SLURM_SUBMIT_DIR:-<unset>}" >&2
    echo "  BASH_SOURCE dir=${bash_dir}" >&2
    exit 1
}

_resolve_script_dir
DEFAULT_HLA_SUMMARY="$(cd "${SCRIPT_DIR}/../.." && pwd)/data/rna/hla_summary.csv"
if [[ -z "${HLA_TSV}" && -f "${DEFAULT_HLA_SUMMARY}" ]]; then
    HLA_TSV="${DEFAULT_HLA_SUMMARY}"
    EXTRA_ARGS="${EXTRA_ARGS} --hla-tsv ${HLA_TSV} --hla-sample-column sample"
fi

# ---------------------------------------------------------------------------
# Dry-run
# ---------------------------------------------------------------------------
if [[ "${DRY_RUN}" == "true" ]]; then
    echo "[dry-run] python3 ${SCRIPT_DIR}/build_mutated_protein_db.py \\"
    echo "    --orf-fasta  ${ORF_FASTA} \\"
    echo "    --mutations  ${MUTATIONS} \\"
    echo "    --outdir     ${OUTDIR} \\"
    echo "    --per-sample ${EXTRA_ARGS}"
    [[ -n "${SAMPLE_MANIFEST}" ]] && echo "    --sample-manifest ${SAMPLE_MANIFEST}"
    exit 0
fi

mkdir -p "${OUTDIR}" "$(dirname "${OUTDIR}")/logs"

# ---------------------------------------------------------------------------
# Validate inputs
# ---------------------------------------------------------------------------
errors=0
[[ ! -f "${ORF_FASTA}"  ]] && echo "[01_build_mutated_db.sh] ERROR: ORF FASTA not found: ${ORF_FASTA}"  >&2 && (( errors++ )) || true
[[ ! -f "${MUTATIONS}"  ]] && echo "[01_build_mutated_db.sh] ERROR: Mutations file not found: ${MUTATIONS}"   >&2 && (( errors++ )) || true
[[ -n "${SAMPLE_MANIFEST}" && ! -f "${SAMPLE_MANIFEST}" ]] \
    && echo "[01_build_mutated_db.sh] ERROR: Manifest not found: ${SAMPLE_MANIFEST}" >&2 && (( errors++ )) || true
[[ -n "${EXPRESSION_TSV}" && ! -f "${EXPRESSION_TSV}" ]] \
    && echo "[01_build_mutated_db.sh] ERROR: expression TSV not found: ${EXPRESSION_TSV}" >&2 && (( errors++ )) || true
[[ -n "${RNA_VARIANTS_TSV}" && ! -f "${RNA_VARIANTS_TSV}" ]] \
    && echo "[01_build_mutated_db.sh] ERROR: RNA variant TSV not found: ${RNA_VARIANTS_TSV}" >&2 && (( errors++ )) || true
[[ -n "${HLA_TSV}" && ! -f "${HLA_TSV}" ]] \
    && echo "[01_build_mutated_db.sh] ERROR: HLA TSV not found: ${HLA_TSV}" >&2 && (( errors++ )) || true
[[ -n "${BINDING_TSV}" && ! -f "${BINDING_TSV}" ]] \
    && echo "[01_build_mutated_db.sh] ERROR: binding TSV not found: ${BINDING_TSV}" >&2 && (( errors++ )) || true
[[ -n "${LABELS_TSV}" && ! -f "${LABELS_TSV}" ]] \
    && echo "[01_build_mutated_db.sh] ERROR: labels TSV not found: ${LABELS_TSV}" >&2 && (( errors++ )) || true
[[ "${errors}" -gt 0 ]] && exit 1

# ---------------------------------------------------------------------------
# Load environment
# ---------------------------------------------------------------------------
module load python/3.10 2>/dev/null || true
python3 -c "import pandas" 2>/dev/null || {
    echo "[01_build_mutated_db.sh] ERROR: pandas is not available in this environment." >&2
    echo "  Load a Python module/conda environment with pandas before submitting." >&2
    exit 1
}

# ---------------------------------------------------------------------------
# Run
# ---------------------------------------------------------------------------
echo "======================================"
echo " Build variant-centered mutant database"
echo " ORF FASTA : ${ORF_FASTA}"
echo " Mutations : ${MUTATIONS}"
echo " Output    : ${OUTDIR}"
echo " Started   : $(date)"
echo "======================================"

MANIFEST_ARG=""
[[ -n "${SAMPLE_MANIFEST}" ]] && MANIFEST_ARG="--sample-manifest ${SAMPLE_MANIFEST}"

python3 "${SCRIPT_DIR}/build_mutated_protein_db.py" \
    --orf-fasta       "${ORF_FASTA}"  \
    --mutations       "${MUTATIONS}"  \
    --outdir          "${OUTDIR}"     \
    --per-sample      \
    ${MANIFEST_ARG}   \
    ${EXTRA_ARGS}

echo "======================================"
echo " Finished: $(date)"
echo "======================================"

# ---------------------------------------------------------------------------
# Sanity check
# ---------------------------------------------------------------------------
COMBINED="${OUTDIR}/all_samples_mutated_proteins.fasta"
if [[ -f "${COMBINED}" ]]; then
    N=$(grep -c "^>" "${COMBINED}" || true)
    echo "Combined FASTA: ${N} records → ${COMBINED}"
fi

REPORT="${OUTDIR}/mutation_application_report.tsv"
[[ -f "${REPORT}.gz" ]] && REPORT="${REPORT}.gz"
if [[ -f "${REPORT}" ]]; then
    if [[ "${REPORT}" == *.gz ]]; then
        APPLIED=$(gzip -cd "${REPORT}" | awk -F'\t' 'NR>1 && $8=="applied"' | wc -l)
        SKIPPED=$(gzip -cd "${REPORT}" | awk -F'\t' 'NR>1 && $8!="applied"' | wc -l)
    else
        APPLIED=$(awk -F'\t' 'NR>1 && $8=="applied"' "${REPORT}" | wc -l)
        SKIPPED=$(awk -F'\t' 'NR>1 && $8!="applied"' "${REPORT}" | wc -l)
    fi
    echo "Mutations: ${APPLIED} applied, ${SKIPPED} skipped"
fi
