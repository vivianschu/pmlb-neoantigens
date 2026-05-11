#!/bin/bash
# =============================================================================
# 03_merge_fragpipe_db.sh
#
# Merge mutated protein FASTAs (from 01_build_mutated_db.sh) with novel
# transcript protein FASTAs (from 02_novel_orfs.sh) into final per-sample
# and combined FragPipe protein databases.
#
# Cluster paths (edit only the four NEO_BASE lines to change layout):
#   Mutated proteins : neo/mutation_reference/mutated_proteins/per_sample/
#   Novel proteins   : /cluster/projects/livingbank/workspace/vivian/neo/mutation_reference/novel_proteins/
#   Final output     : /cluster/projects/livingbank/workspace/vivian/neo/mutation_reference/final/
#
# Output structure:
#   FINAL_DIR/
#   ├── all_samples_combined.fasta
#   ├── per_sample/
#   │   ├── BPTO0051.TPO.fragpipe_db.fasta
#   │   ├── PPTO0002.TPO.fragpipe_db.fasta
#   │   └── ...
#   └── merge_summary.tsv
#
# Modes:
#   --mode combined      (default) mutated + novel per sample
#   --mode mutated-only  only mutated proteins, no novel
#   --mode novel-only    only novel proteins, no mutation calls
#
# Usage:
#   sbatch 03_merge_fragpipe_db.sh
#   sbatch 03_merge_fragpipe_db.sh --mode mutated-only
#   bash   03_merge_fragpipe_db.sh --dry-run
#
# =============================================================================

#SBATCH --job-name=merge_fragpipe_db
#SBATCH --account=hansengroup
#SBATCH -p all
#SBATCH -N 1
#SBATCH --mail-user=vivian.chu@mail.utoronto.ca
#SBATCH --mail-type=ALL
#SBATCH --time=01:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=2
#SBATCH -o /cluster/projects/livingbank/workspace/vivian/neo/mutation_reference/logs/merge_fragpipe_db_%A_%a.out

set -euo pipefail

# ---------------------------------------------------------------------------
# Paths  (edit ONLY these lines to change cluster/repo layout)
# ---------------------------------------------------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
NEO_BASE="/cluster/projects/livingbank/workspace/vivian/neo"

MUTATED_DIR="${NEO_BASE}/mutation_reference/mutated_proteins/per_sample"
NOVEL_DIR="${NEO_BASE}/mutation_reference/novel_proteins"
FINAL_DIR="${NEO_BASE}/mutation_reference/final"

# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------
MODE="combined"
DRY_RUN=false

while [[ $# -gt 0 ]]; do
    case "$1" in
        --mutated-dir) MUTATED_DIR="$2"; shift 2 ;;
        --novel-dir)   NOVEL_DIR="$2";   shift 2 ;;
        --outdir)      FINAL_DIR="$2";   shift 2 ;;
        --mode)        MODE="$2";        shift 2 ;;
        --dry-run)     DRY_RUN=true;     shift   ;;
        *) echo "[03_merge_fragpipe_db.sh] ERROR: Unknown argument: $1" >&2; exit 1 ;;
    esac
done

case "${MODE}" in
    combined|mutated-only|novel-only) ;;
    *) echo "[03_merge_fragpipe_db.sh] ERROR: --mode must be combined, mutated-only, or novel-only" >&2; exit 1 ;;
esac

if [[ "${DRY_RUN}" == "true" ]]; then
    echo "[dry-run] mode=${MODE}"
    echo "[dry-run] mutated-dir=${MUTATED_DIR}"
    echo "[dry-run] novel-dir=${NOVEL_DIR}"
    echo "[dry-run] final-dir=${FINAL_DIR}"
    echo "[dry-run] Would write per-sample FASTAs to ${FINAL_DIR}/per_sample/"
    echo "[dry-run] Would write ${FINAL_DIR}/all_samples_combined.fasta"
    exit 0
fi

mkdir -p "${FINAL_DIR}/per_sample" \
         "$(dirname "${FINAL_DIR}")/logs"

echo "======================================"
echo " Merge FragPipe protein database"
echo " mode        : ${MODE}"
echo " mutated-dir : ${MUTATED_DIR}"
echo " novel-dir   : ${NOVEL_DIR}"
echo " final-dir   : ${FINAL_DIR}"
echo " started     : $(date)"
echo "======================================"

# ---------------------------------------------------------------------------
# Validate
# ---------------------------------------------------------------------------
if [[ "${MODE}" != "novel-only" && ! -d "${MUTATED_DIR}" ]]; then
    echo "[03_merge_fragpipe_db.sh] ERROR: mutated-dir not found: ${MUTATED_DIR}" >&2
    echo "  → Run 01_build_mutated_db.sh first." >&2
    exit 1
fi
if [[ "${MODE}" != "mutated-only" && ! -d "${NOVEL_DIR}" ]]; then
    echo "[03_merge_fragpipe_db.sh] ERROR: novel-dir not found: ${NOVEL_DIR}" >&2
    echo "  → Run 02_novel_orfs.sh first." >&2
    exit 1
fi

# ---------------------------------------------------------------------------
# Collect per-sample FASTAs
# ---------------------------------------------------------------------------
MUTATED_SAMPLES=()
NOVEL_SAMPLES=()

if [[ "${MODE}" != "novel-only" ]]; then
    for f in "${MUTATED_DIR}"/*.mutated_proteins.fasta; do
        [[ -f "${f}" ]] || continue
        s=$(basename "${f}" .mutated_proteins.fasta)
        MUTATED_SAMPLES+=("${s}")
    done
    echo "Mutated FASTAs : ${#MUTATED_SAMPLES[@]}"
fi

if [[ "${MODE}" != "mutated-only" ]]; then
    for f in "${NOVEL_DIR}"/*.novel_proteins.fasta; do
        [[ -f "${f}" ]] || continue
        s=$(basename "${f}" .novel_proteins.fasta)
        NOVEL_SAMPLES+=("${s}")
    done
    echo "Novel FASTAs   : ${#NOVEL_SAMPLES[@]}"
fi

# Union of all sample IDs
ALL_SAMPLES=()
while IFS= read -r sample; do
    [[ -n "${sample}" ]] || continue
    ALL_SAMPLES+=("${sample}")
done < <(
    {
        if [[ "${#MUTATED_SAMPLES[@]}" -gt 0 ]]; then
            for s in "${MUTATED_SAMPLES[@]}"; do printf "%s\n" "${s}"; done
        fi
        if [[ "${#NOVEL_SAMPLES[@]}" -gt 0 ]]; then
            for s in "${NOVEL_SAMPLES[@]}"; do printf "%s\n" "${s}"; done
        fi
    } | sort -u
)

if [[ "${#ALL_SAMPLES[@]}" -eq 0 ]]; then
    echo "[03_merge_fragpipe_db.sh] ERROR: No FASTA files found in the input directories." >&2
    exit 1
fi
echo "Total samples  : ${#ALL_SAMPLES[@]}"
echo ""

# ---------------------------------------------------------------------------
# Merge per sample
# ---------------------------------------------------------------------------
COMBINED="${FINAL_DIR}/all_samples_combined.fasta"
SUMMARY="${FINAL_DIR}/merge_summary.tsv"
: > "${COMBINED}"

printf "sample_id\tmutated_proteins\tnovel_proteins\ttotal_proteins\toutput_fasta\n" \
    > "${SUMMARY}"

for sample in "${ALL_SAMPLES[@]}"; do
    mut_fasta="${MUTATED_DIR}/${sample}.mutated_proteins.fasta"
    nov_fasta="${NOVEL_DIR}/${sample}.novel_proteins.fasta"
    out_fasta="${FINAL_DIR}/per_sample/${sample}.fragpipe_db.fasta"

    n_mut=0
    n_nov=0
    : > "${out_fasta}"

    if [[ -n "${mut_fasta}" && -f "${mut_fasta}" ]]; then
        cat "${mut_fasta}" >> "${out_fasta}"
        n_mut=$(grep -c "^>" "${mut_fasta}" 2>/dev/null || true)
    elif [[ "${MODE}" == "combined" ]]; then
        echo "  WARNING: no mutated FASTA for ${sample}" >&2
    fi

    if [[ -n "${nov_fasta}" && -f "${nov_fasta}" ]]; then
        cat "${nov_fasta}" >> "${out_fasta}"
        n_nov=$(grep -c "^>" "${nov_fasta}" 2>/dev/null || true)
    elif [[ "${MODE}" == "combined" ]]; then
        echo "  WARNING: no novel FASTA for ${sample}" >&2
    fi

    n_total=$(( n_mut + n_nov ))
    cat "${out_fasta}" >> "${COMBINED}"

    printf "%s\t%d\t%d\t%d\t%s\n" \
        "${sample}" "${n_mut}" "${n_nov}" "${n_total}" \
        "per_sample/${sample}.fragpipe_db.fasta" \
        >> "${SUMMARY}"

    echo "  ${sample}: ${n_mut} mutated + ${n_nov} novel = ${n_total}"
done

# ---------------------------------------------------------------------------
# Final summary
# ---------------------------------------------------------------------------
N_COMBINED=$(grep -c "^>" "${COMBINED}" 2>/dev/null || true)
echo ""
echo "======================================"
echo " Done"
echo "  Combined FASTA : ${COMBINED}"
echo "             (${N_COMBINED} total protein records)"
echo "  Per-sample dir : ${FINAL_DIR}/per_sample/"
echo "  Summary        : ${SUMMARY}"
echo "  Finished       : $(date)"
echo "======================================"
