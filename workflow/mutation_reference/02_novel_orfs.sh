#!/bin/bash
# =============================================================================
# 02_novel_orfs.sh
#
# SLURM array job: predict and reformat novel transcript proteins for FragPipe.
#
# Per-sample novel transcript FASTAs were extracted by extract_novel.sh and
# ORF-predicted by TransDecoder in this script, then converted to a
# FragPipe-compatible protein FASTA per sample.
#
# Input layout:
#   NOVEL_DIR/<SAMPLE>/<SAMPLE>_novel_sequences.fasta   ← gffcompare output
#   or NOVEL_DIR/<SAMPLE>_novel_sequences.fasta         ← flat layout
# Cluster paths:
#   NOVEL_DIR = /cluster/projects/livingbank/workspace/vivian/neo/novel/novel_transcripts
#   OUT_DIR   = /cluster/projects/livingbank/workspace/vivian/neo/mutation_reference/novel_proteins
#
# Behaviour:
#   - If a .pep file already exists under OUT_DIR/transdecoder_pep, the
#     existing .pep is reformatted directly.
#   - If no .pep exists, TransDecoder is run and the output is then reformatted.
#
# Output per sample:
#   OUT_DIR/<SAMPLE>.novel_proteins.fasta    ← FragPipe-ready protein FASTA
#
# Output directory:
#   OUT_DIR = /cluster/projects/livingbank/workspace/vivian/neo/mutation_reference/novel_proteins
#
# FragPipe FASTA header format:
#   >SAMPLE|TRANSCRIPT_ID|novel|source=transdecoder|orf_type=TYPE|strand=STRAND
#
# Submit:
#   N=$(find /cluster/projects/livingbank/workspace/vivian/neo/novel/novel_transcripts \
#          -type f -name '*_novel_sequences.fasta' | wc -l)
#   sbatch --array=0-$((N-1)) 02_novel_orfs.sh
#
#   Dry-run a single task:
#   SLURM_ARRAY_TASK_ID=0 bash 02_novel_orfs.sh --dry-run
#
# =============================================================================

#SBATCH --job-name=novel_orfs
#SBATCH --account=hansengroup
#SBATCH -p himem
#SBATCH -N 1
#SBATCH --mail-user=vivian.chu@mail.utoronto.ca
#SBATCH --mail-type=ALL
#SBATCH --time=24:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=4
#SBATCH -o /cluster/projects/livingbank/workspace/vivian/neo/mutation_reference/logs/novel_orfs_%A_%a.out

set -euo pipefail

echo "======================================"
echo "Job started   : $(date)"
echo "Host          : $(hostname)"
echo "Array task    : ${SLURM_ARRAY_TASK_ID:-0}"
echo "======================================"

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
NEO_BASE="/cluster/projects/livingbank/workspace/vivian/neo"
NOVEL_DIR="${NEO_BASE}/novel/novel_transcripts"
OUT_DIR="${NEO_BASE}/mutation_reference/novel_proteins"
PEP_DIR=""
WORK_BASE=""
MIN_AA=25
DRY_RUN=false

# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------
while [[ $# -gt 0 ]]; do
    case "$1" in
        --novel-dir)      NOVEL_DIR="$2";      shift 2 ;;
        --out-dir)        OUT_DIR="$2";        shift 2 ;;
        --pep-dir)        PEP_DIR="$2";        shift 2 ;;
        --work-dir)       WORK_BASE="$2";      shift 2 ;;
        --min-aa)         MIN_AA="$2";         shift 2 ;;
        --dry-run)        DRY_RUN=true;        shift   ;;
        *) echo "Unknown argument: $1" >&2; exit 1 ;;
    esac
done

LOG_DIR="${NOVEL_DIR}/logs"
STATUS_DIR="${OUT_DIR}/status"
PEP_DIR="${PEP_DIR:-${OUT_DIR}/transdecoder_pep}"
WORK_BASE="${WORK_BASE:-${OUT_DIR}/transdecoder_work}"
mkdir -p "${OUT_DIR}" "${LOG_DIR}" "${STATUS_DIR}" "${PEP_DIR}" "${WORK_BASE}"

# ---------------------------------------------------------------------------
# Build ordered sample list from actual novel FASTA files.
#
# Supports both:
#   NOVEL_DIR/<SAMPLE>_novel_sequences.fasta
#   NOVEL_DIR/<SAMPLE>/<SAMPLE>_novel_sequences.fasta
# ---------------------------------------------------------------------------
SAMPLE_ROWS=()
while IFS= read -r row; do
    SAMPLE_ROWS+=("${row}")
done < <(
    find "${NOVEL_DIR}" -type f -name "*_novel_sequences.fasta" \
        ! -path "${LOG_DIR}/*" \
        | sort \
        | while IFS= read -r fasta; do
            sample="$(basename "${fasta}" _novel_sequences.fasta)"
            printf "%s\t%s\n" "${sample}" "${fasta}"
        done
)

if [[ "${#SAMPLE_ROWS[@]}" -eq 0 ]]; then
    echo "ERROR: No per-sample *_novel_sequences.fasta files found under ${NOVEL_DIR}" >&2
    exit 1
fi
echo "Novel FASTAs found : ${#SAMPLE_ROWS[@]}"

# ---------------------------------------------------------------------------
# Select this array task
# ---------------------------------------------------------------------------
TASK_ID="${SLURM_ARRAY_TASK_ID:-0}"
if [[ "${TASK_ID}" -ge "${#SAMPLE_ROWS[@]}" ]]; then
    echo "ERROR: TASK_ID ${TASK_ID} out of range (${#SAMPLE_ROWS[@]} novel FASTAs)" >&2
    exit 1
fi

IFS=$'\t' read -r SAMPLE FASTA <<< "${SAMPLE_ROWS[${TASK_ID}]}"
FRAGPIPE_FASTA="${OUT_DIR}/${SAMPLE}.novel_proteins.fasta"
STATUS_TSV="${STATUS_DIR}/${SAMPLE}.novel_orfs.status.tsv"

echo "Sample        : ${SAMPLE}"
echo "Novel FASTA   : ${FASTA}"
echo "Out dir       : ${OUT_DIR}"
echo "TransDecoder peptides: ${PEP_DIR}"
echo "TransDecoder work dir: ${WORK_BASE}"
echo ""

SEQ_COUNT=$(grep -c "^>" "${FASTA}" 2>/dev/null || true)
if [[ "${SEQ_COUNT}" -eq 0 ]]; then
    echo "[SKIP] ${SAMPLE}: empty FASTA"
    : > "${FRAGPIPE_FASTA}"
    printf "sample_id\tstatus\treason\tinput_fasta\tpep_file\toutput_fasta\tinput_sequences\tpredicted_orfs\twritten_records\n" > "${STATUS_TSV}"
    printf "%s\tskipped\tempty_input_fasta\t%s\t\t%s\t0\t0\t0\n" "${SAMPLE}" "${FASTA}" "${FRAGPIPE_FASTA}" >> "${STATUS_TSV}"
    exit 0
fi
echo "Input sequences: ${SEQ_COUNT}"

# ---------------------------------------------------------------------------
# Dry-run
# ---------------------------------------------------------------------------
FASTA_BASENAME=$(basename "${FASTA}")
PEP_FILE="${PEP_DIR}/${FASTA_BASENAME}.transdecoder.pep"

if [[ "${DRY_RUN}" == "true" ]]; then
    if [[ -f "${PEP_FILE}" ]]; then
        echo "[dry-run] Would reformat existing ${PEP_FILE}"
    else
        echo "[dry-run] Would run TransDecoder on ${FASTA}"
        echo "[dry-run] Would write TransDecoder .pep to ${PEP_FILE}"
    fi
    echo "[dry-run] Would write ${FRAGPIPE_FASTA}"
    exit 0
fi

# ---------------------------------------------------------------------------
# Load modules
# ---------------------------------------------------------------------------
module load TransDecoder 2>/dev/null || true
module load python/3.10   2>/dev/null || true

# ---------------------------------------------------------------------------
# Step 1 — TransDecoder
# ---------------------------------------------------------------------------
if [[ -f "${PEP_FILE}" ]]; then
    PEP_N=$(grep -c "^>" "${PEP_FILE}" 2>/dev/null || true)
    echo "[1/3] Using existing TransDecoder output: ${PEP_FILE} (${PEP_N} ORFs)"
else
    echo "[1/3] Running TransDecoder (min ${MIN_AA} aa) ..."

    ABS_FASTA=$(realpath "${FASTA}")
    WORKDIR="${WORK_BASE}/${SAMPLE}_wd"
    mkdir -p "${WORKDIR}"

    (
        cd "${WORKDIR}"

        TransDecoder.LongOrfs \
            -t "${ABS_FASTA}" \
            -m "${MIN_AA}"    \
            --output_dir .    \
            2>&1

        TransDecoder.Predict \
            -t "${ABS_FASTA}"   \
            --single_best_only  \
            --output_dir .      \
            2>&1

        LOCAL_PEP="${FASTA_BASENAME}.transdecoder.pep"
        if [[ -f "${LOCAL_PEP}" ]]; then
            cp "${LOCAL_PEP}" "${PEP_FILE}"
        else
            echo "WARNING: TransDecoder produced no .pep for ${SAMPLE}" >&2
            touch "${PEP_FILE}"
        fi
    )

    PEP_N=$(grep -c "^>" "${PEP_FILE}" 2>/dev/null || true)
    echo "  TransDecoder predicted ${PEP_N} ORFs"
fi

if [[ ! -s "${PEP_FILE}" ]]; then
    echo "[SKIP] ${SAMPLE}: no ORFs in ${PEP_FILE}"
    : > "${FRAGPIPE_FASTA}"
    printf "sample_id\tstatus\treason\tinput_fasta\tpep_file\toutput_fasta\tinput_sequences\tpredicted_orfs\twritten_records\n" > "${STATUS_TSV}"
    printf "%s\tskipped\tno_predicted_orfs\t%s\t%s\t%s\t%s\t0\t0\n" "${SAMPLE}" "${FASTA}" "${PEP_FILE}" "${FRAGPIPE_FASTA}" "${SEQ_COUNT}" >> "${STATUS_TSV}"
    exit 0
fi

# ---------------------------------------------------------------------------
# Step 2 — Reformat headers to FragPipe format
#
# TransDecoder header:
#   >MSTRG.2.1.p1 GENE.MSTRG.2.1~~MSTRG.2.1.p1  ORF type:complete len:127 (+) ...
#
# Output header:
#   >SAMPLE|MSTRG.2.1.p1|novel|source=transdecoder|orf_type=complete|strand=+
# ---------------------------------------------------------------------------
echo "[2/3] Reformatting headers → ${FRAGPIPE_FASTA} ..."

SAMPLE_ID="${SAMPLE}" \
PEP_FILE="${PEP_FILE}" \
FRAGPIPE_FASTA="${FRAGPIPE_FASTA}" \
python3 - <<'PYEOF'
import os, re, sys

sample_id = os.environ["SAMPLE_ID"]
pep_file  = os.environ["PEP_FILE"]
out_file  = os.environ["FRAGPIPE_FASTA"]
written = 0

with open(pep_file) as fh_in, open(out_file, "w") as fh_out:
    header = None
    seq_parts = []

    def flush(hdr, parts):
        global written
        seq = "".join(parts)
        if not seq:
            return
        transcript_id = hdr.split()[0]
        orf_type = "unknown"
        m = re.search(r"ORF type:(\S+)", hdr)
        if m:
            orf_type = m.group(1)
        strand = "."
        m = re.search(r"\(([+\-.])\)", hdr)
        if m:
            strand = m.group(1)
        new_hdr = (
            f">{sample_id}|{transcript_id}|novel"
            f"|source=transdecoder|orf_type={orf_type}|strand={strand}"
        )
        fh_out.write(new_hdr + "\n")
        for i in range(0, len(seq), 60):
            fh_out.write(seq[i:i+60] + "\n")
        written += 1

    for line in fh_in:
        line = line.rstrip()
        if line.startswith(">"):
            if header is not None:
                flush(header, seq_parts)
            header = line[1:]
            seq_parts = []
        elif line:
            seq_parts.append(line)
    if header is not None:
        flush(header, seq_parts)

print(f"  {written} novel protein records → {out_file}")
PYEOF

# ---------------------------------------------------------------------------
# Step 3 — Summary
# ---------------------------------------------------------------------------
N_OUT=$(grep -c "^>" "${FRAGPIPE_FASTA}" 2>/dev/null || true)
printf "sample_id\tstatus\treason\tinput_fasta\tpep_file\toutput_fasta\tinput_sequences\tpredicted_orfs\twritten_records\n" > "${STATUS_TSV}"
printf "%s\tcompleted\t\t%s\t%s\t%s\t%s\t%s\t%s\n" "${SAMPLE}" "${FASTA}" "${PEP_FILE}" "${FRAGPIPE_FASTA}" "${SEQ_COUNT}" "${PEP_N}" "${N_OUT}" >> "${STATUS_TSV}"
echo "[3/3] ${SAMPLE}: ${N_OUT} novel proteins → ${FRAGPIPE_FASTA}"

echo ""
echo "======================================"
echo "Job finished: $(date)"
echo "======================================"
