#!/bin/bash
# =============================================================================
# 00_fetch_cdna_orfs.sh
#
# Extract CDS sequences from the GENCODE v47 protein-coding transcript FASTA
# and write a cDNA ORF FASTA for use in build_mutated_protein_db.py.
#
# Source file (read-only):
#   /cluster/projects/livingbank/workspace/references/gencode.v47.pc_transcripts.fa.gz
#
# The pc_transcripts FASTA contains full mRNA sequences (5'UTR + CDS + 3'UTR).
# HGVSc coordinates (c.1, c.2, …) are numbered from the start of the CDS, so
# we must extract only the CDS portion of each transcript.
#
# CDS coordinates are encoded in the GENCODE header:
#   >ENST…|ENSG…|…|GENE|…|protein_coding|CDS:203-1384
#
# Extraction strategy (in priority order):
#   1. Parse CDS:start-end from the header and slice the sequence.
#   2. If no CDS field: find the first ATG and use that as the CDS start.
#   3. If neither works: skip the transcript and report it.
#
# Output:
#   <OUTDIR>/cdna_orfs.fasta    headers: >GENE|TRANSCRIPT_ID|ORF
#
# Usage:
#   sbatch 00_fetch_cdna_orfs.sh
#   bash   00_fetch_cdna_orfs.sh --dry-run
#   bash   00_fetch_cdna_orfs.sh --outdir /custom/path
#
# =============================================================================

#SBATCH --job-name=fetch_cdna_orfs
#SBATCH --account=hansengroup
#SBATCH -p all
#SBATCH -N 1
#SBATCH --mail-user=vivian.chu@mail.utoronto.ca
#SBATCH --mail-type=ALL
#SBATCH --time=01:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=2
#SBATCH -o /cluster/projects/livingbank/workspace/vivian/neo/mutation_reference/logs/fetch_cdna_orfs_%A_%a.out


set -euo pipefail

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
NEO_BASE="/cluster/projects/livingbank/workspace/vivian/neo"

SOURCE_GZ="/cluster/projects/livingbank/workspace/references/gencode.v47.pc_transcripts.fa.gz"
OUTDIR="${NEO_BASE}/mutation_reference/orf"
DRY_RUN=false

# ---------------------------------------------------------------------------
# Parse arguments
# ---------------------------------------------------------------------------
while [[ $# -gt 0 ]]; do
    case "$1" in
        --source-gz) SOURCE_GZ="$2"; shift 2 ;;
        --outdir)    OUTDIR="$2";    shift 2 ;;
        --dry-run)   DRY_RUN=true;   shift   ;;
        *) echo "Unknown arg: $1" >&2; exit 1 ;;
    esac
done

mkdir -p "${OUTDIR}" logs

OUT_FASTA="${OUTDIR}/cdna_orfs.fasta"

echo "Source  : ${SOURCE_GZ}"
echo "Output  : ${OUT_FASTA}"

if [[ ! -f "${SOURCE_GZ}" ]]; then
    echo "ERROR: Source file not found: ${SOURCE_GZ}" >&2
    exit 1
fi

if [[ "${DRY_RUN}" == "true" ]]; then
    echo "[dry-run] Would extract CDS from ${SOURCE_GZ} → ${OUT_FASTA}"
    exit 0
fi

# ---------------------------------------------------------------------------
# Load modules
# ---------------------------------------------------------------------------
module load python/3.10 2>/dev/null || true

# ---------------------------------------------------------------------------
# Extract CDS and reformat headers
#
# GENCODE pipe-delimited header fields:
#   0: transcript_id   4: gene_name   7+: may include CDS:start-end (1-based)
#
# CDS coordinates in the header are 1-based, inclusive on both ends.
# Python slicing is 0-based, so: seq[cds_start-1 : cds_end]
# ---------------------------------------------------------------------------
echo "Extracting CDS sequences ..."

SOURCE_GZ="${SOURCE_GZ}" OUT_FASTA="${OUT_FASTA}" python3 - <<'PYEOF'
import gzip, os, re

source = os.environ["SOURCE_GZ"]
out    = os.environ["OUT_FASTA"]

n_cds_coords  = 0   # extracted using CDS:start-end from header
n_first_atg   = 0   # fell back to first ATG
n_skipped     = 0   # no CDS coords and no ATG found

with gzip.open(source, "rt") as fh_in, open(out, "w") as fh_out:
    header    = None
    seq_parts = []

    def flush(hdr, parts):
        global n_cds_coords, n_first_atg, n_skipped

        seq    = "".join(parts).upper()
        fields = hdr.split("|")
        tid    = fields[0].strip()
        gene   = fields[4].strip() if len(fields) > 4 else tid

        # --- Strategy 1: CDS coordinates in header ---
        cds_seq = None
        m = re.search(r"CDS:(\d+)-(\d+)", hdr)
        if m:
            # 1-based inclusive → 0-based slice
            cds_start = int(m.group(1)) - 1
            cds_end   = int(m.group(2))
            cds_seq   = seq[cds_start:cds_end]
            if cds_seq.startswith("ATG"):
                n_cds_coords += 1
            else:
                # Header coordinates point to a non-ATG start — fall through
                cds_seq = None

        # --- Strategy 2: first ATG in the transcript ---
        if cds_seq is None:
            atg_pos = seq.find("ATG")
            if atg_pos >= 0:
                cds_seq = seq[atg_pos:]
                n_first_atg += 1
            else:
                n_skipped += 1
                return

        fh_out.write(f">{gene}|{tid}|ORF\n")
        for i in range(0, len(cds_seq), 60):
            fh_out.write(cds_seq[i:i+60] + "\n")

    for line in fh_in:
        line = line.rstrip()
        if line.startswith(">"):
            if header is not None:
                flush(header, seq_parts)
            header    = line[1:]
            seq_parts = []
        elif line:
            seq_parts.append(line)
    if header is not None:
        flush(header, seq_parts)

total = n_cds_coords + n_first_atg + n_skipped
print(f"  CDS coords from header : {n_cds_coords:,}")
print(f"  Fell back to first ATG : {n_first_atg:,}")
print(f"  Skipped (no ATG found) : {n_skipped:,}")
print(f"  Total written          : {n_cds_coords + n_first_atg:,} / {total:,} transcripts")
PYEOF

# ---------------------------------------------------------------------------
# Validate
# ---------------------------------------------------------------------------
N=$(grep -c "^>" "${OUT_FASTA}" || true)
echo "Total ORF records: ${N}"
[[ "${N}" -eq 0 ]] && echo "ERROR: No records written." >&2 && exit 1

echo "Done → ${OUT_FASTA}"
