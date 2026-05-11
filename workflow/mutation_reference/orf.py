"""
ORF loading and validation.

Parses a cDNA ORF FASTA and validates each record for use in mutation
application and translation. Issues are written to a report rather than
silently discarding problematic ORFs.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Tuple

VALID_NUCLEOTIDES = set("ACGTNacgtn")
STOP_CODONS = {"TAA", "TAG", "TGA"}


@dataclass
class ORFRecord:
    gene: str
    transcript_id: str
    sequence: str  # uppercase nucleotides
    raw_header: str

    @property
    def length_nt(self) -> int:
        return len(self.sequence)

    @property
    def length_aa_expected(self) -> int:
        return self.length_nt // 3


@dataclass
class ValidationIssue:
    level: str  # 'WARNING' or 'ERROR'
    gene: str
    transcript_id: str
    issue_type: str
    message: str


def _parse_header(header: str) -> Tuple[str, str]:
    """
    Parse ORF FASTA header into (gene, transcript_id).

    Supports formats:
      >GENE|TRANSCRIPT_ID|...
      >GENE|TRANSCRIPT_ID
      >GENE
    """
    parts = header.split("|")
    gene = parts[0].strip()
    transcript_id = parts[1].strip() if len(parts) >= 2 else gene
    return gene, transcript_id


def load_orfs(
    fasta_path: str,
) -> Tuple[Dict[str, ORFRecord], List[ValidationIssue]]:
    """
    Parse a cDNA ORF FASTA file.

    Returns:
        records: dict keyed by transcript_id (or 'gene|transcript_id' if duplicate)
        issues:  list of ValidationIssue for any problems found
    """
    records: Dict[str, ORFRecord] = {}
    issues: List[ValidationIssue] = []

    current_header: str | None = None
    seq_parts: List[str] = []

    def _flush(header: str, parts: List[str]) -> None:
        seq = "".join(parts).upper()
        gene, transcript_id = _parse_header(header)

        def warn(issue_type: str, msg: str) -> None:
            issues.append(ValidationIssue(
                level="WARNING", gene=gene, transcript_id=transcript_id,
                issue_type=issue_type, message=msg,
            ))

        invalid = set(seq) - VALID_NUCLEOTIDES
        if invalid:
            warn("orf_invalid_nucleotides", f"Invalid characters: {sorted(invalid)}")

        if len(seq) % 3 != 0:
            warn("orf_length_not_divisible_by_three",
                 f"Length {len(seq)} not divisible by 3")

        if not seq.startswith("ATG"):
            warn("orf_missing_start_codon",
                 f"Does not start with ATG; starts with {seq[:3]!r}")

        # Check internal stops (all codons except the last)
        internal_stops = [
            (i * 3 + 1, seq[i * 3: i * 3 + 3])
            for i in range(len(seq) // 3 - 1)
            if seq[i * 3: i * 3 + 3] in STOP_CODONS
        ]
        if internal_stops:
            warn("orf_internal_stop",
                 f"Internal stop codons at cDNA positions: {internal_stops}")

        key = transcript_id
        if key in records:
            key = f"{gene}|{transcript_id}"
        records[key] = ORFRecord(
            gene=gene, transcript_id=transcript_id, sequence=seq, raw_header=header
        )

    with open(fasta_path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if current_header is not None:
                    _flush(current_header, seq_parts)
                current_header = line[1:]
                seq_parts = []
            elif line:
                seq_parts.append(line)

    if current_header is not None:
        _flush(current_header, seq_parts)

    return records, issues
