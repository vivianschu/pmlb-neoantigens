"""
FragPipe-compatible FASTA writing utilities.

Header format:
  >SAMPLE_ID|GENE|TRANSCRIPT_ID|mutated|variants=c.X,c.Y|protein_change=p.Z

Reference proteins (when --include-reference is used):
  >SAMPLE_ID|GENE|TRANSCRIPT_ID|reference
"""
from __future__ import annotations

import os
from typing import IO, List, Optional, Sequence, Tuple

_LINE_WIDTH = 60

FastaRecord = Tuple[str, str]   # (header_line_with_gt, sequence)


def format_header(
    sample_id: str,
    gene: str,
    transcript_id: str,
    variants: Sequence[str],
    protein_change: Optional[str] = None,
    sequence_id: Optional[str] = None,
    mutation_type: Optional[str] = None,
    peptide_coordinates: Optional[str] = None,
    is_reference: bool = False,
) -> str:
    label = "reference" if is_reference else "mutated"
    parts = [sequence_id or sample_id, gene, transcript_id, label]
    if sequence_id:
        parts.append(f"sample={sample_id}")
    if not is_reference and variants:
        parts.append(f"variants={','.join(variants)}")
    if not is_reference and protein_change:
        parts.append(f"protein_change={protein_change}")
    if not is_reference and mutation_type:
        parts.append(f"mutation_type={mutation_type}")
    if not is_reference and peptide_coordinates:
        parts.append(f"coords={peptide_coordinates}")
    return ">" + "|".join(parts)


def _write_record(fh: IO[str], header: str, sequence: str) -> None:
    fh.write(header + "\n")
    for i in range(0, len(sequence), _LINE_WIDTH):
        fh.write(sequence[i: i + _LINE_WIDTH] + "\n")


def write_fasta_file(path: str, records: List[FastaRecord]) -> None:
    """Write a list of (header, sequence) tuples to a FASTA file."""
    parent = os.path.dirname(path)
    if parent:
        os.makedirs(parent, exist_ok=True)
    with open(path, "w") as fh:
        for header, sequence in records:
            _write_record(fh, header, sequence)


def append_fasta_record(fh: IO[str], header: str, sequence: str) -> None:
    """Append a single record to an already-open file handle."""
    _write_record(fh, header, sequence)
