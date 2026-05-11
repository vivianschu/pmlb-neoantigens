"""
Mutation loading and application to cDNA ORF sequences.

Supports MAF input (primary format). Mutation coordinates are taken from HGVSc
notation and applied directly to the cDNA ORF nucleotide sequence.

HGVSc uses 'N' as a reference-base placeholder (common in tumor-only MAFs).
When the ref allele is 'N' or absent, the actual reference is read from the
ORF sequence and used — strict checking still applies to non-N alleles.
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

import pandas as pd


# ---------------------------------------------------------------------------
# Data classes
# ---------------------------------------------------------------------------

@dataclass
class Mutation:
    sample_id: str
    gene: str
    transcript_id: Optional[str]
    variant_id: str
    cdna_start: int          # 1-based start position in the ORF
    cdna_end: int            # 1-based end position (== cdna_start for SNVs)
    ref_allele: str          # '' means unknown / read from ORF
    alt_allele: str
    variant_type: str        # SNV | MNV | INS | DEL | INDEL
    variant_classification: Optional[str] = None
    tumor_vaf: Optional[float] = None
    t_alt_count: Optional[int] = None
    t_depth: Optional[int] = None
    hgvsc: Optional[str] = None
    hgvsp: Optional[str] = None


@dataclass
class MutationApplicationResult:
    sample_id: str
    gene: str
    transcript_id: Optional[str]
    variant_id: str
    cdna_position: int
    ref_allele: str
    alt_allele: str
    status: str          # applied | skipped_*
    status_reason: str
    reference_context: str
    mutated_context: str
    coordinate_shift: int = 0


# ---------------------------------------------------------------------------
# MAF column defaults
# ---------------------------------------------------------------------------

MAF_DEFAULTS: Dict[str, str] = {
    "sample_col":       "Tumor_Sample_Barcode",
    "gene_col":         "Hugo_Symbol",
    "transcript_col":   "Transcript_ID",
    "variant_class_col": "Variant_Classification",
    "variant_type_col": "Variant_Type",
    "ref_col":          "Reference_Allele",
    "alt_col":          "Tumor_Seq_Allele2",
    "hgvsc_col":        "HGVSc",
    "hgvsp_col":        "HGVSp_Short",
    "alt_count_col":    "t_alt_count",
    "depth_col":        "t_depth",
}

# MAF Variant_Classification → normalised type
_VC_TO_TYPE: Dict[str, str] = {
    "SNP": "SNV", "DNP": "MNV", "TNP": "MNV", "ONP": "MNV", "MNP": "MNV",
    "INS": "INS", "DEL": "DEL",
    "Frame_Shift_Del": "DEL", "Frame_Shift_Ins": "INS",
    "In_Frame_Del":    "DEL", "In_Frame_Ins":    "INS",
    "Missense_Mutation":  "SNV", "Nonsense_Mutation": "SNV",
    "Nonstop_Mutation":   "SNV", "Silent":             "SNV",
    "Splice_Site":        "SNV", "Start_Codon_SNP":    "SNV",
    "Stop_Codon_Del": "DEL",     "Stop_Codon_Ins": "INS",
}

# MAF Variant_Type column (simpler)
_VT_TO_TYPE: Dict[str, str] = {
    "SNP": "SNV", "DNP": "MNV", "TNP": "MNV", "ONP": "MNV",
    "INS": "INS", "DEL": "DEL",
}

_SILENT_CLASSIFICATIONS = {
    "Silent", "Synonymous", "Synonymous_Mutation", "RNA",
}


# ---------------------------------------------------------------------------
# HGVSc parsing
# ---------------------------------------------------------------------------

def _is_n_placeholder(seq: str) -> bool:
    """Return True if seq is entirely N placeholders (or empty)."""
    return not seq or all(c in "Nn" for c in seq)


def parse_hgvsc(
    hgvsc: str,
) -> Optional[Tuple[int, int, str, str, str]]:
    """
    Parse an HGVSc string into (cdna_start, cdna_end, ref, alt, variant_type).
    All positions are 1-based.  Returns None if outside CDS or unparseable.

    Handles:
      c.743G>A          SNV
      c.743N>A          SNV with N ref placeholder
      c.743_744GC>TA    MNV / INDEL
      c.743_744NN>TA    MNV with N ref
      c.743delG         DEL (single base, explicit ref)
      c.743delN         DEL (single base, N placeholder)
      c.743_745delGCC   DEL (multi-base, explicit)
      c.743_745delNNN   DEL (multi-base, N placeholder)
      c.743_744insAA    INS
      c.743_745delinsAA INDEL (delins)
      c.743dupA         DUP treated as INS
    """
    if not hgvsc or (isinstance(hgvsc, float)):
        return None

    hgvsc = str(hgvsc).strip()

    # Strip transcript prefix: "ENST00000269305.9:c.743G>A" → "c.743G>A"
    if ":" in hgvsc:
        hgvsc = hgvsc.split(":", 1)[1]

    if not hgvsc.startswith("c."):
        return None

    notation = hgvsc[2:]

    # Skip UTR / intronic (c.*, c.-, c.+)
    if not notation or notation[0] not in "0123456789":
        return None

    # --- SNV: 743G>A or 743N>A ---
    m = re.match(r"^(\d+)([ACGTNacgtn])>([ACGTNacgtn])$", notation)
    if m:
        pos = int(m.group(1))
        ref = m.group(2).upper()
        alt = m.group(3).upper()
        if _is_n_placeholder(ref):
            ref = ""
        return (pos, pos, ref, alt, "SNV")

    # --- Range substitution / MNV / INDEL: 743_744GC>TA ---
    m = re.match(
        r"^(\d+)_(\d+)([ACGTNacgtn]+)>([ACGTNacgtn]+)$", notation
    )
    if m:
        start, end = int(m.group(1)), int(m.group(2))
        ref = m.group(3).upper()
        alt = m.group(4).upper()
        if _is_n_placeholder(ref):
            ref = ""
        vtype = (
            "MNV"   if len(alt) == len(ref) and len(ref) > 1
            else "INDEL" if ref and alt
            else "SNV"
        )
        return (start, end, ref, alt, vtype)

    # --- DEL with bases: 743delG, 743delN, 743_745delGCC, 743_745delNNN ---
    m = re.match(r"^(\d+)(?:_(\d+))?del([ACGTNacgtn]*)$", notation)
    if m:
        start = int(m.group(1))
        end   = int(m.group(2)) if m.group(2) else start
        ref   = m.group(3).upper() if m.group(3) else ""
        if _is_n_placeholder(ref):
            ref = ""  # will be read from ORF; length = end - start + 1
        return (start, end, ref, "", "DEL")

    # --- INS: 743_744insAA ---
    m = re.match(r"^(\d+)_(\d+)ins([ACGTNacgtn]+)$", notation)
    if m:
        start = int(m.group(1))
        end   = int(m.group(2))
        alt   = m.group(3).upper()
        return (start, end, "", alt, "INS")

    # --- delins / INDEL: 743delinsAA, 743_745delinsAA ---
    m = re.match(r"^(\d+)(?:_(\d+))?delins([ACGTNacgtn]+)$", notation)
    if m:
        start = int(m.group(1))
        end   = int(m.group(2)) if m.group(2) else start
        alt   = m.group(3).upper()
        return (start, end, "", alt, "INDEL")

    # --- DUP: 743dupA, 743_745dup ---
    m = re.match(r"^(\d+)(?:_(\d+))?dup([ACGTNacgtn]*)$", notation)
    if m:
        start = int(m.group(1))
        end   = int(m.group(2)) if m.group(2) else start
        alt   = m.group(3).upper()
        return (start, end, "", alt, "INS")

    return None


# ---------------------------------------------------------------------------
# MAF loader
# ---------------------------------------------------------------------------

def load_mutations_maf(
    path: str,
    col_map: Optional[Dict[str, str]] = None,
    sample_manifest: Optional[List[str]] = None,
) -> List[Mutation]:
    """
    Load mutations from a MAF file.

    col_map overrides any of the keys in MAF_DEFAULTS, e.g.
        {'sample_col': 'patient_id', 'gene_col': 'gene_name'}
    """
    cols = {**MAF_DEFAULTS, **(col_map or {})}

    df = pd.read_csv(path, sep="\t", comment="#", low_memory=False)
    df.columns = df.columns.str.strip()

    manifest_set = set(sample_manifest) if sample_manifest else None

    mutations: List[Mutation] = []

    for _, row in df.iterrows():
        sample_id = str(row.get(cols["sample_col"], "")).strip()
        if not sample_id or sample_id == "nan":
            continue
        if manifest_set and sample_id not in manifest_set:
            continue

        gene           = str(row.get(cols["gene_col"],         "")).strip()
        transcript_id  = str(row.get(cols["transcript_col"],   "")).strip() or None
        variant_class  = str(row.get(cols["variant_class_col"],"")).strip()
        variant_type_raw = str(row.get(cols["variant_type_col"],"")).strip()
        hgvsc          = str(row.get(cols["hgvsc_col"],        "")).strip()
        hgvsp          = str(row.get(cols["hgvsp_col"],        "")).strip()

        if variant_class in _SILENT_CLASSIFICATIONS:
            continue

        if hgvsc in ("", "nan"):
            continue

        parsed = parse_hgvsc(hgvsc)
        if parsed is None:
            continue

        cdna_start, cdna_end, ref_allele, alt_allele, vtype = parsed

        # Resolve variant type from MAF columns when HGVS parse is ambiguous
        if variant_type_raw in _VT_TO_TYPE:
            vtype = _VT_TO_TYPE[variant_type_raw]
        elif variant_class in _VC_TO_TYPE:
            vtype = _VC_TO_TYPE[variant_class]

        try:
            alt_count = int(row.get(cols["alt_count_col"], 0))
        except (ValueError, TypeError):
            alt_count = None

        try:
            depth = int(row.get(cols["depth_col"], 0))
        except (ValueError, TypeError):
            depth = None

        vaf = (alt_count / depth) if (alt_count is not None and depth) else None

        variant_id = f"{gene}:{hgvsc}"
        if transcript_id and transcript_id != "nan":
            variant_id = f"{transcript_id}:{hgvsc}"

        mutations.append(Mutation(
            sample_id=sample_id,
            gene=gene,
            transcript_id=transcript_id if transcript_id != "nan" else None,
            variant_id=variant_id,
            cdna_start=cdna_start,
            cdna_end=cdna_end,
            ref_allele=ref_allele,
            alt_allele=alt_allele,
            variant_type=vtype,
            variant_classification=variant_class,
            tumor_vaf=vaf,
            t_alt_count=alt_count,
            t_depth=depth,
            hgvsc=hgvsc,
            hgvsp=hgvsp if hgvsp != "nan" else None,
        ))

    return mutations


# ---------------------------------------------------------------------------
# Mutation application
# ---------------------------------------------------------------------------

def apply_mutations_to_sequence(
    sequence: str,
    mutations: List[Mutation],
    strict_ref_check: bool = True,
) -> Tuple[str, List[MutationApplicationResult]]:
    """
    Apply a list of mutations to a cDNA ORF sequence.

    Mutations are sorted by cdna_start and applied in order. A running
    cumulative_offset tracks how earlier indels shift the positions of
    later mutations.

    When ref_allele is empty ('N' placeholder), the actual reference is read
    from the ORF at the mapped position and substituted in before checking.
    strict_ref_check then only fires on non-empty, explicitly provided refs.

    Returns:
        mutated_sequence: str
        results: list of MutationApplicationResult (one per input mutation)
    """
    results: List[MutationApplicationResult] = []
    seq = list(sequence.upper())
    cumulative_offset = 0

    for mut in sorted(mutations, key=lambda m: m.cdna_start):
        # Convert to 0-based with offset from prior indels
        pos0 = mut.cdna_start - 1 + cumulative_offset

        # Context window for reporting
        ctx_lo = max(0, pos0 - 5)
        ref_ctx_hi = min(len(seq), pos0 + max(len(mut.ref_allele), len(mut.alt_allele), 1) + 5)
        ref_context = "".join(seq[ctx_lo:ref_ctx_hi])

        def _skip(status: str, reason: str) -> MutationApplicationResult:
            return MutationApplicationResult(
                sample_id=mut.sample_id, gene=mut.gene,
                transcript_id=mut.transcript_id, variant_id=mut.variant_id,
                cdna_position=mut.cdna_start, ref_allele=mut.ref_allele,
                alt_allele=mut.alt_allele, status=status, status_reason=reason,
                reference_context=ref_context,
                mutated_context="".join(seq[ctx_lo:ref_ctx_hi]),
            )

        def _applied(shift: int = 0) -> MutationApplicationResult:
            mut_ctx_hi = min(len(seq), ctx_lo + len(ref_context) + shift)
            return MutationApplicationResult(
                sample_id=mut.sample_id, gene=mut.gene,
                transcript_id=mut.transcript_id, variant_id=mut.variant_id,
                cdna_position=mut.cdna_start, ref_allele=mut.ref_allele,
                alt_allele=mut.alt_allele, status="applied", status_reason="",
                reference_context=ref_context,
                mutated_context="".join(seq[ctx_lo:mut_ctx_hi]),
                coordinate_shift=shift,
            )

        # Bounds check
        if pos0 < 0 or pos0 >= len(seq):
            results.append(_skip(
                "skipped_invalid_cdna_position",
                f"Position {mut.cdna_start} (offset {cumulative_offset}) "
                f"out of bounds for sequence length {len(sequence)}",
            ))
            continue

        vtype = mut.variant_type

        # --- SNV ---
        if vtype == "SNV":
            actual_ref = seq[pos0]
            expected_ref = mut.ref_allele or actual_ref  # '' → don't check
            if strict_ref_check and mut.ref_allele and actual_ref != expected_ref:
                results.append(_skip(
                    "skipped_ref_mismatch",
                    f"Expected {expected_ref!r} at c.{mut.cdna_start}, found {actual_ref!r}",
                ))
                continue
            seq[pos0] = mut.alt_allele
            results.append(_applied(0))

        # --- MNV / multi-nucleotide substitution ---
        elif vtype in ("MNV", "DNP", "TNP", "ONP"):
            ref_len = len(mut.ref_allele) if mut.ref_allele else (mut.cdna_end - mut.cdna_start + 1)
            actual_ref = "".join(seq[pos0: pos0 + ref_len])
            if strict_ref_check and mut.ref_allele and actual_ref != mut.ref_allele:
                results.append(_skip(
                    "skipped_ref_mismatch",
                    f"Expected {mut.ref_allele!r} at c.{mut.cdna_start}, found {actual_ref!r}",
                ))
                continue
            seq[pos0: pos0 + ref_len] = list(mut.alt_allele)
            shift = len(mut.alt_allele) - ref_len
            cumulative_offset += shift
            results.append(_applied(shift))

        # --- DEL ---
        elif vtype == "DEL":
            # Length from explicit ref or from the coordinate span
            if mut.ref_allele:
                ref_len = len(mut.ref_allele)
            else:
                ref_len = mut.cdna_end - mut.cdna_start + 1
            actual_ref = "".join(seq[pos0: pos0 + ref_len])
            if strict_ref_check and mut.ref_allele and actual_ref != mut.ref_allele:
                results.append(_skip(
                    "skipped_ref_mismatch",
                    f"Expected {mut.ref_allele!r} at c.{mut.cdna_start}, found {actual_ref!r}",
                ))
                continue
            del seq[pos0: pos0 + ref_len]
            shift = -ref_len
            cumulative_offset += shift
            results.append(_applied(shift))

        # --- INS ---
        elif vtype == "INS":
            # Insert alt_allele after pos0 (between cdna_start and cdna_start+1)
            if not mut.alt_allele:
                results.append(_skip(
                    "skipped_unsupported_variant_type",
                    "INS with empty alt_allele",
                ))
                continue
            ins_pos = pos0 + 1
            seq[ins_pos:ins_pos] = list(mut.alt_allele)
            shift = len(mut.alt_allele)
            cumulative_offset += shift
            results.append(_applied(shift))

        # --- INDEL (delins) ---
        elif vtype == "INDEL":
            if mut.ref_allele:
                ref_len = len(mut.ref_allele)
            else:
                ref_len = mut.cdna_end - mut.cdna_start + 1
            actual_ref = "".join(seq[pos0: pos0 + ref_len])
            if strict_ref_check and mut.ref_allele and actual_ref != mut.ref_allele:
                results.append(_skip(
                    "skipped_ref_mismatch",
                    f"Expected {mut.ref_allele!r} at c.{mut.cdna_start}, found {actual_ref!r}",
                ))
                continue
            seq[pos0: pos0 + ref_len] = list(mut.alt_allele)
            shift = len(mut.alt_allele) - ref_len
            cumulative_offset += shift
            results.append(_applied(shift))

        else:
            results.append(_skip(
                "skipped_unsupported_variant_type",
                f"Variant type {mut.variant_type!r} is not supported",
            ))

    return "".join(seq), results
