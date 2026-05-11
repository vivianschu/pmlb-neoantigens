#!/usr/bin/env python3
"""
Build sample-specific variant-centered mutant FASTA databases for FragPipe
immunopeptidomics workflows.

Pipeline:
  cDNA ORF FASTA + mutation calls  →  map each coding variant to an ORF  →
  translate local mutant protein context  →  write compact FASTA entries

Usage:
  python build_mutated_protein_db.py \\
      --orf-fasta    data/orf/cdna_orfs.fasta \\
      --mutations    data/mutation/Org_exome.data_mutations_extended.gt4.202507.txt \\
      --outdir       results/mutation_reference \\
      --per-sample

Optional:
  --sample-manifest  samples.txt        one sample_id per line, or a TSV with a 'sample' column
  --mutation-format  maf|tsv            (default: maf)
  --include-reference                   also write unmodified reference proteins
  --allow-ref-mismatch                  apply mutations even when HGVSc ref != ORF
  --translate-through-stop              do not halt at first stop codon
  --sample-column    COL                override default MAF column names
  --gene-column      COL
  --transcript-column COL
  --ref-column       COL
  --alt-column       COL
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import os
import re
import sys
from collections import defaultdict
from typing import Dict, Iterable, List, Optional

# Ensure the package directory is on sys.path when invoked directly
sys.path.insert(0, os.path.dirname(__file__))

# Strips RNA-seq batch/replicate suffixes to recover the canonical MAF barcode.
# e.g. CSC0073.TXO_novo → CSC0073.TXO, PPTO0002.TPO_pmlb → PPTO0002.TPO
_SUFFIX_RE = re.compile(r"(_pmlb|_novo\d*|\.r\d+|_r\d+)+$", re.IGNORECASE)


def _normalize_sample_id(sample_id: str) -> str:
    return _SUFFIX_RE.sub("", sample_id)

from fasta import format_header, write_fasta_file, append_fasta_record
from mutations import load_mutations_maf, apply_mutations_to_sequence, Mutation
from orf import load_orfs, ORFRecord, ValidationIssue
from translate import translate, TranslationResult


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Build mutated protein FASTA for FragPipe immunopeptidomics",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--orf-fasta",   required=True,
                   help="FASTA containing cDNA ORF sequences (nucleotide)")
    p.add_argument("--mutations",   required=True,
                   help="Mutation file (MAF, VCF, or TSV)")
    p.add_argument("--outdir",      default="results/mutation_reference",
                   help="Output directory")
    p.add_argument("--mutation-format", default="maf",
                   choices=["maf", "vcf", "tsv"],
                   help="Mutation file format")
    p.add_argument("--sample-manifest",
                   help="Plain text (one sample_id per line) or TSV with a 'sample' column; "
                        "batch suffixes (_pmlb, _novo, etc.) are stripped to match MAF barcodes")
    p.add_argument("--per-sample", action="store_true",
                   help="Write a separate FASTA per sample under outdir/per_sample/")
    p.add_argument("--include-reference", action="store_true",
                   help="Also write unmodified reference protein for each mutated gene")
    p.add_argument("--allow-ref-mismatch", action="store_true",
                   help="Apply mutations even when the ORF ref does not match HGVSc ref")
    p.add_argument("--translate-through-stop", action="store_true",
                   help="Continue translation past stop codons")
    p.add_argument("--deduplicate", action=argparse.BooleanOptionalAction,
                   default=True,
                   help="Deduplicate identical protein sequences per sample/cohort FASTA")
    p.add_argument("--compress-tables", action="store_true",
                   help="Write TSV reports as .tsv.gz")
    p.add_argument("--peptide-min-length", type=int, default=8,
                   help="Minimum mutant peptide candidate length")
    p.add_argument("--peptide-max-length", type=int, default=15,
                   help="Maximum mutant peptide candidate length")
    p.add_argument("--frameshift-tail-aa", type=int, default=120,
                   help="Maximum altered frameshift/stop-loss tail length used for peptide candidates")
    p.add_argument("--variant-context-flank-aa", type=int, default=30,
                   help="Reference flanking amino acids to include on each side of the altered interval in mutant FASTA records")
    p.add_argument("--expression-tsv",
                   help="Optional TSV with TPM expression values")
    p.add_argument("--rna-variants-tsv",
                   help="Optional TSV with RNA VAF/depth/alternate-read counts")
    p.add_argument("--hla-tsv",
                   help="Optional TSV mapping sample IDs to one or more HLA alleles")
    p.add_argument("--binding-tsv",
                   help="Optional TSV with binding predictions keyed by sample/peptide/HLA")
    p.add_argument("--labels-tsv",
                   help="Optional TSV with binder/presented/immunogenic labels")
    p.add_argument("--expression-sample-column", default="sample_id")
    p.add_argument("--expression-gene-column", default="gene")
    p.add_argument("--expression-transcript-column", default="transcript_id")
    p.add_argument("--expression-tpm-column", default="TPM")
    p.add_argument("--rna-sample-column", default="sample_id")
    p.add_argument("--rna-gene-column", default="gene")
    p.add_argument("--rna-transcript-column", default="transcript_id")
    p.add_argument("--rna-variant-column", default="variant_id")
    p.add_argument("--rna-hgvsc-column", default="HGVSc")
    p.add_argument("--rna-vaf-column", default="rna_vaf")
    p.add_argument("--rna-depth-column", default="rna_depth")
    p.add_argument("--rna-alt-count-column", default="rna_alt_count")
    p.add_argument("--hla-sample-column", default="sample_id")
    p.add_argument("--hla-allele-column", default="hla_allele")
    p.add_argument("--binding-sample-column", default="sample_id")
    p.add_argument("--binding-peptide-column", default="mutant_peptide_sequence")
    p.add_argument("--binding-hla-column", default="hla_allele")
    p.add_argument("--binding-ic50-column", default="predicted_binding_ic50")
    p.add_argument("--binding-rank-column", default="predicted_binding_rank")
    p.add_argument("--label-sample-column", default="sample_id")
    p.add_argument("--label-peptide-column", default="mutant_peptide_sequence")
    p.add_argument("--label-hla-column", default="hla_allele")
    p.add_argument("--binder-label-column", default="binder_label")
    p.add_argument("--presented-label-column", default="presented_label")
    p.add_argument("--immunogenic-label-column", default="immunogenic_label")
    p.add_argument("--binder-ic50-threshold", type=float,
                   help="If set, fill missing binder_label from IC50 <= threshold")
    p.add_argument("--binder-rank-threshold", type=float,
                   help="If set, fill missing binder_label from rank <= threshold")
    p.add_argument("--min-tpm", type=float,
                   help="Optional TPM filter; missing TPM passes unless --require-rna-filters is set")
    p.add_argument("--min-rna-vaf", type=float,
                   help="Optional RNA VAF filter; missing RNA VAF passes unless --require-rna-filters is set")
    p.add_argument("--min-rna-depth", type=int,
                   help="Optional RNA depth filter; missing RNA depth passes unless --require-rna-filters is set")
    p.add_argument("--require-rna-filters", action="store_true",
                   help="Fail/filter variants with missing RNA data when RNA filters are requested")
    # Column name overrides
    p.add_argument("--sample-column")
    p.add_argument("--gene-column")
    p.add_argument("--transcript-column")
    p.add_argument("--ref-column")
    p.add_argument("--alt-column")
    return p


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _load_manifest(path: str) -> List[str]:
    """
    Load sample IDs from a manifest file.

    Accepts either:
      - A plain text file with one sample ID per line, OR
      - A TSV with a 'sample' header column (e.g. 260313_manifest_final.tsv)

    In both cases, RNA-seq batch suffixes (_pmlb, _novo, _novo1/_novo2, etc.)
    are stripped so the resulting IDs match MAF Tumor_Sample_Barcode values.
    Duplicate canonical IDs are collapsed.
    """
    with open(path) as fh:
        first = fh.readline()
        fields = first.strip().split("\t")

        if "sample" in fields:
            sample_col = fields.index("sample")
            raw_ids = [ln.strip().split("\t")[sample_col] for ln in fh if ln.strip()]
        else:
            # plain one-column file; first line is a sample ID, not a header
            raw_ids = [first.strip()] + [ln.strip() for ln in fh if ln.strip()]

    canonical = dict.fromkeys(
        _normalize_sample_id(s) for s in raw_ids
        if s and not s.startswith("#")
    )
    return list(canonical)


def _write_tsv(path: str, rows: list, columns: list) -> None:
    parent = os.path.dirname(path)
    if parent:
        os.makedirs(parent, exist_ok=True)
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "wt", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=columns, delimiter="\t",
                           extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def _table_path(outdir: str, name: str, compress: bool) -> str:
    return os.path.join(outdir, name + (".gz" if compress else ""))


def _open_text(path: str):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def _read_tsv(path: Optional[str]) -> list[dict]:
    if not path:
        return []
    with _open_text(path) as fh:
        sample = fh.read(4096)
        fh.seek(0)
        delimiter = "," if sample.count(",") > sample.count("\t") else "\t"
        return list(csv.DictReader(fh, delimiter=delimiter))


def _clean(value) -> Optional[str]:
    if value is None:
        return None
    value = str(value).strip()
    if not value or value.lower() in {"nan", "none", "null", "na"}:
        return None
    return value


def _to_float(value) -> Optional[float]:
    value = _clean(value)
    if value is None:
        return None
    try:
        return float(value)
    except ValueError:
        return None


def _to_int(value) -> Optional[int]:
    value = _clean(value)
    if value is None:
        return None
    try:
        return int(float(value))
    except ValueError:
        return None


def _load_expression(args: argparse.Namespace) -> dict[tuple[str, str, str], Optional[float]]:
    values: dict[tuple[str, str, str], Optional[float]] = {}
    rows = _read_tsv(args.expression_tsv)
    if not rows:
        return values

    fieldnames = list(rows[0].keys())
    is_long = (
        args.expression_sample_column in fieldnames
        and args.expression_tpm_column in fieldnames
    )

    if not is_long:
        gene_col = args.expression_gene_column
        if gene_col not in fieldnames:
            gene_col = fieldnames[0]
        metadata_cols = {gene_col, args.expression_transcript_column, "gene_name", "symbol"}
        for row in rows:
            gene = _clean(row.get(gene_col))
            tx = _clean(row.get(args.expression_transcript_column)) or ""
            if not gene:
                continue
            for sample_col in fieldnames:
                if sample_col in metadata_cols:
                    continue
                tpm = _to_float(row.get(sample_col))
                if tpm is None:
                    continue
                sample = _normalize_sample_id(sample_col)
                for key in ((sample, gene, tx), (sample, gene, "")):
                    current = values.get(key)
                    values[key] = tpm if current is None else max(current, tpm)
        return values

    for row in rows:
        sample = _clean(row.get(args.expression_sample_column))
        gene = _clean(row.get(args.expression_gene_column))
        tx = _clean(row.get(args.expression_transcript_column)) or ""
        if not sample or not gene:
            continue
        sample = _normalize_sample_id(sample)
        tpm = _to_float(row.get(args.expression_tpm_column))
        values[(sample, gene, tx)] = tpm
        values.setdefault((sample, gene, ""), tpm)
    return values


def _load_rna_variants(args: argparse.Namespace) -> dict[tuple[str, str, str, str], dict]:
    values: dict[tuple[str, str, str, str], dict] = {}
    for row in _read_tsv(args.rna_variants_tsv):
        sample = _clean(row.get(args.rna_sample_column))
        gene = _clean(row.get(args.rna_gene_column)) or ""
        tx = _clean(row.get(args.rna_transcript_column)) or ""
        variant = _clean(row.get(args.rna_variant_column)) or _clean(row.get(args.rna_hgvsc_column))
        if not sample or not variant:
            continue
        payload = {
            "rna_vaf": _to_float(row.get(args.rna_vaf_column)),
            "rna_depth": _to_int(row.get(args.rna_depth_column)),
            "rna_alt_count": _to_int(row.get(args.rna_alt_count_column)),
        }
        for key in {
            (sample, gene, tx, variant),
            (sample, gene, "", variant),
            (sample, "", "", variant),
        }:
            values.setdefault(key, payload)
    return values


def _load_hla(args: argparse.Namespace) -> dict[str, list[str]]:
    def _add(sample: str, allele: str) -> None:
        sample = _normalize_sample_id(sample)
        allele = _normalize_hla_allele(allele)
        if allele and allele not in alleles[sample]:
            alleles[sample].append(allele)

    alleles: dict[str, list[str]] = defaultdict(list)
    for row in _read_tsv(args.hla_tsv):
        status = _clean(row.get("status"))
        if status and status.lower() not in {"success", "ok", "completed", "complete"}:
            continue

        sample = _clean(row.get(args.hla_sample_column)) or _clean(row.get("sample"))
        if not sample:
            continue

        allele_field = _clean(row.get(args.hla_allele_column))
        if allele_field:
            for allele in re.split(r"[,; ]+", allele_field):
                _add(sample, allele)
            continue

        for col in ("hla_a1", "hla_a2", "hla_b1", "hla_b2", "hla_c1", "hla_c2"):
            allele = _clean(row.get(col))
            if allele:
                _add(sample, allele)
    return alleles


def _normalize_hla_allele(allele: Optional[str]) -> Optional[str]:
    allele = _clean(allele)
    if allele is None:
        return None
    if allele.upper() == "NA":
        return None
    allele = allele.replace("HLA-", "")
    locus = allele.split("*", 1)[0].upper()
    if "*" not in allele or locus not in {"A", "B", "C"}:
        return allele
    fields = allele.split("*", 1)[1].split(":")
    if len(fields) < 2:
        return f"HLA-{locus}*{fields[0]}"
    return f"HLA-{locus}*{fields[0]}:{fields[1]}"


def _load_binding(args: argparse.Namespace) -> dict[tuple[str, str, str], dict]:
    values = {}
    for row in _read_tsv(args.binding_tsv):
        sample = _clean(row.get(args.binding_sample_column)) or ""
        peptide = _clean(row.get(args.binding_peptide_column))
        hla = _clean(row.get(args.binding_hla_column)) or ""
        if not peptide:
            continue
        values[(sample, peptide, hla)] = {
            "predicted_binding_ic50": _to_float(row.get(args.binding_ic50_column)),
            "predicted_binding_rank": _to_float(row.get(args.binding_rank_column)),
        }
    return values


def _load_labels(args: argparse.Namespace) -> dict[tuple[str, str, str], dict]:
    values = {}
    for row in _read_tsv(args.labels_tsv):
        sample = _clean(row.get(args.label_sample_column)) or ""
        peptide = _clean(row.get(args.label_peptide_column))
        hla = _clean(row.get(args.label_hla_column)) or ""
        if not peptide:
            continue
        values[(sample, peptide, hla)] = {
            "binder_label": _clean(row.get(args.binder_label_column)),
            "presented_label": _clean(row.get(args.presented_label_column)),
            "immunogenic_label": _clean(row.get(args.immunogenic_label_column)),
        }
    return values


def _lookup_expression(expr: dict, sample: str, gene: str, tx: str) -> Optional[float]:
    if (sample, gene, tx) in expr:
        return expr[(sample, gene, tx)]
    return expr.get((sample, gene, ""))


def _lookup_rna_variant(rna: dict, mut: Mutation) -> dict:
    tx = mut.transcript_id or ""
    for variant in (mut.variant_id, mut.hgvsc):
        if not variant:
            continue
        for key in (
            (mut.sample_id, mut.gene, tx, variant),
            (mut.sample_id, mut.gene, "", variant),
            (mut.sample_id, "", "", variant),
        ):
            if key in rna:
                return rna[key]
    return {"rna_vaf": None, "rna_depth": None, "rna_alt_count": None}


def _passes_rna_filters(args: argparse.Namespace, tpm: Optional[float], rna: dict) -> tuple[bool, str]:
    checks = [
        ("TPM", tpm, args.min_tpm),
        ("RNA_VAF", rna.get("rna_vaf"), args.min_rna_vaf),
        ("RNA_depth", rna.get("rna_depth"), args.min_rna_depth),
    ]
    reasons = []
    for name, observed, threshold in checks:
        if threshold is None:
            continue
        if observed is None:
            if args.require_rna_filters:
                reasons.append(f"missing_{name}")
            continue
        if observed < threshold:
            reasons.append(f"{name}_below_{threshold}")
    return (not reasons, ";".join(reasons))


def _sequence_id(sample_id: str, sequence: str) -> str:
    digest = hashlib.sha1(sequence.encode()).hexdigest()[:12]
    return f"mutprot_{digest}"


def _first_difference(ref: str, mut: str) -> int:
    for i, (a, b) in enumerate(zip(ref, mut)):
        if a != b:
            return i
    return min(len(ref), len(mut))


def _altered_interval(ref: str, mut: str, max_tail: int) -> tuple[int, int]:
    start = _first_difference(ref, mut)
    ref_suffix = len(ref)
    mut_suffix = len(mut)
    while ref_suffix > start and mut_suffix > start and ref[ref_suffix - 1] == mut[mut_suffix - 1]:
        ref_suffix -= 1
        mut_suffix -= 1
    end = max(start + 1, mut_suffix)
    if end - start > max_tail:
        end = start + max_tail
    return start, min(end, len(mut))


def _variant_context_sequence(
    mut_protein: str,
    ref_protein: str,
    flank_aa: int,
    max_tail: int,
) -> tuple[str, int, int, int, int]:
    """
    Return a compact mutant FASTA segment around the altered amino-acid interval.

    Coordinates are 0-based half-open on the full mutant protein:
      segment_start, segment_end, altered_start, altered_end
    """
    altered_start, altered_end = _altered_interval(ref_protein, mut_protein, max_tail)
    segment_start = max(0, altered_start - flank_aa)
    segment_end = min(len(mut_protein), altered_end + flank_aa)
    return (
        mut_protein[segment_start:segment_end],
        segment_start,
        segment_end,
        altered_start,
        altered_end,
    )


def _peptide_candidates(
    mut_protein: str,
    ref_protein: str,
    min_len: int,
    max_len: int,
    max_tail: int,
) -> Iterable[dict]:
    start, end = _altered_interval(ref_protein, mut_protein, max_tail)
    seen: set[tuple[str, int, int]] = set()
    for length in range(min_len, max_len + 1):
        lo = max(0, start - length + 1)
        hi = min((end - 1), len(mut_protein) - length)
        for pep_start in range(lo, hi + 1):
            pep_end = pep_start + length
            if pep_end > len(mut_protein):
                continue
            mut_pep = mut_protein[pep_start:pep_end]
            if "*" in mut_pep or "X" in mut_pep:
                continue
            wt_pep = ref_protein[pep_start:pep_end] if pep_end <= len(ref_protein) else ""
            if wt_pep == mut_pep:
                continue
            key = (mut_pep, pep_start + 1, pep_end)
            if key in seen:
                continue
            seen.add(key)
            yield {
                "mutant_peptide_sequence": mut_pep,
                "wildtype_peptide_sequence": wt_pep,
                "peptide_start_aa": pep_start + 1,
                "peptide_end_aa": pep_end,
                "peptide_length": length,
            }


def _orf_lookup(
    orfs: Dict[str, ORFRecord],
) -> tuple[Dict[str, list[ORFRecord]], Dict[str, ORFRecord]]:
    """Build gene→[ORF] and transcript_id→ORF lookup dicts."""
    by_gene: Dict[str, list[ORFRecord]] = defaultdict(list)
    by_transcript: Dict[str, ORFRecord] = {}
    for rec in orfs.values():
        by_gene[rec.gene].append(rec)
        # Also index by bare ENST (without version suffix)
        by_transcript[rec.transcript_id] = rec
        bare = rec.transcript_id.split(".")[0]
        if bare != rec.transcript_id:
            by_transcript.setdefault(bare, rec)
    return by_gene, by_transcript


def _find_orf(
    gene: str,
    transcript_id: Optional[str],
    by_gene: Dict[str, list[ORFRecord]],
    by_transcript: Dict[str, ORFRecord],
) -> Optional[ORFRecord]:
    """Resolve the best matching ORF record for a given gene/transcript."""
    if transcript_id:
        if transcript_id in by_transcript:
            return by_transcript[transcript_id]
        bare = transcript_id.split(".")[0]
        if bare in by_transcript:
            return by_transcript[bare]
    candidates = by_gene.get(gene, [])
    if not candidates:
        return None
    if len(candidates) == 1:
        return candidates[0]
    # Prefer the longest ORF when multiple transcripts exist for one gene
    return max(candidates, key=lambda r: r.length_nt)


def _variant_summary(muts: List[Mutation], max_show: int = 3) -> tuple[list[str], Optional[str]]:
    """Return (variant_strings, first_protein_change)."""
    hgvscs = [m.hgvsc for m in muts if m.hgvsc]
    if not hgvscs:
        hgvscs = [m.variant_id for m in muts]
    shown = hgvscs[:max_show]
    if len(hgvscs) > max_show:
        shown.append(f"...+{len(hgvscs) - max_show}")
    protein_change = next((m.hgvsp for m in muts if m.hgvsp), None)
    return shown, protein_change


# ---------------------------------------------------------------------------
# Core processing
# ---------------------------------------------------------------------------

def run(args: argparse.Namespace) -> None:
    if args.peptide_min_length < 1 or args.peptide_max_length < args.peptide_min_length:
        sys.exit("ERROR: peptide length bounds are invalid")
    if args.variant_context_flank_aa < 0:
        sys.exit("ERROR: --variant-context-flank-aa must be >= 0")
    if args.require_rna_filters and not (args.expression_tsv or args.rna_variants_tsv):
        sys.exit("ERROR: --require-rna-filters needs --expression-tsv and/or --rna-variants-tsv")

    os.makedirs(args.outdir, exist_ok=True)
    per_sample_dir = os.path.join(args.outdir, "per_sample")
    if args.per_sample:
        os.makedirs(per_sample_dir, exist_ok=True)

    strict_ref = not args.allow_ref_mismatch
    stop_at_stop = not args.translate_through_stop

    expression = _load_expression(args)
    rna_variants = _load_rna_variants(args)
    hla_by_sample = _load_hla(args)
    binding = _load_binding(args)
    labels = _load_labels(args)
    print(
        f"Optional feature inputs: expression={len(expression):,} keys, "
        f"RNA variants={len(rna_variants):,} keys, "
        f"HLA samples={len(hla_by_sample):,}, binding={len(binding):,}, labels={len(labels):,}"
    )

    # ── Load ORFs ────────────────────────────────────────────────────────────
    print(f"Loading ORFs from {args.orf_fasta} ...")
    orfs, orf_issues = load_orfs(args.orf_fasta)
    print(f"  {len(orfs):,} ORF records  |  {len(orf_issues):,} validation issues")

    by_gene, by_transcript = _orf_lookup(orfs)

    # ── Load sample manifest ─────────────────────────────────────────────────
    manifest: Optional[List[str]] = None
    if args.sample_manifest:
        manifest = _load_manifest(args.sample_manifest)
        print(f"  Manifest: {len(manifest):,} canonical sample IDs (after suffix normalization)")

    # ── Load mutations ───────────────────────────────────────────────────────
    print(f"Loading mutations from {args.mutations} ...")
    col_map: dict = {}
    if args.sample_column:      col_map["sample_col"]      = args.sample_column
    if args.gene_column:        col_map["gene_col"]        = args.gene_column
    if args.transcript_column:  col_map["transcript_col"]  = args.transcript_column
    if args.ref_column:         col_map["ref_col"]         = args.ref_column
    if args.alt_column:         col_map["alt_col"]         = args.alt_column

    if args.mutation_format == "maf":
        mutations = load_mutations_maf(args.mutations, col_map=col_map,
                                       sample_manifest=manifest)
    else:
        sys.exit(f"ERROR: format '{args.mutation_format}' not yet implemented")

    print(f"  {len(mutations):,} parseable cDNA mutations loaded")

    # Group: sample → (gene, transcript) → [Mutation]
    by_sample: Dict[str, Dict[tuple, List[Mutation]]] = defaultdict(
        lambda: defaultdict(list)
    )
    for m in mutations:
        by_sample[m.sample_id][(m.gene, m.transcript_id or "")].append(m)

    # ── Process each sample ───────────────────────────────────────────────────
    mut_report_rows: list = []
    tx_report_rows: list = []
    protein_metadata_rows: list = []
    rna_feature_rows: list = []
    peptide_rows: list = []
    fasta_seen: set[str] = set()

    samples = sorted(by_sample)
    print(f"Processing {len(samples):,} samples ...")

    combined_path = os.path.join(args.outdir, "all_samples_mutated_proteins.fasta")
    combined_fh = open(combined_path, "w")

    for sample_id in samples:
        sample_records: list = []
        sample_seen: set[str] = set()

        for (gene, transcript_key), muts in by_sample[sample_id].items():
            orf = _find_orf(gene, transcript_key or None, by_gene, by_transcript)

            if orf is None:
                for m in muts:
                    mut_report_rows.append({
                        "sample_id": sample_id, "gene": gene,
                        "transcript_id": transcript_key,
                        "variant_id": m.variant_id,
                        "cdna_position": m.cdna_start,
                        "ref_allele": m.ref_allele, "alt_allele": m.alt_allele,
                        "status": "skipped_no_matching_orf",
                        "status_reason": "No ORF found for gene/transcript",
                        "reference_context": "", "mutated_context": "",
                        "coordinate_shift": 0,
                })
                continue

            tpm = _lookup_expression(expression, sample_id, gene, orf.transcript_id)
            eligible_muts: list[Mutation] = []
            for m in muts:
                rna = _lookup_rna_variant(rna_variants, m)
                passes_filter, filter_reason = _passes_rna_filters(args, tpm, rna)
                rna_feature_rows.append({
                    "sample_id": sample_id,
                    "gene": gene,
                    "transcript_id": orf.transcript_id,
                    "variant_id": m.variant_id,
                    "hgvsc": m.hgvsc,
                    "TPM": tpm,
                    "rna_vaf": rna.get("rna_vaf"),
                    "rna_depth": rna.get("rna_depth"),
                    "rna_alt_count": rna.get("rna_alt_count"),
                    "rna_filter_pass": passes_filter,
                    "rna_filter_reason": filter_reason,
                })
                if not passes_filter:
                    mut_report_rows.append({
                        "sample_id": sample_id, "gene": gene,
                        "transcript_id": orf.transcript_id,
                        "variant_id": m.variant_id,
                        "cdna_position": m.cdna_start,
                        "ref_allele": m.ref_allele, "alt_allele": m.alt_allele,
                        "status": "skipped_rna_filter",
                        "status_reason": filter_reason,
                        "reference_context": "", "mutated_context": "",
                        "coordinate_shift": 0,
                    })
                    continue
                eligible_muts.append(m)

            if not eligible_muts:
                continue

            ref_tx = translate(orf.sequence, stop_at_first_stop=stop_at_stop)
            hla_alleles = hla_by_sample.get(sample_id, [""])

            for mut in eligible_muts:
                mutated_seq, app_results = apply_mutations_to_sequence(
                    orf.sequence, [mut], strict_ref_check=strict_ref
                )

                for r in app_results:
                    mut_report_rows.append({
                        "sample_id": r.sample_id, "gene": r.gene,
                        "transcript_id": r.transcript_id,
                        "variant_id": r.variant_id,
                        "cdna_position": r.cdna_position,
                        "ref_allele": r.ref_allele, "alt_allele": r.alt_allele,
                        "status": r.status, "status_reason": r.status_reason,
                        "reference_context": r.reference_context,
                        "mutated_context": r.mutated_context,
                        "coordinate_shift": r.coordinate_shift,
                    })

                applied = [r for r in app_results if r.status == "applied"]
                if not applied:
                    continue

                mut_tx = translate(
                    mutated_seq,
                    reference_protein=ref_tx.protein_sequence,
                    ref_orf_length_nt=len(orf.sequence),
                    stop_at_first_stop=stop_at_stop,
                    ref_stop_codon_position=ref_tx.stop_codon_position,
                )

                variant_strs, protein_change = _variant_summary([mut])
                variant_ids = mut.variant_id
                hgvscs = mut.hgvsc or ""

                if mut_tx.is_synonymous:
                    tx_report_rows.append({
                        "sample_id": sample_id, "gene": gene,
                        "transcript_id": orf.transcript_id,
                        "reference_orf_length_nt": len(orf.sequence),
                        "mutated_orf_length_nt": len(mutated_seq),
                        "reference_protein_length_aa": len(ref_tx.protein_sequence),
                        "mutated_protein_length_aa": len(mut_tx.protein_sequence),
                        "mutation_count": len(applied),
                        "has_frameshift": mut_tx.has_frameshift,
                        "has_premature_stop": mut_tx.has_premature_stop,
                        "has_start_loss": mut_tx.has_start_loss,
                        "has_stop_loss": mut_tx.has_stop_loss,
                        "protein_change_summary": mut_tx.mutation_effect,
                        "sequence_status": "skipped_synonymous",
                        "fasta_id": "",
                    })
                    continue

                context_seq, ctx_start, ctx_end, alt_start, alt_end = _variant_context_sequence(
                    mut_tx.protein_sequence,
                    ref_tx.protein_sequence,
                    args.variant_context_flank_aa,
                    args.frameshift_tail_aa,
                )
                if not context_seq:
                    continue

                fasta_id = _sequence_id(sample_id, f"{mut.variant_id}|{context_seq}")
                coords = (
                    f"aa{ctx_start + 1}-{ctx_end};"
                    f"altered=aa{alt_start + 1}-{alt_end}"
                )
                header = format_header(
                    sample_id, gene, orf.transcript_id,
                    variant_strs, protein_change,
                    sequence_id=fasta_id,
                    mutation_type=mut_tx.mutation_effect,
                    peptide_coordinates=coords,
                )

                tx_report_rows.append({
                    "sample_id": sample_id, "gene": gene,
                    "transcript_id": orf.transcript_id,
                    "reference_orf_length_nt": len(orf.sequence),
                    "mutated_orf_length_nt": len(mutated_seq),
                    "reference_protein_length_aa": len(ref_tx.protein_sequence),
                    "mutated_protein_length_aa": len(mut_tx.protein_sequence),
                    "mutation_count": len(applied),
                    "has_frameshift": mut_tx.has_frameshift,
                    "has_premature_stop": mut_tx.has_premature_stop,
                    "has_start_loss": mut_tx.has_start_loss,
                    "has_stop_loss": mut_tx.has_stop_loss,
                    "protein_change_summary": mut_tx.mutation_effect,
                    "sequence_status": "written_variant_context",
                    "fasta_id": fasta_id,
                })

                protein_metadata_rows.append({
                    "fasta_id": fasta_id,
                    "sample_id": sample_id,
                    "gene": gene,
                    "transcript_id": orf.transcript_id,
                    "variant_ids": variant_ids,
                    "hgvsc": hgvscs,
                    "amino_acid_change": protein_change,
                    "mutation_type": mut_tx.mutation_effect,
                    "protein_length_aa": len(context_seq),
                    "sequence_sha1": hashlib.sha1(context_seq.encode()).hexdigest(),
                    "fasta_header": header[1:],
                })

                rna = _lookup_rna_variant(rna_variants, mut)
                write_record = (not args.deduplicate) or fasta_id not in sample_seen
                if write_record:
                    sample_records.append((header, context_seq))
                    sample_seen.add(fasta_id)
                if (not args.deduplicate) or fasta_id not in fasta_seen:
                    append_fasta_record(combined_fh, header, context_seq)
                    fasta_seen.add(fasta_id)

                for pep in _peptide_candidates(
                    mut_tx.protein_sequence,
                    ref_tx.protein_sequence,
                    args.peptide_min_length,
                    args.peptide_max_length,
                    args.frameshift_tail_aa,
                ):
                    for hla in hla_alleles:
                        bind = (
                            binding.get((sample_id, pep["mutant_peptide_sequence"], hla))
                            or binding.get(("", pep["mutant_peptide_sequence"], hla))
                            or binding.get((sample_id, pep["mutant_peptide_sequence"], ""))
                            or {}
                        )
                        label = (
                            labels.get((sample_id, pep["mutant_peptide_sequence"], hla))
                            or labels.get(("", pep["mutant_peptide_sequence"], hla))
                            or labels.get((sample_id, pep["mutant_peptide_sequence"], ""))
                            or {}
                        )
                        binder_label = label.get("binder_label")
                        ic50 = bind.get("predicted_binding_ic50")
                        rank = bind.get("predicted_binding_rank")
                        if binder_label is None and args.binder_ic50_threshold is not None and ic50 is not None:
                            binder_label = int(ic50 <= args.binder_ic50_threshold)
                        if binder_label is None and args.binder_rank_threshold is not None and rank is not None:
                            binder_label = int(rank <= args.binder_rank_threshold)
                        peptide_rows.append({
                            **pep,
                            "sample_id": sample_id,
                            "gene": gene,
                            "transcript_id": orf.transcript_id,
                            "protein_fasta_id": fasta_id,
                            "variant_id": variant_ids,
                            "hgvsc": hgvscs,
                            "amino_acid_change": protein_change,
                            "mutation_type": mut_tx.mutation_effect,
                            "hla_allele": hla,
                            "hla_pseudosequence_key": hla,
                            "TPM": tpm,
                            "rna_vaf": rna.get("rna_vaf"),
                            "rna_depth": rna.get("rna_depth"),
                            "rna_alt_count": rna.get("rna_alt_count"),
                            "dna_vaf": mut.tumor_vaf,
                            "dna_depth": mut.t_depth,
                            "dna_alt_count": mut.t_alt_count,
                            "predicted_binding_ic50": ic50,
                            "predicted_binding_rank": rank,
                            "binder_label": binder_label,
                            "presented_label": label.get("presented_label"),
                            "immunogenic_label": label.get("immunogenic_label"),
                        })

            if args.include_reference:
                ref_header = format_header(
                    sample_id, gene, orf.transcript_id, [], is_reference=True
                )
                append_fasta_record(combined_fh, ref_header, ref_tx.protein_sequence)
                sample_records.append((ref_header, ref_tx.protein_sequence))

        if args.per_sample and sample_records:
            sample_path = os.path.join(per_sample_dir,
                                       f"{sample_id}.mutated_proteins.fasta")
            write_fasta_file(sample_path, sample_records)

    combined_fh.close()
    print(f"Wrote combined FASTA → {combined_path}")

    # ── Reports ───────────────────────────────────────────────────────────────
    MUT_COLS = [
        "sample_id", "gene", "transcript_id", "variant_id", "cdna_position",
        "ref_allele", "alt_allele", "status", "status_reason",
        "reference_context", "mutated_context", "coordinate_shift",
    ]
    TX_COLS = [
        "sample_id", "gene", "transcript_id",
        "reference_orf_length_nt", "mutated_orf_length_nt",
        "reference_protein_length_aa", "mutated_protein_length_aa",
        "mutation_count", "has_frameshift", "has_premature_stop",
        "has_start_loss", "has_stop_loss", "protein_change_summary",
        "sequence_status", "fasta_id",
    ]
    PROTEIN_META_COLS = [
        "fasta_id", "sample_id", "gene", "transcript_id", "variant_ids",
        "hgvsc", "amino_acid_change", "mutation_type", "protein_length_aa",
        "sequence_sha1", "fasta_header",
    ]
    RNA_COLS = [
        "sample_id", "gene", "transcript_id", "variant_id", "hgvsc", "TPM",
        "rna_vaf", "rna_depth", "rna_alt_count", "rna_filter_pass",
        "rna_filter_reason",
    ]
    PEPTIDE_COLS = [
        "sample_id", "gene", "transcript_id", "protein_fasta_id", "variant_id",
        "hgvsc", "amino_acid_change", "mutation_type",
        "mutant_peptide_sequence", "wildtype_peptide_sequence",
        "peptide_start_aa", "peptide_end_aa", "peptide_length",
        "hla_allele", "hla_pseudosequence_key", "TPM", "rna_vaf",
        "dna_vaf", "rna_depth", "rna_alt_count", "dna_depth",
        "dna_alt_count", "predicted_binding_ic50", "predicted_binding_rank",
        "binder_label", "presented_label", "immunogenic_label",
    ]
    VAL_COLS = [
        "level", "sample_id", "gene", "transcript_id", "issue_type", "message",
    ]

    _write_tsv(_table_path(args.outdir, "mutation_application_report.tsv", args.compress_tables),
               mut_report_rows, MUT_COLS)
    _write_tsv(_table_path(args.outdir, "translation_report.tsv", args.compress_tables),
               tx_report_rows, TX_COLS)
    _write_tsv(_table_path(args.outdir, "validation_report.tsv", args.compress_tables),
               [{"level": i.level, "sample_id": "", "gene": i.gene,
                 "transcript_id": i.transcript_id,
                 "issue_type": i.issue_type, "message": i.message}
                for i in orf_issues],
               VAL_COLS)
    _write_tsv(_table_path(args.outdir, "variant_to_protein_metadata.tsv", args.compress_tables),
               protein_metadata_rows, PROTEIN_META_COLS)
    _write_tsv(_table_path(args.outdir, "rna_feature_table.tsv", args.compress_tables),
               rna_feature_rows, RNA_COLS)
    _write_tsv(_table_path(args.outdir, "ml_training_table.tsv", args.compress_tables),
               peptide_rows, PEPTIDE_COLS)

    n_applied  = sum(1 for r in mut_report_rows if r["status"] == "applied")
    n_skipped  = len(mut_report_rows) - n_applied
    n_proteins = len(protein_metadata_rows)

    print(
        f"\nSummary\n"
        f"  Mutations applied  : {n_applied:,}\n"
        f"  Mutations skipped  : {n_skipped:,}\n"
        f"  Protein records    : {n_proteins:,}\n"
        f"  Peptide candidates : {len(peptide_rows):,}\n"
        f"  Reports → {args.outdir}/"
    )


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    run(_build_parser().parse_args())
