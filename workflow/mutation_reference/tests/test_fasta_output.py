"""Tests for FASTA header formatting and file writing."""
import os
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))

from fasta import format_header, write_fasta_file, _LINE_WIDTH
from argparse import Namespace
from build_mutated_protein_db import (
    _load_expression,
    _load_hla,
    _peptide_candidates,
    _variant_context_sequence,
)


def test_header_mutated_with_protein_change():
    h = format_header(
        "PPTO0002.TPO", "TP53", "ENST00000269305",
        ["c.743G>A"], "p.R248Q",
    )
    assert h == ">PPTO0002.TPO|TP53|ENST00000269305|mutated|variants=c.743G>A|protein_change=p.R248Q"


def test_header_mutated_multiple_variants():
    h = format_header(
        "S1", "KRAS", "ENST00000256078",
        ["c.35G>A", "c.34G>T"],
    )
    assert "variants=c.35G>A,c.34G>T" in h
    assert "|mutated|" in h


def test_header_reference():
    h = format_header(
        "PPTO0002.TPO", "TP53", "ENST00000269305", [], is_reference=True
    )
    assert h == ">PPTO0002.TPO|TP53|ENST00000269305|reference"


def test_header_reference_ignores_variants():
    h = format_header(
        "S1", "BRCA1", "ENST1", ["c.1A>T"], "p.M1I", is_reference=True
    )
    assert "variants" not in h
    assert "protein_change" not in h
    assert "|reference" in h


def test_header_mutated_with_traceability_metadata():
    h = format_header(
        "S1", "TP53", "ENST1", ["c.4G>A"], "p.A2T",
        sequence_id="mutprot_abc123",
        mutation_type="missense",
        peptide_coordinates="aa1-4",
    )
    assert h.startswith(">mutprot_abc123|TP53|ENST1|mutated|sample=S1")
    assert "variants=c.4G>A" in h
    assert "protein_change=p.A2T" in h
    assert "mutation_type=missense" in h
    assert "coords=aa1-4" in h


def test_write_fasta_basic():
    with tempfile.TemporaryDirectory() as d:
        path = os.path.join(d, "out.fasta")
        records = [
            (">SEQ1|GENE1", "MALPQA"),
            (">SEQ2|GENE2", "MWRFKL"),
        ]
        write_fasta_file(path, records)
        content = open(path).read()
        assert ">SEQ1|GENE1\n" in content
        assert "MALPQA\n" in content
        assert ">SEQ2|GENE2\n" in content
        assert "MWRFKL\n" in content


def test_write_fasta_line_wrapping():
    with tempfile.TemporaryDirectory() as d:
        path = os.path.join(d, "out.fasta")
        seq = "M" * (_LINE_WIDTH * 2 + 10)
        write_fasta_file(path, [(">LONGSEQ", seq)])
        lines = open(path).readlines()
        assert lines[0].strip() == ">LONGSEQ"
        assert len(lines[1].strip()) == _LINE_WIDTH
        assert len(lines[2].strip()) == _LINE_WIDTH
        assert len(lines[3].strip()) == 10


def test_write_fasta_creates_parent_dirs():
    with tempfile.TemporaryDirectory() as d:
        path = os.path.join(d, "subdir", "out.fasta")
        write_fasta_file(path, [(">SEQ1", "MAPL")])
        assert os.path.exists(path)


def test_write_fasta_empty():
    with tempfile.TemporaryDirectory() as d:
        path = os.path.join(d, "empty.fasta")
        write_fasta_file(path, [])
        assert os.path.exists(path)
        assert open(path).read() == ""


def test_load_hla_wide_csv_normalizes_samples_and_alleles():
    with tempfile.TemporaryDirectory() as d:
        path = os.path.join(d, "hla_summary.csv")
        with open(path, "w") as fh:
            fh.write("sample,status,reason,hla_a1,hla_a2,hla_b1,hla_b2,hla_c1,hla_c2\n")
            fh.write("CSC0073.TXO_novo,success,completed,A*11:01:01,A*11:01:01,B*55:01:01,NA,C*03:03:01,C*03:03:01\n")
        alleles = _load_hla(Namespace(
            hla_tsv=path,
            hla_sample_column="sample_id",
            hla_allele_column="hla_allele",
        ))
        assert alleles["CSC0073.TXO"] == [
            "HLA-A*11:01",
            "HLA-B*55:01",
            "HLA-C*03:03",
        ]


def test_load_expression_wide_csv_normalizes_duplicate_samples():
    with tempfile.TemporaryDirectory() as d:
        path = os.path.join(d, "expression.csv")
        with open(path, "w") as fh:
            fh.write("gene_id,S1_pmlb,S1_novo,S2\n")
            fh.write("TP53,1.5,3.0,4.0\n")
        expr = _load_expression(Namespace(
            expression_tsv=path,
            expression_sample_column="sample_id",
            expression_gene_column="gene_id",
            expression_transcript_column="transcript_id",
            expression_tpm_column="TPM",
        ))
        assert expr[("S1", "TP53", "")] == 3.0
        assert expr[("S2", "TP53", "")] == 4.0


def test_variant_context_sequence_is_not_full_protein():
    ref = "M" + "A" * 40 + "R" + "G" * 40
    mut = "M" + "A" * 40 + "Q" + "G" * 40
    context, start, end, altered_start, altered_end = _variant_context_sequence(
        mut, ref, flank_aa=5, max_tail=120
    )
    assert context == mut[36:47]
    assert len(context) == 11
    assert start == 36
    assert end == 47
    assert altered_start == 41
    assert altered_end == 42


def test_peptide_candidates_overlap_mutated_residue_only():
    ref = "MABCDE"
    mut = "MABQDE"
    peptides = list(_peptide_candidates(mut, ref, min_len=2, max_len=3, max_tail=120))
    assert peptides
    assert all(
        row["peptide_start_aa"] <= 4 <= row["peptide_end_aa"]
        for row in peptides
    )
    assert all(
        row["mutant_peptide_sequence"] != row["wildtype_peptide_sequence"]
        for row in peptides
    )


if __name__ == "__main__":
    test_header_mutated_with_protein_change()
    test_header_mutated_multiple_variants()
    test_header_reference()
    test_header_reference_ignores_variants()
    test_header_mutated_with_traceability_metadata()
    test_write_fasta_basic()
    test_write_fasta_line_wrapping()
    test_write_fasta_creates_parent_dirs()
    test_write_fasta_empty()
    test_load_hla_wide_csv_normalizes_samples_and_alleles()
    test_load_expression_wide_csv_normalizes_duplicate_samples()
    test_variant_context_sequence_is_not_full_protein()
    test_peptide_candidates_overlap_mutated_residue_only()
    print("All FASTA output tests passed.")
