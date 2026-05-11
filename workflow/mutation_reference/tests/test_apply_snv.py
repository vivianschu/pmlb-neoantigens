"""Tests for SNV application logic."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))

from mutations import Mutation, apply_mutations_to_sequence
from mutations import load_mutations_maf


def _snv(pos, ref, alt, gene="G1", sample="S1", tid="T1"):
    return Mutation(
        sample_id=sample, gene=gene, transcript_id=tid,
        variant_id=f"{gene}:c.{pos}{ref}>{alt}",
        cdna_start=pos, cdna_end=pos,
        ref_allele=ref, alt_allele=alt, variant_type="SNV",
    )


def test_snv_basic():
    seq = "ATGGCCTTA"
    result, reports = apply_mutations_to_sequence(seq, [_snv(4, "G", "A")])
    assert result == "ATGACCTTA"
    assert reports[0].status == "applied"


def test_snv_last_base():
    seq = "ATGGCC"
    result, reports = apply_mutations_to_sequence(seq, [_snv(6, "C", "T")])
    assert result == "ATGGCT"
    assert reports[0].status == "applied"


def test_snv_n_placeholder_ref_is_accepted():
    """HGVSc with N ref (empty ref_allele) should always apply."""
    seq = "ATGGCCTTA"
    mut = _snv(4, "", "A")   # empty ref = N placeholder
    result, reports = apply_mutations_to_sequence(seq, [mut], strict_ref_check=True)
    assert result == "ATGACCTTA"
    assert reports[0].status == "applied"


def test_snv_ref_mismatch_strict():
    seq = "ATGGCCTTA"
    result, reports = apply_mutations_to_sequence(
        seq, [_snv(4, "T", "A")], strict_ref_check=True
    )
    assert result == "ATGGCCTTA"
    assert reports[0].status == "skipped_ref_mismatch"


def test_snv_ref_mismatch_lenient():
    seq = "ATGGCCTTA"
    result, reports = apply_mutations_to_sequence(
        seq, [_snv(4, "T", "A")], strict_ref_check=False
    )
    assert result == "ATGACCTTA"
    assert reports[0].status == "applied"


def test_snv_out_of_bounds():
    seq = "ATGGCC"
    result, reports = apply_mutations_to_sequence(seq, [_snv(100, "A", "T")])
    assert result == "ATGGCC"
    assert reports[0].status == "skipped_invalid_cdna_position"


def test_multiple_snvs_no_offset():
    seq = "ATGGCCTTA"
    muts = [_snv(4, "G", "A"), _snv(7, "T", "C")]
    result, reports = apply_mutations_to_sequence(seq, muts)
    assert result == "ATGACCCTA"
    assert all(r.status == "applied" for r in reports)


def test_snv_first_position():
    seq = "ATGGCC"
    result, reports = apply_mutations_to_sequence(seq, [_snv(1, "A", "C")])
    assert result == "CTGGCC"
    assert reports[0].status == "applied"


def test_load_mutations_maf_drops_silent(tmp_path):
    maf = tmp_path / "muts.maf"
    maf.write_text(
        "Tumor_Sample_Barcode\tHugo_Symbol\tTranscript_ID\t"
        "Variant_Classification\tVariant_Type\tReference_Allele\t"
        "Tumor_Seq_Allele2\tHGVSc\tHGVSp_Short\tt_alt_count\tt_depth\n"
        "S1\tTP53\tENST1\tSilent\tSNP\tC\tT\tENST1:c.6C>T\tp.A2A\t4\t20\n"
        "S1\tTP53\tENST1\tMissense_Mutation\tSNP\tG\tA\tENST1:c.4G>A\tp.A2T\t5\t20\n"
    )
    muts = load_mutations_maf(str(maf))
    assert len(muts) == 1
    assert muts[0].hgvsp == "p.A2T"


if __name__ == "__main__":
    test_snv_basic()
    test_snv_last_base()
    test_snv_n_placeholder_ref_is_accepted()
    test_snv_ref_mismatch_strict()
    test_snv_ref_mismatch_lenient()
    test_snv_out_of_bounds()
    test_multiple_snvs_no_offset()
    test_snv_first_position()
    import tempfile
    from pathlib import Path
    with tempfile.TemporaryDirectory() as d:
        test_load_mutations_maf_drops_silent(Path(d))
    print("All SNV tests passed.")
