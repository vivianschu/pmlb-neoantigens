"""Tests for DEL, INS, MNV, and INDEL application with coordinate tracking."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))

from mutations import Mutation, apply_mutations_to_sequence


def _mut(pos, ref, alt, vtype, end=None, gene="G1", sample="S1", tid="T1"):
    return Mutation(
        sample_id=sample, gene=gene, transcript_id=tid,
        variant_id=f"{gene}:c.{pos}{ref}>{alt}",
        cdna_start=pos, cdna_end=end if end is not None else pos + max(len(ref), 1) - 1,
        ref_allele=ref, alt_allele=alt, variant_type=vtype,
    )


# ── DEL ────────────────────────────────────────────────────────────────────

def test_del_single_base():
    seq = "ATGGCCTTA"   # A=1 T=2 G=3 G=4 C=5 C=6 T=7 T=8 A=9
    result, reports = apply_mutations_to_sequence(
        seq, [_mut(4, "G", "", "DEL", end=4)]
    )
    assert result == "ATGCCTTA"
    assert reports[0].status == "applied"


def test_del_two_bases():
    seq = "ATGGCCTTA"
    result, reports = apply_mutations_to_sequence(
        seq, [_mut(4, "GC", "", "DEL", end=5)]
    )
    assert result == "ATGCTTA"
    assert reports[0].status == "applied"


def test_del_n_placeholder():
    """Deletion with empty ref (N placeholder) uses coordinate span for length."""
    seq = "ATGGCCTTA"
    # cdna_start=4, cdna_end=5 → delete 2 bases
    result, reports = apply_mutations_to_sequence(
        seq, [_mut(4, "", "", "DEL", end=5)]
    )
    assert len(result) == len(seq) - 2
    assert reports[0].status == "applied"


def test_del_ref_mismatch():
    seq = "ATGGCCTTA"
    result, reports = apply_mutations_to_sequence(
        seq, [_mut(4, "TT", "", "DEL", end=5)], strict_ref_check=True
    )
    assert result == "ATGGCCTTA"
    assert reports[0].status == "skipped_ref_mismatch"


# ── INS ────────────────────────────────────────────────────────────────────

def test_ins_single_base():
    # Insert 'A' after position 3 → between c.3 and c.4
    seq = "ATGGCC"
    result, reports = apply_mutations_to_sequence(
        seq, [_mut(3, "", "A", "INS", end=4)]
    )
    assert result == "ATGAGCC"
    assert reports[0].status == "applied"


def test_ins_three_bases():
    seq = "ATGGCC"
    result, reports = apply_mutations_to_sequence(
        seq, [_mut(3, "", "GAG", "INS", end=4)]
    )
    assert result == "ATGGAGGCC"
    assert len(result) == len(seq) + 3
    assert reports[0].status == "applied"


def test_ins_frameshift():
    seq = "ATGGCCTTA"
    result, reports = apply_mutations_to_sequence(
        seq, [_mut(3, "", "T", "INS", end=4)]
    )
    assert len(result) == len(seq) + 1
    assert reports[0].status == "applied"


# ── MNV ────────────────────────────────────────────────────────────────────

def test_mnv_two_bases():
    seq = "ATGGCC"
    result, reports = apply_mutations_to_sequence(
        seq, [_mut(4, "GC", "TT", "MNV", end=5)]
    )
    assert result == "ATGTTC"
    assert reports[0].status == "applied"


# ── INDEL ──────────────────────────────────────────────────────────────────

def test_indel_shrink():
    seq = "ATGGCC"
    result, reports = apply_mutations_to_sequence(
        seq, [_mut(4, "GC", "T", "INDEL", end=5)]
    )
    assert result == "ATGTC"
    assert reports[0].status == "applied"


def test_indel_grow():
    seq = "ATGGCC"
    result, reports = apply_mutations_to_sequence(
        seq, [_mut(4, "G", "TAA", "INDEL", end=4)]
    )
    assert result == "ATGTAACC"
    assert reports[0].status == "applied"


# ── Coordinate offset propagation ─────────────────────────────────────────

def test_del_then_snv_offset():
    """Deleting 2 bases should shift all subsequent mutations by -2."""
    seq = "ATGGCCTTAG"   # len=10
    del_mut = _mut(4, "GC", "", "DEL", end=5)    # delete pos 4-5
    snv_mut = _mut(9, "A", "C", "SNV", end=9)    # original pos 9 = 'A'
    result, reports = apply_mutations_to_sequence(seq, [del_mut, snv_mut])
    # After del: ATGCTTAG (pos 9 → offset 7 → index 6 = 'A' → 'C')
    assert reports[0].status == "applied"
    assert reports[1].status == "applied"
    assert len(result) == 8   # 10 - 2
    assert result[6] == "C"


def test_ins_then_snv_offset():
    """Inserting 2 bases should shift all subsequent mutations by +2."""
    seq = "ATGGCC"       # len=6
    ins_mut = _mut(3, "", "TT", "INS", end=4)   # insert TT after pos 3
    snv_mut = _mut(5, "C", "A", "SNV", end=5)   # original pos 5 = 'C'
    result, reports = apply_mutations_to_sequence(seq, [ins_mut, snv_mut])
    # After ins: ATGTТGCC → then SNV at original pos 5 (offset +2) → index 6 = 'C' → 'A'
    assert reports[0].status == "applied"
    assert reports[1].status == "applied"
    assert len(result) == 8   # 6 + 2


if __name__ == "__main__":
    test_del_single_base()
    test_del_two_bases()
    test_del_n_placeholder()
    test_del_ref_mismatch()
    test_ins_single_base()
    test_ins_three_bases()
    test_ins_frameshift()
    test_mnv_two_bases()
    test_indel_shrink()
    test_indel_grow()
    test_del_then_snv_offset()
    test_ins_then_snv_offset()
    print("All indel tests passed.")
