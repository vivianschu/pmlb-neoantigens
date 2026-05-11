"""Tests for translation correctness and mutation effect classification."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))

from translate import translate


def test_basic_translation():
    result = translate("ATGGCCTTA")
    assert result.protein_sequence == "MAL"


def test_stop_codon_halts_translation():
    # ATG GCC TAA TTA → M A stop
    result = translate("ATGGCCTAATTA", stop_at_first_stop=True)
    assert result.protein_sequence == "MA"
    assert result.stop_codon_position == 3


def test_translate_through_stop():
    # ATG GCC TAA TTT → M A * F
    result = translate("ATGGCCTAATTT", stop_at_first_stop=False)
    assert "*" in result.protein_sequence
    assert result.protein_sequence.startswith("MA*")


def test_frameshift_insertion():
    ref_seq = "ATGGCCTTA"       # MAL (9 nt)
    mut_seq  = "ATGGCCTTTA"     # +1 → not divisible by 3 relative to ref
    ref_result = translate(ref_seq)
    mut_result = translate(mut_seq,
                           reference_protein=ref_result.protein_sequence,
                           ref_orf_length_nt=len(ref_seq))
    assert mut_result.has_frameshift
    assert mut_result.mutation_effect == "frameshift"


def test_frameshift_deletion():
    ref_seq = "ATGGCCTTAAAG"    # MALK (12 nt)
    mut_seq  = "ATGGCTTAAAG"    # del one base at pos 4 → frameshift
    ref_result = translate(ref_seq)
    mut_result = translate(mut_seq,
                           reference_protein=ref_result.protein_sequence,
                           ref_orf_length_nt=len(ref_seq))
    assert mut_result.has_frameshift
    assert mut_result.mutation_effect == "frameshift"


def test_inframe_deletion():
    ref_seq = "ATGGCCTTAAAG"    # MALK (12 nt)
    mut_seq  = "ATGTTAAAG"      # del 3 bases (GCC) → in-frame, 9 nt
    ref_result = translate(ref_seq)
    mut_result = translate(mut_seq,
                           reference_protein=ref_result.protein_sequence,
                           ref_orf_length_nt=len(ref_seq))
    assert not mut_result.has_frameshift
    assert mut_result.mutation_effect == "in_frame_del"
    assert mut_result.protein_sequence == "MLK"


def test_inframe_insertion():
    ref_seq = "ATGGCCAAG"       # MAK (9 nt)
    mut_seq  = "ATGGCCTTAAAG"   # +3 nt → MALK
    ref_result = translate(ref_seq)
    mut_result = translate(mut_seq,
                           reference_protein=ref_result.protein_sequence,
                           ref_orf_length_nt=len(ref_seq))
    assert not mut_result.has_frameshift
    assert mut_result.mutation_effect == "in_frame_ins"
    assert mut_result.protein_sequence == "MALK"


def test_premature_stop():
    ref_seq = "ATGGCCTTAAAG"    # MALK (12 nt)
    mut_seq  = "ATGTAATTAAAG"   # ATG TAA → M stop
    ref_result = translate(ref_seq)
    mut_result = translate(mut_seq,
                           reference_protein=ref_result.protein_sequence,
                           ref_orf_length_nt=len(ref_seq))
    assert mut_result.has_premature_stop
    assert mut_result.mutation_effect == "nonsense"


def test_synonymous():
    ref_seq = "ATGGCC"   # MA
    mut_seq  = "ATGGCT"  # MA (GCC → GCT, both Ala)
    ref_result = translate(ref_seq)
    mut_result = translate(mut_seq,
                           reference_protein=ref_result.protein_sequence,
                           ref_orf_length_nt=len(ref_seq))
    assert mut_result.protein_sequence == ref_result.protein_sequence
    assert mut_result.is_synonymous
    assert mut_result.mutation_effect == "synonymous"


def test_missense():
    ref_seq = "ATGGCC"   # MA
    mut_seq  = "ATGACC"  # MT (GCC → ACC, Ala → Thr)
    ref_result = translate(ref_seq)
    mut_result = translate(mut_seq,
                           reference_protein=ref_result.protein_sequence,
                           ref_orf_length_nt=len(ref_seq))
    assert mut_result.protein_sequence == "MT"
    assert not mut_result.is_synonymous
    assert mut_result.mutation_effect == "missense"


def test_start_loss():
    result = translate("TTGGCC")  # no ATG
    assert result.has_start_loss
    assert result.mutation_effect == "start_loss"


def test_stop_loss():
    # Reference ORF ends with a stop codon; mutation converts stop → sense codon
    ref_seq = "ATGGCCTAA"   # MA + stop → protein = "MA", stop_codon_position = 3
    mut_seq  = "ATGGCCTTT"  # MA + F → no stop codon → stop loss (protein = "MAF")
    ref_result = translate(ref_seq)
    assert ref_result.stop_codon_position == 3
    mut_result = translate(mut_seq,
                           reference_protein=ref_result.protein_sequence,
                           ref_orf_length_nt=len(ref_seq),
                           ref_stop_codon_position=ref_result.stop_codon_position)
    assert mut_result.has_stop_loss
    assert mut_result.mutation_effect == "stop_loss"


def test_unknown_codons_become_X():
    result = translate("ATGNNN")  # NNN is not in codon table → 'X'
    assert "X" in result.protein_sequence


if __name__ == "__main__":
    test_basic_translation()
    test_stop_codon_halts_translation()
    test_translate_through_stop()
    test_frameshift_insertion()
    test_frameshift_deletion()
    test_inframe_deletion()
    test_inframe_insertion()
    test_premature_stop()
    test_synonymous()
    test_missense()
    test_start_loss()
    test_stop_loss()
    test_unknown_codons_become_X()
    print("All translation tests passed.")
