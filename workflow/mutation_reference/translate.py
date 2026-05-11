"""
Nucleotide-to-protein translation with effect classification.

Translation starts from the first base of the ORF (index 0) and proceeds
in triplets. Stop codons are represented as '*'. By default translation
halts at the first stop codon; pass stop_at_first_stop=False to read through.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

# Standard genetic code
CODON_TABLE: dict[str, str] = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

STOP_CODONS = {"TAA", "TAG", "TGA"}


@dataclass
class TranslationResult:
    protein_sequence: str
    stop_codon_position: Optional[int]   # 1-based AA position where stop was hit
    has_start_loss: bool
    has_premature_stop: bool
    has_stop_loss: bool
    has_frameshift: bool
    is_synonymous: bool
    mutation_effect: str
    # missense | synonymous | nonsense | frameshift |
    # in_frame_ins | in_frame_del | stop_loss | start_loss | unknown


def translate(
    sequence: str,
    reference_protein: Optional[str] = None,
    ref_orf_length_nt: Optional[int] = None,
    stop_at_first_stop: bool = True,
    ref_stop_codon_position: Optional[int] = None,
) -> TranslationResult:
    """
    Translate a nucleotide sequence to protein.

    Args:
        sequence:          Nucleotide string (reading frame starts at index 0).
        reference_protein: Translated reference protein for effect comparison.
        ref_orf_length_nt: Length of the original (pre-mutation) ORF in nt,
                           used to detect frameshifts when reference_protein
                           is not supplied.
        stop_at_first_stop: If True (default), stop at the first stop codon.
    """
    sequence = sequence.upper()
    aas: list[str] = []
    stop_aa_pos: Optional[int] = None

    for i in range(0, len(sequence) - 2, 3):
        codon = sequence[i: i + 3]
        if len(codon) < 3:
            break
        aa = CODON_TABLE.get(codon, "X")
        if aa == "*":
            stop_aa_pos = len(aas) + 1   # 1-based
            if stop_at_first_stop:
                break
            aas.append(aa)
        else:
            aas.append(aa)

    protein = "".join(aas)

    # Start-loss: first codon is not ATG or first AA is not M
    has_start_loss = not sequence.startswith("ATG") or (bool(protein) and protein[0] != "M")

    # Frameshift: mutated ORF length differs from reference by a non-multiple of 3
    has_frameshift = False
    if ref_orf_length_nt is not None:
        delta = len(sequence) - ref_orf_length_nt
        has_frameshift = delta % 3 != 0

    # Premature stop: stop codon appears before the reference protein ends
    ref_aa_len: Optional[int] = None
    if reference_protein is not None:
        ref_aa_len = len(reference_protein)
    elif ref_orf_length_nt is not None:
        ref_aa_len = ref_orf_length_nt // 3

    has_premature_stop = (
        stop_aa_pos is not None
        and ref_aa_len is not None
        and len(protein) < ref_aa_len
    )

    # Stop-loss: reference had a stop codon and mutant lost it, extending the protein.
    # Requires ref_stop_codon_position to be set (from translating the reference
    # ORF first) so we don't confuse an in-frame insertion with a stop-loss.
    has_stop_loss = (
        stop_aa_pos is None
        and ref_stop_codon_position is not None
        and ref_aa_len is not None
        and len(protein) > ref_aa_len
    )

    # Effect classification (priority order matters)
    if has_start_loss:
        effect = "start_loss"
    elif has_frameshift:
        effect = "frameshift"
    elif has_stop_loss:
        effect = "stop_loss"
    elif has_premature_stop:
        effect = "nonsense"
    elif reference_protein is not None:
        if protein == reference_protein:
            effect = "synonymous"
        elif len(protein) > len(reference_protein):
            effect = "in_frame_ins"
        elif len(protein) < len(reference_protein):
            effect = "in_frame_del"
        else:
            effect = "missense"
    else:
        effect = "unknown"

    return TranslationResult(
        protein_sequence=protein,
        stop_codon_position=stop_aa_pos,
        has_start_loss=has_start_loss,
        has_premature_stop=has_premature_stop,
        has_stop_loss=has_stop_loss,
        has_frameshift=has_frameshift,
        is_synonymous=(effect == "synonymous"),
        mutation_effect=effect,
    )
