from __future__ import annotations

from dataclasses import dataclass

DNA_ALPHABET = frozenset("ACGTRYSWKMBDHVN")
PROTEIN_ALPHABET = frozenset("ACDEFGHIKLMNPQRSTVWYBXZUOJ*")
BASES = "TCAG"
AMINO_ACIDS = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODON_TABLE = {
    first + second + third: AMINO_ACIDS[index]
    for index, (first, second, third) in enumerate(
        (a, b, c) for a in BASES for b in BASES for c in BASES
    )
}


@dataclass(frozen=True, slots=True)
class TranslationCheck:
    valid: bool
    reason: str
    exact: bool
    mismatches: int
    translated: str


def translate_cds(sequence: str) -> str:
    sequence = sequence.upper()
    return "".join(
        CODON_TABLE.get(sequence[index : index + 3], "X")
        for index in range(0, len(sequence) - 2, 3)
    )


def _strip_terminal_stop(sequence: str) -> str:
    return sequence[:-1] if sequence.endswith("*") else sequence


def validate_cds_pep(cds: str, pep: str) -> TranslationCheck:
    cds = cds.strip().upper().replace("-", "")
    pep = pep.strip().upper().replace("-", "")
    if not cds:
        return TranslationCheck(False, "empty_cds", False, 0, "")
    if not pep:
        return TranslationCheck(False, "empty_protein", False, 0, "")
    if not set(cds) <= DNA_ALPHABET:
        return TranslationCheck(False, "invalid_cds_alphabet", False, 0, "")
    if not set(pep) <= PROTEIN_ALPHABET:
        return TranslationCheck(False, "invalid_protein_alphabet", False, 0, "")
    if len(cds) % 3:
        return TranslationCheck(False, "cds_length_not_multiple_of_3", False, 0, "")

    translated = translate_cds(cds)
    translated_clean = _strip_terminal_stop(translated)
    protein_clean = _strip_terminal_stop(pep)
    if "*" in translated_clean:
        return TranslationCheck(False, "internal_stop", False, 0, translated)
    if "*" in protein_clean:
        return TranslationCheck(False, "internal_stop", False, 0, translated)
    if len(translated_clean) != len(protein_clean):
        return TranslationCheck(False, "translation_length_mismatch", False, 0, translated)

    mismatches = sum(
        1
        for translated_aa, protein_aa in zip(translated_clean, protein_clean)
        if not (
            translated_aa == protein_aa
            or translated_aa == "X"
            or protein_aa in {"X", "B", "Z", "J"}
        )
    )
    if mismatches:
        return TranslationCheck(False, "translation_mismatch", False, mismatches, translated)
    return TranslationCheck(
        True,
        "",
        translated_clean == protein_clean,
        mismatches,
        translated,
    )

