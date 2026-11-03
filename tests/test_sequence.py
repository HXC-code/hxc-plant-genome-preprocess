from __future__ import annotations

import pytest

from pgci_preprocess.sequence import translate_cds, validate_cds_pep


def test_translate_standard_code_and_terminal_stop() -> None:
    assert translate_cds("ATGGCTTAA") == "MA*"


def test_translation_accepts_terminal_stop_and_ambiguity() -> None:
    assert validate_cds_pep("ATGGCTTAA", "MA").valid
    assert validate_cds_pep("ATGGCTTAA", "MA*").valid
    assert validate_cds_pep("ATGGCTAAN", "MAX").valid
    assert validate_cds_pep("ATGAAG", "MK").valid
    assert validate_cds_pep("ATGGCTTAA", "MJ").valid


def test_translation_rejects_non_multiple_of_three() -> None:
    result = validate_cds_pep("ATGG", "M")
    assert not result.valid
    assert result.reason == "cds_length_not_multiple_of_3"


def test_translation_rejects_internal_stop_even_with_wildcard_protein() -> None:
    result = validate_cds_pep("ATGTAAGCT", "MXA")
    assert not result.valid
    assert result.reason == "internal_stop"


def test_translation_rejects_invalid_alphabet_and_mismatch() -> None:
    assert validate_cds_pep("ATGGCT1AA", "MA").reason == "invalid_cds_alphabet"
    assert validate_cds_pep("ATGGCTTAA", "MA1").reason == "invalid_protein_alphabet"

    result = validate_cds_pep("ATGGCTTAA", "MW")
    assert not result.valid
    assert result.reason == "translation_mismatch"


@pytest.mark.parametrize("cds,pep", [("", "M"), ("ATG", "")])
def test_translation_rejects_empty_sequences(cds: str, pep: str) -> None:
    assert not validate_cds_pep(cds, pep).valid

