from __future__ import annotations

import gzip
from pathlib import Path

from pgci_preprocess.fasta import read_cds_fasta, read_pep_fasta


def test_read_ncbi_cds_uses_protein_id(tmp_path: Path) -> None:
    path = tmp_path / "sample.cds.fna.gz"
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        handle.write(
            ">lcl|chr1_cds_XP_1.1_1 [gene=G1] [protein_id=XP_1.1] [location=1..9]\n"
            "ATGGCTTAA\n"
            ">lcl|chr1_cds_XP_2.1_1 [gene=G2] [protein_id=XP_2.1]\n"
            "ATGAAATAA\n"
        )

    result = read_cds_fasta(path)
    assert result.sequences == {"XP_1.1": "ATGGCTTAA", "XP_2.1": "ATGAAATAA"}
    assert result.counts["records"] == 2


def test_read_pep_detects_duplicate_conflict_and_empty_record(tmp_path: Path) -> None:
    path = tmp_path / "sample.pep.faa"
    path.write_text(
        ">P1 description one\nMA\n>P1 description two\nMW\n>P2 empty\n>P3\nMAA\n",
        encoding="utf-8",
    )

    result = read_pep_fasta(path)
    assert result.sequences["P1"] == "MA"
    assert result.sequences["P3"] == "MAA"
    assert result.counts["duplicate_conflict"] == 1
    assert result.counts["empty_sequence"] == 1

