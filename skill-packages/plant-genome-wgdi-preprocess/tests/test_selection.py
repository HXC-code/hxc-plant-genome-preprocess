from __future__ import annotations

from collections import Counter

from pgci_preprocess.models import Candidate, Region
from pgci_preprocess.selection import choose_representative, select_analysis_sequences


def _candidate(transcript: str, protein_length: int, cds_length: int) -> Candidate:
    return Candidate(
        source_gene_id="G",
        source_gene_name="G",
        gene_biotype="protein_coding",
        source_transcript_id=transcript,
        source_protein_id=transcript + "-P",
        seqid="chr1",
        start=1,
        end=100,
        strand="+",
        cds_length=cds_length,
        protein_length=protein_length,
        translation_exact=True,
        translation_mismatches=0,
    )


def test_representative_is_deterministic_and_prefers_longest_protein() -> None:
    chosen = choose_representative(
        [
            _candidate("T3", 10, 33),
            _candidate("T1", 20, 60),
            _candidate("T2", 20, 63),
        ]
    )
    assert chosen.source_transcript_id == "T2"


def _region(seqid: str, *, primary: bool, organelle: bool = False) -> Region:
    return Region(
        seqid=seqid,
        length=1000,
        genome_type="mitochondrion" if organelle else ("chromosome" if primary else "genomic"),
        chromosome="1" if primary else "",
        taxon_id="1",
        source_order=0,
        is_primary=primary,
    )


def test_sequences_use_primary_chromosomes_and_ignore_scaffolds() -> None:
    regions = {
        "chr1": _region("chr1", primary=True),
        "scaffold": _region("scaffold", primary=False),
        "mt": _region("mt", primary=False, organelle=True),
    }
    selection = select_analysis_sequences(
        regions,
        Counter({"chr1": 100, "scaffold": 1000, "mt": 20}),
        min_genes=50,
        target_coverage=0.9,
        max_sequences=200,
    )
    assert selection.keep == {"chr1"}
    assert selection.mode == "primary_chromosome"
    assert selection.ratio == 0.0909
    assert selection.ready == "no"


def test_fragmented_assembly_falls_back_to_gene_coverage() -> None:
    regions = {
        "s1": _region("s1", primary=False),
        "s2": _region("s2", primary=False),
        "s3": _region("s3", primary=False),
        "mt": _region("mt", primary=False, organelle=True),
    }
    selection = select_analysis_sequences(
        regions,
        Counter({"s1": 45, "s2": 35, "s3": 20, "mt": 10}),
        min_genes=50,
        target_coverage=0.9,
        max_sequences=200,
    )
    assert selection.keep == {"s1", "s2", "s3"}
    assert selection.mode == "coverage"
    assert selection.ratio == 1.0
    assert selection.ready == "yes"

