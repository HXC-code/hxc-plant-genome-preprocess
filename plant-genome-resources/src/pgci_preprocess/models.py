from __future__ import annotations

from collections import Counter
from dataclasses import dataclass, field


@dataclass(slots=True)
class Anomaly:
    level: str
    code: str
    detail: str


@dataclass(slots=True)
class Region:
    seqid: str
    length: int
    genome_type: str
    chromosome: str
    taxon_id: str
    source_order: int
    is_primary: bool

    @property
    def is_organelle(self) -> bool:
        return self.genome_type in {"chloroplast", "mitochondrion", "plastid", "apicoplast"}


@dataclass(slots=True)
class Gene:
    source_id: str
    seqid: str
    start: int
    end: int
    strand: str
    name: str
    biotype: str
    feature_type: str


@dataclass(slots=True)
class Transcript:
    source_id: str
    gene_id: str
    seqid: str
    start: int
    end: int
    strand: str
    feature_type: str
    transcript_id: str = ""
    protein_ids: set[str] = field(default_factory=set)
    cds_count: int = 0
    exon_count: int = 0


@dataclass(slots=True)
class GffData:
    regions: dict[str, Region]
    genes: dict[str, Gene]
    transcripts: dict[str, Transcript]
    feature_counts: Counter[str]
    anomalies: list[Anomaly]
    assembly_build: str = ""
    annotation_source: str = ""


@dataclass(slots=True)
class FastaData:
    sequences: dict[str, str]
    counts: Counter[str]
    anomalies: list[Anomaly]


@dataclass(slots=True)
class Candidate:
    source_gene_id: str
    source_gene_name: str
    gene_biotype: str
    source_transcript_id: str
    source_protein_id: str
    seqid: str
    start: int
    end: int
    strand: str
    cds_length: int
    protein_length: int
    translation_exact: bool
    translation_mismatches: int


@dataclass(slots=True)
class SequenceSelection:
    keep: set[str]
    mode: str
    selected_genes: int
    total_genes: int
    ratio: float
    ready: str

