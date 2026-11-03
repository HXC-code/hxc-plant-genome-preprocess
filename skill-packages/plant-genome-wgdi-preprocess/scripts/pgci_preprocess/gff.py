from __future__ import annotations

import re
from collections import Counter, defaultdict
from pathlib import Path
from urllib.parse import unquote

from pgci_preprocess.fasta import open_text
from pgci_preprocess.models import Anomaly, Gene, GffData, Region, Transcript

REGION_TYPES = {"region"}
TRANSCRIPT_TYPES = {
    "mRNA",
    "transcript",
    "pseudogenic_transcript",
    "lnc_RNA",
    "tRNA",
    "rRNA",
    "snRNA",
    "snoRNA",
    "ncRNA",
    "miRNA",
}
ORGANELLE_TYPES = {"chloroplast", "mitochondrion", "plastid", "apicoplast"}
UNKNOWN_VALUES = {"", ".", "na", "n/a", "nan", "none", "null", "unknown"}
TAXON_RE = re.compile(r"(?:^|,)taxon:(\d+)")


def parse_attributes(field: str) -> dict[str, str]:
    attributes: dict[str, str] = {}
    for item in field.split(";"):
        if not item:
            continue
        key, separator, value = item.partition("=")
        key = key.strip()
        if not key:
            continue
        attributes[key] = unquote(value.strip()) if separator else ""
    return attributes


def extract_taxon(attributes: dict[str, str]) -> str:
    match = TAXON_RE.search(attributes.get("Dbxref", ""))
    return match.group(1) if match else ""


def strip_prefix(value: str, prefix: str) -> str:
    return value[len(prefix) :] if value.startswith(prefix) else value


def _signed_int(value: str, line_number: int, anomalies: list[Anomaly]) -> int | None:
    try:
        return int(value)
    except ValueError:
        anomalies.append(
            Anomaly("warning", "invalid_coordinate", f"line={line_number} value={value!r}")
        )
        return None


def _register_region(
    region: Region,
    regions: dict[str, Region],
    anomalies: list[Anomaly],
) -> None:
    previous = regions.get(region.seqid)
    if previous and previous.length != region.length:
        anomalies.append(
            Anomaly(
                "fatal",
                "conflicting_region_length",
                f"{region.seqid}: {previous.length} != {region.length}",
            )
        )
        return
    if previous:
        if (
            previous.genome_type == "genomic"
            and not previous.chromosome
            and not previous.taxon_id
            and not previous.is_primary
        ):
            previous.genome_type = region.genome_type
            previous.chromosome = region.chromosome
            previous.taxon_id = region.taxon_id
            previous.is_primary = region.is_primary
            return
        anomalies.append(
            Anomaly("warning", "duplicate_region", f"duplicate region row for {region.seqid}")
        )
        return
    regions[region.seqid] = region


def parse_gff(path: Path | str) -> GffData:
    path = Path(path)
    regions: dict[str, Region] = {}
    genes: dict[str, Gene] = {}
    transcripts: dict[str, Transcript] = {}
    transcript_links: dict[str, set[str]] = defaultdict(set)
    transcript_cds_counts: Counter[str] = Counter()
    transcript_exon_counts: Counter[str] = Counter()
    counts: Counter[str] = Counter()
    anomalies: list[Anomaly] = []
    assembly_build = ""
    annotation_source = ""
    source_order = 0

    with open_text(path) as handle:
        for line_number, raw_line in enumerate(handle, 1):
            if raw_line.startswith("#"):
                if raw_line.startswith("##sequence-region"):
                    parts = raw_line.rstrip("\n").split()
                    if len(parts) >= 4:
                        try:
                            length = int(parts[3])
                        except ValueError:
                            anomalies.append(
                                Anomaly(
                                    "warning",
                                    "invalid_sequence_region",
                                    f"line={line_number} {raw_line.strip()}",
                                )
                            )
                            continue
                        _register_region(
                            Region(
                                seqid=parts[1],
                                length=length,
                                genome_type="genomic",
                                chromosome="",
                                taxon_id="",
                                source_order=source_order,
                                is_primary=False,
                            ),
                            regions,
                            anomalies,
                        )
                elif raw_line.startswith("#!genome-build "):
                    assembly_build = raw_line.split(None, 1)[1].strip()
                elif raw_line.startswith("#!annotation-source "):
                    annotation_source = raw_line.split(None, 1)[1].strip()
                continue

            if not raw_line.strip():
                continue
            fields = raw_line.rstrip("\r\n").split("\t")
            if len(fields) != 9:
                counts["malformed_line"] += 1
                anomalies.append(
                    Anomaly("warning", "gff_column_count", f"line={line_number} columns={len(fields)}")
                )
                continue

            seqid, _source, feature_type, start_text, end_text, _score, strand, _phase, attr_text = fields
            counts[feature_type] += 1
            attributes = parse_attributes(attr_text)

            if feature_type in REGION_TYPES:
                start = _signed_int(start_text, line_number, anomalies)
                end = _signed_int(end_text, line_number, anomalies)
                if start is None or end is None or start < 1 or end < start:
                    continue
                genome_type = attributes.get("genome", "genomic").strip().lower() or "genomic"
                chromosome = attributes.get("chromosome", "").strip()
                chromosome_known = chromosome.casefold() not in UNKNOWN_VALUES
                is_primary = (
                    genome_type not in ORGANELLE_TYPES
                    and (genome_type == "chromosome" or chromosome_known)
                )
                region = Region(
                    seqid=seqid,
                    length=end - start + 1,
                    genome_type=genome_type,
                    chromosome=chromosome,
                    taxon_id=extract_taxon(attributes),
                    source_order=source_order,
                    is_primary=is_primary,
                )
                _register_region(region, regions, anomalies)
                source_order += 1
                continue

            if feature_type in {"gene", "pseudogene"}:
                gene_id = strip_prefix(attributes.get("ID", ""), "gene-")
                if not gene_id:
                    anomalies.append(
                        Anomaly("warning", "gene_without_id", f"{seqid}:{start_text}-{end_text}")
                    )
                    continue
                start = _signed_int(start_text, line_number, anomalies)
                end = _signed_int(end_text, line_number, anomalies)
                if start is None or end is None or start < 1 or end < start:
                    continue
                biotype = (
                    attributes.get("gene_biotype")
                    or attributes.get("gbkey")
                    or ("pseudogene" if feature_type == "pseudogene" else "unknown")
                )
                name = attributes.get("Name") or attributes.get("gene") or gene_id
                gene = Gene(
                    source_id=gene_id,
                    seqid=seqid,
                    start=start,
                    end=end,
                    strand=strand,
                    name=name,
                    biotype=biotype,
                    feature_type=feature_type,
                )
                previous_gene = genes.get(gene_id)
                if previous_gene:
                    previous_region = regions.get(previous_gene.seqid)
                    current_region = regions.get(gene.seqid)
                    both_organelle = bool(
                        previous_region
                        and current_region
                        and previous_region.is_organelle
                        and current_region.is_organelle
                    )
                    level = (
                        "warning"
                        if (previous_gene.seqid, previous_gene.start, previous_gene.end)
                        == (gene.seqid, gene.start, gene.end)
                        or both_organelle
                        else "fatal"
                    )
                    anomalies.append(
                        Anomaly(
                            level,
                            "duplicate_organelle_gene" if both_organelle else "duplicate_gene",
                            f"{gene_id}:{gene.seqid}:{gene.start}-{gene.end} line={line_number}",
                        )
                    )
                    continue
                genes[gene_id] = gene
                continue

            if feature_type in TRANSCRIPT_TYPES:
                transcript_id = strip_prefix(attributes.get("ID", ""), "rna-")
                parent = attributes.get("Parent", "").split(",", 1)[0].strip()
                gene_id = strip_prefix(parent, "gene-")
                if not transcript_id or not gene_id:
                    anomalies.append(
                        Anomaly(
                            "warning",
                            "transcript_without_id_or_parent",
                            f"{seqid}:{start_text}-{end_text}",
                        )
                    )
                    continue
                start = _signed_int(start_text, line_number, anomalies)
                end = _signed_int(end_text, line_number, anomalies)
                if start is None or end is None or start < 1 or end < start:
                    continue
                transcript_record = Transcript(
                    source_id=transcript_id,
                    gene_id=gene_id,
                    seqid=seqid,
                    start=start,
                    end=end,
                    strand=strand,
                    feature_type=feature_type,
                    transcript_id=attributes.get("transcript_id", ""),
                )
                previous_transcript = transcripts.get(transcript_id)
                if previous_transcript:
                    anomalies.append(
                        Anomaly("warning", "duplicate_transcript", f"{transcript_id} line={line_number}")
                    )
                    continue
                transcripts[transcript_id] = transcript_record
                continue

            if feature_type == "CDS":
                parents = [strip_prefix(item.strip(), "rna-") for item in attributes.get("Parent", "").split(",") if item.strip()]
                protein_id = attributes.get("protein_id", "")
                if not protein_id:
                    protein_id = strip_prefix(attributes.get("ID", ""), "cds-")
                for transcript_id in parents:
                    transcript_links[transcript_id].add(protein_id)
                    transcript_cds_counts[transcript_id] += 1
                if not parents:
                    anomalies.append(
                        Anomaly("warning", "cds_without_parent", f"{seqid}:{start_text}-{end_text}")
                    )
                continue

            if feature_type == "exon":
                for parent in attributes.get("Parent", "").split(","):
                    parent = strip_prefix(parent.strip(), "rna-")
                    if parent:
                        transcript_exon_counts[parent] += 1

    for transcript_id, protein_ids in transcript_links.items():
        linked_transcript = transcripts.get(transcript_id)
        if linked_transcript is None:
            anomalies.append(
                Anomaly(
                    "warning",
                    "cds_parent_missing_transcript",
                    f"transcript={transcript_id}",
                )
            )
            continue
        linked_transcript.protein_ids = {
            protein_id for protein_id in protein_ids if protein_id
        }
        linked_transcript.cds_count = transcript_cds_counts[transcript_id]
        linked_transcript.exon_count = transcript_exon_counts[transcript_id]

    return GffData(
        regions=regions,
        genes=genes,
        transcripts=transcripts,
        feature_counts=counts,
        anomalies=anomalies,
        assembly_build=assembly_build,
        annotation_source=annotation_source,
    )

