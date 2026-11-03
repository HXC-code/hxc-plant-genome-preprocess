from __future__ import annotations

import csv
import hashlib
import json
import os
import re
import shutil
import uuid
from collections import Counter, defaultdict
from collections.abc import Mapping, Sequence
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

from pgci_preprocess.fasta import read_cds_fasta, read_pep_fasta
from pgci_preprocess.gff import parse_gff
from pgci_preprocess.models import Candidate
from pgci_preprocess.selection import (
    choose_representative,
    natural_key,
    region_sort_key,
    select_analysis_sequences,
)
from pgci_preprocess.sequence import validate_cds_pep

ACCESSION_RE = re.compile(r"^(GC[AF]_\d+\.\d+)_")
PIPELINE_VERSION = "2.0.0"


@dataclass(slots=True)
class SpeciesInput:
    show_id: str
    code: str
    species: str
    order: str
    family: str
    directory: Path
    gff_name: str
    cds_name: str
    pep_name: str
    accession: str


@dataclass(slots=True)
class PipelineOptions:
    project_root: Path
    clean_dirname: str = "clean_v2"
    analysis_dirname: str = "analysis_v2"
    min_genes_per_sequence: int = 50
    target_gene_coverage: float = 0.90
    max_sequences: int = 200
    overwrite: bool = True


class PipelineFailure(RuntimeError):
    def __init__(self, reason: str, detail: str = "") -> None:
        super().__init__(detail or reason)
        self.reason = reason
        self.detail = detail


def _accession_from_filename(name: str) -> str:
    match = ACCESSION_RE.match(name)
    return match.group(1) if match else ""


def discover_species(project_root: Path | str) -> list[SpeciesInput]:
    project_root = Path(project_root)
    inventory = project_root / "resources" / "download_inventory.tsv"
    if not inventory.is_file():
        raise FileNotFoundError(f"download inventory not found: {inventory}")

    downloads = project_root / "data" / "downloads"
    directories = {
        path.name.split("_", 1)[0]: path
        for path in downloads.iterdir()
        if path.is_dir() and "_" in path.name
    }
    species: list[SpeciesInput] = []
    with inventory.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            show_id = row.get("showid", "").strip()
            if not show_id or show_id not in directories:
                continue
            files = [row.get("gff", ""), row.get("cds", ""), row.get("pep", "")]
            accessions = {_accession_from_filename(name) for name in files if name}
            accessions.discard("")
            accession = next(iter(accessions)) if len(accessions) == 1 else "NA"
            species.append(
                SpeciesInput(
                    show_id=show_id,
                    code=f"s{show_id}",
                    species=row.get("Species", "").strip(),
                    order=row.get("Order", "").strip(),
                    family=row.get("Family", "").strip(),
                    directory=directories[show_id],
                    gff_name=row.get("gff", "").strip(),
                    cds_name=row.get("cds", "").strip(),
                    pep_name=row.get("pep", "").strip(),
                    accession=accession,
                )
            )
    species.sort(key=lambda item: int(item.show_id) if item.show_id.isdigit() else 10**12)
    return species


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _remove_staging(path: Path, parent: Path) -> None:
    path = path.resolve()
    parent = parent.resolve()
    if path.parent != parent or not path.name.startswith("."):
        raise ValueError(f"refusing to remove unexpected staging path: {path}")
    if path.exists():
        shutil.rmtree(path)


def _replace_directory(staging: Path, destination: Path, overwrite: bool) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists() and not overwrite:
        raise FileExistsError(f"output already exists: {destination}")
    backup: Path | None = None
    if destination.exists():
        backup = destination.with_name(f".{destination.name}.previous-{uuid.uuid4().hex}")
        os.replace(destination, backup)
    try:
        os.replace(staging, destination)
    except Exception:
        if backup is not None and backup.exists() and not destination.exists():
            os.replace(backup, destination)
        raise
    if backup is not None and backup.exists():
        shutil.rmtree(backup)


def _write_fasta(path: Path, records: list[dict[str, object]]) -> None:
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        for record in records:
            canonical_id = str(record["canonical_gene_id"])
            source_gene = str(record["source_gene_id"])
            source_protein = str(record["source_protein_id"])
            handle.write(f">{canonical_id} source_gene={source_gene} source_protein={source_protein}\n")
            sequence = str(record["sequence"])
            for index in range(0, len(sequence), 80):
                handle.write(sequence[index : index + 80] + "\n")


def _write_gff(path: Path, entries: list[dict[str, object]]) -> None:
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        for entry in entries:
            handle.write(
                "\t".join(
                    [
                        str(entry["clean_sequence_id"]),
                        str(entry["canonical_gene_id"]),
                        str(entry["start"]),
                        str(entry["end"]),
                        str(entry["strand"]),
                        str(entry["order"]),
                        str(entry["source_gene_id"]),
                    ]
                )
                + "\n"
            )


def _write_lens(path: Path, region_rows: list[dict[str, object]]) -> None:
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        for row in region_rows:
            handle.write(
                f"{row['clean_sequence_id']}\t{row['length']}\t{row['gene_count']}\n"
            )


def _read_fasta_ids(path: Path) -> list[str]:
    with path.open(encoding="utf-8") as handle:
        return [
            line[1:].split(None, 1)[0]
            for line in handle
            if line.startswith(">") and len(line) > 1
        ]


def _read_lens(path: Path) -> dict[str, tuple[int, int]]:
    result: dict[str, tuple[int, int]] = {}
    with path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, 1):
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) != 3:
                raise PipelineFailure("invalid_lens", f"{path.name}:{line_number}")
            sequence_id, length, gene_count = fields
            result[sequence_id] = (int(length), int(gene_count))
    return result


def _validate_bundle(
    gff_path: Path,
    lens_path: Path,
    cds_path: Path,
    pep_path: Path,
) -> dict[str, object]:
    pep_ids = _read_fasta_ids(pep_path)
    cds_ids = _read_fasta_ids(cds_path)
    if len(pep_ids) != len(set(pep_ids)):
        raise PipelineFailure("duplicate_pep_id", str(pep_path))
    if len(cds_ids) != len(set(cds_ids)):
        raise PipelineFailure("duplicate_cds_id", str(cds_path))

    lens = _read_lens(lens_path)
    gff_ids: list[str] = []
    orders: dict[str, list[int]] = defaultdict(list)
    per_sequence_counts: Counter[str] = Counter()
    with gff_path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, 1):
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) != 7:
                raise PipelineFailure("invalid_gff_column_count", f"{gff_path.name}:{line_number}")
            sequence_id, gene_id, start, end, strand, order_text, _source_gene = fields
            if sequence_id not in lens:
                raise PipelineFailure("gff_sequence_missing_from_lens", sequence_id)
            start_int, end_int, order = int(start), int(end), int(order_text)
            if start_int < 1 or end_int < start_int or end_int > lens[sequence_id][0]:
                raise PipelineFailure("invalid_gff_coordinates", gene_id)
            if strand not in {"+", "-"}:
                raise PipelineFailure("invalid_strand", gene_id)
            orders[sequence_id].append(order)
            per_sequence_counts[sequence_id] += 1
            gff_ids.append(gene_id)

    if len(gff_ids) != len(set(gff_ids)):
        raise PipelineFailure("duplicate_gff_id", gff_path.name)
    if set(gff_ids) != set(cds_ids) or set(gff_ids) != set(pep_ids):
        raise PipelineFailure(
            "bundle_id_mismatch",
            f"gff={len(gff_ids)} cds={len(cds_ids)} pep={len(pep_ids)}",
        )
    for sequence_id, sequence_orders in orders.items():
        if sequence_orders != list(range(1, len(sequence_orders) + 1)):
            raise PipelineFailure("non_contiguous_order", sequence_id)
    for sequence_id, (_length, expected_count) in lens.items():
        if per_sequence_counts[sequence_id] != expected_count:
            raise PipelineFailure("lens_gene_count_mismatch", sequence_id)
    return {
        "n_genes": len(gff_ids),
        "n_sequences": len(lens),
    }


def _add_rejection(
    rows: list[dict[str, str]],
    scope: str,
    identifier: str,
    reason: str,
    detail: str = "",
) -> None:
    rows.append(
        {
            "scope": scope,
            "id": identifier,
            "reason": reason,
            "detail": detail,
        }
    )


def _discover_accession(species: SpeciesInput) -> str:
    accessions = {
        accession
        for accession in (
            _accession_from_filename(species.gff_name),
            _accession_from_filename(species.cds_name),
            _accession_from_filename(species.pep_name),
        )
        if accession
    }
    if len(accessions) != 1:
        raise PipelineFailure(
            "accession_mismatch",
            f"expected one accession, found {sorted(accessions)}",
        )
    return next(iter(accessions))


def _prepare_staging(
    options: PipelineOptions,
    code: str,
) -> tuple[Path, Path, Path, Path]:
    clean_root = options.project_root / options.clean_dirname
    analysis_root = options.project_root / options.analysis_dirname
    clean_root.mkdir(parents=True, exist_ok=True)
    analysis_root.mkdir(parents=True, exist_ok=True)
    suffix = uuid.uuid4().hex
    clean_staging = clean_root / f".{code}.tmp-{suffix}"
    analysis_staging = analysis_root / f".{code}.tmp-{suffix}"
    clean_staging.mkdir()
    analysis_staging.mkdir()
    return clean_root, analysis_root, clean_staging, analysis_staging


def process_species(species: SpeciesInput, options: PipelineOptions) -> dict[str, object]:
    clean_root, analysis_root, clean_staging, analysis_staging = _prepare_staging(
        options, species.code
    )
    clean_final = clean_root / species.code
    analysis_final = analysis_root / species.code
    try:
        result = _process_species_into(
            species,
            options,
            clean_staging,
            analysis_staging,
        )
        _replace_directory(clean_staging, clean_final, options.overwrite)
        _replace_directory(analysis_staging, analysis_final, options.overwrite)
        return result
    except PipelineFailure as exc:
        _remove_staging(clean_staging, clean_root)
        _remove_staging(analysis_staging, analysis_root)
        return {
            "show_id": species.show_id,
            "code": species.code,
            "species": species.species,
            "accession": species.accession,
            "status": "FAIL",
            "reason": exc.reason,
            "detail": exc.detail,
        }
    except Exception:
        _remove_staging(clean_staging, clean_root)
        _remove_staging(analysis_staging, analysis_root)
        raise


def _process_species_into(
    species: SpeciesInput,
    options: PipelineOptions,
    clean_dir: Path,
    analysis_dir: Path,
) -> dict[str, object]:
    accession = _discover_accession(species)
    gff_path = species.directory / species.gff_name
    cds_path = species.directory / species.cds_name
    pep_path = species.directory / species.pep_name
    for path in (gff_path, cds_path, pep_path):
        if not path.is_file() or path.stat().st_size == 0:
            raise PipelineFailure("missing_or_empty_input", str(path))

    parsed = parse_gff(gff_path)
    fatal_gff = [item for item in parsed.anomalies if item.level == "fatal"]
    if fatal_gff:
        raise PipelineFailure(fatal_gff[0].code, fatal_gff[0].detail)

    nuclear_taxa = {
        region.taxon_id
        for region in parsed.regions.values()
        if region.taxon_id and not region.is_organelle
    }
    if len(nuclear_taxa) > 1:
        raise PipelineFailure(
            "multiple_nuclear_taxa",
            f"{species.species}: taxon IDs {sorted(nuclear_taxa)}",
        )

    cds = read_cds_fasta(cds_path)
    pep = read_pep_fasta(pep_path)
    fatal_fasta = [
        item for item in [*cds.anomalies, *pep.anomalies] if item.level == "fatal"
    ]
    if fatal_fasta:
        raise PipelineFailure(fatal_fasta[0].code, fatal_fasta[0].detail)

    rejected: list[dict[str, str]] = []
    candidates: list[Candidate] = []
    reason_counts: Counter[str] = Counter()
    for transcript_id, transcript in parsed.transcripts.items():
        gene = parsed.genes.get(transcript.gene_id)
        if gene is None:
            _add_rejection(
                rejected,
                "transcript",
                transcript_id,
                "gene_missing",
                transcript.gene_id,
            )
            reason_counts["gene_missing"] += 1
            continue
        region = parsed.regions.get(transcript.seqid)
        if region is None:
            _add_rejection(
                rejected,
                "transcript",
                transcript_id,
                "region_missing",
                transcript.seqid,
            )
            reason_counts["region_missing"] += 1
            continue
        if region.is_organelle:
            _add_rejection(
                rejected,
                "transcript",
                transcript_id,
                "organelle_sequence",
                f"{transcript.seqid}:{region.genome_type}",
            )
            reason_counts["organelle_sequence"] += 1
            continue
        if transcript.strand not in {"+", "-"}:
            _add_rejection(
                rejected,
                "transcript",
                transcript_id,
                "invalid_strand",
                transcript.strand,
            )
            reason_counts["invalid_strand"] += 1
            continue
        if len(transcript.protein_ids) != 1:
            reason = "protein_id_missing" if not transcript.protein_ids else "multiple_protein_ids"
            _add_rejection(
                rejected,
                "transcript",
                transcript_id,
                reason,
                ",".join(sorted(transcript.protein_ids)),
            )
            reason_counts[reason] += 1
            continue
        protein_id = next(iter(transcript.protein_ids))
        cds_sequence = cds.sequences.get(protein_id)
        pep_sequence = pep.sequences.get(protein_id)
        if cds_sequence is None or pep_sequence is None:
            _add_rejection(
                rejected,
                "transcript",
                transcript_id,
                "sequence_missing",
                f"cds={protein_id in cds.sequences} pep={protein_id in pep.sequences}",
            )
            reason_counts["sequence_missing"] += 1
            continue
        check = validate_cds_pep(cds_sequence, pep_sequence)
        if not check.valid:
            _add_rejection(
                rejected,
                "transcript",
                transcript_id,
                check.reason,
                protein_id,
            )
            reason_counts[check.reason] += 1
            continue
        if (
            transcript.start < 1
            or transcript.end < transcript.start
            or transcript.end > region.length
        ):
            _add_rejection(
                rejected,
                "transcript",
                transcript_id,
                "coordinate_out_of_region",
                f"{transcript.start}-{transcript.end}/{region.length}",
            )
            reason_counts["coordinate_out_of_region"] += 1
            continue
        candidates.append(
            Candidate(
                source_gene_id=gene.source_id,
                source_gene_name=gene.name,
                gene_biotype=gene.biotype,
                source_transcript_id=transcript_id,
                source_protein_id=protein_id,
                seqid=transcript.seqid,
                start=transcript.start,
                end=transcript.end,
                strand=transcript.strand,
                cds_length=len(cds_sequence),
                protein_length=len(pep_sequence.rstrip("*")),
                translation_exact=check.exact,
                translation_mismatches=check.mismatches,
            )
        )

    protein_to_genes: dict[str, set[str]] = defaultdict(set)
    for candidate in candidates:
        protein_to_genes[candidate.source_protein_id].add(candidate.source_gene_id)
    ambiguous_proteins = {
        protein_id for protein_id, gene_ids in protein_to_genes.items() if len(gene_ids) > 1
    }
    if ambiguous_proteins:
        retained: list[Candidate] = []
        for candidate in candidates:
            if candidate.source_protein_id in ambiguous_proteins:
                _add_rejection(
                    rejected,
                    "transcript",
                    candidate.source_transcript_id,
                    "protein_id_reused_across_genes",
                    candidate.source_protein_id,
                )
                reason_counts["protein_id_reused_across_genes"] += 1
            else:
                retained.append(candidate)
        candidates = retained

    by_gene: dict[str, list[Candidate]] = defaultdict(list)
    for candidate in candidates:
        by_gene[candidate.source_gene_id].append(candidate)
    representatives: list[Candidate] = []
    for gene_id in sorted(by_gene):
        group = by_gene[gene_id]
        representative = choose_representative(group)
        representatives.append(representative)
        for candidate in group:
            if candidate is not representative:
                _add_rejection(
                    rejected,
                    "transcript",
                    candidate.source_transcript_id,
                    "alternative_isoform",
                    f"representative={representative.source_transcript_id}",
                )
                reason_counts["alternative_isoform"] += 1

    if not representatives:
        raise PipelineFailure("no_valid_representative_genes")

    gene_counts: Counter[str] = Counter(candidate.seqid for candidate in representatives)
    sequence_selection = select_analysis_sequences(
        parsed.regions,
        gene_counts,
        options.min_genes_per_sequence,
        options.target_gene_coverage,
        options.max_sequences,
    )
    if not sequence_selection.keep:
        raise PipelineFailure("empty_wgdi_sequence_set")

    nuclear_with_genes = [
        parsed.regions[seqid]
        for seqid, count in gene_counts.items()
        if count > 0 and seqid in parsed.regions and not parsed.regions[seqid].is_organelle
    ]
    nuclear_with_genes.sort(key=region_sort_key)
    clean_sequence_id = {
        region.seqid: index for index, region in enumerate(nuclear_with_genes, 1)
    }

    by_sequence: dict[str, list[Candidate]] = defaultdict(list)
    for candidate in representatives:
        by_sequence[candidate.seqid].append(candidate)

    entries: list[dict[str, object]] = []
    for region in nuclear_with_genes:
        rows = sorted(
            by_sequence[region.seqid],
            key=lambda candidate: (
                candidate.start,
                candidate.end,
                candidate.source_gene_id,
            ),
        )
        for order, candidate in enumerate(rows, 1):
            canonical_id = (
                f"{species.code}{clean_sequence_id[candidate.seqid]}"
                f"g{order:04d}"
            )
            entries.append(
                {
                    **asdict(candidate),
                    "canonical_gene_id": canonical_id,
                    "clean_sequence_id": clean_sequence_id[candidate.seqid],
                    "order": order,
                    "selected_for_wgdi": candidate.seqid in sequence_selection.keep,
                    "cds_sequence": cds.sequences[candidate.source_protein_id],
                    "pep_sequence": pep.sequences[candidate.source_protein_id].rstrip("*"),
                }
            )

    selected_entries = [
        entry for entry in entries if bool(entry["selected_for_wgdi"])
    ]
    if not selected_entries:
        raise PipelineFailure("empty_wgdi_gene_set")

    all_region_rows: list[dict[str, object]] = []
    selected_region_rows: list[dict[str, object]] = []
    for region in nuclear_with_genes:
        row = {
            "clean_sequence_id": clean_sequence_id[region.seqid],
            "source_seqid": region.seqid,
            "length": region.length,
            "gene_count": gene_counts[region.seqid],
            "selected": region.seqid in sequence_selection.keep,
        }
        all_region_rows.append(row)
        if row["selected"]:
            selected_region_rows.append(row)

    all_cds_records = [
        {
            "canonical_gene_id": entry["canonical_gene_id"],
            "source_gene_id": entry["source_gene_id"],
            "source_protein_id": entry["source_protein_id"],
            "sequence": entry["cds_sequence"],
        }
        for entry in entries
    ]
    all_pep_records = [
        {
            "canonical_gene_id": entry["canonical_gene_id"],
            "source_gene_id": entry["source_gene_id"],
            "source_protein_id": entry["source_protein_id"],
            "sequence": entry["pep_sequence"],
        }
        for entry in entries
    ]
    selected_cds_records = [
        record
        for entry, record in zip(entries, all_cds_records)
        if bool(entry["selected_for_wgdi"])
    ]
    selected_pep_records = [
        record
        for entry, record in zip(entries, all_pep_records)
        if bool(entry["selected_for_wgdi"])
    ]

    code = species.code
    all_cds_path = analysis_dir / f"{code}.all.cds.fa"
    all_pep_path = analysis_dir / f"{code}.all.pep.fa"
    all_gff_path = analysis_dir / f"{code}.all.gff"
    all_lens_path = analysis_dir / f"{code}.all.lens"
    cds_out = analysis_dir / f"{code}.cds"
    pep_out = analysis_dir / f"{code}.pep"
    gff_out = analysis_dir / f"{code}.new.gff"
    lens_out = analysis_dir / f"{code}.lens"
    _write_fasta(all_cds_path, all_cds_records)
    _write_fasta(all_pep_path, all_pep_records)
    _write_gff(all_gff_path, entries)
    _write_lens(all_lens_path, all_region_rows)
    _write_fasta(cds_out, selected_cds_records)
    _write_fasta(pep_out, selected_pep_records)
    _write_gff(gff_out, selected_entries)
    _write_lens(lens_out, selected_region_rows)

    all_validation = _validate_bundle(
        all_gff_path, all_lens_path, all_cds_path, all_pep_path
    )
    selected_validation = _validate_bundle(gff_out, lens_out, cds_out, pep_out)

    gene_map_path = clean_dir / f"{code}.gene_map.tsv"
    gene_map_fields = [
        "canonical_gene_id",
        "clean_sequence_id",
        "canonical_order",
        "source_gene_id",
        "source_gene_name",
        "source_transcript_id",
        "source_protein_id",
        "source_seqid",
        "start",
        "end",
        "strand",
        "gene_biotype",
        "cds_length",
        "protein_length",
        "translation_exact",
        "translation_mismatches",
        "selected_for_wgdi",
    ]
    with gene_map_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=gene_map_fields, delimiter="\t")
        writer.writeheader()
        for entry in entries:
            writer.writerow(
                {
                    "canonical_gene_id": entry["canonical_gene_id"],
                    "clean_sequence_id": entry["clean_sequence_id"],
                    "canonical_order": entry["order"],
                    "source_gene_id": entry["source_gene_id"],
                    "source_gene_name": entry["source_gene_name"],
                    "source_transcript_id": entry["source_transcript_id"],
                    "source_protein_id": entry["source_protein_id"],
                    "source_seqid": entry["seqid"],
                    "start": entry["start"],
                    "end": entry["end"],
                    "strand": entry["strand"],
                    "gene_biotype": entry["gene_biotype"],
                    "cds_length": entry["cds_length"],
                    "protein_length": entry["protein_length"],
                    "translation_exact": "Y" if entry["translation_exact"] else "N",
                    "translation_mismatches": entry["translation_mismatches"],
                    "selected_for_wgdi": "Y" if entry["selected_for_wgdi"] else "N",
                }
            )

    chr_map_path = clean_dir / f"{code}.chr_map.tsv"
    with chr_map_path.open("w", encoding="utf-8", newline="") as handle:
        fields = [
            "source_seqid",
            "clean_sequence_id",
            "length",
            "genome_type",
            "taxon_id",
            "is_primary",
            "is_organelle",
            "gene_count",
            "in_wgdi_set",
            "exclusion_reason",
        ]
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for region in sorted(
            parsed.regions.values(),
            key=lambda item: (
                item.is_organelle,
                *region_sort_key(item),
            ),
        ):
            count = gene_counts.get(region.seqid, 0)
            clean_id = clean_sequence_id.get(region.seqid, "")
            reason = ""
            if region.is_organelle:
                reason = "organelle_sequence"
            elif count == 0:
                reason = "no_valid_protein_gene"
            elif region.seqid not in sequence_selection.keep:
                reason = "outside_wgdi_sequence_set"
            writer.writerow(
                {
                    "source_seqid": region.seqid,
                    "clean_sequence_id": clean_id,
                    "length": region.length,
                    "genome_type": region.genome_type,
                    "taxon_id": region.taxon_id,
                    "is_primary": "Y" if region.is_primary else "N",
                    "is_organelle": "Y" if region.is_organelle else "N",
                    "gene_count": count,
                    "in_wgdi_set": "Y" if region.seqid in sequence_selection.keep else "N",
                    "exclusion_reason": reason,
                }
            )

    rejected_path = clean_dir / f"{code}.rejected.tsv"
    with rejected_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=["scope", "id", "reason", "detail"], delimiter="\t"
        )
        writer.writeheader()
        writer.writerows(rejected)

    input_info = {}
    for label, path in (("gff", gff_path), ("cds", cds_path), ("pep", pep_path)):
        input_info[label] = {
            "path": str(path),
            "bytes": path.stat().st_size,
            "sha256": _sha256(path),
        }
    qc = {
        "show_id": species.show_id,
        "code": species.code,
        "species": species.species,
        "order": species.order,
        "family": species.family,
        "accession": accession,
        "assembly_build": parsed.assembly_build,
        "annotation_source": parsed.annotation_source,
        "pipeline_version": PIPELINE_VERSION,
        "inputs": input_info,
        "n_regions": len(parsed.regions),
        "n_genes": len(parsed.genes),
        "n_transcripts": len(parsed.transcripts),
        "n_candidates": len(candidates),
        "n_representative_genes": len(representatives),
        "n_organelle_genes": len(
            {
                transcript.gene_id
                for transcript in parsed.transcripts.values()
                if transcript.seqid in parsed.regions
                and parsed.regions[transcript.seqid].is_organelle
            }
        ),
        "n_cds_records": cds.counts["records"],
        "n_pep_records": pep.counts["records"],
        "n_cds_duplicate_conflicts": cds.counts["duplicate_conflict"],
        "n_pep_duplicate_conflicts": pep.counts["duplicate_conflict"],
        "rejection_counts": dict(sorted(reason_counts.items())),
        "nuclear_taxa": sorted(nuclear_taxa),
        "sequence_selection": {
            "mode": sequence_selection.mode,
            "selected_genes": sequence_selection.selected_genes,
            "total_genes": sequence_selection.total_genes,
            "ratio": sequence_selection.ratio,
            "wgdi_ready": sequence_selection.ready,
            "min_genes_per_sequence": options.min_genes_per_sequence,
            "target_gene_coverage": options.target_gene_coverage,
            "max_sequences": options.max_sequences,
        },
        "all_bundle": all_validation,
        "wgdi_bundle": selected_validation,
        "gff_feature_counts": dict(sorted(parsed.feature_counts.items())),
        "anomalies": [asdict(item) for item in [*parsed.anomalies, *cds.anomalies, *pep.anomalies]],
        "status": "OK",
    }
    with (analysis_dir / f"{code}.qc.json").open("w", encoding="utf-8") as handle:
        json.dump(qc, handle, ensure_ascii=False, indent=2)
    shutil.copy2(gene_map_path, analysis_dir / f"{code}.gene_map.tsv")
    shutil.copy2(chr_map_path, analysis_dir / f"{code}.chr_map.tsv")

    return {
        "show_id": species.show_id,
        "code": species.code,
        "species": species.species,
        "accession": accession,
        "status": "OK",
        "assembly_build": parsed.assembly_build,
        "n_regions": len(parsed.regions),
        "n_representative_genes": len(representatives),
        "n_all_genes": all_validation["n_genes"],
        "n_wgdi_genes": selected_validation["n_genes"],
        "n_wgdi_sequences": selected_validation["n_sequences"],
        "sequence_mode": sequence_selection.mode,
        "retain_ratio": sequence_selection.ratio,
        "wgdi_ready": sequence_selection.ready,
        "rejection_counts": dict(sorted(reason_counts.items())),
        "nuclear_taxa": ",".join(sorted(nuclear_taxa)),
        "n_anomalies": len(parsed.anomalies) + len(cds.anomalies) + len(pep.anomalies),
    }


def _write_tsv(
    path: Path,
    rows: Sequence[Mapping[str, Any]],
    fields: list[str],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def _read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        return []
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def run_project(
    project_root: Path | str,
    *,
    show_ids: set[str] | None = None,
    options: PipelineOptions | None = None,
) -> list[dict[str, object]]:
    project_root = Path(project_root)
    options = options or PipelineOptions(project_root=project_root)
    all_species = discover_species(project_root)
    species = all_species
    if show_ids is not None:
        species = [item for item in species if item.show_id in show_ids]
    selected_codes = {item.code for item in species}

    results: list[dict[str, object]] = []
    for item in species:
        try:
            result = process_species(item, options)
        except Exception as exc:
            result = {
                "show_id": item.show_id,
                "code": item.code,
                "species": item.species,
                "accession": item.accession,
                "status": "ERROR",
                "reason": type(exc).__name__,
                "detail": str(exc),
            }
        results.append(result)
        print(
            f"[{result['status']}] {item.show_id} {item.species} "
            f"genes={result.get('n_wgdi_genes', '')} "
            f"seqs={result.get('n_wgdi_sequences', '')} "
            f"ready={result.get('wgdi_ready', '')} "
            f"reason={result.get('reason', '')}",
            flush=True,
        )

    resources = project_root / "resources"
    registry_rows = [
        {
            "show_id": item.show_id,
            "code": item.code,
            "species": item.species,
            "order": item.order,
            "family": item.family,
            "accession": item.accession,
            "gff": item.gff_name,
            "cds": item.cds_name,
            "pep": item.pep_name,
            "directory": item.directory.name,
        }
        for item in all_species
    ]
    _write_tsv(
        resources / "species_registry_v2.tsv",
        registry_rows,
        [
            "show_id",
            "code",
            "species",
            "order",
            "family",
            "accession",
            "gff",
            "cds",
            "pep",
            "directory",
        ],
    )
    summary_fields = [
        "show_id",
        "code",
        "species",
        "accession",
        "status",
        "assembly_build",
        "n_regions",
        "n_representative_genes",
        "n_all_genes",
        "n_wgdi_genes",
        "n_wgdi_sequences",
        "sequence_mode",
        "retain_ratio",
        "wgdi_ready",
        "nuclear_taxa",
        "n_anomalies",
        "rejection_counts",
        "reason",
        "detail",
    ]
    summary_path = resources / "stage_preprocess_v2_summary.tsv"
    existing_summary = _read_tsv(summary_path)
    merged_summary: list[dict[str, Any]] = (
        [row for row in existing_summary if row.get("code") not in selected_codes]
        if show_ids is not None
        else []
    )
    merged_summary.extend(results)
    merged_summary.sort(key=lambda row: natural_key(str(row.get("code", ""))))
    _write_tsv(summary_path, merged_summary, summary_fields)

    index_rows: list[dict[str, object]] = []
    for result in results:
        code = str(result.get("code", ""))
        if result.get("status") == "OK":
            base = f"{options.analysis_dirname}/{code}"
            index_rows.append(
                {
                    **result,
                    "pep": f"{base}/{code}.pep",
                    "cds": f"{base}/{code}.cds",
                    "gff": f"{base}/{code}.new.gff",
                    "lens": f"{base}/{code}.lens",
                    "all_pep": f"{base}/{code}.all.pep.fa",
                    "all_cds": f"{base}/{code}.all.cds.fa",
                    "all_gff": f"{base}/{code}.all.gff",
                    "all_lens": f"{base}/{code}.all.lens",
                    "gene_map": f"{options.clean_dirname}/{code}/{code}.gene_map.tsv",
                    "chr_map": f"{options.clean_dirname}/{code}/{code}.chr_map.tsv",
                    "rejected": f"{options.clean_dirname}/{code}/{code}.rejected.tsv",
                    "qc": f"{base}/{code}.qc.json",
                }
            )
    index_fields = [
        "code",
        "species",
        "accession",
        "wgdi_ready",
        "pep",
        "cds",
        "gff",
        "lens",
        "all_pep",
        "all_cds",
        "all_gff",
        "all_lens",
        "gene_map",
        "chr_map",
        "rejected",
        "qc",
    ]
    index_path = resources / "wgdi_v2_index.tsv"
    existing_index = _read_tsv(index_path)
    merged_index: list[dict[str, Any]] = (
        [row for row in existing_index if row.get("code") not in selected_codes]
        if show_ids is not None
        else []
    )
    merged_index.extend(index_rows)
    merged_index.sort(key=lambda row: natural_key(str(row.get("code", ""))))
    _write_tsv(index_path, merged_index, index_fields)

    manifest_rows: list[dict[str, object]] = []
    successful = {str(row.get("code", "")) for row in merged_index if row.get("code")}
    for layer in (options.clean_dirname, options.analysis_dirname):
        layer_root = project_root / layer
        for code in sorted(successful, key=natural_key):
            code_dir = layer_root / code
            if not code_dir.is_dir():
                continue
            for path in sorted(code_dir.rglob("*")):
                if not path.is_file():
                    continue
                manifest_rows.append(
                    {
                        "code": code,
                        "layer": layer,
                        "relative_path": path.relative_to(project_root).as_posix(),
                        "bytes": path.stat().st_size,
                        "sha256": _sha256(path),
                    }
                )
    _write_tsv(
        resources / "wgdi_v2_manifest_sha256.tsv",
        manifest_rows,
        ["code", "layer", "relative_path", "bytes", "sha256"],
    )
    return results

