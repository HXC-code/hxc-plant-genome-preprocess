from __future__ import annotations

import re
from collections import Counter

from pgci_preprocess.models import Candidate, Region, SequenceSelection

DIGIT_RE = re.compile(r"(\d+)")


def natural_key(value: str) -> tuple[object, ...]:
    return tuple(int(part) if part.isdigit() else part.casefold() for part in DIGIT_RE.split(value))


def region_sort_key(region: Region) -> tuple[object, ...]:
    chromosome = region.chromosome
    if chromosome.casefold() not in {"", ".", "unknown", "na", "n/a"}:
        return (0, *natural_key(chromosome), *natural_key(region.seqid))
    return (1, *natural_key(region.seqid))


def choose_representative(candidates: list[Candidate]) -> Candidate:
    if not candidates:
        raise ValueError("candidates must not be empty")
    return min(
        candidates,
        key=lambda candidate: (
            -candidate.protein_length,
            -candidate.cds_length,
            candidate.source_transcript_id,
            candidate.source_protein_id,
        ),
    )


def select_analysis_sequences(
    regions: dict[str, Region],
    gene_counts: Counter[str],
    min_genes: int,
    target_coverage: float,
    max_sequences: int,
) -> SequenceSelection:
    total = sum(
        count
        for seqid, count in gene_counts.items()
        if seqid in regions and not regions[seqid].is_organelle
    )
    if total == 0:
        return SequenceSelection(set(), "empty", 0, 0, 0.0, "no")

    primary = {
        seqid
        for seqid, count in gene_counts.items()
        if count > 0 and seqid in regions and regions[seqid].is_primary and not regions[seqid].is_organelle
    }
    if primary:
        selected = sum(gene_counts[seqid] for seqid in primary)
        ratio = selected / total
        ready = "yes" if ratio >= 0.7 and len(primary) <= max_sequences else "marginal"
        if ratio < 0.5:
            ready = "no"
        return SequenceSelection(primary, "primary_chromosome", selected, total, round(ratio, 4), ready)

    eligible = {
        seqid: count
        for seqid, count in gene_counts.items()
        if count > 0 and seqid in regions and not regions[seqid].is_organelle
    }
    threshold_keep = {
        seqid for seqid, count in eligible.items() if count >= min_genes
    }
    threshold_genes = sum(eligible[seqid] for seqid in threshold_keep)
    if (
        threshold_keep
        and len(threshold_keep) <= max_sequences
        and threshold_genes / total >= 0.5
    ):
        ratio = threshold_genes / total
        ready = "yes" if ratio >= 0.7 else "marginal"
        return SequenceSelection(
            threshold_keep,
            "threshold",
            threshold_genes,
            total,
            round(ratio, 4),
            ready,
        )

    keep: list[str] = []
    selected = 0
    for seqid, count in sorted(
        eligible.items(),
        key=lambda item: (-item[1], *natural_key(item[0])),
    ):
        if len(keep) >= max_sequences:
            break
        keep.append(seqid)
        selected += count
        if selected / total >= target_coverage:
            break
    ratio = selected / total
    if ratio >= 0.7 and len(keep) <= max_sequences:
        ready = "yes"
    elif ratio >= 0.5:
        ready = "marginal"
    else:
        ready = "no"
    return SequenceSelection(
        set(keep),
        "coverage",
        selected,
        total,
        round(ratio, 4),
        ready,
    )

