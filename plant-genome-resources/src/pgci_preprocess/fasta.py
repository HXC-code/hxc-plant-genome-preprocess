from __future__ import annotations

import gzip
import re
from collections import Counter
from pathlib import Path
from typing import TextIO

from pgci_preprocess.models import Anomaly, FastaData

PROTEIN_ID_RE = re.compile(r"\[protein_id=([^\]]+)\]")
CDS_FALLBACK_RE = re.compile(r"(?:^|_)cds_([^_\\s]+?)(?:_\\d+)?(?:\\s|$)", re.IGNORECASE)


def open_text(path: Path | str) -> TextIO:
    path = Path(path)
    if path.suffix.casefold() == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace", newline="")
    return path.open("rt", encoding="utf-8", errors="replace", newline="")


def _read_fasta(path: Path | str, header_id) -> FastaData:
    path = Path(path)
    sequences: dict[str, str] = {}
    counts: Counter[str] = Counter()
    anomalies: list[Anomaly] = []
    current_id: str | None = None
    buffer: list[str] = []

    def flush() -> None:
        nonlocal current_id, buffer
        if current_id is None:
            return
        sequence = "".join(buffer).replace(" ", "").upper()
        counts["records"] += 1
        if not sequence:
            counts["empty_sequence"] += 1
            anomalies.append(Anomaly("warning", "empty_fasta_sequence", current_id))
        if current_id in sequences:
            if sequences[current_id] == sequence:
                counts["duplicate_identical"] += 1
                anomalies.append(Anomaly("warning", "duplicate_fasta_id_identical", current_id))
            else:
                counts["duplicate_conflict"] += 1
                anomalies.append(Anomaly("fatal", "duplicate_fasta_id_conflict", current_id))
        else:
            sequences[current_id] = sequence
        current_id = None
        buffer = []

    with open_text(path) as handle:
        for line_number, raw_line in enumerate(handle, 1):
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                flush()
                header = line[1:].strip()
                current_id = header_id(header)
                if not current_id:
                    counts["invalid_header"] += 1
                    anomalies.append(
                        Anomaly("warning", "fasta_header_without_id", f"line={line_number}")
                    )
                continue
            if current_id is None:
                counts["sequence_before_header"] += 1
                anomalies.append(
                    Anomaly("warning", "sequence_before_header", f"line={line_number}")
                )
                continue
            buffer.append("".join(line.split()))
    flush()
    return FastaData(sequences=sequences, counts=counts, anomalies=anomalies)


def extract_cds_id(header: str) -> str:
    match = PROTEIN_ID_RE.search(header)
    if match:
        return match.group(1).strip()
    match = CDS_FALLBACK_RE.search(header)
    if match:
        return match.group(1).strip()
    return header.split(None, 1)[0].strip()


def read_cds_fasta(path: Path | str) -> FastaData:
    return _read_fasta(path, extract_cds_id)


def read_pep_fasta(path: Path | str) -> FastaData:
    return _read_fasta(path, lambda header: header.split(None, 1)[0].strip())

