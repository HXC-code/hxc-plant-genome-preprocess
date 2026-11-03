from __future__ import annotations

import csv
import hashlib
import json
import os
import re
import shutil
import uuid
from collections import Counter
from dataclasses import dataclass
from pathlib import Path

INDEX_RELATIVE = Path("resources") / "wgdi_v2_index.tsv"
MANIFEST_RELATIVE = Path("resources") / "wgdi_v2_manifest_sha256.tsv"
BUNDLE_KEYS = ("pep", "cds", "gff", "lens")
SPECIES_DIR_RE = re.compile(r"^s\d+$")
DNA_ALPHABET = frozenset("ACGTRYSWKMBDHVN")
PROTEIN_ALPHABET = frozenset("ACDEFGHIKLMNPQRSTVWYBXZUOJ*")


@dataclass(slots=True)
class ValidationOptions:
    project_root: Path
    expected_count: int | None = None
    require_compat_links: bool = True
    verify_hashes: bool = True
    report_path: Path | None = None


@dataclass(slots=True)
class ValidationReport:
    ok: bool
    rows: list[dict[str, object]]
    report_path: Path


def resolve_layout(project_root: Path | str) -> tuple[Path, Path]:
    project_root = Path(project_root).resolve()
    if (project_root / "plant-genome-resources" / INDEX_RELATIVE).is_file():
        return project_root, project_root / "plant-genome-resources"
    if (project_root / INDEX_RELATIVE).is_file():
        workspace = project_root.parent if project_root.name == "plant-genome-resources" else project_root
        return workspace, project_root
    raise FileNotFoundError(f"wgdi_v2_index.tsv not found under {project_root}")


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _safe_project_path(project_root: Path, relative: str) -> Path:
    path = (project_root / relative).resolve()
    if path != project_root and project_root not in path.parents:
        raise ValueError(f"path escapes project root: {relative}")
    return path


def _read_index(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _read_hash_manifest(path: Path) -> dict[str, str]:
    if not path.is_file():
        return {}
    with path.open(encoding="utf-8", newline="") as handle:
        return {
            row["relative_path"]: row["sha256"]
            for row in csv.DictReader(handle, delimiter="\t")
            if row.get("relative_path") and row.get("sha256")
        }


def _check_fasta(path: Path, kind: str) -> tuple[set[str], list[str]]:
    ids: set[str] = set()
    problems: list[str] = []
    current_id: str | None = None
    has_sequence = False
    invalid_chars: set[str] = set()
    allowed = DNA_ALPHABET if kind == "cds" else PROTEIN_ALPHABET

    def finish_record() -> None:
        nonlocal current_id, has_sequence, invalid_chars
        if current_id is not None and not has_sequence:
            problems.append(f"empty_sequence_{kind}:{current_id}")
        if invalid_chars:
            problems.append(f"invalid_alphabet_{kind}:{current_id}:{''.join(sorted(invalid_chars))}")
        current_id = None
        has_sequence = False
        invalid_chars = set()

    with path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped:
                continue
            if stripped.startswith(">"):
                finish_record()
                current_id = stripped[1:].split(None, 1)[0]
                if not current_id:
                    problems.append(f"empty_header_{kind}")
                elif current_id in ids:
                    problems.append(f"duplicate_header_{kind}:{current_id}")
                else:
                    ids.add(current_id)
                continue
            if current_id is None:
                problems.append(f"sequence_before_header_{kind}")
                continue
            sequence = stripped.replace(" ", "").upper()
            if sequence:
                has_sequence = True
                invalid_chars.update(set(sequence) - allowed)
    finish_record()
    return ids, problems


def _check_species(
    row: dict[str, str],
    canonical_root: Path,
) -> tuple[dict[str, object], set[str]]:
    code = row.get("code", "").strip()
    problems: list[str] = []
    paths: dict[str, Path] = {}
    for key in BUNDLE_KEYS:
        relative = row.get(key, "").strip()
        if not relative:
            problems.append(f"missing_index_path_{key}")
            continue
        try:
            path = _safe_project_path(canonical_root, relative)
        except ValueError as exc:
            problems.append(f"unsafe_index_path_{key}:{exc}")
            continue
        paths[key] = path
        if not path.is_file() or path.stat().st_size == 0:
            problems.append(f"missing_or_empty_{key}")

    if problems:
        return (
            {
                "code": code,
                "status": "FAIL",
                "genes": "",
                "sequences": "",
                **{key: str(paths.get(key, "")) for key in BUNDLE_KEYS},
                "problems": ";".join(problems),
            },
            set(),
        )

    gff_ids: set[str] = set()
    sequence_rows: dict[str, list[tuple[int, int, int]]] = {}
    genes = 0
    with paths["gff"].open(encoding="utf-8", errors="replace") as handle:
        for line_number, line in enumerate(handle, 1):
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) != 7:
                problems.append(f"gff_column_count:{line_number}")
                continue
            sequence_id, gene_id, start_text, end_text, strand, order_text, _source_id = fields
            if not sequence_id.isdigit():
                problems.append(f"non_numeric_sequence_id:{sequence_id}")
            if not gene_id:
                problems.append(f"empty_gene_id:{line_number}")
            elif gene_id in gff_ids:
                problems.append(f"duplicate_gene_id:{gene_id}")
            else:
                gff_ids.add(gene_id)
            try:
                start = int(start_text)
                end = int(end_text)
                order = int(order_text)
            except ValueError:
                problems.append(f"non_integer_coordinate_or_order:{line_number}")
                continue
            if start < 1 or end < start:
                problems.append(f"invalid_coordinate:{gene_id}")
            if strand not in {"+", "-"}:
                problems.append(f"invalid_strand:{gene_id}")
            if order < 1:
                problems.append(f"invalid_order:{gene_id}")
            sequence_rows.setdefault(sequence_id, []).append((order, start, end))
            genes += 1

    lens: dict[str, tuple[int, int]] = {}
    with paths["lens"].open(encoding="utf-8", errors="replace") as handle:
        for line_number, line in enumerate(handle, 1):
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) != 3:
                problems.append(f"lens_column_count:{line_number}")
                continue
            sequence_id, length_text, count_text = fields
            if sequence_id in lens:
                problems.append(f"duplicate_lens_sequence:{sequence_id}")
                continue
            try:
                length = int(length_text)
            except ValueError:
                problems.append(f"invalid_lens_length:{sequence_id}")
                length = 0
            try:
                count = int(count_text)
            except ValueError:
                problems.append(f"invalid_lens_gene_count:{sequence_id}")
                count = 0
            if length == 0 and count == 0 and not length_text.isdigit() and not count_text.isdigit():
                continue
            if length <= 0:
                problems.append(f"invalid_lens_length:{sequence_id}")
            if count <= 0:
                problems.append(f"invalid_lens_gene_count:{sequence_id}")
            lens[sequence_id] = (length, count)

    if not gff_ids:
        problems.append("empty_gff")
    if not lens:
        problems.append("empty_lens")
    if set(sequence_rows) != set(lens):
        problems.append("sequence_set_mismatch")

    for sequence_id, rows in sequence_rows.items():
        rows.sort()
        expected_orders = list(range(1, len(rows) + 1))
        actual_orders = [order for order, _start, _end in rows]
        if actual_orders != expected_orders:
            problems.append(f"non_contiguous_order:{sequence_id}")
        previous_start = -1
        for _order, start, end in rows:
            if start < previous_start:
                problems.append(f"unsorted_gene_coordinates:{sequence_id}")
                break
            previous_start = start
            if sequence_id in lens and end > lens[sequence_id][0]:
                problems.append(f"coordinate_out_of_range:{sequence_id}")
                break
        if sequence_id in lens and lens[sequence_id][1] != len(rows):
            problems.append(f"lens_gene_count_mismatch:{sequence_id}")

    cds_ids, cds_problems = _check_fasta(paths["cds"], "cds")
    pep_ids, pep_problems = _check_fasta(paths["pep"], "pep")
    problems.extend(cds_problems)
    problems.extend(pep_problems)
    if gff_ids != cds_ids:
        problems.append("gff_cds_id_mismatch")
    if gff_ids != pep_ids:
        problems.append("gff_pep_id_mismatch")

    qc_relative = row.get("qc", "").strip()
    if qc_relative:
        try:
            qc_path = _safe_project_path(canonical_root, qc_relative)
            qc = json.loads(qc_path.read_text(encoding="utf-8"))
            if qc.get("status") != "OK":
                problems.append("qc_status_not_ok")
            if qc.get("code") != code:
                problems.append("qc_code_mismatch")
        except (OSError, ValueError, json.JSONDecodeError):
            problems.append("qc_unreadable")
    else:
        problems.append("qc_missing_from_index")

    return (
        {
            "code": code,
            "status": "PASS" if not problems else "FAIL",
            "genes": genes,
            "sequences": len(lens),
            **{key: str(paths[key]) for key in BUNDLE_KEYS},
            "problems": ";".join(problems),
        },
        gff_ids,
    )


def _atomic_write_report(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp-{uuid.uuid4().hex}")
    try:
        with temporary.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=[
                    "code",
                    "status",
                    "genes",
                    "sequences",
                    "pep",
                    "cds",
                    "gff",
                    "lens",
                    "problems",
                ],
                delimiter="\t",
            )
            writer.writeheader()
            writer.writerows(rows)
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def _add_species_problem(rows: list[dict[str, object]], code: str, problem: str) -> None:
    for row in reversed(rows):
        if row.get("code") == code and row.get("code") not in {"_layout", "_hashes", "_compat"}:
            existing = str(row.get("problems", ""))
            row["problems"] = f"{existing};{problem}".strip(";")
            row["status"] = "FAIL"
            return
    rows.append(
        {
            "code": code,
            "status": "FAIL",
            "genes": "",
            "sequences": "",
            "pep": "",
            "cds": "",
            "gff": "",
            "lens": "",
            "problems": problem,
        }
    )


def validate_project(options: ValidationOptions) -> ValidationReport:
    workspace, canonical_root = resolve_layout(options.project_root)
    index_path = canonical_root / INDEX_RELATIVE
    index_rows = _read_index(index_path)
    report_path = options.report_path or workspace / "data-validation-report.tsv"
    rows: list[dict[str, object]] = []

    code_counts = Counter(row.get("code", "") for row in index_rows)
    duplicate_codes = {code for code, count in code_counts.items() if count > 1}
    if not index_rows:
        rows.append(
            {
                "code": "_layout",
                "status": "FAIL",
                "genes": "",
                "sequences": "",
                "pep": "",
                "cds": "",
                "gff": "",
                "lens": "",
                "problems": "empty_manifest",
            }
        )
    expected_count = options.expected_count
    if expected_count is not None and len(index_rows) != expected_count:
        rows.append(
            {
                "code": "_layout",
                "status": "FAIL",
                "genes": "",
                "sequences": "",
                "pep": "",
                "cds": "",
                "gff": "",
                "lens": "",
                "problems": f"expected_count_mismatch:expected={expected_count},actual={len(index_rows)}",
            }
        )
    if duplicate_codes:
        rows.append(
            {
                "code": "_layout",
                "status": "FAIL",
                "genes": "",
                "sequences": "",
                "pep": "",
                "cds": "",
                "gff": "",
                "lens": "",
                "problems": "duplicate_manifest_codes:" + ",".join(sorted(duplicate_codes)),
            }
        )

    global_gene_ids: dict[str, str] = {}
    for row in index_rows:
        result, gene_ids = _check_species(row, canonical_root)
        code = str(result["code"])
        for gene_id in gene_ids:
            previous_code = global_gene_ids.get(gene_id)
            if previous_code and previous_code != code:
                result["status"] = "FAIL"
                result["problems"] = (
                    f"{result['problems']};duplicate_gene_id_across_species:{gene_id}"
                ).strip(";")
            else:
                global_gene_ids[gene_id] = code
        rows.append(result)

    if options.verify_hashes:
        hash_manifest = _read_hash_manifest(canonical_root / MANIFEST_RELATIVE)
        if not hash_manifest:
            rows.append(
                {
                    "code": "_hashes",
                    "status": "FAIL",
                    "genes": "",
                    "sequences": "",
                    "pep": "",
                    "cds": "",
                    "gff": "",
                    "lens": "",
                    "problems": "hash_manifest_missing_or_empty",
                }
            )
        else:
            for row in index_rows:
                code = row.get("code", "")
                problems = []
                for key in BUNDLE_KEYS:
                    relative = row.get(key, "")
                    expected_hash = hash_manifest.get(relative)
                    try:
                        path = _safe_project_path(canonical_root, relative)
                    except ValueError:
                        problems.append(f"unsafe_hash_path:{key}")
                        continue
                    if not expected_hash:
                        problems.append(f"hash_missing:{key}")
                    elif not path.is_file() or _sha256(path) != expected_hash:
                        problems.append(f"hash_mismatch:{key}")
                if problems:
                    for problem in problems:
                        _add_species_problem(rows, code, problem)

    if options.require_compat_links:
        expected_codes = {row.get("code", "") for row in index_rows}
        for row in index_rows:
            code = row.get("code", "")
            problems = []
            compatibility_dir = workspace / code
            for key in BUNDLE_KEYS:
                source = _safe_project_path(canonical_root, row.get(key, ""))
                destination = compatibility_dir / source.name
                if not destination.is_file():
                    problems.append(f"compat_missing_{key}")
                else:
                    try:
                        same = os.path.samefile(source, destination)
                    except OSError:
                        same = False
                    if not same and _sha256(source) != _sha256(destination):
                        problems.append(f"compat_content_mismatch_{key}")
            if problems:
                for problem in problems:
                    _add_species_problem(rows, code, problem)
        extra_dirs = sorted(
            path.name
            for path in workspace.iterdir()
            if path.is_dir()
            and SPECIES_DIR_RE.fullmatch(path.name)
            and path.name not in expected_codes
        )
        if extra_dirs:
            rows.append(
                {
                    "code": "_compat",
                    "status": "FAIL",
                    "genes": "",
                    "sequences": "",
                    "pep": "",
                    "cds": "",
                    "gff": "",
                    "lens": "",
                    "problems": "extra_compat_directories:" + ",".join(extra_dirs),
                }
            )

    ok = bool(rows) and all(row["status"] == "PASS" for row in rows)
    _atomic_write_report(report_path, rows)
    return ValidationReport(ok=ok, rows=rows, report_path=report_path)


def sync_compat_links(canonical_root: Path | str, workspace: Path | str) -> list[str]:
    canonical_root = Path(canonical_root).resolve()
    workspace = Path(workspace).resolve()
    index_rows = _read_index(canonical_root / INDEX_RELATIVE)
    messages: list[str] = []
    for row in index_rows:
        code = row.get("code", "").strip()
        if not code:
            continue
        destination_dir = workspace / code
        destination_dir.mkdir(parents=True, exist_ok=True)
        for key in BUNDLE_KEYS:
            source = _safe_project_path(canonical_root, row.get(key, ""))
            if not source.is_file():
                messages.append(f"missing source: {source}")
                continue
            destination = destination_dir / source.name
            try:
                if destination.exists() and os.path.samefile(source, destination):
                    continue
            except OSError:
                pass
            destination.unlink(missing_ok=True)
            try:
                os.link(source, destination)
                messages.append(f"linked: {destination}")
            except OSError:
                shutil.copy2(source, destination)
                messages.append(f"copied: {destination}")
    return messages

