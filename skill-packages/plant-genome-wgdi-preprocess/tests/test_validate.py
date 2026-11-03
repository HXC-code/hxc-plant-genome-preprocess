from __future__ import annotations

from pathlib import Path

from pgci_preprocess.validate import ValidationOptions, sync_compat_links, validate_project
from tests.support import build_validation_project


def test_valid_project_passes_manifest_and_compatibility_checks(tmp_path: Path) -> None:
    paths = build_validation_project(tmp_path)

    result = validate_project(
        ValidationOptions(
            project_root=paths["workspace"],
            require_compat_links=True,
        )
    )

    assert result.ok
    assert result.rows[0]["status"] == "PASS"
    assert result.report_path.is_file()
    assert result.report_path.name == "data-validation-report.tsv"


def test_empty_manifest_fails_instead_of_returning_zero_species_pass(tmp_path: Path) -> None:
    paths = build_validation_project(tmp_path, species_count=0)

    result = validate_project(ValidationOptions(project_root=paths["workspace"]))

    assert not result.ok
    assert result.rows[0]["code"] == "_layout"
    assert "empty_manifest" in result.rows[0]["problems"]


def test_expected_species_count_mismatch_fails(tmp_path: Path) -> None:
    paths = build_validation_project(tmp_path, species_count=1)

    result = validate_project(
        ValidationOptions(project_root=paths["workspace"], expected_count=2)
    )

    assert not result.ok
    assert "expected_count_mismatch" in result.rows[0]["problems"]


def test_species_with_missing_pep_is_reported(tmp_path: Path) -> None:
    paths = build_validation_project(tmp_path)
    (paths["analysis"] / "s1" / "s1.pep").unlink()

    result = validate_project(
        ValidationOptions(project_root=paths["workspace"], require_compat_links=False)
    )

    assert not result.ok
    assert "missing_or_empty_pep" in result.rows[0]["problems"]


def test_invalid_strand_lens_length_and_coordinate_are_reported(tmp_path: Path) -> None:
    paths = build_validation_project(tmp_path)
    species = paths["analysis"] / "s1"
    (species / "s1.new.gff").write_text(
        "1\ts11g0001\t10\t2000\t.\t1\tG1\n", encoding="utf-8"
    )
    (species / "s1.lens").write_text("1\tabc\t1\n", encoding="utf-8")

    result = validate_project(
        ValidationOptions(project_root=paths["workspace"], require_compat_links=False)
    )

    problems = result.rows[0]["problems"]
    assert not result.ok
    assert "invalid_strand" in problems
    assert "invalid_lens_length" in problems
    assert "coordinate_out_of_range" in problems


def test_cross_species_duplicate_gene_id_is_reported(tmp_path: Path) -> None:
    paths = build_validation_project(tmp_path, species_count=2, compatible_links=False)
    for code in ("s1", "s2"):
        species = paths["analysis"] / code
        (species / f"{code}.new.gff").write_text(
            "1\tduplicate_gene\t10\t100\t+\t1\tG\n", encoding="utf-8"
        )
        (species / f"{code}.cds").write_text(">duplicate_gene\nATGGCTTAA\n", encoding="utf-8")
        (species / f"{code}.pep").write_text(">duplicate_gene\nMA\n", encoding="utf-8")

    result = validate_project(
        ValidationOptions(project_root=paths["workspace"], require_compat_links=False)
    )

    assert not result.ok
    assert any("duplicate_gene_id_across_species" in row["problems"] for row in result.rows)


def test_stale_compatibility_file_is_rejected_and_can_be_repaired(tmp_path: Path) -> None:
    paths = build_validation_project(tmp_path)
    compatibility_pep = paths["workspace"] / "s1" / "s1.pep"
    compatibility_pep.unlink()
    compatibility_pep.write_text(">s11g0001\nMAA\n", encoding="utf-8")

    failed = validate_project(ValidationOptions(project_root=paths["workspace"]))
    assert not failed.ok
    assert "compat_content_mismatch_pep" in failed.rows[0]["problems"]

    sync_compat_links(paths["project"], paths["workspace"])
    repaired = validate_project(ValidationOptions(project_root=paths["workspace"]))
    assert repaired.ok

