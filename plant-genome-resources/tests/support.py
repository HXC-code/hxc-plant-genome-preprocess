from __future__ import annotations

import csv
import hashlib
import json
import os
from pathlib import Path


def build_validation_project(
    tmp_path: Path,
    *,
    species_count: int = 1,
    compatible_links: bool = True,
) -> dict[str, Path]:
    workspace = tmp_path / "workspace"
    project = workspace / "plant-genome-resources"
    resources = project / "resources"
    analysis = project / "analysis_v2"
    clean = project / "clean_v2"
    resources.mkdir(parents=True)
    analysis.mkdir(parents=True)
    clean.mkdir(parents=True)

    index_path = resources / "wgdi_v2_index.tsv"
    fields = [
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
    rows: list[dict[str, str]] = []
    for index in range(1, species_count + 1):
        code = f"s{index}"
        accession = f"GCF_00000000{index}.1"
        species_dir = analysis / code
        clean_dir = clean / code
        species_dir.mkdir()
        clean_dir.mkdir()

        gff = f"1\t{code}1g0001\t10\t100\t+\t1\tG{index}\n"
        lens = "1\t1000\t1\n"
        cds = f">{code}1g0001\nATGGCTTAA\n"
        pep = f">{code}1g0001\nMA\n"

        files = {
            "pep": species_dir / f"{code}.pep",
            "cds": species_dir / f"{code}.cds",
            "gff": species_dir / f"{code}.new.gff",
            "lens": species_dir / f"{code}.lens",
            "gene_map": clean_dir / f"{code}.gene_map.tsv",
            "chr_map": clean_dir / f"{code}.chr_map.tsv",
            "rejected": clean_dir / f"{code}.rejected.tsv",
            "qc": species_dir / f"{code}.qc.json",
        }
        files["pep"].write_text(pep, encoding="utf-8")
        files["cds"].write_text(cds, encoding="utf-8")
        files["gff"].write_text(gff, encoding="utf-8")
        files["lens"].write_text(lens, encoding="utf-8")
        files["gene_map"].write_text("canonical_gene_id\n", encoding="utf-8")
        files["chr_map"].write_text("source_seqid\n", encoding="utf-8")
        files["rejected"].write_text("scope\tid\treason\tdetail\n", encoding="utf-8")
        files["qc"].write_text(
            json.dumps({"code": code, "status": "OK"}, ensure_ascii=False),
            encoding="utf-8",
        )

        if compatible_links:
            compatibility_dir = workspace / code
            compatibility_dir.mkdir()
            for kind in ("pep", "cds", "gff", "lens"):
                source = files[kind]
                destination = compatibility_dir / source.name
                try:
                    os.link(source, destination)
                except OSError:
                    destination.write_bytes(source.read_bytes())

        rows.append(
            {
                "code": code,
                "species": f"Testus plantus {index}",
                "accession": accession,
                "wgdi_ready": "yes",
                "pep": files["pep"].relative_to(project).as_posix(),
                "cds": files["cds"].relative_to(project).as_posix(),
                "gff": files["gff"].relative_to(project).as_posix(),
                "lens": files["lens"].relative_to(project).as_posix(),
                "all_pep": files["pep"].relative_to(project).as_posix(),
                "all_cds": files["cds"].relative_to(project).as_posix(),
                "all_gff": files["gff"].relative_to(project).as_posix(),
                "all_lens": files["lens"].relative_to(project).as_posix(),
                "gene_map": files["gene_map"].relative_to(project).as_posix(),
                "chr_map": files["chr_map"].relative_to(project).as_posix(),
                "rejected": files["rejected"].relative_to(project).as_posix(),
                "qc": files["qc"].relative_to(project).as_posix(),
            }
        )

    with index_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    manifest_path = resources / "wgdi_v2_manifest_sha256.tsv"
    manifest_rows: list[dict[str, object]] = []
    for row in rows:
        for kind in ("pep", "cds", "gff", "lens"):
            path = project / row[kind]
            manifest_rows.append(
                {
                    "code": row["code"],
                    "layer": "analysis_v2",
                    "relative_path": row[kind],
                    "bytes": path.stat().st_size,
                    "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
                }
            )
    with manifest_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["code", "layer", "relative_path", "bytes", "sha256"],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(manifest_rows)
    return {
        "workspace": workspace,
        "project": project,
        "analysis": analysis,
        "clean": clean,
        "index": index_path,
        "manifest": manifest_path,
    }

