from __future__ import annotations

from pathlib import Path

from pgci_preprocess.gff import parse_attributes, parse_gff


def test_parse_attributes_decodes_gff3_values() -> None:
    attrs = parse_attributes("ID=gene-A;Name=alpha%2Cbeta;Dbxref=taxon:123")
    assert attrs["ID"] == "gene-A"
    assert attrs["Name"] == "alpha,beta"
    assert attrs["Dbxref"] == "taxon:123"


def test_parse_gff_links_gene_mrna_cds_and_regions(tmp_path: Path) -> None:
    gff = tmp_path / "sample.gff"
    gff.write_text(
        "\n".join(
            [
                "##sequence-region chr1 1 1000",
                (
                    "chr1\tsrc\tregion\t1\t1000\t.\t+\t.\t"
                    "ID=chr1:1..1000;Dbxref=taxon:123;chromosome=1;genome=chromosome"
                ),
                (
                    "chr1\tsrc\tgene\t10\t200\t.\t+\t.\t"
                    "ID=gene-A;Name=A;gbkey=Gene;gene_biotype=protein_coding"
                ),
                "chr1\tsrc\tmRNA\t10\t200\t.\t+\t.\tID=rna-A1;Parent=gene-A;transcript_id=A1",
                "chr1\tsrc\tCDS\t10\t90\t.\t+\t0\tID=cds-P1;Parent=rna-A1;protein_id=P1;gene=A",
                "chr1\tsrc\texon\t10\t90\t.\t+\t.\tID=exon-A1;Parent=rna-A1",
                (
                    "chr1\tsrc\tgene\t300\t500\t.\t-\t.\t"
                    "ID=gene-B;Name=B;gbkey=Gene;gene_biotype=protein_coding"
                ),
                "chr1\tsrc\tmRNA\t300\t500\t.\t-\t.\tID=rna-B1;Parent=gene-B;transcript_id=B1",
                "chr1\tsrc\tmRNA\t300\t500\t.\t-\t.\tID=rna-B2;Parent=gene-B;transcript_id=B2",
                "chr1\tsrc\tCDS\t300\t380\t.\t-\t0\tID=cds-P2;Parent=rna-B1;protein_id=P2;gene=B",
                "chr1\tsrc\tCDS\t400\t470\t.\t-\t0\tID=cds-P3;Parent=rna-B2;protein_id=P3;gene=B",
                (
                    "chr1\tsrc\tgene\t600\t700\t.\t.\t.\t"
                    "ID=gene-NC;Name=NC;gbkey=Gene;gene_biotype=lncRNA"
                ),
                "chr1\tsrc\tlnc_RNA\t600\t700\t.\t.\t.\tID=rna-NC;Parent=gene-NC",
                (
                    "chr2\tsrc\tregion\t1\t500\t.\t+\t.\t"
                    "ID=chr2:1..500;Dbxref=taxon:123;chromosome=2;genome=chromosome"
                ),
                (
                    "organelle\tsrc\tregion\t1\t200\t.\t+\t.\t"
                    "ID=organelle:1..200;Dbxref=taxon:999;genome=mitochondrion"
                ),
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    parsed = parse_gff(gff)

    assert parsed.regions["chr1"].length == 1000
    assert parsed.regions["chr1"].is_primary
    assert parsed.regions["chr1"].taxon_id == "123"
    assert parsed.regions["organelle"].genome_type == "mitochondrion"
    assert not parsed.regions["organelle"].is_primary
    assert parsed.genes["A"].name == "A"
    assert parsed.transcripts["A1"].gene_id == "A"
    assert parsed.transcripts["A1"].protein_ids == {"P1"}
    assert parsed.transcripts["B1"].cds_count == 1
    assert parsed.feature_counts["region"] == 3
    assert parsed.feature_counts["CDS"] == 3


def test_parse_gff_reports_conflicting_region_length(tmp_path: Path) -> None:
    gff = tmp_path / "duplicate_region.gff"
    gff.write_text(
        (
            "chr1\tsrc\tregion\t1\t100\t.\t+\t.\tID=r1;Dbxref=taxon:1;genome=chromosome\n"
            "chr1\tsrc\tregion\t1\t200\t.\t+\t.\tID=r2;Dbxref=taxon:1;genome=chromosome\n"
        ),
        encoding="utf-8",
    )

    parsed = parse_gff(gff)
    assert any(item.code == "conflicting_region_length" for item in parsed.anomalies)


def test_repeated_organelle_gene_id_is_not_a_fatal_duplicate(tmp_path: Path) -> None:
    gff = tmp_path / "organelle_duplicate.gff"
    gff.write_text(
        (
            "mt\tsrc\tregion\t1\t1000\t.\t+\t.\tID=mt;genome=mitochondrion;Dbxref=taxon:1\n"
            "mt\tsrc\tgene\t10\t100\t.\t+\t.\tID=gene-M1;Name=nad1;gene_biotype=protein_coding\n"
            "mt\tsrc\tgene\t300\t400\t.\t+\t.\tID=gene-M1;Name=nad1;gene_biotype=protein_coding\n"
        ),
        encoding="utf-8",
    )

    parsed = parse_gff(gff)
    assert any(item.code == "duplicate_organelle_gene" for item in parsed.anomalies)
    assert not any(item.code == "duplicate_gene" and item.level == "fatal" for item in parsed.anomalies)

