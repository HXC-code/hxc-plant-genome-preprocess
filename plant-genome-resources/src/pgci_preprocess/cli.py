from __future__ import annotations

import argparse
import sys
from pathlib import Path

from pgci_preprocess.pipeline import PipelineOptions, run_project


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Build strict one-gene-one-transcript WGDI inputs from NCBI downloads."
    )
    parser.add_argument(
        "--project-root",
        type=Path,
        default=Path(__file__).resolve().parents[2],
        help="Directory containing resources/ and data/downloads/.",
    )
    parser.add_argument(
        "--species",
        default="",
        help="Comma-separated show_id values. Empty means all species.",
    )
    parser.add_argument("--min-genes-per-sequence", type=int, default=50)
    parser.add_argument("--coverage", type=float, default=0.90)
    parser.add_argument("--max-sequences", type=int, default=200)
    parser.add_argument("--no-overwrite", action="store_true")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    show_ids = {
        item.strip() for item in args.species.split(",") if item.strip()
    } or None
    options = PipelineOptions(
        project_root=args.project_root,
        min_genes_per_sequence=args.min_genes_per_sequence,
        target_gene_coverage=args.coverage,
        max_sequences=args.max_sequences,
        overwrite=not args.no_overwrite,
    )
    results = run_project(args.project_root, show_ids=show_ids, options=options)
    failures = sum(1 for result in results if result.get("status") != "OK")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())

