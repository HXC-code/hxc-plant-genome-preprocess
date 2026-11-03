#!/usr/bin/env python3
from __future__ import annotations

import argparse
import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from pgci_preprocess.validate import (  # noqa: E402
    ValidationOptions,
    resolve_layout,
    sync_compat_links,
    validate_project,
)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Validate WGDI preprocessing output from the project manifest."
    )
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--expected-count", type=int)
    parser.add_argument("--no-compat-links", action="store_true")
    parser.add_argument("--skip-hashes", action="store_true")
    parser.add_argument("--refresh-links", action="store_true")
    parser.add_argument("--report", type=Path)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    workspace, canonical_root = resolve_layout(args.project_root)
    if args.refresh_links:
        for message in sync_compat_links(canonical_root, workspace):
            print(message)
    result = validate_project(
        ValidationOptions(
            project_root=workspace,
            expected_count=args.expected_count,
            require_compat_links=not args.no_compat_links,
            verify_hashes=not args.skip_hashes,
            report_path=args.report,
        )
    )
    for row in result.rows:
        print(f"[{row['status']}] {row['code']} {row['problems']}")
    passed = sum(1 for row in result.rows if row["status"] == "PASS")
    failed = sum(1 for row in result.rows if row["status"] == "FAIL")
    print(f"Total: {passed} PASS / {failed} FAIL")
    print(f"Report: {result.report_path}")
    return 0 if result.ok else 1


if __name__ == "__main__":
    raise SystemExit(main())

