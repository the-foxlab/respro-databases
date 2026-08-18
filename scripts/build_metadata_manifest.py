#!/usr/bin/env python3
"""Generate a repository-wide manifest of database metadata files."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

REQUIRED_METADATA_KEYS = (
    "maintainers",
    "contact",
    "publication_pmid",
    "website",
    "description",
    "maintainer_update",
    "license",
    "tsv_checksum",
)

# Optional metadata key supported by ResPro init but intentionally excluded
# from repository-level manifest entries.
EXCLUDED_METADATA_KEYS = {
    "interpretation_algorithms",
}


def load_metadata(metadata_path: Path) -> dict[str, Any]:
    with metadata_path.open("r", encoding="utf-8") as fh:
        metadata = json.load(fh)

    if not isinstance(metadata, dict):
        raise ValueError(f"Metadata must be a JSON object: {metadata_path}")

    keys = set(metadata.keys())
    required = set(REQUIRED_METADATA_KEYS)

    missing = sorted(required - keys)
    unknown = sorted(keys - required - EXCLUDED_METADATA_KEYS)
    if missing or unknown:
        parts: list[str] = []
        if missing:
            parts.append(f"missing keys={missing}")
        if unknown:
            parts.append(f"unknown keys={unknown}")
        details = "; ".join(parts)
        raise ValueError(f"Invalid metadata schema in {metadata_path}: {details}")

    return {k: metadata[k] for k in REQUIRED_METADATA_KEYS}


def build_manifest(repo_root: Path) -> dict[str, Any]:
    databases_dir = repo_root / "databases"
    entries: list[dict[str, Any]] = []

    for metadata_path in sorted(databases_dir.glob("*/output/metadata.json")):
        source_dir = metadata_path.parent.parent
        source_name = source_dir.name
        rules_path = metadata_path.parent / "rules.tsv"
        formula_rules_path = metadata_path.parent / "formula-rules.tsv"
        example_fasta_path = source_dir / "example" / "example.fasta"

        if not rules_path.exists():
            raise FileNotFoundError(f"Missing rules.tsv for {source_name}: {rules_path}")

        metadata = load_metadata(metadata_path)

        entry = {
            "source_name": source_name,
            "metadata_path": str(metadata_path.relative_to(repo_root)),
            "rules_path": str(rules_path.relative_to(repo_root)),
            "formula_rules_path": (
                str(formula_rules_path.relative_to(repo_root))
                if formula_rules_path.exists()
                else ""
            ),
            "example_fasta_path": (
                str(example_fasta_path.relative_to(repo_root))
                if example_fasta_path.exists()
                else ""
            ),
            "metadata": metadata,
        }
        entries.append(entry)

    return {
        "manifest_version": 1,
        "databases": entries,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--repo-root",
        default=".",
        help="Path to repository root (default: current directory)",
    )
    parser.add_argument(
        "--output",
        default="databases/manifest.json",
        help="Manifest output path relative to repo root",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    repo_root = Path(args.repo_root).resolve()
    output_path = repo_root / args.output

    manifest = build_manifest(repo_root)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8") as fh:
        json.dump(manifest, fh, ensure_ascii=False, indent=2)
        fh.write("\n")

    print(f"Written {output_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
