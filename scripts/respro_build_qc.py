#!/usr/bin/env python3
"""Build a ResPro project database from a source's converted artifacts.

This is the autobump quality-control step: it takes a database's generated
``rules.tsv`` (plus optional ``formula-rules.tsv`` and ``metadata.json``),
fetches the GenBank reference record(s) named by the ``reference_identifier``
column from NCBI, runs ``respro init`` to build a SQLite project database, and
emits a structured summary of the build outcome (success / hard-failure /
skipped rules). The summary is consumed by the autobump workflow and rendered
into the pull-request body so a failed or degraded build is explicit.

The script intentionally does NOT modify any repository artifacts. It writes
the GenBank files, the built ``.db``, and the summary into a working directory
(default: a fresh temp directory) that the caller may discard after the run.

Network access is required only for the NCBI GenBank fetch; ``respro init`` is
invoked with ``--no-additional-info`` so no PubChem/PubMed enrichment happens
and the build is deterministic given fixed inputs.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Any

NCBI_EFETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
# Matches a single NCBI nucleotide accession with an optional version suffix.
# Used to validate path safety before writing <accession>.gb files.
ACCESSION_RE = re.compile(r"^[A-Za-z0-9_]+(\.\d+)?$")

# respro log signals we parse to classify the build outcome. These are matched
# against the combined stdout+stderr of ``respro init``.
HARD_FAIL_RE = re.compile(r"Error:\s*Rules validation failed", re.MULTILINE)
REF_AA_SKIP_RE = re.compile(
    r"(\d+)\s+position\(s\) with rules have reference AA mismatches", re.MULTILINE
)
FEATURE_SKIP_RE = re.compile(
    r"(\d+)\s+rule\(s\) skipped\s*\u2014\s*feature\(s\) not found", re.MULTILINE
)
FORMULA_SKIP_RE = re.compile(
    r"(\d+)\s+formula rule\(s\) skipped", re.MULTILINE
)
SUCCESS_RE = re.compile(r"Project initialised:\s*(\S+)", re.MULTILINE)


def eprint(message: str) -> None:
    print(message, file=sys.stderr)


def norm(value: object) -> str:
    if value is None:
        return ""
    return str(value).strip()


def collect_reference_identifiers(rules_path: Path) -> list[str]:
    """Return the unique, non-empty ``reference_identifier`` values in order."""
    seen: set[str] = set()
    accessions: list[str] = []
    with rules_path.open("r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if reader.fieldnames is None or "reference_identifier" not in reader.fieldnames:
            raise ValueError(
                f"{rules_path} is missing a 'reference_identifier' column"
            )
        for row in reader:
            acc = norm(row.get("reference_identifier"))
            if not acc or acc in seen:
                continue
            if not ACCESSION_RE.match(acc):
                raise ValueError(
                    f"{rules_path} contains an invalid reference_identifier: {acc!r}"
                )
            seen.add(acc)
            accessions.append(acc)
    if not accessions:
        raise ValueError(f"{rules_path} contains no reference_identifier values")
    return accessions


def fetch_genbank(
    accession: str,
    dest_dir: Path,
    api_key: str = "",
    retries: int = 3,
    backoff: float = 2.0,
) -> Path:
    """Fetch one GenBank record from NCBI eutils into ``dest_dir``.

    Uses ``rettype=gbwithparts`` so segmented/large records are returned in
    full. Retries with exponential backoff on transient HTTP errors.
    """
    dest_dir.mkdir(parents=True, exist_ok=True)
    out_path = dest_dir / f"{accession}.gb"
    params = {
        "db": "nuccore",
        "id": accession,
        "rettype": "gbwithparts",
        "retmode": "text",
    }
    if api_key:
        params["api_key"] = api_key
    url = f"{NCBI_EFETCH_URL}?{urllib.parse.urlencode(params)}"

    last_error: str = ""
    for attempt in range(1, retries + 1):
        try:
            request = urllib.request.Request(
                url,
                headers={"User-Agent": "respro-databases-qc/1.0"},
            )
            with urllib.request.urlopen(request, timeout=120) as response:
                payload = response.read()
            if not payload.strip():
                raise RuntimeError("empty response from NCBI efetch")
            out_path.write_bytes(payload)
            return out_path
        except Exception as exc:  # noqa: BLE001 - retry any fetch failure
            last_error = str(exc)
            if attempt < retries:
                wait = backoff ** attempt
                eprint(
                    f"WARNING: NCBI fetch for {accession} failed (attempt "
                    f"{attempt}/{retries}): {last_error}; retrying in {wait:.0f}s"
                )
                time.sleep(wait)
    raise RuntimeError(f"Could not fetch GenBank for {accession}: {last_error}")


def build_respro_command(
    source_name: str,
    rules_path: Path,
    genbank_paths: list[Path],
    output_db: Path,
    formula_rules_path: Path | None,
    metadata_path: Path | None,
    respro_bin: str,
) -> list[str]:
    cmd: list[str] = [
        respro_bin,
        "init",
        "--name",
        source_name,
        "--rules",
        str(rules_path),
        "--output",
        str(output_db),
        "--no-additional-info",
        "--overwrite",
    ]
    for gb in genbank_paths:
        cmd.extend(["--genbank", str(gb)])
    if formula_rules_path and formula_rules_path.is_file():
        cmd.extend(["--formula-rules", str(formula_rules_path)])
    if metadata_path and metadata_path.is_file():
        cmd.extend(["--metadata", str(metadata_path)])
    return cmd


def parse_build_log(log: str) -> dict[str, Any]:
    """Classify a ``respro init`` log + exit code into a structured result.

    respro logs through ``rich``, which wraps long lines at the terminal width.
    We collapse all whitespace (including newlines) into single spaces before
    regex matching so wrapped WARNING/Error messages still match.
    """
    # Collapse runs of whitespace (incl. newlines) into single spaces for
    # robust matching against rich-wrapped log output.
    flat = re.sub(r"\s+", " ", log)

    ref_aa = 0
    feature_skip = 0
    formula_skip = 0
    m = REF_AA_SKIP_RE.search(flat)
    if m:
        ref_aa = int(m.group(1))
    m = FEATURE_SKIP_RE.search(flat)
    if m:
        feature_skip = int(m.group(1))
    m = FORMULA_SKIP_RE.search(flat)
    if m:
        formula_skip = int(m.group(1))

    hard_fail = HARD_FAIL_RE.search(flat) is not None
    success = SUCCESS_RE.search(flat) is not None

    # Extract a concise error excerpt for hard failures.
    error_excerpt = ""
    if hard_fail:
        lines = log.splitlines()
        for i, line in enumerate(lines):
            if "Error:" in line:
                # Grab the error line plus up to 3 following indented detail lines.
                excerpt = [line.strip()]
                for follow in lines[i + 1 : i + 4]:
                    if not follow.strip():
                        break
                    excerpt.append(follow.strip())
                error_excerpt = "\n".join(excerpt)
                break

    return {
        "ref_aa_skipped": ref_aa,
        "feature_skipped": feature_skip,
        "formula_skipped": formula_skip,
        "hard_fail": hard_fail,
        "success": success,
        "error_excerpt": error_excerpt,
    }


def render_summary_markdown(
    source_name: str,
    exit_code: int,
    parsed: dict[str, Any],
    log_path: Path,
    imported_rules: int | None,
) -> str:
    """Render the ``## respro build QC`` markdown block for the PR body."""
    total_skipped = (
        parsed["ref_aa_skipped"]
        + parsed["feature_skipped"]
        + parsed["formula_skipped"]
    )

    if parsed["hard_fail"] or exit_code != 0:
        status_line = "\u26a0\ufe0f **`respro init` FAILED** \u2014 the regenerated TSVs do not build and must be fixed before merge."
    elif total_skipped > 0:
        status_line = (
            f"\u26a0\ufe0f `respro init` succeeded with {total_skipped} skipped rule(s)."
        )
    else:
        status_line = "\u2705 `respro init` succeeded (all rules imported)."

    lines = [
        "## respro build QC",
        "",
        f"**Source:** `{source_name}`",
        f"**Status:** {status_line}",
        f"**Exit code:** `{exit_code}`",
    ]

    if imported_rules is not None:
        lines.append(f"**Imported single rules:** {imported_rules}")

    lines.extend(
        [
            f"**Skipped \u2014 reference AA mismatch:** {parsed['ref_aa_skipped']}",
            f"**Skipped \u2014 feature not in GenBank:** {parsed['feature_skipped']}",
            f"**Skipped \u2014 formula member missing:** {parsed['formula_skipped']}",
            f"**Full log:** `{log_path}`",
        ]
    )

    if parsed["error_excerpt"]:
        lines.extend(["", "```", parsed["error_excerpt"], "```"])

    return "\n".join(lines)


def count_imported_rules(db_path: Path, respro_bin: str) -> int | None:
    """Best-effort count of imported single rules via ``respro manage``."""
    if not db_path.is_file():
        return None
    try:
        result = subprocess.run(
            [respro_bin, "manage", "database", str(db_path), "--list-single"],
            capture_output=True,
            text=True,
            timeout=60,
            check=False,
        )
        if result.returncode != 0:
            return None
        # Count non-empty output lines (the table body); this is an approximate
        # count sufficient for a QC summary, not a precise audit.
        return sum(1 for line in result.stdout.splitlines() if line.strip())
    except Exception:  # noqa: BLE001 - count is best-effort
        return None


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--source-name",
        required=True,
        help="Source database name (e.g. stanford-hiv). Used as the respro project name.",
    )
    parser.add_argument(
        "--output-dir",
        required=True,
        help="Directory containing the source's rules.tsv / metadata.json / formula-rules.tsv.",
    )
    parser.add_argument(
        "--respro-bin",
        default="respro",
        help="Path to the respro executable (default: looks up `respro` on PATH).",
    )
    parser.add_argument(
        "--work-dir",
        default=None,
        help="Working directory for GenBank fetch + built .db + logs (default: fresh temp dir).",
    )
    parser.add_argument(
        "--ncbi-api-key",
        default="",
        help=(
            "Optional NCBI eutils API key (raises the rate limit). "
            "Falls back to the NCBI_EUTILS_API_KEY environment variable."
        ),
    )
    parser.add_argument(
        "--summary-md",
        default=None,
        help="Path to write the markdown summary (default: <work-dir>/qc_summary.md).",
    )
    parser.add_argument(
        "--summary-json",
        default=None,
        help="Path to write the machine-readable summary (default: <work-dir>/qc_summary.json).",
    )
    return parser.parse_args()


def write_failure_summary(
    source_name: str,
    message: str,
    summary_md_path: Path,
    summary_json_path: Path,
) -> None:
    """Emit a failed QC summary for infrastructure errors (exit 2 path).

    Ensures the workflow's summary-reading step always has a markdown file to
    interpolate into the PR body, even when the build never ran.
    """
    summary_md = (
        "## respro build QC\n\n"
        f"**Source:** `{source_name}`\n"
        f"**Status:** \u26a0\ufe0f **`respro init` FAILED** \u2014 the QC step "
        f"could not complete.\n"
        f"**Exit code:** `2`\n"
        f"**Error:** {message}\n"
    )
    summary_md_path.write_text(summary_md + "\n", encoding="utf-8")
    summary_json = {
        "source_name": source_name,
        "status": "failed",
        "exit_code": 2,
        "imported_rules": None,
        "ref_aa_skipped": 0,
        "feature_skipped": 0,
        "formula_skipped": 0,
        "error_excerpt": message,
        "log_path": "",
        "db_path": "",
        "summary_md_path": str(summary_md_path),
    }
    summary_json_path.write_text(
        json.dumps(summary_json, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )


def main() -> int:
    args = parse_args()
    output_dir = Path(args.output_dir).resolve()
    rules_path = output_dir / "rules.tsv"
    formula_rules_path = output_dir / "formula-rules.tsv"
    metadata_path = output_dir / "metadata.json"

    ncbi_api_key = args.ncbi_api_key or os.environ.get("NCBI_EUTILS_API_KEY", "")

    work_dir = Path(args.work_dir).resolve() if args.work_dir else Path(
        tempfile.mkdtemp(prefix="respro_qc_")
    )
    work_dir.mkdir(parents=True, exist_ok=True)
    genbank_dir = work_dir / "genbank"
    db_path = work_dir / f"{args.source_name}.db"
    log_path = work_dir / "respro_init.log"
    summary_md_path = Path(args.summary_md) if args.summary_md else work_dir / "qc_summary.md"
    summary_json_path = (
        Path(args.summary_json) if args.summary_json else work_dir / "qc_summary.json"
    )

    if not rules_path.is_file():
        msg = f"rules.tsv not found: {rules_path}"
        eprint(f"ERROR: {msg}")
        write_failure_summary(args.source_name, msg, summary_md_path, summary_json_path)
        return 2

    try:
        accessions = collect_reference_identifiers(rules_path)
    except ValueError as exc:
        msg = str(exc)
        eprint(f"ERROR: {msg}")
        write_failure_summary(args.source_name, msg, summary_md_path, summary_json_path)
        return 2

    eprint(f"Source: {args.source_name}")
    eprint(f"Reference accessions: {', '.join(accessions)}")

    # Fetch GenBank references from NCBI.
    genbank_paths: list[Path] = []
    try:
        for acc in accessions:
            eprint(f"Fetching GenBank {acc} \u2026")
            gb = fetch_genbank(acc, genbank_dir, api_key=ncbi_api_key)
            genbank_paths.append(gb)
    except RuntimeError as exc:
        msg = str(exc)
        eprint(f"ERROR: {msg}")
        write_failure_summary(args.source_name, msg, summary_md_path, summary_json_path)
        return 2

    # Build the respro command and run it.
    cmd = build_respro_command(
        source_name=args.source_name,
        rules_path=rules_path,
        genbank_paths=genbank_paths,
        output_db=db_path,
        formula_rules_path=formula_rules_path,
        metadata_path=metadata_path,
        respro_bin=args.respro_bin,
    )
    eprint(f"Running: {' '.join(cmd)}")
    result = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        timeout=600,
        check=False,
    )
    log = result.stdout + result.stderr
    log_path.write_text(log, encoding="utf-8")

    parsed = parse_build_log(log)
    imported_rules = count_imported_rules(db_path, args.respro_bin)

    summary_md = render_summary_markdown(
        source_name=args.source_name,
        exit_code=result.returncode,
        parsed=parsed,
        log_path=log_path,
        imported_rules=imported_rules,
    )
    summary_md_path.write_text(summary_md + "\n", encoding="utf-8")

    summary_json: dict[str, Any] = {
        "source_name": args.source_name,
        "status": (
            "failed" if parsed["hard_fail"] or result.returncode != 0
            else "skipped" if (
                parsed["ref_aa_skipped"]
                + parsed["feature_skipped"]
                + parsed["formula_skipped"]
            )
            else "ok"
        ),
        "exit_code": result.returncode,
        "imported_rules": imported_rules,
        "ref_aa_skipped": parsed["ref_aa_skipped"],
        "feature_skipped": parsed["feature_skipped"],
        "formula_skipped": parsed["formula_skipped"],
        "error_excerpt": parsed["error_excerpt"],
        "log_path": str(log_path),
        "db_path": str(db_path),
        "summary_md_path": str(summary_md_path),
    }
    summary_json_path.write_text(
        json.dumps(summary_json, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )

    eprint(f"Written {summary_md_path}")
    eprint(f"Written {summary_json_path}")

    # Exit non-zero ONLY on hard respro failure. Soft skips are exit 0 so the
    # PR step still runs; the workflow separately marks the run red on failure.
    if parsed["hard_fail"] or result.returncode != 0:
        eprint(f"ERROR: respro init failed (exit {result.returncode}). See {log_path}")
        return 1

    eprint("Done.")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    finally:
        # The workflow passes an explicit --work-dir it manages; for ad-hoc
        # local runs we leave the temp dir in place for inspection.
        pass
