#!/usr/bin/env python3
"""Convert a CHARMD TSV export into ResPro-compatible artifacts.

The converter only parses and normalizes the CHARMD export; it does not
validate mutations against a reference protein sequence or infer a
fold-IC50-to-phenotype interpretation. Both are respro's responsibility
downstream. Values are passed through as close to the source as possible.

The rule reference is HCMV strain Merlin, GenBank accession X17403.1.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import re
import sys
import urllib.request
from collections.abc import Iterable
from datetime import datetime, timezone
from pathlib import Path

SOURCE_URL = "https://gitlab.com/vtilloy/charmd-data-db/-/raw/master/charmd-data.tsv"
REFERENCE_IDENTIFIER = "X17403.1"
EXPECTED_TAXON_ID = "10359"
SOURCE_LABEL = "CHARMD"

RULES_COLUMNS = [
    "feature",
    "reference_identifier",
    "position",
    "reference",
    "mutation",
    "antiviral",
    "phenotype",
    "fold_ic50",
    "publication",
    "source",
    "comment",
]

NON_MIGRATED_COLUMNS = [
    "reason",
    "source_line",
    "feature",
    "source_mutation",
    "antiviral",
    "details",
]

REQUIRED_INPUT_COLUMNS = {
    "taxon_id",
    "name_mutation",
    "gene",
    "pubmed_id",
    "marker_transfert",
    "source_origin",
    "short_name_antiviral",
    "EC50_fold_increase",
    "phenotype",
    "short_name_virus",
    "assay_details",
}

DRUG_BY_CODE = {
    "GCV": "ganciclovir",
    "CGCV": "ganciclovir",
    "FOS": "foscarnet",
    "CDV": "cidofovir",
    "MBV": "maribavir",
    "LMV": "letermovir",
    "CPV": "cyclopropavir",
    "BCV": "brincidofovir",
}
ALL_ANTIVIRALS = sorted(set(DRUG_BY_CODE.values()))

# Official CNR fold-IC50 (EC50 ratio) breakpoints for HCMV antiviral
# resistance interpretation:
# https://www.unilim.fr/cnr-herpesvirus/outils/codexmv/database/notice/
#
#   EC50 ratio < 0.5 : possible hypersensitivity, no resistance
#   0.5 < ratio < 2.0: no resistance            (threshold 2.0, Chou 2024
#                                                 PMID:37506264; published
#                                                 value 1.9 rounded to 2.0
#                                                 by consortium expertise)
#   2.0 < ratio < 3  : low-level resistance
#   3 < ratio < 5    : intermediate level resistance
#   5 < ratio        : high level resistance
#
# Mapping to the respro rank vocabulary (rules-format.md#phenotype-normalization):
# hypersensitivity has no rank and folds into susceptible (rank 1); the tier
# boundaries are the thresholds at which each higher rank is reached. Boundary
# values use >= semantics, matching respro's breakpoint logic: a fold of 2.0 is
# low-level resistance, 3.0 intermediate, 5.0 high-level.
#
# Used for rule-level phenotype derivation when the source phenotype column
# is empty but a numeric fold-IC50 exists. The fold-IC50 values themselves are
# kept in rules.tsv as data (they feed the report's IC50 distribution plots),
# but respro interpretation is phenotype-only: no by_fold_ic50 algorithm is
# configured, so the CNR breakpoints are applied here, once, at curation time.
CNR_FOLD_IC50_THRESHOLDS = {
    "susceptible": 0.0,
    "low-level resistance": 2.0,
    "intermediate": 3.0,
    "high-level resistance": 5.0,
}


def derive_phenotype_from_fold(value: float) -> str:
    """Classify a numeric fold-IC50 into the respro rank vocabulary via the
    official CNR breakpoints (see CNR_FOLD_IC50_THRESHOLDS)."""
    if value >= 5.0:
        return "high-level resistance"
    if value >= 3.0:
        return "intermediate"
    if value >= 2.0:
        return "low-level resistance"
    return "susceptible"


class ConversionError(ValueError):
    """A source row cannot be represented as a ResPro rule."""

    def __init__(self, reason: str, details: str = "") -> None:
        super().__init__(reason)
        self.reason = reason
        self.details = details


def eprint(message: str) -> None:
    print(message, file=sys.stderr)


def norm(value: object) -> str:
    if value is None:
        return ""
    return str(value).strip()


def clean_text(value: object) -> str:
    text = norm(value)
    text = re.sub(r"[\t\r\n]+", " ", text)
    text = re.sub(r"\s+", " ", text)
    return text.strip(" ;|")


def join_unique(values: Iterable[str], separator: str = ",") -> str:
    seen: set[str] = set()
    result: list[str] = []
    for value in values:
        cleaned = clean_text(value)
        if not cleaned or cleaned in seen:
            continue
        seen.add(cleaned)
        result.append(cleaned)
    return separator.join(result)


def fetch_source_commit_date() -> str:
    """Query the GitLab API for the latest commit date of the source file.

    Returns the commit date as a YYYY-MM-DD string, or today's date as
    fallback. Used for metadata ``maintainer_update`` so it tracks the
    upstream data version rather than the conversion run date.
    """
    # URL-encoded project path (vtilloy/charmd-data-db) for the GitLab API.
    api_url = (
        "https://gitlab.com/api/v4/projects/vtilloy%2Fcharmd-data-db"
        "/repository/commits"
        "?path=charmd-data.tsv&ref_name=master&per_page=1"
    )
    try:
        with urllib.request.urlopen(api_url, timeout=60) as response:
            commits = json.loads(response.read().decode("utf-8"))
        return commits[0]["committed_date"][:10]
    except Exception as exc:
        eprint(
            f"WARNING: Could not fetch commit date from GitLab API ({exc}); using today."
        )
        return datetime.now(timezone.utc).strftime("%Y-%m-%d")


def download_source(source_url: str) -> str:
    request = urllib.request.Request(
        source_url,
        headers={"User-Agent": "respro-databases-charmd-converter/1.0"},
    )
    with urllib.request.urlopen(request, timeout=60) as response:
        payload = response.read()
    try:
        return payload.decode("utf-8-sig")
    except UnicodeDecodeError as exc:
        raise ValueError("Source TSV is not valid UTF-8/UTF-8-SIG") from exc


def validate_input_header(reader: csv.DictReader) -> None:
    header = set(reader.fieldnames or [])
    missing = sorted(REQUIRED_INPUT_COLUMNS - header)
    if missing:
        raise ValueError(f"Missing required source columns: {', '.join(missing)}")


def parse_mutation(source_mutation: str) -> tuple[int, str, str]:
    """Parse a CHARMD mutation token into (position, reference, mutation).

    No reference sequence is consulted; the reference allele is taken
    directly from the source token, and respro validates it downstream.

    Indel policy (see docs/docs/rules-format.md for the respro grammar):

    - ``X{n}del`` (single-residue deletion with anchor AA) is emitted as the
      respro helper form ``X{n}del``; respro resolves the upstream anchor at
      init. Kept.
    - ``{a}-{b}del{SEQ}`` (range deletion with explicit deleted sequence) is
      emitted with ``reference`` = deleted block and ``mutation`` = ``del{a}_{b}``;
      respro resolves the surviving anchor at init. Kept.
    - ``X{n}del{L}``, ``{a}-{b}del`` (no deleted sequence), and ``in{pos}{AA}``
      (no anchor AA) cannot be made deterministic without the reference
      sequence and are rejected as ``indel_undeterministic``.
    """
    text = re.sub(r"\s+", "", source_mutation).upper()

    substitution = re.fullmatch(r"([A-Z])(\d+)([A-Z*]|STOP)", text)
    if substitution:
        ref = substitution.group(1)
        position = int(substitution.group(2))
        alternate = substitution.group(3)
        alternate = "*" if alternate == "STOP" else alternate
        return position, ref, alternate

    # in{pos}{AA}: insertion without an anchor AA — not deterministic.
    if re.fullmatch(r"IN(\d+)([A-Z]+)", text):
        raise ConversionError(
            "indel_undeterministic",
            f"Insertion {source_mutation} has no anchor AA; cannot derive a deterministic reference",
        )

    # X{n}del{L}: single-residue deletion with explicit length — the second
    # deleted residue is unspecified, so the deleted block is unknown.
    if re.fullmatch(r"([A-Z])(\d+)DEL\d+", text):
        raise ConversionError(
            "indel_undeterministic",
            f"Deletion {source_mutation} specifies a length but only the first anchor AA; "
            "deleted block cannot be derived without the reference sequence",
        )

    # X{n}del: single-residue deletion with anchor AA — respro helper form.
    single_deletion = re.fullmatch(r"([A-Z])(\d+)DEL", text)
    if single_deletion:
        ref = single_deletion.group(1)
        position = int(single_deletion.group(2))
        return position, ref, f"{ref}{position}del"

    # {a}-{b}del with no deleted sequence — cannot populate reference/mutation.
    if re.fullmatch(r"(\d+)-(\d+)DEL", text):
        raise ConversionError(
            "indel_undeterministic",
            f"Range deletion {source_mutation} supplies no deleted sequence; "
            "deleted block cannot be derived without the reference sequence",
        )

    # {a}-{b}del{SEQ}: range deletion with explicit deleted sequence.
    # Emitted as the respro anchor-less helper form {SEQ}{a}del (e.g. NSS449del),
    # matching respro's _RE_ANCHORLESS_DEL = ^([A-Za-z]+)\d+del$. respro resolves
    # the upstream anchor at init. The reference column holds the first deleted
    # residue (the AA at position `a`), which respro validates against the
    # GenBank translation before anchor resolution overwrites it.
    range_deletion = re.fullmatch(r"(\d+)-(\d+)DEL([A-Z]+)", text)
    if range_deletion:
        start = int(range_deletion.group(1))
        end = int(range_deletion.group(2))
        deleted = range_deletion.group(3)
        if end < start:
            raise ConversionError("invalid_deletion_range", f"Deletion range {start}-{end} is reversed")
        if len(deleted) != end - start + 1:
            raise ConversionError(
                "deletion_sequence_length_mismatch",
                f"Deleted sequence {deleted} does not match range {start}-{end} length {end - start + 1}",
            )
        return start, deleted[0], f"{deleted}{start}del"

    raise ConversionError("unparseable_mutation", f"Unsupported mutation token: {source_mutation}")


def normalize_source_phenotype(value: str) -> str:
    token = clean_text(value).lower()
    if not token:
        return ""
    if "polymorph" in token:
        return "sensitive"
    if token in {"resistant", "resistance", "res", "r"}:
        return "resistant"
    if token in {"intermediate", "interm", "i"}:
        return "intermediate"
    if token in {"sensitive", "susceptible", "sensi", "sens", "s"}:
        return "sensitive"
    return ""


def normalize_fold(raw: str) -> str:
    # Only normalize the French decimal comma; the value is otherwise passed
    # through unchanged, including ranges or inequalities, with no bounds assumed.
    return clean_text(raw).replace(",", ".")


def parse_fold_ic50(raw: str) -> tuple[float | None, str]:
    """Parse a fold-IC50 token into (numeric_upper_bound, original_string).

    Plain numbers are returned as ``(value, "")`` — the original string is
    empty because no annotation is needed.

    Ranges (``a-b``) and inequalities (``>n``, ``<n``) are reduced to their
    *upper boundary* so a deterministic numeric value feeds the mean. The
    verbatim source token is returned in the second element so it can be
    recorded in the comment, preserving the original measurement semantics.
    """
    text = normalize_fold(raw)
    if not text:
        return None, ""

    # Plain number.
    try:
        return float(text), ""
    except ValueError:
        pass

    # Range: a-b  → upper boundary is b.
    range_match = re.fullmatch(r"([0-9]*\.?[0-9]+)\s*-\s*([0-9]*\.?[0-9]+)", text)
    if range_match:
        upper = float(range_match.group(2))
        return upper, text

    # Inequality: >n, >=n, <n, <=n  → boundary is n.
    ineq_match = re.fullmatch(r"([<>]=?)\s*([0-9]*\.?[0-9]+)", text)
    if ineq_match:
        boundary = float(ineq_match.group(2))
        return boundary, text

    # Unrecognised non-numeric token: no numeric value, kept verbatim.
    return None, text


def publication_values(row: dict[str, str]) -> list[str]:
    # PMIDs from the dedicated pubmed_id column plus any cross-referenced in
    # the observation field ("Autre publi ID : <PMID> ..."). Both sources are
    # merged so a rule carries every supporting publication.
    tokens = re.findall(r"\b\d{7,9}\b", norm(row.get("pubmed_id")))
    tokens.extend(re.findall(r"\b\d{7,9}\b", norm(row.get("observation"))))
    return [f"PMID:{token}" for token in tokens]


def row_comment(row: dict[str, str]) -> str:
    parts: list[str] = []
    for field, label in (
        ("source_origin", "Source origin"),
        ("marker_transfert", "Marker transfer"),
        ("assay_details", "Assay details"),
    ):
        value = clean_text(row.get(field))
        if value:
            parts.append(f"{label}: {value}")

    # The observation column is unstructured free text that often carries
    # cross-referenced PMIDs ("Autre publi ID"). Those are harvested into the
    # publication column by publication_values(); the residual text (clinical
    # context, alternative fold values, curation notes) is recorded here.
    observation = clean_text(row.get("observation"))
    if observation:
        residual = re.sub(r"Autre publi ID\s*:", "", observation)
        residual = re.sub(r"\b\d{7,9}\b", "", residual)
        residual = clean_text(residual)
        if residual:
            parts.append(f"Observation: {residual}")

    viability = clean_text(row.get("viability"))
    if viability:
        parts.append(f"Viability: {viability}")

    return join_unique(parts, separator=" | ")


def add_non_migrated(
    bucket: list[dict[str, str]],
    row: dict[str, str],
    reason: str,
    details: str = "",
) -> None:
    # Carry the source-row comment (observation, viability, etc.) into the
    # audit trail so non-migrated rows do not silently lose their context.
    comment = row_comment(row)
    combined = clean_text(details)
    if comment:
        combined = f"{combined} | {comment}" if combined else comment
    bucket.append(
        {
            "reason": reason,
            "source_line": norm(row.get("_source_line")),
            "feature": norm(row.get("gene")),
            "source_mutation": norm(row.get("name_mutation")),
            "antiviral": norm(row.get("short_name_antiviral")),
            "details": combined,
        }
    )


def convert(source_rows: list[dict[str, str]]) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    # Aggregate by (feature, position, reference, mutation, antiviral). Source
    # rows that map to the same rule are merged so that each ResPro rule is
    # emitted once:
    #   - polymorphisms (same feature/pos/ref/alt) collapse to a single row;
    #   - fold-IC50 values are averaged and the contributing values are noted;
    #   - distinct PMIDs are appended to the publication column.
    # The aggregation key includes the antiviral so a mutation tested against
    # two drugs stays as two rules.
    aggregated: dict[tuple, dict[str, object]] = {}
    non_migrated: list[dict[str, str]] = []

    for row in source_rows:
        taxon_id = norm(row.get("taxon_id"))
        if taxon_id and taxon_id != EXPECTED_TAXON_ID:
            add_non_migrated(non_migrated, row, "unsupported_taxon", f"Expected taxon {EXPECTED_TAXON_ID}")
            continue

        feature = f"HCMV{norm(row.get('gene')).upper()}"
        source_mutation = norm(row.get("name_mutation"))
        drug_code = norm(row.get("short_name_antiviral")).upper()

        if not feature:
            add_non_migrated(non_migrated, row, "missing_feature")
            continue
        if not source_mutation:
            add_non_migrated(non_migrated, row, "missing_mutation")
            continue

        try:
            position, reference, mutation = parse_mutation(source_mutation)
        except ConversionError as exc:
            add_non_migrated(non_migrated, row, exc.reason, exc.details)
            continue

        phenotype = normalize_source_phenotype(row.get("phenotype"))
        fold_ic50 = normalize_fold(row.get("EC50_fold_increase"))
        publication = publication_values(row)
        comment = row_comment(row)

        if not drug_code:
            if phenotype != "sensitive":
                add_non_migrated(
                    non_migrated,
                    row,
                    "missing_antiviral",
                    "ResPro requires antiviral; CHARMD row was retained in the audit only",
                )
                continue
            # A natural polymorphism without a drug association is a sensitive
            # rule against every antiviral in the database. Each fan-out target
            # is its own aggregation key, so duplicate polymorphism rows for the
            # same mutation collapse into one.
            for antiviral in ALL_ANTIVIRALS:
                _aggregate(
                    aggregated,
                    feature,
                    position,
                    reference,
                    mutation,
                    antiviral,
                    "sensitive",
                    fold_ic50,
                    publication,
                    comment,
                )
            continue

        antiviral = DRUG_BY_CODE.get(drug_code)
        if not antiviral:
            add_non_migrated(non_migrated, row, "unsupported_antiviral", f"Unknown CHARMD drug code {drug_code}")
            continue

        _aggregate(
            aggregated,
            feature,
            position,
            reference,
            mutation,
            antiviral,
            phenotype,
            fold_ic50,
            publication,
            comment,
        )

    rules = [_finalize(key, data) for key, data in aggregated.items()]
    rules.sort(
        key=lambda row: (
            row["reference_identifier"],
            row["feature"],
            int(row["position"]),
            row["reference"],
            row["mutation"],
            row["antiviral"],
            row["publication"],
        )
    )
    non_migrated.sort(
        key=lambda row: (
            row["reason"],
            row["feature"],
            row["source_mutation"],
            row["antiviral"],
            int(row["source_line"] or "0"),
        )
    )
    return rules, non_migrated


def _aggregate(
    bucket: dict[tuple, dict[str, object]],
    feature: str,
    position: int,
    reference: str,
    mutation: str,
    antiviral: str,
    phenotype: str,
    fold_ic50: str,
    publication: list[str],
    comment: str,
) -> None:
    key = (feature, position, reference, mutation, antiviral)
    entry = bucket.get(key)
    if entry is None:
        entry = {
            "feature": feature,
            "position": position,
            "reference": reference,
            "mutation": mutation,
            "antiviral": antiviral,
            "phenotypes": [],
            "fold_values": [],
            "fold_tokens": [],
            "publications": [],
            "comments": [],
        }
        bucket[key] = entry

    if phenotype:
        entry["phenotypes"].append(phenotype)
    if fold_ic50:
        numeric, original = parse_fold_ic50(fold_ic50)
        if numeric is not None:
            entry["fold_values"].append(numeric)
        if original:
            # Non-numeric token (range/inequality) reduced to its upper
            # boundary for the mean; the verbatim string is kept for the
            # comment so the original measurement semantics are preserved.
            entry["fold_tokens"].append(original)
    if publication:
        entry["publications"].extend(publication)
    if comment:
        entry["comments"].append(comment)


def _finalize(key: tuple, data: dict[str, object]) -> dict[str, str]:
    feature, position, reference, mutation, antiviral = key

    phenotypes = [p for p in data["phenotypes"] if p]
    if phenotypes:
        unique = set(phenotypes)
        if "resistant" in unique and "sensitive" in unique:
            phenotype = "contradictory"
        else:
            phenotype = phenotypes[0]
    else:
        # No source phenotype label: derive one from the mean fold-IC50 that
        # is displayed in the rule, using the official CNR breakpoints. The
        # mean is the value a fold-IC50-based classification would act on, so
        # the derived phenotype stays consistent with the displayed
        # measurement. Rules without any fold value stay empty (respro stores
        # empty as rank 0, unknown).
        fold_values = data["fold_values"]
        if fold_values:
            mean = sum(fold_values) / len(fold_values)
            phenotype = derive_phenotype_from_fold(mean)
        else:
            phenotype = ""

    fold_values = data["fold_values"]
    fold_tokens = data["fold_tokens"]
    if fold_values:
        mean = sum(fold_values) / len(fold_values)
        fold_ic50 = f"{mean:g}"
    else:
        fold_ic50 = ""

    publication = join_unique(data["publications"])

    comment_parts: list[str] = []
    if len(fold_values) > 1 or fold_tokens:
        # Annotate the contributing fold values. Numeric values are listed
        # directly; non-numeric tokens (ranges/inequalities, reduced to their
        # upper boundary for the mean) are shown verbatim so the original
        # measurement semantics are preserved.
        parts: list[str] = [f"{v:g}" for v in fold_values if isinstance(v, float)]
        parts.extend(fold_tokens)
        values_str = ", ".join(parts)
        comment_parts.append(f"Fold IC50 values: {values_str}; mean displayed")
    source_comment = join_unique(data["comments"], separator=" | ")
    if source_comment:
        comment_parts.append(source_comment)
    comment = " | ".join(comment_parts)

    return {
        "feature": feature,
        "reference_identifier": REFERENCE_IDENTIFIER,
        "position": str(position),
        "reference": reference,
        "mutation": mutation,
        "antiviral": antiviral,
        "phenotype": phenotype,
        "fold_ic50": fold_ic50,
        "publication": publication,
        "source": SOURCE_LABEL,
        "comment": comment,
    }


def tsv_from_rows(rows: list[dict[str, str]], columns: list[str]) -> str:
    buffer = io.StringIO(newline="")
    writer = csv.DictWriter(buffer, fieldnames=columns, delimiter="\t", lineterminator="\n", extrasaction="ignore")
    writer.writeheader()
    writer.writerows(rows)
    return buffer.getvalue()


def checksum(content: str) -> str:
    return "sha256:" + hashlib.sha256(content.encode("utf-8")).hexdigest()


def non_migrated_text(rows: list[dict[str, str]]) -> str:
    lines = [
        "# Non-migrated CHARMD rows",
        f"# Reference: {REFERENCE_IDENTIFIER}",
        "# Columns: " + "\t".join(NON_MIGRATED_COLUMNS),
        "",
    ]
    lines.append(tsv_from_rows(rows, NON_MIGRATED_COLUMNS).rstrip("\n"))
    return "\n".join(lines) + "\n"


def build_metadata(source_date: str, rules_content: str) -> dict[str, object]:
    return {
        "maintainers": ["Valentin Tilloy", "Daniel Diaz-Gonzalez"],
        "contact": "valentin.tilloy@unilim.fr",
        "publication_pmid": "39349222",
        "website": "https://gitlab.com/vtilloy/charmd-data-db",
        "description": (
            "ResPro-conversion of the Comprehensive Herpesviruses Antiviral drug Resistance Mutation Database (CHARMD): https://www.unilim.fr/cnr-herpesvirus/outils/codexmv/)."
        ),
        "maintainer_update": source_date,
        "license": "CC-BY-NC-4.0",
        "tsv_checksum": checksum(rules_content),
        "interpretation_algorithms": [
            {
                "name": "drug_groups",
                "groups": {
                    "Nucleoside analogs": ["ganciclovir"],
                    "Nucleotide analogs": ["cidofovir", "brincidofovir"],
                    "Methylenecyclopropane nucleoside analog": ["cyclopropavir"],
                    "Pyrophosphate analog": ["foscarnet"],
                    "Terminase inhibitor": ["letermovir"],
                    "Viral kinase inhibitor": ["maribavir"],
                },
            },
            {
                "name": "drug_interpretation",
                "method": "by_phenotype",
            },
            {
                "name": "drug_alias",
                "groups": {
                    "ganciclovir": "GCV",
                    "foscarnet": "FOS",
                    "cidofovir": "CDV",
                    "maribavir": "MBV",
                    "letermovir": "LMV",
                    "cyclopropavir": "CPV",
                    "brincidofovir": "BCV",
                },
            },
        ],
    }


def parse_source_rows(text: str) -> list[dict[str, str]]:
    reader = csv.DictReader(io.StringIO(text), delimiter="\t")
    validate_input_header(reader)
    rows: list[dict[str, str]] = []
    for source_line, row in enumerate(reader, start=2):
        if not any(norm(value) for value in row.values()):
            continue
        normalized = {key: norm(value) for key, value in row.items()}
        normalized["_source_line"] = str(source_line)
        rows.append(normalized)
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description="Convert CHARMD TSV to ResPro TSV artifacts")
    parser.add_argument("--source-url", default=SOURCE_URL, help="Upstream TSV URL")
    parser.add_argument(
        "--output-dir",
        default=str(Path(__file__).resolve().parent.parent / "output"),
        help="Output directory for generated artifacts",
    )
    args = parser.parse_args()

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    maintainer_update = fetch_source_commit_date()
    eprint(f"Source date: {maintainer_update}")
    eprint(f"Downloading {args.source_url} …")
    source_text = download_source(args.source_url)

    source_rows = parse_source_rows(source_text)
    eprint(f"Parsed {len(source_rows)} source rows")
    rules_rows, non_migrated_rows = convert(source_rows)

    rules_content = tsv_from_rows(rules_rows, RULES_COLUMNS)
    rules_path = out_dir / "rules.tsv"
    rules_path.write_text(rules_content, encoding="utf-8")
    eprint(f"Written {rules_path} ({len(rules_rows)} rows).")

    formula_path = out_dir / "formula-rules.tsv"
    formula_path.unlink(missing_ok=True)

    metadata = build_metadata(maintainer_update, rules_content)
    metadata_path = out_dir / "metadata.json"
    metadata_path.write_text(json.dumps(metadata, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    eprint(f"Written {metadata_path}.")

    non_migrated_path = out_dir / "non-migrated-rules.txt"
    non_migrated_path.write_text(non_migrated_text(non_migrated_rows), encoding="utf-8")
    eprint(f"Written {non_migrated_path} ({len(non_migrated_rows)} rows).")
    eprint("Done.")


if __name__ == "__main__":
    main()
