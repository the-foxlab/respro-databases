#!/usr/bin/env python3
"""
Convert HerpesDRG TSV into ResPro-compatible TSV artifacts.

Strategy for handling multiple IC50 values and phenotypes from the same publication:
- When a mutation-antiviral pair appears in multiple sources (by publication),
  the data is aggregated:
    1. Fold-change values: Calculate mean. Add comment about count only if N>1.
  2. Phenotypes: Detect conflicts (resistant vs sensitive -> contradictory).
  3. Publications: Join with commas.

Outputs:
- rules.tsv (required, contains aggregated atomic rules)
- formula-rules.tsv (optional, only when grouped co-mutation rows are emitted)
- metadata.json (required)
- non-migrated-rules.txt (audit trail)
"""

import argparse
import csv
import hashlib
import json
import re
import sys
import urllib.request
from datetime import datetime, timezone
from pathlib import Path


def eprint(msg: str) -> None:
    """Print a message to stderr."""
    print(msg, file=sys.stderr)

SOURCE_URL = "https://raw.githubusercontent.com/ojcharles/herpesdrg-db/main/herpesdrg-db.tsv"

RULES_COLUMNS = [
    "feature",
    "reference_identifier",
    "position",
    "reference",
    "mutation",
    "antiviral",
    "member_id",
    "phenotype",
    "fold_ic50",
    "publication",
    "source",
    "comment",
]

FORMULA_COLUMNS = [
    "group_id",
    "antiviral",
    "expression",
    "phenotype",
    "fold_ic50",
    "publication",
    "source",
    "comment",
]

NON_MIGRATED_COLUMNS = [
    "reason",
    "mutation_id",
    "virus",
    "feature",
    "aa_change",
    "co_feature",
    "co_aa",
    "status",
    "note",
    "details",
]

ANTIVIRAL_COLUMNS = [
    "Ganciclovir",
    "Aciclovir",
    "Cidofovir",
    "Foscarnet",
    "Brincidofovir",
    "Letermovir",
    "Brivudine",
    "Penciclovir",
    "Tomeglovir",
    "Maribavir",
    "Cyclopropavir",
    "Amenamevir",
    "Pritelivir",
]

REQUIRED_INPUT_COLUMNS = {
    "mutation_id",
    "virus",
    "gene",
    "aa_change",
    "ref_title",
    "ref_link",
    "ref_doi",
    "co_gene",
    "co_aa",
    "note",
    "status",
    *ANTIVIRAL_COLUMNS,
}

# User-provided fixed mappings.
REFERENCE_BY_VIRUS = {
    "hcmv": "NC_006273",
    "vzv": "NC_001348",
    "adeno": "AC_000008.1",
    "hsv1": "NC_001806",
    "hsv2": "NC_001798",
    "hhv6b": "NC_001664",
}


def norm(v: object) -> str:
    if v is None:
        return ""
    return str(v).strip()


def normalize_virus(raw_virus: str) -> str:
    return raw_virus.strip().lower()


def should_exclude(status: str) -> tuple[bool, str]:
    if status.strip().upper() != "A":
        return True, "status_not_active"
    return False, ""


def parse_mutation(aa_change: str) -> tuple[str, str, int]:
    text = aa_change.strip().replace(" ", "")

    insertion = re.fullmatch(r"([A-Za-z])(\d+)(?:insert|ins)([A-Za-z]+)", text, flags=re.IGNORECASE)
    if insertion:
        ref = insertion.group(1).upper()
        pos = int(insertion.group(2))
        inserted = insertion.group(3).upper()
        return ref, f"{ref}{pos}{ref}{inserted}", pos

    sub = re.fullmatch(r"([A-Za-z*])(\d+)([A-Za-z*])", text)
    if sub:
        ref = sub.group(1).upper()
        pos = int(sub.group(2))
        alt = sub.group(3).upper()
        return ref, "*" if alt == "*" else alt, pos

    stop_word = re.fullmatch(r"([A-Za-z])(\d+)(stop)", text, flags=re.IGNORECASE)
    if stop_word:
        return stop_word.group(1).upper(), "*", int(stop_word.group(2))

    fs = re.fullmatch(r"([A-Za-z])(\d+)(frameshift\*?|fs.*)", text, flags=re.IGNORECASE)
    if fs:
        ref = fs.group(1).upper()
        pos = int(fs.group(2))
        return ref, f"{ref}{pos}fsX", pos

    pref_del = re.fullmatch(r"del([A-Za-z])(\d+)", text, flags=re.IGNORECASE)
    if pref_del:
        ref = pref_del.group(1).upper()
        pos = int(pref_del.group(2))
        return ref, f"{ref}{pos}del", pos

    simple_del = re.fullmatch(r"([A-Za-z])(\d+)del", text, flags=re.IGNORECASE)
    if simple_del:
        ref = simple_del.group(1).upper()
        pos = int(simple_del.group(2))
        return ref, f"{ref}{pos}del", pos

    range_del = re.fullmatch(r"([A-Za-z])?(\d+)-(\d+)del", text, flags=re.IGNORECASE)
    if range_del:
        ref = (range_del.group(1) or "").upper()
        start = int(range_del.group(2))
        end = int(range_del.group(3))
        mut = f"del{start}_{end}" if not ref else f"{ref}{start}-{end}del"
        return ref, mut, start

    alt_pref_del = re.fullmatch(r"([A-Za-z])del(\d+)", text, flags=re.IGNORECASE)
    if alt_pref_del:
        ref = alt_pref_del.group(1).upper()
        pos = int(alt_pref_del.group(2))
        return ref, f"{ref}{pos}del", pos

    # Insertion wildcard: patterns like N301+, E686+, G684+, P1112+
    # The "+" suffix indicates an insertion at this position with unknown sequence.
    # Emit INS_any as the mutation token per ResPro formatting spec.
    ins_wildcard = re.fullmatch(r"([A-Za-z])(\d+)\+", text, flags=re.IGNORECASE)
    if ins_wildcard:
        ref = ins_wildcard.group(1).upper()
        pos = int(ins_wildcard.group(2))
        return ref, "INS_any", pos

    # Insertion wildcard without reference AA: patterns like 301+
    ins_wildcard_noref = re.fullmatch(r"(\d+)\+", text, flags=re.IGNORECASE)
    if ins_wildcard_noref:
        pos = int(ins_wildcard_noref.group(1))
        return "", "INS_any", pos

    numeric_only = re.fullmatch(r"\d+", text)
    if numeric_only:
        raise ValueError("numeric_only_mutation_not_supported")

    raise ValueError("unparseable_mutation")


def parse_fold_ic50_and_phenotype(value: str) -> tuple[str, str]:
    token = value.strip()
    if not token:
        return "", ""

    lowered = token.lower()
    if lowered == "resistant":
        return "", "resistant"
    if lowered == "polymorphism":
        return "", "sensitive"

    try:
        numeric = float(token)
    except ValueError:
        return "", ""

    if numeric < 0:
        return "", ""
    return f"{numeric:g}", ""


def tsv_from_rows(rows: list[dict], columns: list[str]) -> str:
    lines = ["\t".join(columns)]
    for row in rows:
        lines.append("\t".join(norm(row.get(col, "")) for col in columns))
    return "\n".join(lines) + "\n"


def checksum(content: str) -> str:
    return "sha256:" + hashlib.sha256(content.encode("utf-8")).hexdigest()


def non_migrated_text(rows: list[dict]) -> str:
    lines = [
        "# Non-migrated HerpesDRG rows",
        "# Columns: " + "\t".join(NON_MIGRATED_COLUMNS),
        "",
        "\t".join(NON_MIGRATED_COLUMNS),
    ]
    for row in sorted(
        rows,
        key=lambda r: (norm(r.get("reason")), norm(r.get("virus")), norm(r.get("feature")), norm(r.get("aa_change"))),
    ):
        lines.append("\t".join(norm(row.get(col, "")) for col in NON_MIGRATED_COLUMNS))
    return "\n".join(lines) + "\n"


def add_non_migrated(bucket: list[dict], source_row: dict, reason: str, details: str = "") -> None:
    bucket.append(
        {
            "reason": reason,
            "mutation_id": norm(source_row.get("mutation_id")),
            "virus": norm(source_row.get("virus")),
            "feature": norm(source_row.get("gene")),
            "aa_change": norm(source_row.get("aa_change")),
            "co_feature": norm(source_row.get("co_gene")),
            "co_aa": norm(source_row.get("co_aa")),
            "status": norm(source_row.get("status")),
            "note": norm(source_row.get("note")),
            "details": details,
        }
    )


def build_publication(doi: str, link: str) -> str:
    if doi:
        return doi if doi.lower().startswith("doi:") else f"doi:{doi}"
    if link:
        # Fall back to extracting a PMID from a PubMed URL, e.g.
        # https://pubmed.ncbi.nlm.nih.gov/15230476/
        m = re.search(r"pubmed\.ncbi\.nlm\.nih\.gov/(\d+)", link)
        if m:
            return f"PMID:{m.group(1)}"
    return ""


def split_combo_mutations(aa_change: str) -> list[str]:
    return [part.strip() for part in aa_change.split(";") if part.strip()]


def make_member_id(group_id: str, feature: str, position: int, mutation: str, idx: int) -> str:
    token = re.sub(r"[^A-Za-z0-9_]+", "_", mutation).strip("_") or "mut"
    return f"{group_id}_{feature}_{position}_{token}_{idx}"


def join_unique(values: list[str]) -> str:
    seen = set()
    out = []
    for value in values:
        cleaned = norm(value)
        if not cleaned or cleaned in seen:
            continue
        seen.add(cleaned)
        out.append(cleaned)
    return ",".join(out)


def convert(source_rows: list[dict]) -> tuple[list[dict], list[dict], list[dict]]:
    rules_rows = []
    formula_rows = []
    non_migrated = []
    # Dictionary to aggregate single-mutation rules: (feature, ref_id, position, mutation, antiviral) -> aggregation
    aggregated_rules = {}
    # Dictionary to aggregate combo rules: (feature, ref_id, antiviral, parsed_mutations) -> aggregation
    aggregated_combos = {}
    group_counter = 0

    def next_group_id() -> str:
        nonlocal group_counter
        group_counter += 1
        return f"G{group_counter:05d}"

    for row in source_rows:
        virus_raw = norm(row.get("virus"))
        feature = norm(row.get("gene"))
        aa_change = norm(row.get("aa_change"))
        note = norm(row.get("note"))
        status = norm(row.get("status"))

        virus_key = normalize_virus(virus_raw)
        if not virus_key:
            add_non_migrated(non_migrated, row, "unknown_virus", "No supported virus mapping")
            continue

        if virus_key not in REFERENCE_BY_VIRUS:
            add_non_migrated(non_migrated, row, "missing_reference_mapping", "Supported virus but no fixed reference")
            continue

        exclude, reason = should_exclude(status)
        if exclude:
            add_non_migrated(non_migrated, row, reason)
            continue

        if not feature:
            add_non_migrated(non_migrated, row, "missing_feature")
            continue
        if not aa_change:
            add_non_migrated(non_migrated, row, "missing_aa_change")
            continue

        mutation_parts = split_combo_mutations(aa_change)
        parsed_mutations = []
        for part in mutation_parts:
            try:
                ref_aa, mutation_token, pos = parse_mutation(part)
            except ValueError as exc:
                add_non_migrated(non_migrated, row, str(exc), details=f"failed_part={part}")
                parsed_mutations = []
                break
            parsed_mutations.append((ref_aa, mutation_token, pos))
        if not parsed_mutations:
            continue

        publication = build_publication(norm(row.get("ref_doi")), norm(row.get("ref_link")))
        source_label = "HerpesDRG"

        has_antiviral_data = False
        emitted_for_row = 0
        for antiviral_column in ANTIVIRAL_COLUMNS:
            raw_value = norm(row.get(antiviral_column))
            if not raw_value:
                continue

            fold_ic50, phenotype = parse_fold_ic50_and_phenotype(raw_value)
            if not fold_ic50 and not phenotype:
                continue

            has_antiviral_data = True
            antiviral_norm = antiviral_column.lower()

            if len(parsed_mutations) == 1:
                ref_aa, mutation_token, pos = parsed_mutations[0]
                ref_id = REFERENCE_BY_VIRUS[virus_key]

                # Aggregation key: (feature, ref_id, position, mutation, antiviral) without publication
                aggr_key = (feature, ref_id, pos, mutation_token, antiviral_norm)
                
                if aggr_key not in aggregated_rules:
                    aggregated_rules[aggr_key] = {
                        "feature": feature,
                        "reference_identifier": ref_id,
                        "position": pos,
                        "reference": ref_aa,
                        "mutation": mutation_token,
                        "antiviral": antiviral_norm,
                        "member_id": "",
                        "source": source_label,
                        "publications": [],
                        "phenotypes": [],
                        "ic50_values": [],
                        "test_methods": [],
                        "created_dates": [],
                        "co_mutations": [],
                        "notes": [],
                    }
                
                # Aggregate data from this row
                if publication:
                    aggregated_rules[aggr_key]["publications"].append(publication)
                if phenotype:
                    aggregated_rules[aggr_key]["phenotypes"].append(phenotype)
                if fold_ic50:
                    try:
                        aggregated_rules[aggr_key]["ic50_values"].append(float(fold_ic50))
                    except ValueError:
                        pass
                test_method = norm(row.get("test_method"))
                if test_method:
                    aggregated_rules[aggr_key]["test_methods"].append(test_method)
                created_date = norm(row.get("created_date"))
                if created_date:
                    aggregated_rules[aggr_key]["created_dates"].append(created_date)
                co_gene_val = norm(row.get("co_gene"))
                co_aa_val = norm(row.get("co_aa"))
                if co_gene_val or co_aa_val:
                    co_str = f"{co_gene_val}:{co_aa_val}" if co_gene_val else co_aa_val
                    aggregated_rules[aggr_key]["co_mutations"].append(co_str)
                if note:
                    aggregated_rules[aggr_key]["notes"].append(note)
                emitted_for_row += 1
                continue

            ref_id = REFERENCE_BY_VIRUS[virus_key]
            combo_parts = tuple((ref_aa, mutation_token, pos) for ref_aa, mutation_token, pos in parsed_mutations)
            combo_key = (feature, ref_id, antiviral_norm, combo_parts)
            if combo_key not in aggregated_combos:
                aggregated_combos[combo_key] = {
                    "feature": feature,
                    "reference_identifier": ref_id,
                    "antiviral": antiviral_norm,
                    "combo_parts": combo_parts,
                    "source": source_label,
                    "publications": [],
                    "phenotypes": [],
                    "ic50_values": [],
                    "comments": [],
                    "test_methods": [],
                    "created_dates": [],
                }
            if publication:
                aggregated_combos[combo_key]["publications"].append(publication)
            if phenotype:
                aggregated_combos[combo_key]["phenotypes"].append(phenotype)
            if fold_ic50:
                try:
                    aggregated_combos[combo_key]["ic50_values"].append(float(fold_ic50))
                except ValueError:
                    pass
            if note:
                aggregated_combos[combo_key]["comments"].append(note)
            test_method = norm(row.get("test_method"))
            if test_method:
                aggregated_combos[combo_key]["test_methods"].append(test_method)
            created_date = norm(row.get("created_date"))
            if created_date:
                aggregated_combos[combo_key]["created_dates"].append(created_date)
            emitted_for_row += 1

        if not has_antiviral_data:
            add_non_migrated(non_migrated, row, "no_antiviral_signal")

    # Finalize aggregated single-mutation rules
    for aggr_key, aggr_data in aggregated_rules.items():
        # Calculate mean fold-change if present
        ic50_values = aggr_data["ic50_values"]
        if ic50_values:
            mean_ic50 = sum(ic50_values) / len(ic50_values)
            fold_ic50 = f"{mean_ic50:g}"
        else:
            fold_ic50 = ""

        # Determine phenotype: check for conflicts
        phenotypes = aggr_data["phenotypes"]
        phenotype = ""
        if phenotypes:
            unique_phenotypes = set(phenotypes)
            # Check for contradiction: resistant vs sensitive
            if "resistant" in unique_phenotypes and "sensitive" in unique_phenotypes:
                phenotype = "contradictory"
            else:
                # Use the first phenotype (or most common if needed)
                phenotype = phenotypes[0]

        # Join publications with comma and keep first-seen uniqueness.
        publication = join_unique(aggr_data["publications"])

        # Add comment about fold-change values if multiple
        comment = ""
        if len(ic50_values) > 1:
            values_str = ", ".join(f"{v:g}" for v in ic50_values)
            comment = f"Fold IC50 values: {values_str}; mean displayed"
        
        # Append test_method, created_date, and co-mutation context to comment
        test_methods = join_unique(aggr_data["test_methods"])
        created_dates = join_unique(aggr_data["created_dates"])
        co_mutations = join_unique(aggr_data.get("co_mutations", []))
        notes = join_unique(aggr_data.get("notes", []))
        if test_methods or created_dates or co_mutations or notes:
            comment_parts = [comment] if comment else []
            if test_methods:
                comment_parts.append(f"Test method: {test_methods}")
            if created_dates:
                comment_parts.append(f"Rule created at: {created_dates}")
            if co_mutations:
                comment_parts.append(f"Co-mutation context in source: {co_mutations}")
            if notes:
                comment_parts.append(notes)
            comment = ", ".join(comment_parts)

        rule = {
            "feature": aggr_data["feature"],
            "reference_identifier": aggr_data["reference_identifier"],
            "position": aggr_data["position"],
            "reference": aggr_data["reference"],
            "mutation": aggr_data["mutation"],
            "antiviral": aggr_data["antiviral"],
            "member_id": "",
            "phenotype": phenotype,
            "fold_ic50": fold_ic50,
            "publication": publication,
            "source": aggr_data["source"],
            "comment": comment,
        }
        rules_rows.append(rule)

    # Finalize aggregated combination rules (grouped atomic rows + one formula per combo key).
    for combo_data in aggregated_combos.values():
        group_id = next_group_id()
        member_ids = []
        for idx, (ref_aa, mutation_token, pos) in enumerate(combo_data["combo_parts"], start=1):
            member_id = make_member_id(group_id, combo_data["feature"], pos, mutation_token, idx)
            member_ids.append(member_id)
            rules_rows.append(
                {
                    "feature": combo_data["feature"],
                    "reference_identifier": combo_data["reference_identifier"],
                    "position": pos,
                    "reference": ref_aa,
                    "mutation": mutation_token,
                    "antiviral": "",
                    "member_id": member_id,
                    "phenotype": "",
                    "fold_ic50": "",
                    "publication": join_unique(combo_data["publications"]),
                    "source": combo_data["source"],
                    "comment": join_unique(combo_data["comments"]),
                }
            )

        phenotypes = combo_data["phenotypes"]
        phenotype = ""
        if phenotypes:
            unique_phenotypes = set(phenotypes)
            if "resistant" in unique_phenotypes and "sensitive" in unique_phenotypes:
                phenotype = "contradictory"
            else:
                phenotype = phenotypes[0]

        ic50_values = sorted(combo_data["ic50_values"])
        if ic50_values:
            mean_ic50 = sum(ic50_values) / len(ic50_values)
            fold_ic50 = f"{mean_ic50:g}"
        else:
            fold_ic50 = ""

        formula_comment = ""
        if len(ic50_values) > 1:
            values_str = ", ".join(f"{v:g}" for v in ic50_values)
            formula_comment = f"Fold IC50 values: {values_str}; mean displayed"
        
        # Append test_method and created_date to formula comment
        test_methods = join_unique(combo_data["test_methods"])
        created_dates = join_unique(combo_data["created_dates"])
        if test_methods or created_dates:
            comment_parts = [formula_comment] if formula_comment else []
            if test_methods:
                comment_parts.append(f"Test method: {test_methods}")
            if created_dates:
                comment_parts.append(f"Rule created at: {created_dates}")
            formula_comment = ", ".join(comment_parts)

        formula_rows.append(
            {
                "group_id": group_id,
                "antiviral": combo_data["antiviral"],
                "expression": "(" + " AND ".join(member_ids) + ")",
                "phenotype": phenotype,
                "fold_ic50": fold_ic50,
                "publication": join_unique(combo_data["publications"]),
                "source": combo_data["source"],
                "comment": formula_comment,
            }
        )

    rules_rows.sort(
        key=lambda r: (
            norm(r["reference_identifier"]),
            norm(r["feature"]),
            int(r["position"]),
            norm(r["mutation"]),
            norm(r["antiviral"]),
            norm(r["source"]),
        )
    )
    formula_rows.sort(key=lambda r: (norm(r["group_id"]), norm(r["antiviral"])))
    return rules_rows, formula_rows, non_migrated


def fetch_source_commit_date(source_url: str) -> str:
    """Query the GitHub API for the latest commit date of the source file.

    Returns the committer date as a YYYY-MM-DD string, or today's date as fallback.
    """
    # Parse owner/repo/path from a raw.githubusercontent.com URL
    # e.g. https://raw.githubusercontent.com/owner/repo/branch/path
    import re as _re

    match = _re.match(
        r"https://raw\.githubusercontent\.com/([^/]+)/([^/]+)/([^/]+)/(.+)",
        source_url,
    )
    if not match:
        eprint(
            f"WARNING: Cannot parse GitHub raw URL to fetch commit date: {source_url}"
        )
        return datetime.now(timezone.utc).strftime("%Y-%m-%d")

    owner, repo, branch, path = match.groups()
    api_url = (
        f"https://api.github.com/repos/{owner}/{repo}/commits"
        f"?path={path}&sha={branch}&per_page=1"
    )
    try:
        with urllib.request.urlopen(api_url) as response:
            commits = json.loads(response.read().decode("utf-8"))
        commit_date = commits[0]["commit"]["committer"]["date"][:10]
        return commit_date
    except Exception as exc:
        eprint(
            f"WARNING: Could not fetch commit date from GitHub API ({exc}); using today."
        )
        return datetime.now(timezone.utc).strftime("%Y-%m-%d")


def download_source(source_url: str) -> str:
    with urllib.request.urlopen(source_url) as response:
        return response.read().decode("utf-8")


def validate_input_header(reader: csv.DictReader) -> None:
    header = set(reader.fieldnames or [])
    missing = sorted(REQUIRED_INPUT_COLUMNS - header)
    if missing:
        raise ValueError(f"Missing required source columns: {', '.join(missing)}")


def main() -> None:
    parser = argparse.ArgumentParser(description="Convert HerpesDRG TSV to ResPro TSV artifacts")
    parser.add_argument("--source-url", default=SOURCE_URL, help="Upstream TSV URL")
    parser.add_argument(
        "--output-dir",
        default=str(Path(__file__).resolve().parent.parent / "output"),
        help="Output directory for generated artifacts",
    )
    args = parser.parse_args()

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    maintainer_update = fetch_source_commit_date(args.source_url)
    eprint(f"Source date: {maintainer_update}")
    eprint(f"Downloading {args.source_url} …")
    source_text = download_source(args.source_url)

    reader = csv.DictReader(source_text.splitlines(), delimiter="\t")
    validate_input_header(reader)
    source_rows = list(reader)
    eprint(f"Parsed {len(source_rows)} source rows")

    rules_rows, formula_rows, non_migrated_rows = convert(source_rows)

    rules_content = tsv_from_rows(rules_rows, RULES_COLUMNS)
    rules_path = out_dir / "rules.tsv"
    rules_path.write_text(rules_content, encoding="utf-8")
    eprint(f"Written {rules_path} ({len(rules_rows)} rows).")

    formula_path = out_dir / "formula-rules.tsv"
    if formula_rows:
        formula_content = tsv_from_rows(formula_rows, FORMULA_COLUMNS)
        formula_path.write_text(formula_content, encoding="utf-8")
        eprint(f"Written {formula_path} ({len(formula_rows)} rows).")
    else:
        formula_path.unlink(missing_ok=True)

    metadata = {
        "maintainers": ["Oscar Charles"],
        "contact": "Oscar.Charles@syngenta.com",
        "publication_pmid": "39192205",
        "website": "https://github.com/ojcharles/herpesdrg-db",
        "description": "Comprehensive resource for human herpesvirus antiviral drug resistance genotyping.",
        "maintainer_update": maintainer_update,
        "license": "MIT",
        "tsv_checksum": checksum(rules_content),
        "interpretation_algorithms": [
            {
                "name": "drug_groups",
                "groups": {
                    "Nucleoside analogs": [
                        "aciclovir",
                        "brivudine",
                        "ganciclovir",
                        "penciclovir",
                    ],
                    "Nucleotide analogs": [
                        "cidofovir",
                        "brincidofovir",
                    ],
                    "Methylenecyclopropane nucleoside analog": [
                        "cyclopropavir",
                    ],
                    "Pyrophosphate analog": [
                        "foscarnet",
                    ],
                    "Helicase-Primase inhibitors": [
                        "amenamevir",
                        "pritelivir",
                        "tomeglovir",
                    ],
                    "Terminase inhibitor": [
                        "letermovir",
                    ],
                    "Viral kinase inhibitor": [
                        "maribavir",
                    ],
                },
            },
            {
                "name": "drug_interpretation",
                "method": "by_phenotype",
            },
            {
                "name": "drug_interpretation",
                "method": "by_fold_ic50",
                "thresholds": {
                    "susceptible": 0.0,
                    "low-level resistance": 2.0,
                    "intermediate": 5.0,
                    "high-level resistance": 15.0,
                },
            },
            {
                "name": "drug_alias",
                "groups": {
                    "aciclovir": "ACV",
                    "amenamevir": "AMV",
                    "brincidofovir": "BCV",
                    "brivudine": "BVDU",
                    "cidofovir": "CDV",
                    "cyclopropavir": "CPV",
                    "foscarnet": "PFA",
                    "ganciclovir": "GCV",
                    "letermovir": "LTV",
                    "maribavir": "MBV",
                    "penciclovir": "PCV",
                    "pritelivir": "PTV",
                    "tomeglovir": "TGVR",
                },
            },
            {
                "name": "effect_as_resistant",
                "rules": [
                    {
                        "feature": "UL23",
                        "effect": ["frameshift", "stop_gained"],
                        "reference": "NC_001806",
                        "drug": "aciclovir",
                    },
                    {
                        "feature": "UL23",
                        "effect": ["frameshift", "stop_gained"],
                        "reference": "NC_001806",
                        "drug": "penciclovir",
                    },
                    {
                        "feature": "UL23",
                        "effect": ["frameshift", "stop_gained"],
                        "reference": "NC_001806",
                        "drug": "brivudine",
                    },
                    {
                        "feature": "UL23",
                        "effect": ["frameshift", "stop_gained"],
                        "reference": "NC_001798",
                        "drug": "aciclovir",
                    },
                    {
                        "feature": "UL23",
                        "effect": ["frameshift", "stop_gained"],
                        "reference": "NC_001798",
                        "drug": "penciclovir",
                    },
                    {
                        "feature": "UL23",
                        "effect": ["frameshift", "stop_gained"],
                        "reference": "NC_001798",
                        "drug": "brivudine",
                    },
                    {
                        "feature": "ORF36",
                        "effect": ["frameshift", "stop_gained"],
                        "reference": "NC_001348",
                        "drug": "aciclovir",
                    },
                    {
                        "feature": "ORF36",
                        "effect": ["frameshift", "stop_gained"],
                        "reference": "NC_001348",
                        "drug": "penciclovir",
                    },
                    {
                        "feature": "ORF36",
                        "effect": ["frameshift", "stop_gained"],
                        "reference": "NC_001348",
                        "drug": "brivudine",
                    },
                ],
            },
        ],
    }
    metadata_path = out_dir / "metadata.json"
    metadata_path.write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    eprint(f"Written {metadata_path}.")

    non_migrated_content = non_migrated_text(non_migrated_rows)
    non_migrated_path = out_dir / "non-migrated-rules.txt"
    non_migrated_path.write_text(non_migrated_content, encoding="utf-8")
    eprint(f"Written {non_migrated_path} ({len(non_migrated_rows)} aggregated entries).")

    eprint("Done.")


if __name__ == "__main__":
    main()
