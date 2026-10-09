#!/usr/bin/env python3
"""
Converter for HSV Drug Resistance Database (Dähne et al. 2025, Zenodo 15149867).

Downloads the Excel supplementary table, converts it to ResPro-compatible TSV files,
and writes metadata.json. The downloaded file is removed after extraction.

Usage:
    python convert.py [--output-dir <path>]
"""

import argparse
import hashlib
import json
import re
import sys
import tempfile
import urllib.request
from collections import defaultdict
from pathlib import Path


def eprint(msg: str) -> None:
    """Print a message to stderr."""
    print(msg, file=sys.stderr)

import openpyxl

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

ZENODO_URL = (
    "https://zenodo.org/records/15149867/files/"
    "Supplementary%20table%203.docx.xlsx"
)

FEATURE_MAP = {
    "Thymidine Kinase": "UL23",
    "DNA Polymerase": "UL30",
    "UL5": "UL5",
    "UL52": "UL52",
}

REFERENCE_IDS = {
    "HSV-1": "NC_001806",
    "HSV-2": "NC_001798",
}

# Sheets grouped by virus type
SHEET_VIRUS = {
    "HSV-1 TK": "HSV-1",
    "HSV-1 Pol": "HSV-1",
    "HSV-1 HPC": "HSV-1",
    "HSV-2 TK": "HSV-2",
    "HSV-2 Pol": "HSV-2",
    "HSV-2 HPC": "HSV-2",
}

RULES_COLUMNS = [
    "feature",
    "reference_identifier",
    "position",
    "reference",
    "mutation",
    "antiviral",
    "member_id",
    "phenotype",
    "clinical_phenotype",
    "publication",
    "source",
    "comment",
]

FORMULA_COLUMNS = [
    "group_id",
    "antiviral",
    "expression",
    "phenotype",
    "clinical_phenotype",
    "publication",
    "source",
]

NON_MIGRATED_COLUMNS = [
    "reason",
    "sheet",
    "feature",
    "drug",
    "aa_change",
    "aa_position",
    "resistance",
    "pmids",
    "details",
]


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def norm_str(v) -> str:
    if v is None:
        return ""
    return str(v).strip()


def count_pmids(pmid_str: str) -> int:
    if not pmid_str:
        return 0
    return len([p for p in re.split(r"[\s,;]+", pmid_str) if p.strip()])


def parse_pmids(pmid_str: str) -> str:
    """Normalise PMIDs to a comma-separated list with PMID: prefix."""
    if not pmid_str:
        return ""
    # Source rows may separate PMIDs with commas, semicolons, or whitespace.
    raw_parts = [p.strip() for p in re.split(r"[\s,;]+", pmid_str) if p.strip()]
    norm_parts = []
    for part in raw_parts:
        cleaned = re.sub(r"^PMID:\s*", "", part, flags=re.IGNORECASE)
        norm_parts.append(f"PMID:{cleaned}")
    return ",".join(norm_parts)


def source_website_from_url(source_url: str) -> str:
    match = re.search(r"/records/(\d+)/", source_url)
    if match:
        return f"https://zenodo.org/records/{match.group(1)}"
    return "https://zenodo.org/records/15149867"


def zenodo_record_id_from_url(source_url: str) -> str:
    match = re.search(r"/records/(\d+)/", source_url)
    if not match:
        raise ValueError(f"Could not determine Zenodo record id from source URL: {source_url}")
    return match.group(1)


def zenodo_record_metadata(record_id: str) -> dict:
    api_url = f"https://zenodo.org/api/records/{record_id}"
    with urllib.request.urlopen(api_url) as response:
        return json.load(response)


def zenodo_record_update_date(record_metadata: dict) -> str:
    for key in ("updated", "created"):
        value = str(record_metadata.get(key, "")).strip()
        if value:
            return value.split("T", 1)[0]
    raise ValueError("Zenodo record metadata did not contain created/updated timestamp")


def parse_aa_change(aa_change: str):
    """
    Return (reference_aa, position, mutation_token) from an AA change string.

    Supported forms:
      A336V          — substitution
      A336*          — stop
      N23Stop        — stop with ref AA
      8Stop          — stop without ref AA
      Y101fsx        — frameshift with ref AA
      144fsx         — frameshift without ref AA
      I194Del        — single-residue deletion
      1-248Del       — range deletion (no leading AAs)
      DD676-677Del   — multi-residue range deletion
      V813M*         — substitution with trailing asterisk annotation (strip *)

    Returns a list of (ref_aa, position_int, mutation_token_str) tuples.
    Returns an empty list when the mutation cannot be represented.
    Dual-allele patterns like "E43A/D" are expanded into two separate entries.
    """
    aa_change = aa_change.strip().replace("\n", "")

    # Dual allele like "E43A/D" — split into two separate substitutions
    dual_allele = re.match(r"^([A-Z])(\d+)([A-Z])/([A-Z])$", aa_change, re.IGNORECASE)
    if dual_allele:
        ref_aa = dual_allele.group(1).upper()
        pos = int(dual_allele.group(2))
        alt1 = dual_allele.group(3).upper()
        alt2 = dual_allele.group(4).upper()
        return [(ref_aa, pos, alt1), (ref_aa, pos, alt2)]

    # Strip spurious trailing asterisk that is not itself a stop notation
    # e.g. V813M* where * is an annotation marker
    trailing_star = re.match(r"^([A-Z])(\d+)([A-Z])\*$", aa_change, re.IGNORECASE)
    if trailing_star:
        ref_aa = trailing_star.group(1).upper()
        pos = int(trailing_star.group(2))
        alt_aa = trailing_star.group(3).upper()
        return [(ref_aa, pos, alt_aa)]

    # frameshift with ref AA: e.g. Y101fsx, K715fsATFF*, Y101fs
    fs_match = re.match(r"^([A-Z])(\d+)(fs.*)$", aa_change, re.IGNORECASE)
    if fs_match:
        ref_aa = fs_match.group(1).upper()
        pos = int(fs_match.group(2))
        return [(ref_aa, pos, f"{ref_aa}fsX")]

    # frameshift without ref AA: e.g. 144fsx, 145fs
    fs_nore = re.match(r"^(\d+)(fs.*)$", aa_change, re.IGNORECASE)
    if fs_nore:
        pos = int(fs_nore.group(1))
        return [(None, pos, "fsX")]

    # standard substitution / stop: A336V, A336*
    sub_match = re.match(r"^([A-Z])(\d+)([A-Z\*])$", aa_change, re.IGNORECASE)
    if sub_match:
        ref_aa = sub_match.group(1).upper()
        pos = int(sub_match.group(2))
        alt_aa = sub_match.group(3).upper()
        if alt_aa == "*":
            return [(ref_aa, pos, "*")]
        return [(ref_aa, pos, alt_aa)]

    # stop with ref AA: N23Stop, M182Stop
    stop_re = re.match(r"^([A-Z])(\d+)[Ss]top$", aa_change, re.IGNORECASE)
    if stop_re:
        ref_aa = stop_re.group(1).upper()
        pos = int(stop_re.group(2))
        return [(ref_aa, pos, "*")]

    # stop without ref AA: 8Stop, 44Stop
    stop_nore = re.match(r"^(\d+)[Ss]top$", aa_change, re.IGNORECASE)
    if stop_nore:
        pos = int(stop_nore.group(1))
        return [(None, pos, "*")]

    # single-residue deletion: e.g. "I194Del"
    single_del = re.match(r"^([A-Z])(\d+)[Dd]el$", aa_change, re.IGNORECASE)
    if single_del:
        ref_aa = single_del.group(1).upper()
        pos = int(single_del.group(2))
        return [(ref_aa, pos, f"{ref_aa}{pos}del")]

    # Single-position insertion without range: e.g. "A301Ins", "684Ins"
    # Cannot determine inserted sequence — emit INS_any wildcard per ResPro formatting spec.
    ins_wildcard = re.match(r"^([A-Z])(\d+)[Ii]ns$", aa_change, re.IGNORECASE)
    if ins_wildcard:
        ref_aa = ins_wildcard.group(1).upper()
        pos = int(ins_wildcard.group(2))
        return [(ref_aa, pos, "INS_any")]

    ins_wildcard_noref = re.match(r"^(\d+)[Ii]ns$", aa_change, re.IGNORECASE)
    if ins_wildcard_noref:
        pos = int(ins_wildcard_noref.group(1))
        return [(None, pos, "INS_any")]

    # Multi-residue range deletion/insertion with leading AAs:
    # e.g. DD676-677Del, PGDEPA1106-1111Del, ED684-685Ins
    multi_indel = re.match(
        r"^([A-Z]+)(\d+)-(\d+)(Del|Ins)$", aa_change, re.IGNORECASE
    )
    if multi_indel:
        ref_aas = multi_indel.group(1).upper()
        start = int(multi_indel.group(2))
        kind = multi_indel.group(4).lower()
        ref_aa = ref_aas[0]
        if kind == "del":
            # ResPro deletion notation: reference = full deleted span (anchor + deleted residues),
            # mutation = anchor only (first residue survives, the rest are deleted).
            # e.g. DD676-677Del → reference="DD", mutation="D", position=676
            return [(ref_aas, start, ref_aa)]
        else:
            # ResPro insertion notation: reference = anchor residue,
            # mutation = anchor + inserted payload (e.g. "ED" means E is anchor, D is inserted).
            return [(ref_aa, start, ref_aas)]

    # bare range deletion without leading AAs: e.g. 1-248Del
    bare_del = re.match(r"^(\d+)-(\d+)[Dd]el$", aa_change)
    if bare_del:
        start = int(bare_del.group(1))
        return [(None, start, f"del{start}_{bare_del.group(2)}")]

    # Insertion wildcard with "+" suffix: e.g. N301+, E686+
    # The "+" suffix indicates an insertion at this position with unknown sequence.
    # Emit INS_any as the mutation token per ResPro formatting spec.
    ins_plus = re.match(r"^([A-Z])(\d+)\+$", aa_change, re.IGNORECASE)
    if ins_plus:
        ref_aa = ins_plus.group(1).upper()
        pos = int(ins_plus.group(2))
        return [(ref_aa, pos, "INS_any")]

    ins_plus_noref = re.match(r"^(\d+)\+$", aa_change, re.IGNORECASE)
    if ins_plus_noref:
        pos = int(ins_plus_noref.group(1))
        return [(None, pos, "INS_any")]

    # unrecognised
    return []


def determine_phenotypes(resistance: str, cell_culture: str, clinical: str):
    """
    Returns (phenotype, clinical_phenotype) according to mapping rules.
    """

    def define_phenotype(resistance, cc_cl):
        if resistance == "sensitive":
            if cc_cl == "yes":
                return "sensitive"
            elif cc_cl == "contradiction":
                return "resistant"
            else:
                return "unknown"
        if resistance == "resistant":
            if cc_cl == "yes":
                return "resistant"
            elif cc_cl == "contradiction":
                return "sensitive"
            else:
                return "unknown"
        if resistance == "sensitive/resistant":
            if cc_cl == "yes":
                return "contradictory"
            elif cc_cl == "contradiction":
                return "contradictory"
            else:
                return "unknown"
        
        return "unknown"

    r = resistance.strip().lower()
    cc = cell_culture.strip().lower()
    cl = clinical.strip().lower()

    return define_phenotype(r, cc), define_phenotype(r, cl)


def tsv_checksum(content: str) -> str:
    return "sha256:" + hashlib.sha256(content.encode("utf-8")).hexdigest()


def add_non_migrated(
    non_migrated_rows,
    reason: str,
    sheet: str,
    feature: str,
    drug: str,
    aa_change: str,
    aa_position,
    resistance: str,
    pmids: str,
    details: str = "",
) -> None:
    non_migrated_rows.append(
        {
            "reason": reason,
            "sheet": sheet,
            "feature": feature,
            "drug": drug,
            "aa_change": aa_change,
            "aa_position": "" if aa_position is None else str(aa_position),
            "resistance": resistance,
            "pmids": pmids,
            "details": details,
        }
    )


# ---------------------------------------------------------------------------
# Core extraction
# ---------------------------------------------------------------------------


def extract_rows(wb):
    """
    Extract all relevant rows from the workbook.

    Returns:
        single_rows: list of dicts for rules.tsv (non-combo rows + combo member rows)
        formula_rows: list of dicts for formula-rules.tsv
    """
    single_rows = []
    formula_rows = []
    non_migrated_rows = []
    group_counter = [0]
    # Registry of already-emitted combo members: (feature, ref_id, pos, mut_token) -> member_id.
    # Shared members across groups reuse the same row and member_id.
    combo_member_registry: dict = {}

    def next_group_id():
        group_counter[0] += 1
        return f"G{group_counter[0]:04d}"

    for sheet_name, virus in SHEET_VIRUS.items():
        if sheet_name not in wb.sheetnames:
            eprint(f"WARNING: Sheet '{sheet_name}' not found, skipping.")
            continue

        ws = wb[sheet_name]
        ref_id = REFERENCE_IDS[virus]

        for row in ws.iter_rows(min_row=2, values_only=True):
            feature_raw = norm_str(row[0])
            if not feature_raw:
                continue

            drug = norm_str(row[1]).strip()
            aa_change = norm_str(row[2])
            aa_pos_raw = row[3]
            resistance = norm_str(row[4])
            pmids_raw = norm_str(row[5])
            cell_culture = norm_str(row[6])
            clinical = norm_str(row[7])
            comment = norm_str(row[10]) if len(row) > 10 else ""

            pmids = parse_pmids(pmids_raw)

            # Skip rows where AA change is absent (nucleotide-only entries)
            if not aa_change:
                eprint(
                    f"SKIPPED nucleotide-only row: sheet={sheet_name} "
                    f"pos_col={aa_pos_raw!r} drug={drug} "
                    f"(no amino-acid change recorded)")
                add_non_migrated(
                    non_migrated_rows,
                    reason="nucleotide_only_row",
                    sheet=sheet_name,
                    feature=feature_raw,
                    drug=drug,
                    aa_change=aa_change,
                    aa_position=aa_pos_raw,
                    resistance=resistance,
                    pmids=pmids,
                    details="No amino-acid change value provided in source row.",
                )
                continue

            # Slash-separated mutations define a combo rule.
            # If features are slash-separated too, members are matched by index.
            # If only one feature is provided, it applies to all mutation members.
            parts_aa = [part.strip() for part in aa_change.split("/") if part.strip()]
            is_dual_allele = bool(re.match(r"^[A-Z]\d+[A-Z]/[A-Z]$", aa_change, re.IGNORECASE))
            is_combo = len(parts_aa) > 1 and not is_dual_allele
            if is_combo:
                raw_parts_feature = [part.strip() for part in feature_raw.split("/") if part.strip()]
                if len(raw_parts_feature) == 1:
                    parts_feature = raw_parts_feature * len(parts_aa)
                else:
                    parts_feature = raw_parts_feature

                # Position column may also be slash-split ("356/222").
                if aa_pos_raw and "/" in str(aa_pos_raw):
                    parts_pos = [part.strip() for part in str(aa_pos_raw).split("/") if part.strip()]
                else:
                    parts_pos = [str(aa_pos_raw)] * len(parts_aa)

                if len(parts_feature) != len(parts_aa) or len(parts_pos) != len(parts_aa):
                    eprint(
                        f"WARNING: Cannot parse combo row feature='{feature_raw}' "
                        f"aa='{aa_change}' in sheet '{sheet_name}' — skipping.")
                    add_non_migrated(
                        non_migrated_rows,
                        reason="invalid_combo_format",
                        sheet=sheet_name,
                        feature=feature_raw,
                        drug=drug,
                        aa_change=aa_change,
                        aa_position=aa_pos_raw,
                        resistance=resistance,
                        pmids=pmids,
                        details=(
                            "Expected the same number of mutation and position members, and "
                            "either one feature or one feature per mutation member."
                        ),
                    )
                    continue

                gid = next_group_id()
                member_ids = []
                combo_skip = False

                for sub_feature, sub_aa, sub_pos_raw in zip(parts_feature, parts_aa, parts_pos):
                    mapped_feature = FEATURE_MAP.get(sub_feature, sub_feature)
                    parsed_mutations = parse_aa_change(sub_aa)
                    if not parsed_mutations:
                        eprint(
                            f"WARNING: Cannot parse mutation '{sub_aa}' in combo row "
                            f"(sheet={sheet_name}, feature={sub_feature}) — skipping entire combo.")
                        add_non_migrated(
                            non_migrated_rows,
                            reason="combo_member_unparseable_or_missing_ref",
                            sheet=sheet_name,
                            feature=sub_feature,
                            drug=drug,
                            aa_change=sub_aa,
                            aa_position=sub_pos_raw,
                            resistance=resistance,
                            pmids=pmids,
                            details="Combo row dropped because one member could not be represented.",
                        )
                        combo_skip = True
                        break

                    for ref_aa, pos, mut_token in parsed_mutations:
                        if ref_aa is None:
                            eprint(
                                f"WARNING: Missing ref AA for mutation '{sub_aa}' in combo row "
                                f"(sheet={sheet_name}, feature={sub_feature}) — skipping entire combo.")
                            add_non_migrated(
                                non_migrated_rows,
                                reason="combo_member_unparseable_or_missing_ref",
                                sheet=sheet_name,
                                feature=sub_feature,
                                drug=drug,
                                aa_change=sub_aa,
                                aa_position=sub_pos_raw,
                                resistance=resistance,
                                pmids=pmids,
                                details="Combo row dropped because one member could not be represented.",
                            )
                            combo_skip = True
                            break
                        if pos is None:
                            try:
                                pos = int(sub_pos_raw)
                            except (ValueError, TypeError):
                                pos = 0

                        member_key = (mapped_feature, ref_id, pos, mut_token)
                        if member_key in combo_member_registry:
                            mid = combo_member_registry[member_key]
                        else:
                            mid = f"{mapped_feature}_{pos}_{mut_token}"
                            combo_member_registry[member_key] = mid
                            # Emit member row once. No antiviral or group — the member
                            # is shared across groups; drug context lives in formula-rules.tsv.
                            single_rows.append({
                                "feature": mapped_feature,
                                "reference_identifier": ref_id,
                                "position": pos,
                                "reference": ref_aa or "",
                                "mutation": mut_token,
                                "antiviral": "",
                                "member_id": mid,
                                "phenotype": "",
                                "clinical_phenotype": "",
                                "publication": pmids,
                                "source": "Dähne et al. 2025 (Zenodo 15149867)",
                                "comment": comment,
                                "sheet": sheet_name,
                            })
                        member_ids.append(mid)

                    if combo_skip:
                        break

                if combo_skip or not member_ids:
                    # One member could not be parsed (logged above); skip this formula group.
                    # Already-registered member rows are kept — they are valid atomic mutations.
                    continue

                # formula row
                phenotype, clinical_phenotype = determine_phenotypes(
                    resistance, cell_culture, clinical
                )
                expr = "(" + " AND ".join(member_ids) + ")"
                formula_rows.append({
                    "group_id": gid,
                    "antiviral": drug.lower(),
                    "expression": expr,
                    "phenotype": phenotype,
                    "clinical_phenotype": clinical_phenotype,
                    "publication": pmids,
                    "source": "Dähne et al. 2025 (Zenodo 15149867)",
                })

            else:
                # Regular single-feature row
                mapped_feature = FEATURE_MAP.get(feature_raw)
                if mapped_feature is None:
                    eprint(
                        f"WARNING: Unknown feature '{feature_raw}' in sheet '{sheet_name}' — skipping.")
                    add_non_migrated(
                        non_migrated_rows,
                        reason="unknown_feature",
                        sheet=sheet_name,
                        feature=feature_raw,
                        drug=drug,
                        aa_change=aa_change,
                        aa_position=aa_pos_raw,
                        resistance=resistance,
                        pmids=pmids,
                        details="Feature not mappable to ResPro feature naming.",
                    )
                    continue

                parsed_mutations = parse_aa_change(aa_change)
                if not parsed_mutations:
                    eprint(
                        f"SKIPPED unparseable mutation: sheet={sheet_name} "
                        f"feature={mapped_feature} aa_change={aa_change!r} drug={drug}")
                    add_non_migrated(
                        non_migrated_rows,
                        reason="unparseable_mutation",
                        sheet=sheet_name,
                        feature=mapped_feature,
                        drug=drug,
                        aa_change=aa_change,
                        aa_position=aa_pos_raw,
                        resistance=resistance,
                        pmids=pmids,
                        details="Mutation syntax could not be converted to supported ResPro notation.",
                    )
                    continue

                # Process each parsed mutation (dual alleles like E43A/D produce two entries)
                for ref_aa, pos, mut_token in parsed_mutations:
                    if ref_aa is None:
                        eprint(
                            f"SKIPPED mutation with no reference AA: sheet={sheet_name} "
                            f"feature={mapped_feature} aa_change={aa_change!r} drug={drug} "
                            f"(reference amino acid required by ResPro)")
                        add_non_migrated(
                            non_migrated_rows,
                            reason="missing_reference_amino_acid",
                            sheet=sheet_name,
                            feature=mapped_feature,
                            drug=drug,
                            aa_change=aa_change,
                            aa_position=aa_pos_raw,
                            resistance=resistance,
                            pmids=pmids,
                            details="Reference amino acid is required by rules.tsv schema.",
                        )
                        continue
                    if pos is None:
                        try:
                            pos = int(aa_pos_raw)
                        except (ValueError, TypeError):
                            eprint(
                                f"SKIPPED row with unresolvable position: sheet={sheet_name} "
                                f"feature={mapped_feature} aa_change={aa_change!r} pos_col={aa_pos_raw!r}")
                            add_non_migrated(
                                non_migrated_rows,
                                reason="unresolvable_position",
                                sheet=sheet_name,
                                feature=mapped_feature,
                                drug=drug,
                                aa_change=aa_change,
                                aa_position=aa_pos_raw,
                                resistance=resistance,
                                pmids=pmids,
                                details="Could not resolve amino-acid position to integer.",
                            )
                            continue

                    phenotype, clinical_phenotype = determine_phenotypes(
                        resistance, cell_culture, clinical
                    )

                    single_rows.append({
                        "feature": mapped_feature,
                        "reference_identifier": ref_id,
                        "position": pos,
                        "reference": ref_aa or "",
                        "mutation": mut_token,
                        "antiviral": drug.lower(),
                        "member_id": "",
                        "phenotype": phenotype,
                        "clinical_phenotype": clinical_phenotype,
                        "publication": pmids,
                        "source": "Dähne et al. 2025 (Zenodo 15149867)",
                        "comment": comment,
                        "sheet": sheet_name,
                    })

    return single_rows, formula_rows, non_migrated_rows


# ---------------------------------------------------------------------------
# Deduplication
# ---------------------------------------------------------------------------


def deduplicate_single_rows(rows, non_migrated_rows):
    """
    For rows with identical (feature, reference_identifier, position, mutation, antiviral)
    keep the one with more PMIDs. Print dropped rows.
    """
    # Only deduplicate non-combo rows (member_id == "")
    combo = [r for r in rows if r["member_id"]]
    non_combo = [r for r in rows if not r["member_id"]]

    key_map = defaultdict(list)
    for r in non_combo:
        key = (r["feature"], r["reference_identifier"], r["position"], r["mutation"], r["antiviral"])
        key_map[key].append(r)

    kept = []
    for key, group in key_map.items():
        if len(group) == 1:
            kept.append(group[0])
        else:
            best = max(group, key=lambda r: count_pmids(r["publication"]))
            for r in group:
                if r is not best:
                    eprint(
                        f"DROPPED duplicate row: feature={r['feature']} pos={r['position']} "
                        f"mutation={r['mutation']} drug={r['antiviral']} "
                        f"pmids={r['publication']!r} "
                        f"(kept pmids={best['publication']!r}, "
                        f"reason: fewer PMIDs than retained row)")
                    add_non_migrated(
                        non_migrated_rows,
                        reason="duplicate_rule_dropped_fewer_pmids",
                        sheet=r.get("sheet", ""),
                        feature=r["feature"],
                        drug=r["antiviral"],
                        aa_change=f"{r['reference']}{r['position']}{r['mutation']}",
                        aa_position=r["position"],
                        resistance="",
                        pmids=r["publication"],
                        details=(
                            f"Kept publication set: {best['publication']}; "
                            "this row had fewer PMID entries."
                        ),
                    )
            kept.append(best)

    return kept + combo


# ---------------------------------------------------------------------------
# Serialisation
# ---------------------------------------------------------------------------


def rows_to_tsv(rows, columns) -> str:
    lines = ["\t".join(columns)]
    for r in rows:
        lines.append("\t".join(str(r.get(c, "")) for c in columns))
    return "\n".join(lines) + "\n"


def sort_single_rows(rows):
    return sorted(
        rows,
        key=lambda r: (
            r["reference_identifier"],
            r["feature"],
            r["antiviral"],
            int(r["position"]) if str(r["position"]).isdigit() else 0,
            r["mutation"],
        ),
    )


def sort_formula_rows(rows):
    return sorted(rows, key=lambda r: (r["group_id"], r["antiviral"]))


def non_migrated_rows_to_text(rows) -> str:
    lines = [
        "# Non-migrated rules from source workbook",
        "# Columns: " + "\t".join(NON_MIGRATED_COLUMNS),
        "",
        "\t".join(NON_MIGRATED_COLUMNS),
    ]
    sorted_rows = sorted(
        rows,
        key=lambda r: (
            r["reason"],
            r["sheet"],
            r["feature"],
            r["drug"],
            r["aa_change"],
        ),
    )
    for row in sorted_rows:
        lines.append("\t".join(str(row.get(col, "")) for col in NON_MIGRATED_COLUMNS))
    return "\n".join(lines) + "\n"


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main():
    parser = argparse.ArgumentParser(description="Convert HSV resistance DB from Zenodo.")
    parser.add_argument(
        "--source-url",
        default=ZENODO_URL,
        help="Direct source file URL (defaults to Zenodo record v3 XLSX).",
    )
    parser.add_argument(
        "--output-dir",
        default=str(Path(__file__).resolve().parent.parent / "output"),
        help="Directory to write output files into.",
    )
    args = parser.parse_args()

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    record_id = zenodo_record_id_from_url(args.source_url)
    record_metadata = zenodo_record_metadata(record_id)
    maintainer_update = zenodo_record_update_date(record_metadata)

    # --- Download ---
    with tempfile.NamedTemporaryFile(suffix=".xlsx", delete=False) as tmp_file:
        tmp_path = Path(tmp_file.name)
    eprint(f"Source date: {maintainer_update}")
    eprint(f"Downloading {args.source_url} …")
    urllib.request.urlretrieve(args.source_url, tmp_path)

    try:
        wb = openpyxl.load_workbook(tmp_path)
        single_rows, formula_rows, non_migrated_rows = extract_rows(wb)
    finally:
        tmp_path.unlink(missing_ok=True)
    # Source is a non-tabular XLSX workbook, so the "Parsed <N> source rows"
    # line is omitted per the converter output contract (README).

    # --- Deduplication ---
    single_rows = deduplicate_single_rows(single_rows, non_migrated_rows)

    # --- Sort ---
    single_rows = sort_single_rows(single_rows)
    formula_rows = sort_formula_rows(formula_rows)

    # --- Write rules.tsv ---
    rules_content = rows_to_tsv(single_rows, RULES_COLUMNS)
    rules_path = out_dir / "rules.tsv"
    rules_path.write_text(rules_content, encoding="utf-8")
    eprint(f"Written {rules_path} ({len(single_rows)} rows).")

    # --- Write formula-rules.tsv (only if combo rows exist) ---
    formula_path = out_dir / "formula-rules.tsv"
    if formula_rows:
        formula_content = rows_to_tsv(formula_rows, FORMULA_COLUMNS)
        formula_path.write_text(formula_content, encoding="utf-8")
        eprint(f"Written {formula_path} ({len(formula_rows)} rows).")
    else:
        formula_path.unlink(missing_ok=True)
        formula_content = ""

    # --- Write metadata.json ---
    checksum = tsv_checksum(rules_content)
    metadata = {
        "maintainers": ["Daehne, Theo", "Gosert, Rainer", "Jaki, Lena"],
        "contact": "marcus.panning@uniklinik-freiburg.de",
        "publication_pmid": "40349973",
        "website": source_website_from_url(args.source_url),
        "description": (
            "HSV-1 and HSV-2 drug resistance mutation database curated from published literature. "
            "Covers TK (UL23), DNA Polymerase (UL30), and Helicase-Primase Complex (UL5/UL52) features. "
        ),
        "maintainer_update": maintainer_update,
        "license": "CC-BY-4.0",
        "tsv_checksum": checksum,
        "interpretation_algorithms": [
            {
                "name": "drug_groups",
                "groups": {
                    "Nucleoside analogs": [
                        "aciclovir",
                        "brivudine",
                        "penciclovir",
                    ],
                    "Nucleotide analogs": [
                        "cidofovir",
                    ],
                    "Pyrophosphate analog": [
                        "foscarnet",
                    ],
                    "Helicase-Primase inhibitors": [
                        "amenamevir",
                        "pritelivir",
                    ],
                },
            },
            {
                "name": "drug_interpretation",
                "method": "by_phenotype",
            },
            {
                "name": "drug_alias",
                "groups": {
                    "aciclovir": "ACV",
                    "amenamevir": "AMV",
                    "brivudine": "BVDU",
                    "cidofovir": "CDV",
                    "foscarnet": "PFA",
                    "penciclovir": "PCV",
                    "pritelivir": "PTV",
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
                ],
            },
        ],
    }
    metadata_path = out_dir / "metadata.json"
    metadata_path.write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    eprint(f"Written {metadata_path}.")

    # --- Write non-migrated-rules.txt ---
    non_migrated_path = out_dir / "non-migrated-rules.txt"
    non_migrated_content = non_migrated_rows_to_text(non_migrated_rows)
    non_migrated_path.write_text(non_migrated_content, encoding="utf-8")
    eprint(f"Written {non_migrated_path} ({len(non_migrated_rows)} aggregated entries).")

    eprint("Done.")


if __name__ == "__main__":
    main()
