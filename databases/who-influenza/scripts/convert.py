#!/usr/bin/env python3
"""Convert the three WHO influenza marker-table PDFs into ResPro artifacts.

Sources (downloaded by the workflow into --source-dir):
- human-nai-marker-table.pdf   WHO human NAI marker table (NA feature)
- avian-nai-marker-table.pdf   WHO avian NAI marker table (NA feature)
- pa-marker-table.pdf          WHO PA/baloxavir marker table (PA feature)

Pipeline:
1. Parse the PDFs with who_parsers.py (no network).
2. Normalize substitutions, phenotypes, and fold-changes.
3. Validate every position against the pinned reference CDS translations in
   scripts/references/ (fail-closed: mismatches go to non-migrated-rules.txt).
4. Deduplicate and merge rules across tables by
   (feature, reference_identifier, position, reference, mutation, antiviral):
   publications join with ",", fold_ic50 keeps the max endpoint, the more
   severe phenotype wins, disagreements and all distinct fold values are kept
   in the comment.
5. Emit combination rules as formula groups (AND over member rows; mixed
   alleles become OR sub-expressions).
6. Write rules.tsv, formula-rules.tsv, metadata.json, non-migrated-rules.txt
   and source-state.json (per-PDF provenance for the autobump workflow).

The converter never performs network I/O; citation resolution is consumed
from the reference-lookup.tsv sidecar produced by resolve_references.py.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
from pathlib import Path


# pdfplumber is imported lazily via who_parsers so the pure conversion logic
# is usable without the dependency (e.g. in unit tests).

from who_parsers import MarkerRow, parse_avian_nai, parse_human_nai, parse_pa

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
    "table",
    "page",
    "subtype",
    "substitution",
    "antiviral",
    "details",
]

# Pinned upstream sources (SHA-256 recorded in source-state.json at build time).
SOURCE_URLS = {
    "human-nai": (
        "https://cdn.who.int/media/docs/default-source/influenza/laboratory---network/"
        "quality-assurance/human-nai-marker-table.pdf"
    ),
    "avian-nai": (
        "https://cdn.who.int/media/docs/default-source/influenza/avwg/"
        "avian-nai-marker-who-table.pdf"
    ),
    "pa": (
        "https://cdn.who.int/media/docs/default-source/influenza/laboratory---network/"
        "quality-assurance/antiviral-susceptibility-influenza/"
        "pa-marker-who-table_28-11-2025_updated.pdf"
    ),
}
SOURCE_FILES = {
    "human-nai": "human-nai-marker-table.pdf",
    "avian-nai": "avian-nai-marker-table.pdf",
    "pa": "pa-marker-table.pdf",
}
SOURCE_LABELS = {
    "human-nai": "WHO human NAI",
    "avian-nai": "WHO avian NAI",
    "pa": "WHO PA",
}

# User-provided fixed subtype -> reference accession mappings.
HUMAN_NA_REFERENCE = {
    "A(H1N1)pdm09": "NC_026434",
    "A(H1N1)": "NC_026434",
    "A(H1N1)v": "NC_026434",
    "A(H5N1)": "EF619973",
    "A(H3N2)": "NC_007368",
    "A(H3N2)v": "NC_007368",
    "A(H7N9)": "NC_026429",
    "Type B": "CY115257",
}
AVIAN_NA_REFERENCE = {
    "H5N1": "EF619973",
    "N2": "NC_007368",
    "N3": "AY207514",
    "N4": "MW287218.1",
    "N5": "CY098534.1",
    "N6": "MW287220.1",
    "N7": "MW325955.1",
    "N8": "KX297920.1",
    "N9": "AB292780",
    "H7N9": "NC_026429",
}
PA_REFERENCE_A = "NC_026437"
PA_REFERENCE_B = "CY115260"

# Subtypes whose rows may carry footnote f (full-length N1 numbering).
H5N1_SUBTYPES = {"A(H5N1)", "H5N1"}
# EF619973 (A/turkey/Turkey/1/2005) lacks the 20-aa stalk deletion.
H5N1_STALK_OFFSET = 20

# Phenotype severity for merge conflicts (higher wins).
PHENOTYPE_SEVERITY = {
    "normal inhibition": 1,
    "normal": 1,
    "susceptible": 1,
    "suspected reduced": 2,
    "reduced susceptibility": 3,
    "reduced inhibition": 4,
    "highly reduced inhibition": 5,
}

NAI_DRUGS = {
    "oseltamivir": "oseltamivir",
    "zanamivir": "zanamivir",
    "peramivir": "peramivir",
    "laninamivir": "laninamivir",
}
PA_DRUG = "baloxavir"

# NAI cell phenotype token -> ResPro phenotype, expanded per the WHO legend
# ("NI, normal inhibition; RI, reduced inhibition; HRI, highly reduced
# inhibition"). All three labels are in ResPro's phenotype vocabulary with
# the same severity ranks as the abbreviations.
NAI_TOKEN_MAP = {
    "ni": "normal inhibition",
    "ri": "reduced inhibition",
    "hri": "highly reduced inhibition",
}
# PA classification -> ResPro phenotype.
PA_CLASS_MAP = {
    "normal": "normal",
    "suspected reduced": "suspected reduced",
    "reduced": "reduced susceptibility",
}

# Footnote letters that may suffix substitution tokens (a-k per WHO footnotes).
FOOTNOTE_SUFFIX_RE = re.compile(r"^(?P<base>.+?)(?P<fn>[a-z]+)$")
DELETION_RE = re.compile(r"^del\s*(\d+)\s*[-–−]\s*(\d+)$", re.IGNORECASE)
SUBSTITUTION_RE = re.compile(r"^([A-Z])(\d+)([A-Z](?:/[A-Z])*)$")
CITATION_LIST_RE = re.compile(r"\((\d+(?:\s*[,;]\s*\d+)*(?:\s*[-–−]\s*\d+)?)\)")
FOLD_RANGE_SPLIT_RE = re.compile(r"[-–−]")


def eprint(msg: str) -> None:
    print(msg, file=sys.stderr)


def norm(v: object) -> str:
    if v is None:
        return ""
    return str(v).strip()


def join_unique(values: list[str]) -> str:
    """Join values with ',' preserving first-occurrence order, deduplicating.

    Elements may themselves be comma-joined strings (e.g. merged publication
    lists), so every element is split before deduplication.
    """
    seen: set[str] = set()
    out: list[str] = []
    for value in values:
        for part in norm(value).split(","):
            part = part.strip()
            if not part or part in seen:
                continue
            seen.add(part)
            out.append(part)
    return ",".join(out)


def load_reference_sequences(ref_dir: Path) -> dict[str, str]:
    """Load pinned reference CDS translations keyed by accession."""
    sequences: dict[str, str] = {}
    for fasta in sorted(ref_dir.glob("*.fasta")):
        parts: list[str] = []
        with fasta.open(encoding="utf-8") as fh:
            for line in fh:
                if line.startswith(">"):
                    continue
                parts.append(line.strip())
        sequences[fasta.stem] = "".join(parts)
    return sequences


def strip_footnote(token: str) -> tuple[str, str]:
    """Split a trailing lowercase footnote suffix off a substitution token."""
    m = FOOTNOTE_SUFFIX_RE.match(token)
    if m and m.group("base"):
        return m.group("base"), m.group("fn")
    return token, ""


def parse_substitution_cell(cell: str) -> tuple[list[dict], str]:
    """Parse a substitution cell into components.

    Returns (components, footnote_letters). Components are dicts:
    {"kind": "sub"|"mixed"|"del", "ref": aa, "pos": int,
     "alts": [aa, ...] (sub/mixed), "start"/"end" (del)}.
    Raises ValueError on unparseable syntax (fail-closed).
    """
    text = " ".join(cell.split())
    components: list[dict] = []
    footnotes = ""
    for part in text.split("+"):
        part = part.strip()
        if not part:
            continue
        base, fn = strip_footnote(part)
        footnotes += fn
        del_match = DELETION_RE.match(base)
        if del_match:
            start, end = int(del_match.group(1)), int(del_match.group(2))
            if end < start:
                raise ValueError(f"deletion range inverted: {part!r}")
            components.append({"kind": "del", "start": start, "end": end})
            continue
        sub_match = SUBSTITUTION_RE.match(base)
        if sub_match:
            ref_aa = sub_match.group(1)
            pos = int(sub_match.group(2))
            alts = sub_match.group(3).split("/")
            components.append({"kind": "mixed" if len(alts) > 1 else "sub",
                               "ref": ref_aa, "pos": pos, "alts": alts})
            continue
        raise ValueError(f"unparseable substitution token: {part!r}")
    if not components:
        raise ValueError("empty substitution cell")
    return components, footnotes


def parse_nai_cell(cell: str) -> tuple[str, str, str]:
    """Parse an NAI drug cell like 'RI (11–18)' or 'HRI (>7692)'.

    Returns (phenotype, fold_ic50, fold_note). phenotype is '' for '?' cells.
    """
    text = " ".join(cell.split())
    if not text or text.startswith("?"):
        return "", "", ""
    # Footnote letters may follow the closing parenthesis (e.g. 'NI (<2)j');
    # phenotype tokens may contain spaces around the slash (e.g. 'NI /HRI').
    compact = text.replace(" ", "")
    m = re.match(r"^([A-Za-z/]+)(?:\((.*)\))?[a-z]*$", compact)
    if not m:
        raise ValueError(f"unparseable NAI cell: {cell!r}")
    raw_pheno = m.group(1).strip().lower()
    fold_text = (m.group(2) or "").strip()

    tokens = [t.strip() for t in raw_pheno.split("/") if t.strip()]
    mapped = [NAI_TOKEN_MAP[t] for t in tokens]
    # Slash combinations take the highest-severity phenotype.
    phenotype = max(mapped, key=lambda p: PHENOTYPE_SEVERITY[p])

    fold_ic50, fold_note = parse_fold_value(fold_text, metric="IC50")
    return phenotype, fold_ic50, fold_note


def parse_pa_classification(cell: str) -> str:
    """Map a PA classification cell to a ResPro phenotype."""
    text = " ".join(cell.split()).lower()
    if not text or text == "?":
        return ""
    if text not in PA_CLASS_MAP:
        raise ValueError(f"unparseable PA classification: {cell!r}")
    return PA_CLASS_MAP[text]


def parse_fold_value(text: str, metric: str) -> tuple[str, str]:
    """Parse a fold-change value or range into (fold_ic50, note).

    Ranges keep the max endpoint in fold_ic50; the original text (including
    '>'/'<' qualifiers and thousands separators) is preserved in the note.
    """
    if not text or text == "?":
        return "", ""
    cleaned = text.replace(",", "")
    # Space-grouped thousands separators (e.g. '>10 000').
    cleaned = re.sub(r"(\d)\s+(\d)", r"\1\2", cleaned)
    parts = [p.strip() for p in FOLD_RANGE_SPLIT_RE.split(cleaned) if p.strip()]
    if not parts:
        return "", ""
    try:
        values = [float(p.lstrip("<>")) for p in parts]
    except ValueError:
        raise ValueError(f"unparseable fold-change value: {text!r}") from None
    max_value = max(values)
    fold_ic50 = f"{max_value:.10g}"
    note = ""
    if len(parts) > 1:
        note = f"{metric} fold-change in source: {text}"
    elif any(p.startswith(("<", ">")) for p in parts):
        note = f"{metric} fold-change in source: {text}"
    return fold_ic50, note


def parse_citation_numbers(text: str) -> list[int]:
    """Extract citation numbers from a references cell like '(6, 14)'."""
    numbers: list[int] = []
    for m in CITATION_LIST_RE.finditer(text):
        for token in m.group(1).split(","):
            token = token.strip()
            if not token:
                continue
            range_match = re.fullmatch(r"(\d+)\s*[-–−]\s*(\d+)", token)
            if range_match:
                lo, hi = int(range_match.group(1)), int(range_match.group(2))
                if hi < lo or hi - lo > 50:
                    raise ValueError(f"implausible citation range: {token!r}")
                numbers.extend(range(lo, hi + 1))
            else:
                numbers.append(int(token))
    return sorted(set(numbers))


def resolve_position(
    pos: int,
    ref_aa: str,
    ref_seq: str,
    footnotes: str,
    subtype: str,
) -> tuple[int | None, str]:
    """Resolve a substitution position against the reference sequence.

    Tries the direct 1-based position first; for H5N1 rows carrying footnote f
    (full-length N1 numbering) retries with the 20-aa stalk offset. Returns
    (position_or_None, note).
    """
    if 1 <= pos <= len(ref_seq) and ref_seq[pos - 1] == ref_aa:
        return pos, ""
    if (
        "f" in footnotes
        and subtype in H5N1_SUBTYPES
        and pos - H5N1_STALK_OFFSET >= 1
        and ref_seq[pos - H5N1_STALK_OFFSET - 1] == ref_aa
    ):
        note = (
            "Source position uses full-length N1 numbering (footnote f); "
            f"adjusted by -{H5N1_STALK_OFFSET} for the stalk-deleted reference"
        )
        return pos - H5N1_STALK_OFFSET, note
    return None, f"reference mismatch: {ref_aa} not at position {pos}"


def canonicalize_deletion(
    start: int, end: int, ref_seq: str
) -> tuple[int, str, str] | None:
    """Canonicalize a deletion to the anchored ResPro form.

    Returns (anchor_position, reference_anchor_plus_block, anchor_aa) or None
    if the block does not fit the reference sequence.
    """
    anchor = start - 1
    if anchor < 1 or end > len(ref_seq):
        return None
    reference = ref_seq[anchor - 1:end]
    anchor_aa = ref_seq[anchor - 1]
    return anchor, reference, anchor_aa


class RuleRegistry:
    """Aggregates atomic rules, formula groups, and member identities."""

    def __init__(self) -> None:
        self.atomic: dict[tuple, dict] = {}
        self.combos: dict[tuple, dict] = {}
        # member key -> member_id
        self.members: dict[tuple, str] = {}
        self.group_counter = 0

    def next_group_id(self) -> str:
        self.group_counter += 1
        return f"G{self.group_counter:05d}"

    @staticmethod
    def atomic_key(feature: str, ref_id: str, position: int, reference: str,
                   mutation: str, antiviral: str) -> tuple:
        return (feature, ref_id, position, reference, mutation, antiviral)

    def add_atomic(
        self,
        feature: str,
        ref_id: str,
        position: int,
        reference: str,
        mutation: str,
        antiviral: str,
        phenotype: str,
        fold_ic50: str,
        fold_note: str,
        publication: str,
        source: str,
        comments: list[str],
    ) -> None:
        key = self.atomic_key(feature, ref_id, position, reference, mutation,
                              antiviral)
        entry = self.atomic.setdefault(key, {
            "feature": feature,
            "reference_identifier": ref_id,
            "position": position,
            "reference": reference,
            "mutation": mutation,
            "antiviral": antiviral,
            "phenotypes": [],
            "folds": [],
            "fold_notes": [],
            "publications": [],
            "sources": [],
            "comments": [],
        })
        if phenotype:
            entry["phenotypes"].append(phenotype)
        if fold_ic50:
            try:
                entry["folds"].append(float(fold_ic50))
            except ValueError:
                pass
        if fold_note:
            entry["fold_notes"].append(fold_note)
        if publication:
            entry["publications"].append(publication)
        if source:
            entry["sources"].append(source)
        entry["comments"].extend(c for c in comments if c)

    def add_combo(
        self,
        feature: str,
        ref_id: str,
        member_keys: tuple[tuple, ...],
        antiviral: str,
        phenotype: str,
        fold_ic50: str,
        fold_note: str,
        publication: str,
        source: str,
        comments: list[str],
        components: list[dict],
    ) -> None:
        key = (feature, ref_id, member_keys, antiviral)
        entry = self.combos.setdefault(key, {
            "feature": feature,
            "reference_identifier": ref_id,
            "member_keys": member_keys,
            "antiviral": antiviral,
            "components": components,
            "phenotypes": [],
            "folds": [],
            "fold_notes": [],
            "publications": [],
            "sources": [],
            "comments": [],
        })
        if phenotype:
            entry["phenotypes"].append(phenotype)
        if fold_ic50:
            try:
                entry["folds"].append(float(fold_ic50))
            except ValueError:
                pass
        if fold_note:
            entry["fold_notes"].append(fold_note)
        if publication:
            entry["publications"].append(publication)
        if source:
            entry["sources"].append(source)
        entry["comments"].extend(c for c in comments if c)

    def member_id_for(self, member_key: tuple) -> str:
        """Return the member_id for a component, creating one if needed."""
        if member_key in self.members:
            return self.members[member_key]
        member_id = f"M{len(self.members) + 1:05d}"
        self.members[member_key] = member_id
        return member_id


def merge_phenotype(phenotypes: list[str]) -> tuple[str, str]:
    """Pick the most severe phenotype; return (phenotype, disagreement_note)."""
    if not phenotypes:
        return "", ""
    unique = sorted(set(phenotypes), key=lambda p: (-PHENOTYPE_SEVERITY.get(p, 0), p))
    chosen = unique[0]
    note = ""
    if len(unique) > 1:
        note = "Phenotypes in source: " + ", ".join(
            sorted(set(phenotypes))) + f" (most severe used: {chosen})"
    return chosen, note


def merge_fold(folds: list[float]) -> tuple[str, str]:
    """Merge fold values: max endpoint wins; all values kept in the note."""
    if not folds:
        return "", ""
    fold_ic50 = f"{max(folds):.10g}"
    note = ""
    if len(set(folds)) > 1:
        note = "Fold-change values in source: " + ", ".join(
            f"{v:.10g}" for v in sorted(set(folds)))
    return fold_ic50, note


def build_comment(parts: list[str]) -> str:
    return " | ".join(p for p in parts if p)


SOURCE_TAG_RE = re.compile(
    r" \((WHO human NAI|WHO avian NAI|WHO PA)\)(?=:)"
)


def finalize_comment_parts(
    comments: list[str], sources: list[str]
) -> list[str]:
    """Deduplicate comments and drop redundant source tags.

    Comments are tagged with their source label at parse time so that rules
    merged from several tables show which table each block came from. When a
    rule draws from a single table the tag is redundant (the source column
    already says so) and is stripped.
    """
    parts = list(dict.fromkeys(comments))
    if len(set(sources)) <= 1:
        parts = [SOURCE_TAG_RE.sub("", p) for p in parts]
    return parts


# WHO origin-column abbreviations, expanded per the tables' own annotation
# footnotes. Clin and No differ between the NAI tables and the PA table.
ORIGIN_EXPANSIONS = {
    "rg": "reverse genetics",
    "recna": "recombinant NA",
    "sur": "surveillance studies",
    "p-p": "plaque purification",
    "ose": "oseltamivir used",
    "zan": "zanamivir used",
    "lan": "laninamivir used",
    "per": "peramivir used",
    "cell": "cell culture",
    "mice": "mouse model",
    "in vitro": "in vitro",
    "in vivo": "in vivo",
    "in ovo": "in ovo",
}
ORIGIN_EXPANSIONS_NAI = {
    **ORIGIN_EXPANSIONS,
    "clin": "clinical detection",
    "no": "no NAI used",
}
ORIGIN_EXPANSIONS_PA = {
    **ORIGIN_EXPANSIONS,
    "clin": "clinical trial",
    "no": "no baloxavir exposure",
    "bxa": "substitution selected under baloxavir selection pressure",
}


def expand_origin(origin: str, table: str) -> str:
    """Expand origin-column abbreviations to their full annotation terms.

    Keeps the source separators (';', ',', '/', '+') so the study-context
    structure stays intact, drops empty segments (trailing/doubled
    separators), and passes unknown tokens through verbatim with a warning
    so a future upstream abbreviation change is visible in the logs.
    """
    mapping = ORIGIN_EXPANSIONS_PA if table == "pa" else ORIGIN_EXPANSIONS_NAI
    unknown: set[str] = set()
    segments = re.split(r"\s*([;,])\s*", origin)
    out: list[str] = []
    for i, segment in enumerate(segments):
        if i % 2 == 1:
            out.append("; " if segment == ";" else ", ")
            continue
        pieces = re.split(r"\s*([/+])\s*", segment)
        expanded: list[str] = []
        for j, piece in enumerate(pieces):
            if j % 2 == 1:
                expanded.append(piece)
                continue
            token = piece.strip()
            if not token:
                continue
            expansion = mapping.get(token.lower())
            if expansion is None:
                unknown.add(token)
                expanded.append(token)
            else:
                expanded.append(expansion)
        if expanded:
            out.append("".join(expanded))
    for token in sorted(unknown):
        eprint(f"WARNING: unexpanded origin token {token!r} in {table} table")
    return "".join(out).strip(" ;,")


def lookup_publications(
    citation_numbers: list[int],
    table: str,
    lookup: dict[tuple[str, str], dict[str, str]],
) -> tuple[str, list[str]]:
    """Map citation numbers to publication tokens via the lookup sidecar.

    Returns (publication_column, unresolved_notes).
    """
    tokens: list[str] = []
    unresolved: list[str] = []
    for num in citation_numbers:
        row = lookup.get((table, str(num)))
        if row and row.get("doi"):
            # DOIs are case-insensitive; normalize for deterministic dedup.
            tokens.append("doi:" + row["doi"].lower())
        elif row and row.get("pmid"):
            tokens.append(f"PMID:{row['pmid']}")
        else:
            unresolved.append(str(num))
    return join_unique(tokens), [
        f"reference {n} unresolved" for n in unresolved
    ]


def process_table(
    table: str,
    parsed: object,
    lookup: dict[tuple[str, str], dict[str, str]],
    sequences: dict[str, str],
    registry: RuleRegistry,
    non_migrated: list[dict],
) -> None:
    """Normalize one parsed table into the registries."""
    for row in parsed.rows:
        subtype = norm(row.subtype)
        substitution = " ".join(row.substitution.split())
        context = {
            "table": table,
            "page": row.page,
            "subtype": subtype,
            "substitution": substitution,
        }

        def drop(reason: str, details: str = "", antiviral: str = "") -> None:
            non_migrated.append({
                "reason": reason,
                "table": table,
                "page": row.page,
                "subtype": subtype,
                "substitution": substitution,
                "antiviral": antiviral,
                "details": details,
            })

        if not subtype:
            drop("missing_subtype")
            continue

        if table == "pa":
            ref_id = PA_REFERENCE_B if subtype == "B" else PA_REFERENCE_A
            feature = "PA"
        else:
            mapping = HUMAN_NA_REFERENCE if table == "human-nai" else AVIAN_NA_REFERENCE
            ref_id = mapping.get(subtype)
            feature = "NA"
        if not ref_id:
            drop("missing_reference_mapping")
            continue
        ref_seq = sequences.get(ref_id)
        if ref_seq is None:
            drop("missing_reference_sequence", f"no FASTA for {ref_id}")
            continue

        try:
            components, footnotes = parse_substitution_cell(substitution)
        except ValueError as exc:
            drop("unparseable_substitution", str(exc))
            continue

        # Resolve every component to canonical atomic fields.
        canonical_components: list[dict] = []
        failed = False
        or_group_counter = 0
        for comp in components:
            if comp["kind"] == "del":
                canon = canonicalize_deletion(comp["start"], comp["end"], ref_seq)
                if canon is None:
                    drop("deletion_out_of_range",
                         f"Del {comp['start']}-{comp['end']} vs {len(ref_seq)} aa reference")
                    failed = True
                    break
                anchor, reference, anchor_aa = canon
                canonical_components.append({
                    "kind": "del", "position": anchor, "reference": reference,
                    "mutation": anchor_aa,
                    "member_token": f"{anchor_aa}{anchor}del",
                    "or_group": None,
                })
                continue
            position, note = resolve_position(
                comp["pos"], comp["ref"], ref_seq, footnotes, subtype
            )
            if position is None:
                drop("reference_mismatch", note)
                failed = True
                break
            or_group = None
            if len(comp["alts"]) > 1:
                or_group = or_group_counter
                or_group_counter += 1
            for alt in comp["alts"]:
                canonical_components.append({
                    "kind": comp["kind"], "position": position,
                    "reference": comp["ref"], "mutation": alt,
                    "member_token": f"{comp['ref']}{position}{alt}",
                    "offset_note": note,
                    "or_group": or_group,
                })
        if failed:
            continue

        # Per-drug phenotype/fold values.
        if table == "pa":
            phenotype = parse_pa_classification(row.classification)
            fold_ic50, fold_note = parse_fold_value(
                " ".join(row.fold_change.split()), metric="EC50"
            )
            drug_values = [("baloxavir", "", phenotype, fold_ic50, fold_note)]
        else:
            drug_values = []
            for col, drug in NAI_DRUGS.items():
                cell = norm(getattr(row, col))
                if not cell:
                    continue
                try:
                    phenotype, fold_ic50, fold_note = parse_nai_cell(cell)
                except ValueError as exc:
                    drop("unparseable_drug_cell", str(exc), antiviral=drug)
                    failed = True
                    break
                drug_values.append((drug, "", phenotype, fold_ic50, fold_note))
        if failed:
            continue

        citation_numbers = parse_citation_numbers(" ".join(row.references.split()))
        publication, unresolved_notes = lookup_publications(
            citation_numbers, table, lookup
        )

        # The free-text comments column is deliberately not imported; only
        # the structured origin column (its original heading) is kept. Every
        # per-source comment carries its source label so that merged rows
        # (rules imported from more than one table) show which table each
        # block came from.
        comments = []
        origin = expand_origin(" ".join(row.origin.split()), table)
        if origin:
            comments.append(f"Origin of virus ({SOURCE_LABELS[table]}): " + origin)
        comments.extend(
            f"reference {n} unresolved ({SOURCE_LABELS[table]})"
            for n in unresolved_notes
        )

        if table == "pa":
            drug, _, phenotype, fold_ic50, fold_note = drug_values[0]
            if not phenotype and not fold_ic50:
                drop("all_phenotypes_unknown", antiviral=drug)
                continue
        else:
            if not any(dv[2] or dv[3] for dv in drug_values):
                drop("all_phenotypes_unknown")
                continue

        # Single-component rows (substitutions and deletions alike) become
        # atomic rules; multi-component rows become formula groups.
        if len(canonical_components) == 1:
            comp = canonical_components[0]
            for drug, _, phenotype, fold_ic50, fold_note in drug_values:
                if not phenotype and not fold_ic50:
                    continue
                registry.add_atomic(
                    feature, ref_id, comp["position"], comp["reference"],
                    comp["mutation"], drug, phenotype, fold_ic50, fold_note,
                    publication, SOURCE_LABELS[table], comments,
                )
            continue

        # Combination row: register a formula group. Each component carries
        # its member identity and canonical atomic fields so members can be
        # created or reused at finalization time.
        combo_components = [
            {
                "member_key": (feature, ref_id, comp["member_token"]),
                "canonical": (feature, ref_id, comp["position"],
                              comp["reference"], comp["mutation"]),
                "or_group": comp["or_group"],
                "position": comp["position"],
                "reference": comp["reference"],
                "mutation": comp["mutation"],
            }
            for comp in canonical_components
        ]
        member_keys = tuple(c["member_key"] for c in combo_components)
        for drug, _, phenotype, fold_ic50, fold_note in drug_values:
            if not phenotype and not fold_ic50:
                continue
            registry.add_combo(
                feature, ref_id, member_keys, drug, phenotype, fold_ic50,
                fold_note, publication, SOURCE_LABELS[table], comments,
                combo_components,
            )


def finalize_rules(
    registry: RuleRegistry,
) -> tuple[list[dict], list[dict]]:
    """Materialize atomic rule rows and formula rows from the registries."""
    rules_rows: list[dict] = []
    formula_rows: list[dict] = []

    # Atomic rules.
    atomic_rows: dict[tuple, dict] = {}
    for key in sorted(registry.atomic):
        entry = registry.atomic[key]
        phenotype, pheno_note = merge_phenotype(entry["phenotypes"])
        fold_ic50, fold_note = merge_fold(entry["folds"])
        comment_parts = finalize_comment_parts(
            entry["comments"], entry["sources"]
        )
        comment_parts.extend(dict.fromkeys(entry["fold_notes"]))
        if pheno_note:
            comment_parts.append(pheno_note)
        if fold_note:
            comment_parts.append(fold_note)
        row = {
            "feature": entry["feature"],
            "reference_identifier": entry["reference_identifier"],
            "position": entry["position"],
            "reference": entry["reference"],
            "mutation": entry["mutation"],
            "antiviral": entry["antiviral"],
            "member_id": "",
            "phenotype": phenotype,
            "fold_ic50": fold_ic50,
            "publication": join_unique(entry["publications"]),
            "source": join_unique(entry["sources"]),
            "comment": build_comment(comment_parts),
        }
        rules_rows.append(row)
        atomic_rows[key] = row

    # Canonical atomic identity (without antiviral) for member reuse.
    atomic_canonicals: dict[tuple, dict] = {}
    for key, row in atomic_rows.items():
        atomic_canonicals[key[:5]] = row

    # Formula groups: deterministic group ids ordered by (feature, ref_id,
    # first member position, antiviral).
    combo_entries = sorted(
        registry.combos.values(),
        key=lambda e: (
            e["feature"], e["reference_identifier"],
            min(c["position"] for c in e["components"]), e["antiviral"],
        ),
    )
    for entry in combo_entries:
        group_id = registry.next_group_id()

        # Resolve members: reuse an existing singleton atomic rule (attach the
        # member_id to it) or create a blank-antiviral support row.
        member_ids: list[str | None] = []
        or_groups: list[int | None] = []
        for comp in entry["components"]:
            member_id = registry.member_id_for(comp["member_key"])
            atomic_row = atomic_canonicals.get(comp["canonical"])
            if atomic_row is not None:
                atomic_row["member_id"] = member_id
            else:
                support_key = comp["canonical"] + ("",)
                if support_key not in registry.atomic:
                    registry.atomic[support_key] = {
                        "feature": comp["canonical"][0],
                        "reference_identifier": comp["canonical"][1],
                        "position": comp["canonical"][2],
                        "reference": comp["canonical"][3],
                        "mutation": comp["canonical"][4],
                        "antiviral": "",
                        "phenotypes": [],
                        "folds": [],
                        "fold_notes": [],
                        "publications": list(entry["publications"]),
                        "sources": list(entry["sources"]),
                        "comments": [],
                        "member_id": member_id,
                    }
            member_ids.append(member_id)
            or_groups.append(comp["or_group"])

        # Build the expression: AND over components; components from the same
        # mixed-allele token form an OR sub-expression.
        expression = build_expression(member_ids, or_groups)

        phenotype, pheno_note = merge_phenotype(entry["phenotypes"])
        fold_ic50, fold_note = merge_fold(entry["folds"])
        comment_parts = finalize_comment_parts(
            entry["comments"], entry["sources"]
        )
        comment_parts.extend(dict.fromkeys(entry["fold_notes"]))
        if pheno_note:
            comment_parts.append(pheno_note)
        if fold_note:
            comment_parts.append(fold_note)

        formula_rows.append({
            "group_id": group_id,
            "antiviral": entry["antiviral"],
            "expression": expression,
            "phenotype": phenotype,
            "fold_ic50": fold_ic50,
            "publication": join_unique(entry["publications"]),
            "source": join_unique(entry["sources"]),
            "comment": build_comment(comment_parts),
        })

    # Support rows were added to the registry after the atomic pass; emit them
    # now (blank antiviral, member_id set, combo provenance).
    for key in sorted(registry.atomic):
        if key in atomic_rows:
            continue
        entry = registry.atomic[key]
        rules_rows.append({
            "feature": entry["feature"],
            "reference_identifier": entry["reference_identifier"],
            "position": entry["position"],
            "reference": entry["reference"],
            "mutation": entry["mutation"],
            "antiviral": "",
            "member_id": entry.get("member_id", ""),
            "phenotype": "",
            "fold_ic50": "",
            "publication": join_unique(entry["publications"]),
            "source": join_unique(entry["sources"]),
            "comment": "",
        })

    return rules_rows, formula_rows


def build_expression(member_ids: list[str | None], or_groups: list[int | None]) -> str:
    """Join members with AND; consecutive members sharing an or_group form an
    OR sub-expression (mixed alleles)."""
    parts: list[str] = []
    i = 0
    while i < len(member_ids):
        group = or_groups[i]
        if group is None:
            parts.append(member_ids[i])
            i += 1
            continue
        j = i
        while j < len(or_groups) and or_groups[j] == group:
            j += 1
        parts.append("(" + " OR ".join(member_ids[i:j]) + ")")
        i = j
    return "(" + " AND ".join(parts) + ")"


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Convert WHO influenza marker tables to ResPro artifacts"
    )
    parser.add_argument("--source-dir", required=True, type=Path,
                        help="directory containing the three WHO PDFs")
    parser.add_argument("--lookup-file", required=True, type=Path,
                        help="reference-lookup.tsv sidecar from resolve_references.py")
    parser.add_argument("--output-dir", required=True, type=Path,
                        help="output directory for ResPro artifacts")
    args = parser.parse_args()

    scripts_dir = Path(__file__).resolve().parent
    sequences = load_reference_sequences(scripts_dir / "references")
    if not sequences:
        eprint("ERROR: no reference sequences found in scripts/references/")
        return 2

    lookup: dict[tuple[str, str], dict[str, str]] = {}
    if args.lookup_file.exists():
        lines = args.lookup_file.read_text(encoding="utf-8").splitlines()
        if len(lines) >= 2:
            header = lines[0].split("\t")
            for line in lines[1:]:
                values = line.split("\t")
                if len(values) != len(header):
                    continue
                row = dict(zip(header, values))
                lookup[(row["table"], row["citation_number"])] = row
        eprint(f"Loaded {len(lookup)} citation lookups from {args.lookup_file}")
    else:
        eprint(
            f"WARNING: lookup file {args.lookup_file} not found; "
            "all publications will be unresolved"
        )

    parsers = {
        "human-nai": parse_human_nai,
        "avian-nai": parse_avian_nai,
        "pa": parse_pa,
    }

    registry = RuleRegistry()
    non_migrated: list[dict] = []
    source_state: dict[str, dict] = {}

    for table, parse_fn in parsers.items():
        pdf_path = args.source_dir / SOURCE_FILES[table]
        if not pdf_path.exists():
            eprint(f"ERROR: missing source PDF: {pdf_path}")
            return 2
        parsed = parse_fn(str(pdf_path))
        eprint(f"Source date: {parsed.footer_date} ({SOURCE_LABELS[table]})")
        eprint(f"Parsed {len(parsed.rows)} source rows ({SOURCE_LABELS[table]})")
        sha = hashlib.sha256(pdf_path.read_bytes()).hexdigest()
        source_state[table] = {
            "url": SOURCE_URLS[table],
            "sha256": f"sha256:{sha}",
            "footer_date": parsed.footer_date,
        }
        process_table(table, parsed, lookup, sequences, registry, non_migrated)

    rules_rows, formula_rows = finalize_rules(registry)

    # Fail-closed validation of every emitted rule against the references.
    validated_rows: list[dict] = []
    for row in rules_rows:
        ref_seq = sequences.get(row["reference_identifier"])
        pos = int(row["position"])
        if ref_seq is None:
            non_migrated.append({
                "reason": "missing_reference_sequence", "table": "",
                "page": "", "subtype": "", "substitution": "",
                "antiviral": row["antiviral"],
                "details": f"no FASTA for {row['reference_identifier']}",
            })
            continue
        reference = row["reference"]
        if len(reference) > 1:  # deletion: anchor + deleted block
            block = ref_seq[pos - 1:pos - 1 + len(reference)]
            if block != reference or ref_seq[pos - 1] != row["mutation"]:
                non_migrated.append({
                    "reason": "deletion_block_mismatch", "table": "",
                    "page": "", "subtype": "", "substitution": "",
                    "antiviral": row["antiviral"],
                    "details": f"expected {reference} at {pos}, found {block}",
                })
                continue
        elif ref_seq[pos - 1] != reference:
            non_migrated.append({
                "reason": "reference_mismatch", "table": "",
                "page": "", "subtype": "", "substitution": "",
                "antiviral": row["antiviral"],
                "details": (
                    f"{row['reference_identifier']} position {pos} is "
                    f"{ref_seq[pos - 1]}, rule expects {reference}"
                ),
            })
            continue
        validated_rows.append(row)
    dropped = len(rules_rows) - len(validated_rows)
    if dropped:
        eprint(f"DROPPED: {dropped} rule rows failed reference validation")
    rules_rows = validated_rows

    # Sort deterministically.
    rules_rows.sort(key=lambda r: (
        norm(r["reference_identifier"]), norm(r["feature"]),
        int(r["position"]), norm(r["mutation"]), norm(r["antiviral"]),
    ))
    formula_rows.sort(key=lambda r: (norm(r["group_id"]), norm(r["antiviral"])))

    out_dir = args.output_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    def tsv_from_rows(rows: list[dict], columns: list[str]) -> str:
        lines = ["\t".join(columns)]
        for row in rows:
            lines.append("\t".join(norm(row.get(col, "")) for col in columns))
        return "\n".join(lines) + "\n"

    rules_content = tsv_from_rows(rules_rows, RULES_COLUMNS)
    rules_path = out_dir / "rules.tsv"
    rules_path.write_text(rules_content, encoding="utf-8")
    eprint(f"Written {rules_path} ({len(rules_rows)} rows).")

    formula_path = out_dir / "formula-rules.tsv"
    if formula_rows:
        formula_path.write_text(tsv_from_rows(formula_rows, FORMULA_COLUMNS),
                                encoding="utf-8")
        eprint(f"Written {formula_path} ({len(formula_rows)} rows).")
    else:
        formula_path.unlink(missing_ok=True)

    non_migrated_path = out_dir / "non-migrated-rules.txt"
    if non_migrated:
        lines = [
            "# Non-migrated WHO Influenza rows",
            "# Columns: " + "\t".join(NON_MIGRATED_COLUMNS),
            "",
            "\t".join(NON_MIGRATED_COLUMNS),
        ]
        for row in sorted(non_migrated, key=lambda r: (
            norm(r["reason"]), norm(r["table"]), int(r["page"] or 0),
            norm(r["subtype"]), norm(r["substitution"]),
        )):
            lines.append("\t".join(norm(row.get(col, "")) for col in NON_MIGRATED_COLUMNS))
        non_migrated_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
        eprint(f"Written {non_migrated_path} ({len(non_migrated)} aggregated entries).")
    else:
        non_migrated_path.unlink(missing_ok=True)

    metadata = {
        "maintainers": ["Expert Working Group on Antiviral Susceptibility"],
        "contact": "gisrs-whohq@who.int",
        "publication_pmid": "",
        "website": (
            "https://www.who.int/teams/global-influenza-programme/"
            "laboratory-network/quality-assurance/antiviral-susceptibility-influenza"
        ),
        "description": (
            "WHO influenza antiviral susceptibility marker tables: "
            "neuraminidase inhibitor resistance markers for human and avian "
            "influenza viruses and polymerase acidic (PA) baloxavir markers, "
            "curated from the WHO marker tables."
        ),
        "maintainer_update": max(
            state["footer_date"] for state in source_state.values()
        ),
        "license": "CC-BY-NC-SA-3.0-IGO",
        "tsv_checksum": "sha256:" + hashlib.sha256(
            rules_content.encode("utf-8")
        ).hexdigest(),
        "interpretation_algorithms": [
            {
                "name": "drug_interpretation",
                "method": "by_phenotype",
            },
        ],
    }
    metadata_path = out_dir / "metadata.json"
    metadata_path.write_text(
        json.dumps(metadata, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    eprint(f"Written {metadata_path}.")

    source_state_path = out_dir / "source-state.json"
    source_state_path.write_text(
        json.dumps({"sources": source_state}, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    eprint(f"Written {source_state_path}.")

    eprint("Done.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
