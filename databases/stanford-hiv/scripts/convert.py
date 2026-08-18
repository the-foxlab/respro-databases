#!/usr/bin/env python3
"""
Convert Stanford HIVDB ASI XML into ResPro-compatible TSV artifacts.

Outputs written to databases/stanford-hiv/output/:
  rules.tsv           – one row per atomic mutation or formula member
  formula-rules.tsv   – boolean combination rules (emitted only when needed)
  metadata.json       – provenance, version, and checksum
  non-migrated-rules.txt – score rules that could not be converted

Data sources fetched at runtime:
  HIVDB ASI XML             – drug-resistance scoring rules (Stanford HIVDB)
  hivfacts drugs.json       – drug abbreviation → full name mapping
  hivfacts genes_hiv1.json  – segment reference sequences and genomic coordinates

Pipeline stages:
  1. Fetch XML + hivfacts auxiliary data
  2. Validate XML metadata (ALGNAME, version, date)
  3. Build drug-to-segment map from DEFINITIONS
  4. Extract COMMENT_STRING hints (reference AAs and annotation text)
  5. For each DRUG, parse CONDITION into scored score items
  6. For each score item:
       a. Non-MAX: parse LHS, resolve references, emit atomic or formula row
       b. MAX scope: parse each branch independently:
            - equal-score, no shared members → combine as single OR formula
            - equal-score, shared members   → emit each branch separately
            - mixed-score                   → emit each branch separately
            (Branches are mutually exclusive in single-infection context, so
             emitting them as independent rules is semantically correct.)
  7. Materialise shared formula-member rows
  8. Deduplicate, sort deterministically, write outputs
  9. Validate schema and referential integrity before exit

Mutation rule types handled:
  41L => 5                        – atomic substitution
  65R AND 151M => 10              – AND combination → formula
  67EGNHST => 15                  – compressed alternatives → expanded
  MAX(151L => 30, 151M => 60)     – equal-score: OR formula; mixed: separate rules
  MAX(210W AND 215FY => 10, ...)  – same logic, branch-by-branch
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import sys
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path

# ---------------------------------------------------------------------------
# URLs
# ---------------------------------------------------------------------------

HIVDB_LATEST_URL = (
    "https://raw.githubusercontent.com/hivdb/hivfacts/main/data/algorithms/HIVDB_latest.xml"
)
HIVDB_VERSIONS_URL = (
    "https://raw.githubusercontent.com/hivdb/hivfacts/main/data/algorithms/versions.json"
)
HIVFACTS_DRUGS_URL = (
    "https://raw.githubusercontent.com/hivdb/hivfacts/main/data/drugs.json"
)
HIVFACTS_GENES_URL = (
    "https://raw.githubusercontent.com/hivdb/hivfacts/main/data/genes_hiv1.json"
)

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

# HIV-1 subtype B reference accession used throughout Stanford HIVDB.
REFERENCE_IDENTIFIER = "NC_001802"

SOURCE_LABEL = "Stanford HIVDB"

# Stanford ASI uses short segment codes (CA, PR, RT, IN).
# Maps each to output ResPro feature names.
SEGMENT_TO_GENE: dict[str, str] = {
    "CA": "capsid",
    "PR": "aspartic peptidase",
    "RT": "p66 subunit",
    "IN": "integrase",
}

GENE_ORDER: dict[str, int] = {
    "capsid": 0,
    "aspartic peptidase": 1,
    "p66 subunit": 2,
    "integrase": 3,
}

CANONICAL_AA: frozenset[str] = frozenset("ACDEFGHIKLMNPQRSTVWY")
RESERVED_EXPR_WORDS: frozenset[str] = frozenset({"AND", "OR", "NOT", "XOR"})

RULES_COLUMNS = [
    "feature", "reference_identifier", "position", "reference", "mutation",
    "antiviral", "member_id", "score",
    "ic50", "publication", "source", "comment",
]

FORMULA_COLUMNS = [
    "group_id", "antiviral", "expression", "score",
    "ic50", "publication", "source", "comment",
]

NON_MIGRATED_COLUMNS = ["reason", "drug", "feature", "score", "raw_rule", "details"]

# ---------------------------------------------------------------------------
# Shared utilities
# ---------------------------------------------------------------------------

def norm(v: object) -> str:
    return "" if v is None else str(v).strip()


def sanitize_comment_text(text: str) -> str:
    """Remove unresolved wildcard helper phrases from comment text."""
    cleaned = re.sub(
        r"\s*\$listMutsIn\{[^}]+\}\s*[^.?!]*(?:[.?!]|$)",
        " ",
        norm(text),
        flags=re.IGNORECASE,
    )
    cleaned = re.sub(r"\s+", " ", cleaned).strip()
    return cleaned


def eprint(msg: str) -> None:
    print(msg, file=sys.stderr)


def http_get_bytes(url: str) -> bytes:
    req = urllib.request.Request(
        url,
        headers={
            "User-Agent": "respro-stanford-hivdb-converter/2.0",
            "Accept": "application/vnd.github+json, text/plain, */*",
        },
    )
    with urllib.request.urlopen(req, timeout=60) as resp:
        return resp.read()


def http_get_text(url: str) -> str:
    return http_get_bytes(url).decode("utf-8-sig")


def checksum_text(text: str) -> str:
    return "sha256:" + hashlib.sha256(text.encode("utf-8")).hexdigest()


def checksum_bytes(data: bytes) -> str:
    return "sha256:" + hashlib.sha256(data).hexdigest()


def finite_float(v: str) -> float | None:
    """Return a parsed finite float, or None for non-finite / non-numeric input."""
    try:
        f = float(v)
        return None if (f != f or f in (float("inf"), float("-inf"))) else f
    except (TypeError, ValueError):
        return None


def text_of(el: ET.Element, tag: str, default: str = "") -> str:
    child = el.find(tag)
    return child.text.strip() if (child is not None and child.text) else default


def split_csv(v: str) -> list[str]:
    return [p.strip() for p in v.split(",") if p.strip()]


# ---------------------------------------------------------------------------
# Source resolution
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class SourceInfo:
    xml_url: str
    source_version: str
    source_date: str
    source_sha256: str


def resolve_source(url: str) -> tuple[SourceInfo, bytes]:
    """Download the HIVDB XML (resolving the latest pointer if needed) and
    validate version/date against versions.json.  Returns (SourceInfo, xml_bytes)."""
    url = url.strip()
    if url.endswith("HIVDB_latest.xml"):
        pointer = http_get_text(url).strip()
        if not re.fullmatch(r"HIVDB_.+\.xml", pointer):
            raise ValueError(f"HIVDB_latest.xml returned unexpected content: {pointer!r}")
        xml_url = urllib.parse.urljoin(url, pointer)
    else:
        xml_url = url

    filename = Path(urllib.parse.urlparse(xml_url).path).name
    m = re.fullmatch(r"HIVDB_(.+)\.xml", filename)
    if not m:
        raise ValueError(f"Cannot extract version from XML filename: {filename!r}")
    version = m.group(1)

    versions = json.loads(http_get_text(HIVDB_VERSIONS_URL))
    source_date = ""
    for v, d, virus in versions.get("HIVDB", []):
        if v == version and virus == "HIV1":
            source_date = d
            break
    if not source_date:
        raise ValueError(f"HIVDB {version} HIV1 not found in versions.json")

    xml_bytes = http_get_bytes(xml_url)
    return SourceInfo(
        xml_url=xml_url,
        source_version=version,
        source_date=source_date,
        source_sha256=checksum_bytes(xml_bytes),
    ), xml_bytes


# ---------------------------------------------------------------------------
# hivfacts data loading
# ---------------------------------------------------------------------------

def load_drug_map(drugs_json: str) -> dict[str, str]:
    """Build (displayAbbr.upper()) -> fullName lookup from hivfacts drugs.json.

    Stanford ASI drug codes match hivfacts displayAbbr (e.g. '3TC', 'ATV/R').
    Synonyms (e.g. 'DTG_QD') are also registered.
    """
    mapping: dict[str, str] = {}
    for entry in json.loads(drugs_json):
        abbr = norm(entry.get("displayAbbr")).upper()
        full = norm(entry.get("fullName"))
        if abbr and full:
            mapping[abbr] = full
        for syn in entry.get("synonyms", []):
            if syn:
                mapping[syn.upper()] = full
    return mapping


def load_gene_data(genes_json: str) -> dict[str, tuple[str, str]]:
    """Build segment -> (gene_name, refSequence) from hivfacts genes_hiv1.json."""
    genes_list = json.loads(genes_json)
    by_abstract: dict[str, dict] = {
        g["abstractGene"]: g for g in genes_list if g.get("abstractGene")
    }

    result: dict[str, tuple[str, str]] = {}
    for segment, gene_name in SEGMENT_TO_GENE.items():
        seg = by_abstract.get(segment)
        if seg is None:
            raise ValueError(
                f"hivfacts genes_hiv1.json missing entry for {segment!r}"
            )
        refseq = seg.get("refSequence", "").replace("\n", "").replace(" ", "")
        result[segment] = (gene_name, refseq)

    return result


# ---------------------------------------------------------------------------
# XML extraction
# ---------------------------------------------------------------------------

def extract_score_thresholds(root: ET.Element) -> dict[str, int]:
    """Parse DEFINITIONS/GLOBALRANGE and LEVEL_DEFINITION to derive score thresholds.

    GLOBALRANGE contains score buckets mapped to level ORDER integers, e.g.:
      (-INF TO 9 => 1,  10 TO 14 => 2,  15 TO 29 => 3,  30 TO 59 => 4,  60 TO INF => 5)
    LEVEL_DEFINITION maps ORDER integers to SIR codes (S/I/R).

    Returns a dict with the minimum score for each non-S level name used in
    the ResPro drug_interpretation block, e.g. {"intermediate": 15, "resistant": 60}.
    The threshold is the lowest finite lower-bound score whose level order maps to R or I.
    """
    defs = root.find("DEFINITIONS")
    if defs is None:
        raise ValueError("XML missing DEFINITIONS element")

    # Build order -> SIR mapping from LEVEL_DEFINITION elements.
    order_to_sir: dict[int, str] = {}
    for ld in defs.findall("LEVEL_DEFINITION"):
        order_text = text_of(ld, "ORDER")
        sir = text_of(ld, "SIR").upper()
        if order_text.isdigit() and sir in {"S", "I", "R"}:
            order_to_sir[int(order_text)] = sir

    if not order_to_sir:
        raise ValueError("No LEVEL_DEFINITION entries found in DEFINITIONS")

    # Parse GLOBALRANGE score-bucket strings: "LB TO UB => ORDER"
    global_range_el = defs.find("GLOBALRANGE")
    if global_range_el is None or not global_range_el.text:
        raise ValueError("XML missing DEFINITIONS/GLOBALRANGE element")

    range_text = global_range_el.text.strip().strip("()")
    bucket_re = re.compile(
        r"(-INF|\d+)\s+TO\s+(\d+|INF)\s*=>\s*(\d+)",
        re.IGNORECASE,
    )

    # Collect lowest finite lower-bound for each SIR category.
    sir_min: dict[str, int] = {}
    for m in bucket_re.finditer(range_text):
        lb_raw, _, order_raw = m.group(1), m.group(2), m.group(3)
        order = int(order_raw)
        sir = order_to_sir.get(order)
        if sir is None or sir == "S":
            continue
        lb = None if lb_raw.upper() == "-INF" else int(lb_raw)
        if lb is not None:
            sir_key = "resistant" if sir == "R" else "intermediate"
            if sir_key not in sir_min or lb < sir_min[sir_key]:
                sir_min[sir_key] = lb

    if not sir_min:
        raise ValueError("Could not derive any non-susceptible thresholds from GLOBALRANGE")

    return sir_min


def build_drug_to_gene(root: ET.Element) -> dict[str, str]:
    """Parse DEFINITIONS element -> drug abbreviation -> segment code."""
    defs = root.find("DEFINITIONS")
    if defs is None:
        raise ValueError("XML missing DEFINITIONS element")

    gene_to_classes: dict[str, list[str]] = {}
    class_to_drugs: dict[str, list[str]] = {}

    for gd in defs.findall("GENE_DEFINITION"):
        g = text_of(gd, "NAME")
        if g:
            gene_to_classes[g] = split_csv(text_of(gd, "DRUGCLASSLIST"))

    for dc in defs.findall("DRUGCLASS"):
        c = text_of(dc, "NAME")
        if c:
            class_to_drugs[c] = split_csv(text_of(dc, "DRUGLIST"))

    out: dict[str, str] = {}
    for gene, classes in gene_to_classes.items():
        for cls in classes:
            for drug in class_to_drugs.get(cls, []):
                out[drug] = gene

    if not out:
        raise ValueError("DEFINITIONS yielded no drug->gene mappings")
    return out


@dataclass
class CommentHint:
    """Annotation extracted from a COMMENT_STRING element."""
    reference: str  # reference amino acid at this position
    text: str       # full annotation text (used as comment on matching atomic rules)


_COMMENT_MUT_RE = re.compile(r"\b([A-Z])(\d{1,4})([A-Za-z*_\-/]+)\b")


def extract_comment_hints(
    root: ET.Element,
) -> dict[tuple[str, int, str], CommentHint]:
    """Parse COMMENT_STRING elements into {(segment, position, mutation): CommentHint}.

    COMMENT_STRING ids follow [GENE][pos][mutations], e.g. 'RT227C'.
    The TEXT element starts with mutation notation, e.g. 'F227C is a nonpolymorphic...'.
    We scan the TEXT field for mutation tokens to extract reference AAs and annotation text.
    """
    hints: dict[tuple[str, int, str], CommentHint] = {}

    for cs in root.findall(".//COMMENT_STRING"):
        cs_id = cs.attrib.get("id", "")
        gene_match = re.match(r"^(CA|PR|RT|IN)", cs_id)
        if not gene_match:
            continue
        segment = gene_match.group(1)
        text = text_of(cs, "TEXT")
        text = sanitize_comment_text(text)
        if not text:
            continue

        # Scan the TEXT for mutation tokens; word-boundary rules prevent
        # the gene prefix in cs_id from accidentally yielding false matches.
        for m in _COMMENT_MUT_RE.finditer(text):
            ref = m.group(1)
            pos = int(m.group(2))
            for mut in _expand_muts(m.group(3)):
                if mut and mut != ref:
                    key = (segment, pos, mut)
                    if key not in hints:
                        hints[key] = CommentHint(reference=ref, text=text)

    return hints


# ---------------------------------------------------------------------------
# Mutation token expansion
# ---------------------------------------------------------------------------

def _expand_muts(raw: str) -> list[str]:
    """Expand a compressed mutation string into a list of single-character AA tokens.

    Handles:
      - Compressed alternatives: 'EGNHST' -> ['E','G','N','H','S','T']
      - Stop codon: '*' -> ['*']
      - Deletion shorthands: '-' or 'd' -> ['-']
      - Insertion shorthands: '_' or 'i' -> ['INS_any']  (generic insertion wildcard per ResPro spec)
      - Slash separator: '/' -> ignored (e.g. 'F/Y' -> ['F','Y'])
      - Whole-word forms: 'insertion', 'del', etc.
    """
    lraw = raw.lower().strip()
    if lraw in {"insertion", "insert", "ins"}:
        return ["INS_any"]
    if lraw in {"deletion", "del"}:
        return ["-"]

    out: list[str] = []
    for ch in raw:
        if ch == "*":
            out.append("*")
        elif ch == "-":
            out.append("-")
        elif ch == "d":
            out.append("-")
        elif ch == "i":
            out.append("INS_any")
        elif ch.isalpha() and ch.upper() in CANONICAL_AA:
            out.append(ch.upper())
        # "/" and other separators are silently skipped.
    return out


# ---------------------------------------------------------------------------
# Condition parsing
# ---------------------------------------------------------------------------

def _strip_parens(text: str) -> str:
    """Remove one layer of balanced outer parentheses, if present."""
    text = text.strip()
    if not (text.startswith("(") and text.endswith(")")):
        return text
    depth = 0
    for i, ch in enumerate(text):
        if ch == "(":
            depth += 1
        elif ch == ")":
            depth -= 1
            if depth == 0 and i < len(text) - 1:
                return text  # parens do not wrap the entire string
    return text[1:-1].strip()


def _split_top(text: str, sep: str = ",") -> list[str]:
    """Split *text* at *sep* characters that are not inside parentheses."""
    parts: list[str] = []
    start = depth = 0
    for i, ch in enumerate(text):
        if ch in "([{":
            depth += 1
        elif ch in ")]}":
            depth = max(0, depth - 1)
        elif ch == sep and depth == 0:
            p = text[start:i].strip()
            if p:
                parts.append(p)
            start = i + 1
    tail = text[start:].strip()
    if tail:
        parts.append(tail)
    return parts


@dataclass(frozen=True)
class ScoreItem:
    """One scored branch extracted from a CONDITION string."""
    lhs: str        # mutation expression left of '=>'
    score: str      # numeric score (as string)
    raw: str        # original full text of this branch
    max_scope: str  # non-empty tag when item is inside MAX(...)


def parse_condition(condition: str) -> list[ScoreItem]:
    """Split a CONDITION string into a flat list of ScoreItems.

    Handles nested MAX(...) scopes, assigning a unique scope tag to each.
    Items outside any MAX have an empty max_scope.
    """
    condition = re.sub(r"\s+", " ", norm(condition))
    if not condition:
        return []
    if condition.upper().startswith("SCORE FROM"):
        condition = condition[len("SCORE FROM"):].strip()
    condition = _strip_parens(condition)

    items: list[ScoreItem] = []
    counter = 0

    def visit(frag: str, scope: str = "") -> None:
        nonlocal counter
        frag = _strip_parens(frag.strip())
        for part in _split_top(frag):
            part = part.strip()
            if not part:
                continue
            if part.upper().startswith("MAX"):
                counter += 1
                inner = _strip_parens(part[3:].strip())
                visit(inner, scope=f"MAX{counter:05d}")
                continue
            if "=>" not in part:
                continue
            lhs, score = part.split("=>", 1)
            items.append(ScoreItem(
                lhs=_strip_parens(lhs.strip()),
                score=score.strip(),
                raw=part,
                max_scope=scope,
            ))

    visit(condition)
    return items


# Regex patterns for mutation tokens.
# PAREN form:  '65(NOT V)'  ->  position=65, aas='V', negated
# SIMPLE form: 'K65R', '65R', '67EGNHST'
_TERM_PAREN_RE = re.compile(
    r"(?<![A-Za-z0-9])(?P<pos>\d{1,4})\(\s*(?P<not>NOT\s+)?(?P<aas>[A-Za-z*_\-/]+)\s*\)",
    re.IGNORECASE,
)
_TERM_SIMPLE_RE = re.compile(
    r"(?<![A-Za-z0-9])(?P<ref>[A-Z*]?)(?P<pos>\d{1,4})(?P<aas>[A-Za-z*_\-/]+)(?![A-Za-z0-9])",
    re.IGNORECASE,
)

# ASI markers that indicate unsupported algorithm constructs.
_UNSUPPORTED_MARKERS = frozenset([
    "$", ">=", "<=", "!=", "==", " TO ", "FROM ",
    "SELECT ", "EXCEPT ", "ATLEAST", "NUMBEROF", "COUNT",
])


@dataclass(frozen=True)
class MutTerm:
    position: int
    ref_hint: str   # reference AA hint from the token itself (may be empty)
    mutation: str   # expanded single-character AA token


def _parse_token(text: str) -> tuple[list[MutTerm], bool]:
    """Parse a single mutation token into MutTerms and a negation flag."""
    text = text.strip()

    m = _TERM_PAREN_RE.fullmatch(text)
    if m:
        pos = int(m.group("pos"))
        neg = bool(m.group("not"))
        return [MutTerm(pos, "", aa) for aa in _expand_muts(m.group("aas"))], neg

    m = _TERM_SIMPLE_RE.fullmatch(text)
    if m:
        ref = (m.group("ref") or "").upper()
        pos = int(m.group("pos"))
        terms = []
        for aa in _expand_muts(m.group("aas")):
            # Skip no-op where reference == mutation.
            if ref and aa == ref:
                continue
            terms.append(MutTerm(pos, ref, aa))
        return terms, False

    return [], False


def _has_insertion(lhs: str) -> bool:
    """Return True if the LHS contains a Stanford ASI insertion shorthand ('i' or 'ins').

    Since insertion events are now emitted as INS_any wildcard rules, this
    function is kept for compatibility but always returns False.
    """
    return False


def _normalize_bool(chunk: str) -> str:
    """Normalise boolean operators to canonical uppercase with single spaces."""
    chunk = chunk.replace("&&", " AND ").replace("||", " OR ").replace("+", " AND ")
    for op in ("AND", "OR", "NOT", "XOR"):
        chunk = re.sub(rf"\b{op}\b", f" {op} ", chunk, flags=re.IGNORECASE)
    return re.sub(r"\s+", " ", chunk).strip()


# ---------------------------------------------------------------------------
# Reference resolution
# ---------------------------------------------------------------------------

def resolve_reference(
    segment: str,
    position: int,
    mutation: str,
    ref_hint: str,
    comment_hints: dict[tuple[str, int, str], CommentHint],
    gene_data: dict[str, tuple[str, str]],
) -> str:
    """Determine the reference amino acid for a mutation.

    Resolution order:
      1. Explicit ref_hint embedded in the token (e.g. 'K' in 'K65R').
      2. COMMENT_STRING hint keyed by (segment, position, mutation).
      3. genes_hiv1.json refSequence at (position - 1).
      4. Empty string -- caller treats this as unresolved and skips the rule.
    """
    if ref_hint and ref_hint in CANONICAL_AA:
        return ref_hint

    hint = comment_hints.get((segment, position, mutation))
    if hint and hint.reference in CANONICAL_AA:
        return hint.reference

    _, refseq = gene_data.get(segment, ("", ""))
    if refseq and 1 <= position <= len(refseq):
        return refseq[position - 1]

    return ""


# ---------------------------------------------------------------------------
# Member / group ID generation
# ---------------------------------------------------------------------------

def _make_member_id(gene: str, position: int, reference: str, mutation: str) -> str:
    """Stable, hash-based member ID shared across drugs for the same mutation.

    The hash seed format is identical to the previous converter so that
    member_ids are preserved for mutations unaffected by the rewrite.
    """
    seed = f"{gene}|{position}|{reference}|{mutation}".encode("utf-8")
    mid = "M" + hashlib.sha256(seed).hexdigest()[:14].upper()
    return mid if mid.upper() not in RESERVED_EXPR_WORDS else "M_" + mid


def _make_group_id(drug: str, gene: str, raw_rule: str, score: str) -> str:
    seed = f"{drug}|{gene}|{raw_rule}|{score}".encode("utf-8")
    return "G" + hashlib.sha256(seed).hexdigest()[:12].upper()


# ---------------------------------------------------------------------------
# LHS expression builder
# ---------------------------------------------------------------------------

@dataclass
class ParsedLHS:
    """Result of parsing one ASI LHS expression."""
    expression: str       # ResPro expression using member_ids as atoms
    members: list[dict]   # member dicts (subset of RULES_COLUMNS fields)
    is_formula: bool      # True when boolean operators or >1 member present


def parse_lhs(
    lhs: str,
    segment: str,
    comment_hints: dict[tuple[str, int, str], CommentHint],
    gene_data: dict[str, tuple[str, str]],
) -> ParsedLHS | None:
    """Convert an ASI LHS mutation expression into a ParsedLHS.

    Returns None when the expression is unsupported (unknown ASI constructs)
    or when a reference amino acid cannot be resolved for any mutation term.
    """
    if any(m.upper() in lhs.upper() for m in _UNSUPPORTED_MARKERS):
        return None

    # Collect mutation token spans without overlap.
    occupied = [False] * len(lhs)
    matches: list[re.Match] = []
    for regex in (_TERM_PAREN_RE, _TERM_SIMPLE_RE):
        for m in regex.finditer(lhs):
            if not any(occupied[i] for i in range(m.start(), m.end())):
                for i in range(m.start(), m.end()):
                    occupied[i] = True
                matches.append(m)
    matches.sort(key=lambda x: x.start())
    if not matches:
        return None

    gene_name, _ = gene_data[segment]
    pieces: list[str] = []
    members: list[dict] = []
    last = 0
    saw_bool = False

    for match in matches:
        # Validate the gap between tokens -- only boolean operators and parens allowed.
        gap = _normalize_bool(lhs[last:match.start()])
        if gap:
            cleaned = re.sub(r"\b(AND|OR|NOT|XOR)\b", "", gap, flags=re.IGNORECASE)
            cleaned = cleaned.replace("(", "").replace(")", "").strip()
            if cleaned:
                return None  # unexpected non-boolean content in between tokens
            pieces.append(gap)
            if any(w in gap.upper().split() for w in RESERVED_EXPR_WORDS):
                saw_bool = True

        terms, negative = _parse_token(match.group(0))
        if not terms:
            return None

        token_ids: list[str] = []
        for term in terms:
            ref = resolve_reference(
                segment=segment,
                position=term.position,
                mutation=term.mutation,
                ref_hint=term.ref_hint,
                comment_hints=comment_hints,
                gene_data=gene_data,
            )
            if not ref:
                return None  # unresolved reference -- skip rule

            abs_pos = term.position
            # Stanford ASI '-' (deletion) -> ResPro canonical {ref}{pos}del form.
            mutation = f"{ref}{abs_pos}del" if term.mutation == "-" else term.mutation
            mid = _make_member_id(gene_name, abs_pos, ref, mutation)

            members.append({
                "feature": gene_name,
                "reference_identifier": REFERENCE_IDENTIFIER,
                "position": str(abs_pos),
                "reference": ref,
                "mutation": mutation,
                "member_id": mid,
            })
            token_ids.append(mid)

        if len(token_ids) == 1:
            expr_part = token_ids[0]
        else:
            expr_part = "(" + " OR ".join(token_ids) + ")"
            saw_bool = True

        if negative:
            expr_part = f"(NOT {expr_part})"
            saw_bool = True

        pieces.append(expr_part)
        last = match.end()

    # Validate trailing content.
    tail = _normalize_bool(lhs[last:])
    if tail:
        cleaned = re.sub(r"\b(AND|OR|NOT|XOR)\b", "", tail, flags=re.IGNORECASE)
        cleaned = cleaned.replace("(", "").replace(")", "").strip()
        if cleaned:
            return None
        pieces.append(tail)
        if any(w in tail.upper().split() for w in RESERVED_EXPR_WORDS):
            saw_bool = True

    expr = re.sub(r"\s+", " ", " ".join(p for p in pieces if p)).strip()
    is_formula = saw_bool or len(members) != 1 or expr != members[0]["member_id"]

    return ParsedLHS(expression=expr, members=members, is_formula=is_formula)


# ---------------------------------------------------------------------------
# Conversion context (shared mutable state for one convert() call)
# ---------------------------------------------------------------------------

@dataclass
class _Ctx:
    drug_map: dict[str, str]
    gene_data: dict[str, tuple[str, str]]
    comment_hints: dict[tuple[str, int, str], CommentHint]
    # Absolute-position lookup: (gene_name, abs_pos, mutation) -> annotation text
    abs_comments: dict[tuple[str, int, str], str] = field(default_factory=dict)
    rules_rows: list[dict] = field(default_factory=list)
    member_registry: dict[str, dict] = field(default_factory=dict)
    formula_rows: list[dict] = field(default_factory=list)
    non_migrated: list[dict] = field(default_factory=list)


def _build_abs_comments(
    comment_hints: dict[tuple[str, int, str], CommentHint],
    gene_data: dict[str, tuple[str, str]],
) -> dict[tuple[str, int, str], str]:
    """Re-key COMMENT_STRING hints by segment feature names without coordinate offsets."""
    out: dict[tuple[str, int, str], str] = {}
    for (segment, pos, mut), ch in comment_hints.items():
        if segment in gene_data:
            gene_name, _ = gene_data[segment]
            out[(gene_name, pos, mut)] = ch.text
    return out


# ---------------------------------------------------------------------------
# Rule emission
# ---------------------------------------------------------------------------

def _drug_name(code: str, drug_map: dict[str, str]) -> str:
    name = drug_map.get(code.upper(), code.lower())
    if name.endswith("/r"):
        return name[:-2]
    return name


def _emit_rule(
    parsed: ParsedLHS,
    *,
    drug_code: str,
    score: str,
    raw_rule: str,
    ctx: _Ctx,
    formula_comment: str = "HIVDB ASI formula score rule",
    atomic_comment: str = "",
) -> None:
    """Emit an atomic rule row, or a formula row plus shared member rows."""
    drug = _drug_name(drug_code, ctx.drug_map)

    if not parsed.is_formula:
        mem = parsed.members[0]
        # Prefer an annotation from COMMENT_STRING over the generic fallback.
        comment = (
            ctx.abs_comments.get((mem["feature"], int(mem["position"]), mem["mutation"]))
            or atomic_comment
            or "HIVDB ASI atomic score rule"
        )
        comment = sanitize_comment_text(comment)
        ctx.rules_rows.append({
            "feature": mem["feature"],
            "reference_identifier": mem["reference_identifier"],
            "position": mem["position"],
            "reference": mem["reference"],
            "mutation": mem["mutation"],
            "antiviral": drug,
            "member_id": "",
            "phenotype": "",
            "score": score,
            "ic50": "",
            "publication": "",
            "source": SOURCE_LABEL,
            "comment": comment,
        })
        return

    # Formula rule: register members (shared across drugs) and emit a formula row.
    group_id = _make_group_id(drug_code, parsed.members[0]["feature"], raw_rule, score)

    for mem in parsed.members:
        mid = mem["member_id"]
        if mid not in ctx.member_registry:
            ctx.member_registry[mid] = {**mem}

    ctx.formula_rows.append({
        "group_id": group_id,
        "antiviral": drug,
        "expression": parsed.expression,
        "phenotype": "",
        "score": score,
        "ic50": "",
        "publication": "",
        "source": SOURCE_LABEL,
        "comment": formula_comment,
    })


# ---------------------------------------------------------------------------
# MAX scope handling
# ---------------------------------------------------------------------------

def _handle_max_scope(
    scope_items: list[ScoreItem],
    *,
    drug_code: str,
    segment: str,
    ctx: _Ctx,
) -> None:
    """Process one MAX(...) scope.

    Strategy:
      - Parse every branch individually.  Unparsable branches (insertion events,
        unresolved references, unsupported syntax) are recorded in non_migrated.
      - If all branches share a single score and no member_id would repeat in the
        combined OR expression: emit one combined OR formula.
      - Otherwise: emit each parsed branch as an independent rule.
        This is correct because MAX branches represent mutually exclusive events
        in a single-infection context (the MAX semantics only matter for
        quasi-species / mixed-population scoring, which is outside ResPro scope).
    """
    parsed_branches: list[tuple[ParsedLHS, str]] = []

    for item in scope_items:
        if _has_insertion(item.lhs):
            _add_non_migrated(
                ctx, "non_portable_insertion_event",
                drug=drug_code, gene=segment, score=item.score, raw_rule=item.raw,
                details="Stanford ASI insertion shorthand (e.g. 69i) has no explicit inserted sequence",
            )
            continue

        parsed = parse_lhs(item.lhs, segment, ctx.comment_hints, ctx.gene_data)
        if parsed is None:
            _add_non_migrated(
                ctx, "unsupported_or_unresolved_expression",
                drug=drug_code, gene=segment, score=item.score, raw_rule=item.raw,
            )
            continue

        parsed_branches.append((parsed, item.score))

    if not parsed_branches:
        return

    # Group branches by score.
    by_score: dict[str, list[ParsedLHS]] = defaultdict(list)
    for pb, sc in parsed_branches:
        by_score[sc].append(pb)

    if len(by_score) == 1:
        # Equal-score MAX: attempt a single combined OR formula.
        score = next(iter(by_score))
        _try_combined_or(
            parts=by_score[score],
            score=score,
            raw_rule="MAX:" + "|".join(i.raw for i in scope_items),
            drug_code=drug_code,
            ctx=ctx,
        )
        return

    # Mixed-score MAX: verify that all branches are mutually exclusive before splitting.
    # Mutual exclusivity holds when every pair of branches shares at least one amino-acid
    # position — two branches that overlap on at least one position cannot both be true at
    # the same time in a single-infection context (only one allele can occupy a position).
    # Branches at completely disjoint positions could co-occur and must not be split, because
    # independent scoring would incorrectly add their scores instead of taking the maximum.
    branch_positions = [
        frozenset(int(mem["position"]) for mem in pb.members)
        for pb, _ in parsed_branches
    ]
    can_co_occur = any(
        p1.isdisjoint(p2)
        for i, p1 in enumerate(branch_positions)
        for p2 in branch_positions[i + 1 :]
        if p1 and p2
    )
    if can_co_occur:
        _add_non_migrated(
            ctx, "mixed_score_max_branches_not_mutually_exclusive",
            drug=drug_code, gene=segment, score="",
            raw_rule="MAX:" + "|".join(i.raw for i in scope_items),
            details=(
                "Mixed-score MAX branches involve disjoint amino-acid positions and could "
                "co-occur; splitting into independent rules would over-score additively."
            ),
        )
        return

    # All branches share at least one position, so they are mutually exclusive: only one
    # branch can match at a time.  Emit each as an independent rule; only the matching
    # branch contributes its score to the drug total.
    branch_cmt = "HIVDB ASI MAX branch (verified mutually exclusive: branches share amino-acid positions)"
    for parsed, score in parsed_branches:
        _emit_rule(
            parsed,
            drug_code=drug_code, score=score,
            raw_rule="MAX branch: " + parsed.expression,
            ctx=ctx,
            formula_comment=branch_cmt,
            atomic_comment=branch_cmt,
        )


def _try_combined_or(
    parts: list[ParsedLHS],
    score: str,
    raw_rule: str,
    drug_code: str,
    ctx: _Ctx,
) -> None:
    """Combine equal-score MAX branches into one OR formula when safe to do so.

    Falls back to emitting each branch separately when the combined expression
    would contain duplicate member_ids (forbidden by the ResPro schema validator).
    """
    # Deduplicate branch expressions before combining.
    seen_exprs: set[str] = set()
    unique_exprs: list[str] = []
    all_members: list[dict] = []
    for parsed in parts:
        if parsed.expression not in seen_exprs:
            seen_exprs.add(parsed.expression)
            unique_exprs.append(f"({parsed.expression})")
            all_members.extend(parsed.members)

    combined_expr = " OR ".join(unique_exprs)

    # Check for duplicate member IDs -- disallowed in a single formula expression.
    ids_in_expr = [
        t for t in re.findall(r"\b[A-Za-z_][A-Za-z0-9_]*\b", combined_expr)
        if t.upper() not in RESERVED_EXPR_WORDS
    ]
    if len(ids_in_expr) != len(set(ids_in_expr)):
        # Shared members detected: fall back to separate branch rules.
        cmt = "HIVDB ASI equal-score MAX branch"
        for parsed in parts:
            _emit_rule(
                parsed,
                drug_code=drug_code, score=score,
                raw_rule="MAX branch: " + parsed.expression,
                ctx=ctx,
                formula_comment=cmt, atomic_comment=cmt,
            )
        return

    # Safe to combine into a single OR formula.
    combined = ParsedLHS(expression=combined_expr, members=all_members, is_formula=True)
    _emit_rule(
        combined,
        drug_code=drug_code, score=score, raw_rule=raw_rule,
        ctx=ctx,
        formula_comment="HIVDB ASI equal-score MAX represented as OR",
    )


def _add_non_migrated(
    ctx: _Ctx, reason: str, *,
    drug: str, feature: str, score: str, raw_rule: str, details: str = "",
) -> None:
    ctx.non_migrated.append({
        "reason": reason, "drug": drug, "feature": feature,
        "score": score, "raw_rule": raw_rule, "details": details,
    })


# ---------------------------------------------------------------------------
# Main conversion pipeline
# ---------------------------------------------------------------------------

def convert(
    source_info: SourceInfo,
    xml_bytes: bytes,
    drug_map: dict[str, str],
    gene_data: dict[str, tuple[str, str]],
    output_dir: Path,
) -> None:
    """Run the full conversion pipeline and write all output artifacts."""
    root = ET.fromstring(xml_bytes)

    # Validate XML metadata against resolved source provenance.
    alg_name = text_of(root, "ALGNAME")
    alg_version = text_of(root, "ALGVERSION")
    alg_date = text_of(root, "ALGDATE")
    if alg_name != "HIVDB":
        raise ValueError(f"Expected ALGNAME=HIVDB, got {alg_name!r}")
    if alg_version != source_info.source_version:
        raise ValueError(
            f"Filename version {source_info.source_version!r} != XML ALGVERSION {alg_version!r}"
        )
    if alg_date != source_info.source_date:
        raise ValueError(
            f"versions.json date {source_info.source_date!r} != XML ALGDATE {alg_date!r}"
        )

    drug_to_gene = build_drug_to_gene(root)
    comment_hints = extract_comment_hints(root)
    score_thresholds = extract_score_thresholds(root)

    ctx = _Ctx(
        drug_map=drug_map,
        gene_data=gene_data,
        comment_hints=comment_hints,
        abs_comments=_build_abs_comments(comment_hints, gene_data),
    )

    for drug_el in root.findall("DRUG"):
        drug_code = text_of(drug_el, "NAME")
        segment = drug_to_gene.get(drug_code, "")

        if not drug_code or not segment:
            _add_non_migrated(
                ctx, "drug_without_gene_mapping",
                drug=drug_code, gene=segment, score="", raw_rule="",
                details="Could not map DRUG/NAME to gene via DEFINITIONS",
            )
            continue

        if segment not in gene_data:
            _add_non_migrated(
                ctx, "drug_without_gene_mapping",
                drug=drug_code, gene=segment, score="", raw_rule="",
                details=f"Segment {segment!r} not in supported SEGMENT_TO_GENE",
            )
            continue

        rule_node = drug_el.find("./RULE")
        condition = text_of(rule_node if rule_node is not None else drug_el, "CONDITION")
        items = parse_condition(condition)

        if not items:
            _add_non_migrated(
                ctx, "no_score_items_extracted",
                drug=drug_code, gene=segment, score="", raw_rule=condition[:500],
            )
            continue

        # Separate non-MAX items from MAX-scoped groups.
        max_scopes: dict[str, list[ScoreItem]] = defaultdict(list)
        for item in items:
            max_scopes[item.max_scope].append(item)

        # Process direct (non-MAX) score items.
        for item in max_scopes.pop("", []):
            if _has_insertion(item.lhs):
                _add_non_migrated(
                    ctx, "non_portable_insertion_event",
                    drug=drug_code, gene=segment, score=item.score, raw_rule=item.raw,
                    details="Stanford ASI insertion shorthand (e.g. 69i) has no explicit inserted sequence",
                )
                continue

            parsed = parse_lhs(item.lhs, segment, ctx.comment_hints, ctx.gene_data)
            if parsed is None:
                _add_non_migrated(
                    ctx, "unsupported_or_unresolved_expression",
                    drug=drug_code, gene=segment, score=item.score, raw_rule=item.raw,
                )
                continue

            _emit_rule(parsed, drug_code=drug_code, score=item.score, raw_rule=item.raw, ctx=ctx)

        # Process MAX scopes.
        for scope_items in max_scopes.values():
            _handle_max_scope(scope_items, drug_code=drug_code, segment=segment, ctx=ctx)

    # Build mutation -> member_id mapping for formula members.
    # member_id will be attached to the first scored occurrence of each mutation if it has one.
    mutation_to_member_id: dict[tuple, str] = {}
    for item in ctx.member_registry.values():
        key = (norm(item["feature"]), norm(item["position"]), norm(item["reference"]), norm(item["mutation"]))
        if key not in mutation_to_member_id:
            mutation_to_member_id[key] = item["member_id"]

    # Deduplicate and sort deterministically.
    rules_rows = _dedupe(ctx.rules_rows, RULES_COLUMNS)
    formula_rows = _dedupe(ctx.formula_rows, FORMULA_COLUMNS)

    # Attach member_id to first occurrence of each formula-member mutation with a score.
    seen_mutations: set[tuple] = set()
    for row in rules_rows:
        key = (norm(row["feature"]), norm(row["position"]), norm(row["reference"]), norm(row["mutation"]))
        if key in mutation_to_member_id and key not in seen_mutations and row["antiviral"]:
            row["member_id"] = mutation_to_member_id[key]
            seen_mutations.add(key)

    # Emit placeholder rows for formula members without any scored atomic rules.
    for item in ctx.member_registry.values():
        key = (norm(item["feature"]), norm(item["position"]), norm(item["reference"]), norm(item["mutation"]))
        if key not in seen_mutations:
            rules_rows.append({
            "feature": item["feature"],
                "reference_identifier": item["reference_identifier"],
                "position": item["position"],
                "reference": item["reference"],
                "mutation": item["mutation"],
                "antiviral": "",
                "member_id": item["member_id"],
                "phenotype": "",
                "score": "",
                "ic50": "",
                "publication": "",
                "source": SOURCE_LABEL,
                "comment": "HIVDB ASI formula member",
            })

    rules_rows.sort(key=lambda r: (
        GENE_ORDER.get(norm(r["feature"]), 99),
        norm(r["feature"]),
        int(norm(r.get("position")) or 0),
        norm(r["reference"]),
        norm(r["mutation"]),
        norm(r["antiviral"]),
        norm(r["member_id"]),
        norm(r["score"]),
    ))
    formula_rows.sort(key=lambda r: (
        norm(r["group_id"]),
        norm(r["antiviral"]),
        norm(r["expression"]),
        norm(r["score"]),
    ))

    # Write output files.
    output_dir.mkdir(parents=True, exist_ok=True)

    rules_content = _tsv(rules_rows, RULES_COLUMNS)
    rules_path = output_dir / "rules.tsv"
    rules_path.write_text(rules_content, encoding="utf-8")
    eprint(f"Written {rules_path} ({len(rules_rows)} rows).")

    formula_path = output_dir / "formula-rules.tsv"
    if formula_rows:
        formula_content = _tsv(formula_rows, FORMULA_COLUMNS)
        formula_path.write_text(formula_content, encoding="utf-8")
        eprint(f"Written {formula_path} ({len(formula_rows)} rows).")
    else:
        formula_path.unlink(missing_ok=True)

    metadata = {
        "maintainers": ["Stanford HIVDB Team"],
        "contact": "hivdbteam@lists.stanford.edu",
        "publication_pmid": "22286876",
        "website": "https://hivdb.stanford.edu/",
        "description": (
            "Stanford HIVDB ASI drug-resistance scoring rules."
        ),
        "maintainer_update": source_info.source_date,
        "license": "GNU General Public License v3.0",
        "tsv_checksum": checksum_text(rules_content),
        "interpretation_algorithms": [
            {
                "name": "drug_groups",
                "groups": {
                    "Nucleoside Reverse Transcriptase Inhibitors": [
                        "abacavir",
                        "didanosine",
                        "emtricitabine",
                        "islatravir",
                        "lamivudine",
                        "stavudine",
                        "tenofovir",
                        "zidovudine",
                    ],
                    "Non-Nucleoside Reverse Transcriptase Inhibitors": [
                        "dapivirine",
                        "doravirine",
                        "efavirenz",
                        "etravirine",
                        "nevirapine",
                        "rilpivirine",
                    ],
                    "Protease Inhibitors": [
                        "atazanavir",
                        "darunavir",
                        "fosamprenavir",
                        "indinavir",
                        "lopinavir",
                        "nelfinavir",
                        "saquinavir",
                        "tipranavir",
                    ],
                    "Integrase Strand Transfer Inhibitors": [
                        "bictegravir",
                        "cabotegravir",
                        "dolutegravir",
                        "elvitegravir",
                        "raltegravir",
                    ],
                    "Capsid Inhibitors": ["lenacapavir"],
                },
            },
            {
                "name": "drug_interpretation",
                "method": "by_score",
                "thresholds": score_thresholds,
            },
            {
                "name": "drug_alias",
                "groups": {
                    "abacavir": "ABC",
                    "atazanavir": "ATV",
                    "bictegravir": "BIC",
                    "cabotegravir": "CAB",
                    "dapivirine": "DPV",
                    "darunavir": "DRV",
                    "didanosine": "DDI",
                    "dolutegravir": "DTG",
                    "doravirine": "DOR",
                    "efavirenz": "EFV",
                    "elvitegravir": "EVG",
                    "emtricitabine": "FTC",
                    "etravirine": "ETR",
                    "fosamprenavir": "FPV",
                    "indinavir": "IDV",
                    "islatravir": "ISL",
                    "lamivudine": "3TC",
                    "lenacapavir": "LEN",
                    "lopinavir": "LPV",
                    "nelfinavir": "NFV",
                    "nevirapine": "NVP",
                    "rilpivirine": "RPV",
                    "raltegravir": "RAL",
                    "saquinavir": "SQV",
                    "stavudine": "D4T",
                    "tenofovir": "TDF",
                    "tipranavir": "TPV",
                    "zidovudine": "AZT",
                },
            },
        ],
    }
    metadata_path = output_dir / "metadata.json"
    metadata_path.write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    eprint(f"Written {metadata_path}.")

    nm_content = _non_migrated_text(ctx.non_migrated)
    nm_path = output_dir / "non-migrated-rules.txt"
    nm_path.write_text(nm_content, encoding="utf-8")
    eprint(f"Written {nm_path} ({len(ctx.non_migrated)} aggregated entries).")

    validate_outputs(rules_path, formula_path if formula_rows else None)


# ---------------------------------------------------------------------------
# Output formatting helpers
# ---------------------------------------------------------------------------

def _tsv(rows: list[dict], columns: list[str]) -> str:
    lines = ["\t".join(columns)]
    for row in rows:
        lines.append("\t".join(norm(row.get(col, "")) for col in columns))
    return "\n".join(lines) + "\n"


def _non_migrated_text(rows: list[dict]) -> str:
    lines = [
        "# Non-migrated Stanford HIVDB rules",
        "# Columns: " + "\t".join(NON_MIGRATED_COLUMNS),
        "",
        "\t".join(NON_MIGRATED_COLUMNS),
    ]
    for row in sorted(rows, key=lambda r: (
        norm(r.get("reason")), norm(r.get("feature")),
        norm(r.get("drug")), norm(r.get("raw_rule")),
    )):
        lines.append("\t".join(norm(row.get(c, "")) for c in NON_MIGRATED_COLUMNS))
    return "\n".join(lines) + "\n"


def _dedupe(rows: list[dict], columns: list[str]) -> list[dict]:
    seen: set = set()
    out: list[dict] = []
    for row in rows:
        key = tuple(norm(row.get(c, "")) for c in columns)
        if key not in seen:
            seen.add(key)
            out.append(row)
    return out


# ---------------------------------------------------------------------------
# Output validation
# ---------------------------------------------------------------------------

def validate_outputs(rules_path: Path, formula_path: Path | None) -> None:
    """Validate artifact schemas and referential integrity.

    Checks:
      - rules.tsv has the correct header and no empty required cells.
            - formula-rules.tsv (if present) references only known member_ids.
      - No member_id appears more than once in a single formula expression.
    """
    with rules_path.open(newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if reader.fieldnames != RULES_COLUMNS:
            raise ValueError(f"rules.tsv header mismatch: {reader.fieldnames}")

        member_ids: set[str] = set()
        for lineno, row in enumerate(reader, start=2):
            if not row["feature"]:
                raise ValueError(f"rules.tsv:{lineno} missing feature")
            if not norm(row["position"]).isdigit():
                raise ValueError(f"rules.tsv:{lineno} non-integer position")
            mid = row["member_id"]
            if mid:
                if mid.upper() in RESERVED_EXPR_WORDS:
                    raise ValueError(f"rules.tsv:{lineno} member_id uses reserved keyword")
                member_ids.add(mid)

    if formula_path is None:
        return

    with formula_path.open(newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if reader.fieldnames != FORMULA_COLUMNS:
            raise ValueError(f"formula-rules.tsv header mismatch: {reader.fieldnames}")

        formula_group_ids: set[str] = set()
        for lineno, row in enumerate(reader, start=2):
            gid = row["group_id"]
            expr = row["expression"]
            if not gid:
                raise ValueError(f"formula-rules.tsv:{lineno} missing group_id")
            if not expr:
                raise ValueError(f"formula-rules.tsv:{lineno} missing expression")
            if gid in formula_group_ids:
                raise ValueError(
                    f"formula-rules.tsv:{lineno} duplicate group_id {gid!r}"
                )
            formula_group_ids.add(gid)

            tokens = [
                t for t in re.findall(r"\b[A-Za-z_][A-Za-z0-9_]*\b", expr)
                if t.upper() not in RESERVED_EXPR_WORDS
            ]
            seen_in_expr: set[str] = set()
            for token in tokens:
                if token not in member_ids:
                    raise ValueError(
                        f"formula-rules.tsv:{lineno} unknown member_id {token!r}"
                    )
                if token in seen_in_expr:
                    raise ValueError(
                        f"formula-rules.tsv:{lineno} member_id {token!r} used more than once"
                    )
                seen_in_expr.add(token)


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main() -> int:
    parser = argparse.ArgumentParser(
        description="Convert Stanford HIVDB ASI XML to ResPro TSV artifacts"
    )
    parser.add_argument(
        "--source-url",
        default=HIVDB_LATEST_URL,
        help="URL to HIVDB_latest.xml or a specific HIVDB_X.Y.xml (default: latest)",
    )
    parser.add_argument(
        "--output-dir",
        default=str(Path(__file__).resolve().parent.parent / "output"),
        help="Directory to write output artifacts (default: ../output relative to script)",
    )
    args = parser.parse_args()

    eprint(f"Downloading {args.source_url} …")
    source_info, xml_bytes = resolve_source(args.source_url)
    eprint(f"Source date: {source_info.source_date}")

    eprint(f"Fetching drug map from {HIVFACTS_DRUGS_URL} …")
    drug_map = load_drug_map(http_get_text(HIVFACTS_DRUGS_URL))

    eprint(f"Fetching gene data from {HIVFACTS_GENES_URL} …")
    gene_data = load_gene_data(http_get_text(HIVFACTS_GENES_URL))

    convert(
        source_info=source_info,
        xml_bytes=xml_bytes,
        drug_map=drug_map,
        gene_data=gene_data,
        output_dir=Path(args.output_dir),
    )

    eprint("Done.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
