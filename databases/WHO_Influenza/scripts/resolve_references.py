#!/usr/bin/env python3
"""Resolve WHO marker-table citations to DOIs and PMIDs.

Extracts the numbered reference lists from the three WHO influenza PDFs and
resolves each citation:

- avian NAI and PA tables: DOIs are printed in the citations and are taken
  verbatim; PMIDs are looked up from DOIs via the NCBI ID converter.
- human NAI table: citations carry no DOIs; each citation is resolved against
  the Crossref bibliographic search API, then PMID via the NCBI ID converter.

Results are cached in a TSV sidecar (``reference-lookup.tsv``) keyed by
(table, citation_number). The converter consumes this cache; it never touches
the network. Re-running this script only looks up citations that are missing
or unresolved in the cache, unless ``--refresh`` is given.

All console output goes to stderr.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
import time
import unicodedata
import urllib.parse
import urllib.request
from pathlib import Path


# pdfplumber is imported lazily so the pure citation-parsing helpers are
# usable without the dependency (e.g. in unit tests).
pdfplumber = None


def _pdfplumber():
    global pdfplumber
    if pdfplumber is None:
        import pdfplumber as _pdfplumber_module
        pdfplumber = _pdfplumber_module
    return pdfplumber

LOOKUP_COLUMNS = [
    "table", "citation_number", "citation", "doi", "pmid",
    "match_method", "status", "resolved_title",
]

# A citation entry starts a line with "N. " where N is the sequential number.
CITATION_START_RE = re.compile(r"^(\d+)\.\s+(.*)$")
# Line-wrapped page ranges ("...816–\n24.") must not be mistaken for entries.
DASH_CHARS = "-–−"
YEAR_RE = re.compile(r"\b(19|20)\d{2}\b")
# DOI URLs may be followed by glued page-footer text, so do not anchor at EOL.
DOI_RE = re.compile(r"https?://(?:dx\.)?doi\.org/(10\.\S+)")
# Page footer text glues onto citation lines at page breaks (sometimes even
# into a DOI URL), so strip it wherever it appears on a line.
FOOTER_TEXT_RE = re.compile(r"Last updated.*?Page \d+ of \d+")

CROSSREF_URL = "https://api.crossref.org/works"
IDCONV_URL = "https://www.ncbi.nlm.nih.gov/pmc/utils/idconv/v1.0/"
REQUEST_DELAY_S = 0.5

# Citations shorter than this are not credible reference entries.
MIN_CITATION_LEN = 40


def eprint(*args: object) -> None:
    print(*args, file=sys.stderr)


def _norm_space(text: str) -> str:
    return " ".join(text.split())


def _clean_title(title: str) -> str:
    """Collapse whitespace and strip HTML tags from a Crossref title."""
    return _norm_space(re.sub(r"<[^>]+>", " ", title))


def extract_page_texts(path: str) -> list[str]:
    """Return the text of every page of the PDF."""
    with _pdfplumber().open(path) as pdf:
        return [(page.extract_text() or "") for page in pdf.pages]


def parse_citation_list(page_texts: list[str]) -> dict[int, str]:
    """Extract the sequentially numbered reference list from page texts.

    A new entry starts at a line ``N. text`` where N is exactly the previous
    entry's number plus one, the line is long enough to be a citation, and it
    contains a year. Lines starting with ``N. `` that violate the sequence
    (e.g. a page range wrapped after a dash: ``...1019–`` / ``23. https://…``)
    are treated as continuations of the current entry.
    """
    entries: dict[int, str] = {}
    current: list[str] = []
    current_num = 0

    def flush() -> None:
        if current and current_num:
            entries[current_num] = _norm_space(" ".join(current))

    for page_text in page_texts:
        for line in page_text.splitlines():
            line = FOOTER_TEXT_RE.sub("", line).strip()
            if not line:
                continue
            m = CITATION_START_RE.match(line)
            if m:
                num = int(m.group(1))
                rest = m.group(2)
                # A real entry line is always substantial; the year may sit on
                # a continuation line, so it is not required here. The
                # sequential rule plus the dash guard do the real work.
                plausible = len(line) >= MIN_CITATION_LEN
                sequential = num == current_num + 1
                after_dash = bool(current) and current[-1].endswith(tuple(DASH_CHARS))
                if plausible and sequential and not after_dash:
                    flush()
                    current_num = num
                    current = [rest]
                    continue
            current.append(line)
    flush()
    return entries


def citation_doi(citation: str) -> str:
    """Extract a DOI from a citation that prints its doi.org URL."""
    m = DOI_RE.search(citation)
    if m:
        doi = m.group(1).rstrip(".,;")
        if not doi.endswith("-"):
            return doi
        # The URL was line-wrapped at a hyphen: the DOI tail continues as the
        # next whitespace-delimited token. Rejoin and validate the result.
        rest = citation[m.end():].split()
        tail = rest[0].rstrip(".,;") if rest else ""
        if tail and re.fullmatch(r"[A-Za-z0-9.\-]+", tail):
            candidate = m.group(1) + tail
            if re.fullmatch(r"10\.\d{4,5}/\S+", candidate):
                return candidate
        return doi
    # Repair DOIs whose "doi.org/" URL lost its scheme or was split across a
    # page break: rejoin the text after "doi.org/" and validate as a DOI.
    idx = citation.find("doi.org/")
    if idx != -1:
        candidate = citation[idx + len("doi.org/"):].replace(" ", "")
        candidate = candidate.rstrip(".")
        if re.fullmatch(r"10\.\d{4,5}/\S+", candidate):
            return candidate
    return ""


def _http_get_json(url: str) -> dict | None:
    req = urllib.request.Request(url, headers={"User-Agent": "respro-db/1.0"})
    try:
        with urllib.request.urlopen(req, timeout=30) as resp:
            return json.loads(resp.read().decode("utf-8"))
    except Exception as exc:  # noqa: BLE001 - network errors are expected
        eprint(f"WARNING: request failed: {url}: {exc}")
        return None


def _title_tokens(text: str) -> set[str]:
    text = unicodedata.normalize("NFKD", text.lower())
    return {t for t in re.findall(r"[a-z0-9]+", text) if len(t) > 2}


def crossref_lookup(citation: str) -> tuple[str, str, str]:
    """Resolve a citation via the Crossref bibliographic search.

    Returns (doi, status, resolved_title). status is 'matched', 'ambiguous'
    or 'not-found'.
    """
    # Heuristic title extraction: text after "YYYY. " up to the journal
    # abbreviation (a short token sequence ending before "Journal Vol:Pages").
    m = YEAR_RE.search(citation)
    query = citation
    if m:
        query = citation[m.end():].strip(" .")
    params = urllib.parse.urlencode({
        "query.bibliographic": query,
        "rows": 5,
        "select": "DOI,title,issued,container-title",
        "mailto": "respro-db@users.noreply.github.com",
    })
    data = _http_get_json(f"{CROSSREF_URL}?{params}")
    if not data or data.get("status") != "ok":
        return "", "not-found", ""
    items = data.get("message", {}).get("items", [])
    if not items:
        return "", "not-found", ""

    cite_tokens = _title_tokens(citation)
    scored: list[tuple[float, str, str]] = []
    for item in items:
        doi = item.get("DOI", "")
        title = " ".join(item.get("title") or [])
        item_tokens = _title_tokens(title)
        if not item_tokens:
            continue
        overlap = len(cite_tokens & item_tokens) / max(1, len(item_tokens))
        # Year agreement bonus.
        year = (item.get("issued", {}).get("date-parts") or [[None]])[0][0]
        year_bonus = 0.15 if year and m and m.group(0) == str(year) else 0.0
        scored.append((overlap + year_bonus, doi, title))
    if not scored:
        return "", "not-found", ""
    scored.sort(reverse=True)
    best, best_doi, best_title = scored[0]
    second = scored[1][0] if len(scored) > 1 else 0.0
    best_title = _clean_title(best_title)
    if best >= 0.75 and best - second >= 0.1:
        return best_doi, "matched", best_title
    if best >= 0.5:
        return best_doi, "ambiguous", best_title
    return "", "not-found", ""


def pmid_from_doi(doi: str) -> str:
    """Look up a PMID for a DOI via the NCBI ID converter."""
    if not doi:
        return ""
    params = urllib.parse.urlencode({
        "ids": doi, "format": "json", "tool": "respro-db",
        "email": "respro-db@users.noreply.github.com",
    })
    data = _http_get_json(f"{IDCONV_URL}?{params}")
    if not data:
        return ""
    records = data.get("records", [])
    if records and records[0].get("pmid"):
        return str(records[0]["pmid"])
    return ""


def load_cache(path: Path) -> dict[tuple[str, str], dict[str, str]]:
    cache: dict[tuple[str, str], dict[str, str]] = {}
    if not path.exists():
        return cache
    lines = path.read_text(encoding="utf-8").splitlines()
    if len(lines) < 2:
        return cache
    header = lines[0].split("\t")
    for line in lines[1:]:
        values = line.split("\t")
        row = dict(zip(header, values))
        cache[(row["table"], row["citation_number"])] = row
    return cache


def write_cache(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as fh:
        fh.write("\t".join(LOOKUP_COLUMNS) + "\n")
        for row in rows:
            fh.write("\t".join(row.get(col, "") for col in LOOKUP_COLUMNS) + "\n")


def resolve_table(
    table: str,
    pdf_path: Path,
    cache: dict[tuple[str, str], dict[str, str]],
    refresh: bool,
    delay: float,
) -> list[dict[str, str]]:
    """Resolve all citations of one source table, using and updating the cache."""
    entries = parse_citation_list(extract_page_texts(str(pdf_path)))
    if not entries:
        raise ValueError(f"{pdf_path}: no citation list found")
    if 1 not in entries:
        raise ValueError(f"{pdf_path}: citation list does not start at 1")

    eprint(f"{table}: extracted {len(entries)} citations from {pdf_path.name}")
    rows: list[dict[str, str]] = []
    for num in sorted(entries):
        key = (table, str(num))
        citation = entries[num]
        cached = cache.get(key)
        # 'manual' rows are human-verified pins and always survive re-runs;
        # automated rows are reused unless --refresh.
        if cached and (
            cached.get("match_method") == "manual"
            or (not refresh
                and cached.get("status") in ("matched", "doi-in-source"))
        ):
            rows.append({**cached, "citation": citation,
                         "resolved_title": _clean_title(cached.get("resolved_title", ""))})
            continue
        doi = citation_doi(citation)
        pmid = ""
        title = ""
        if doi:
            method, status = "doi-in-source", "matched"
            pmid = pmid_from_doi(doi)
            time.sleep(delay)
        else:
            doi, status, title = crossref_lookup(citation)
            method = "crossref-bibliographic"
            time.sleep(delay)
            if doi:
                pmid = pmid_from_doi(doi)
                time.sleep(delay)
        if status == "matched":
            eprint(f"  {table} ref {num}: doi={doi or '-'} pmid={pmid or '-'}")
        else:
            eprint(f"  {table} ref {num}: status={status} citation={citation[:70]}...")
        rows.append({
            "table": table, "citation_number": str(num), "citation": citation,
            "doi": doi, "pmid": pmid, "match_method": method,
            "status": status, "resolved_title": title,
        })
    return rows


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-dir", required=True, type=Path,
                        help="directory containing the three WHO PDFs")
    parser.add_argument("--lookup-file", required=True, type=Path,
                        help="output TSV cache of citation lookups")
    parser.add_argument("--refresh", action="store_true",
                        help="re-resolve citations already cached")
    parser.add_argument("--delay", type=float, default=REQUEST_DELAY_S,
                        help="seconds between network requests")
    args = parser.parse_args()

    pdfs = {
        "human-nai": args.source_dir / "human-nai-marker-table.pdf",
        "avian-nai": args.source_dir / "avian-nai-marker-table.pdf",
        "pa": args.source_dir / "pa-marker-table.pdf",
    }
    missing = [str(p) for p in pdfs.values() if not p.exists()]
    if missing:
        eprint(f"ERROR: missing source PDFs: {missing}")
        return 2

    cache = load_cache(args.lookup_file)
    all_rows: list[dict[str, str]] = []
    for table, pdf_path in pdfs.items():
        all_rows.extend(resolve_table(table, pdf_path, cache, args.refresh, args.delay))

    write_cache(args.lookup_file, all_rows)
    resolved = sum(1 for r in all_rows if r["status"] == "matched")
    eprint(f"Written {args.lookup_file} ({len(all_rows)} citations, "
           f"{resolved} resolved).")
    return 0


if __name__ == "__main__":
    sys.exit(main())
