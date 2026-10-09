"""Text-layer parsers for the three WHO influenza antiviral marker tables.

Parses the WHO PDF marker tables into intermediate row dictionaries that
``convert.py`` normalizes into ResPro rules:

- human NAI marker table (neuraminidase inhibitor markers, human influenza)
- avian NAI marker table (neuraminidase inhibitor markers, avian influenza)
- PA marker table (baloxavir markers, polymerase acidic protein)

Parsing strategy
----------------
All three PDFs are digital (text layer present); scanned PDFs are rejected.
Internal table grid lines vary between pages, so cells are mapped to logical
columns by x-center against per-PDF canonical column bounds.

Subtype/section labels live in a sparse leftmost column. Labels are collected
from two sources and merged:

1. merged grid cells inside the label column (precise y-span),
2. free text lines in the label column (no grid cell; span inferred).

Free-text label blocks span from the previous label's block end to the next
label's block start. This handles both top-aligned labels (human table) and
vertically centered labels spanning many rows (avian ``H5N1`` block).
Rows outside any label block carry the previous label forward, across pages.

Comments (human table only) share the leftmost column with subtype labels.
They are extracted as word-grouped lines and assigned to the grid row that
contains each line's vertical center, then joined per row in reading order.
"""

from __future__ import annotations

import re
from dataclasses import dataclass


# pdfplumber is imported lazily so the pure text-parsing helpers are usable
# without the dependency (e.g. in unit tests).
pdfplumber = None


def _pdfplumber():
    global pdfplumber
    if pdfplumber is None:
        import pdfplumber as _pdfplumber_module
        pdfplumber = _pdfplumber_module
    return pdfplumber

MONTHS = {
    "January": 1, "February": 2, "March": 3, "April": 4, "May": 5, "June": 6,
    "July": 7, "August": 8, "September": 9, "October": 10, "November": 11,
    "December": 12,
}

FOOTER_DATE_RE = re.compile(
    r"Last updated(?: on)?\s+(\d{1,2})\s+(January|February|March|April|May|June|"
    r"July|August|September|October|November|December)\s+(\d{4})"
)

# Canonical logical-column x-bounds, stable across the pages of each PDF.
HUMAN_X_BOUNDS: list[tuple[str, float, float]] = [
    ("leftmost", 36.0, 113.0),
    ("substitution", 113.0, 220.0),
    ("n2_numbering", 220.0, 291.0),
    ("oseltamivir", 291.0, 384.0),
    ("zanamivir", 384.0, 477.0),
    ("peramivir", 477.0, 576.0),
    ("laninamivir", 576.0, 655.0),
    ("origin", 655.0, 733.0),
    ("references", 733.0, 806.0),
]

AVIAN_X_BOUNDS: list[tuple[str, float, float]] = [
    ("subtype", 22.0, 78.0),
    ("substitution", 78.0, 177.6),
    ("n2_numbering", 177.6, 247.4),
    ("oseltamivir", 247.4, 347.4),
    ("zanamivir", 347.4, 447.5),
    ("peramivir", 447.5, 542.8),
    ("laninamivir", 542.8, 621.2),
    ("origin", 621.2, 715.2),
    ("references", 715.2, 780.0),
]

PA_X_BOUNDS: list[tuple[str, float, float]] = [
    ("subtype", 29.4, 105.8),
    ("substitution", 105.8, 219.7),
    ("fold_change", 219.7, 331.4),
    ("classification", 331.4, 414.8),
    ("origin", 414.8, 608.8),
    ("references", 608.8, 810.6),
]

# Leftmost label-column name for each table.
HUMAN_LABEL_COL = ("leftmost", 36.0, 113.0)
AVIAN_LABEL_COL = ("subtype", 22.0, 78.0)
PA_LABEL_COL = ("subtype", 29.4, 105.8)

# Subtype/section label patterns. Everything else in the label column is
# either a section header (matched here too) or a comment (human table).
HUMAN_LABEL_RE = re.compile(r"^(Type [AB]|A\([HN]\d[HN]\d\)(pdm09|v)?)$")
AVIAN_LABEL_RE = re.compile(r"^(Group [12]|H5N1|H7N9|N\d)$")
PA_LABEL_RE = re.compile(r"^(A\([HN]\d[HN]\d\)(pdm09|v)?|B)$")

# Header marker in the substitution column.
HEADER_SUBSTITUTION = "Amino acid substitution"

NAI_DRUG_COLUMNS = ("oseltamivir", "zanamivir", "peramivir", "laninamivir")


@dataclass
class MarkerRow:
    """One intermediate data row of a WHO marker table."""

    page: int
    subtype: str = ""
    substitution: str = ""
    n2_numbering: str = ""
    oseltamivir: str = ""
    zanamivir: str = ""
    peramivir: str = ""
    laninamivir: str = ""
    origin: str = ""
    references: str = ""
    comments: str = ""
    # PA-only columns.
    fold_change: str = ""
    classification: str = ""


@dataclass
class ParsedTable:
    """Parsed marker table: footer date plus data rows."""

    footer_date: str  # YYYY-MM-DD
    rows: list[MarkerRow]


def col_for_x(x: float, bounds: list[tuple[str, float, float]]) -> str:
    """Map an x coordinate to its logical column name."""
    for name, lo, hi in bounds:
        if lo <= x < hi:
            return name
    best, bestd = bounds[0][0], 1e9
    for name, lo, hi in bounds:
        d = min(abs(x - lo), abs(x - hi))
        if d < bestd:
            bestd, best = d, name
    return best


def _check_text_layer(pdf: pdfplumber.Pdf, path: str) -> None:
    empty = [i + 1 for i, p in enumerate(pdf.pages)
             if not (p.extract_text() or "").strip()]
    if empty:
        raise ValueError(
            f"{path}: no extractable text layer on pages {empty}; "
            "scanned PDFs are not supported"
        )


def footer_date(pdf: pdfplumber.Pdf) -> str:
    """Extract the 'Last updated' footer date as YYYY-MM-DD."""
    for page in pdf.pages:
        text = page.extract_text() or ""
        m = FOOTER_DATE_RE.search(text)
        if m:
            day, month, year = m.group(1), m.group(2), m.group(3)
            return f"{year}-{MONTHS[month]:02d}-{int(day):02d}"
    raise ValueError("no 'Last updated' footer date found in PDF")


def _grid_rows(
    page: pdfplumber.Page, bounds: list[tuple[str, float, float]]
) -> tuple[list[tuple[dict, tuple[float, float, float, float] | None]], object]:
    """Extract per-row logical cells and row bboxes from the page's table."""
    tables = page.find_tables(
        table_settings={"vertical_strategy": "lines", "horizontal_strategy": "lines"}
    )
    if not tables:
        return [], None
    table = tables[0]
    x_edges = sorted(
        {round(c[0], 1) for c in table.cells} | {round(c[2], 1) for c in table.cells}
    )
    centers = [(x_edges[i] + x_edges[i + 1]) / 2 for i in range(len(x_edges) - 1)]
    grid_to_logical = [col_for_x(c, bounds) for c in centers]

    rows: list[tuple[dict, tuple | None]] = []
    data = table.extract()
    for ri, row in enumerate(data):
        logical: dict[str, str] = {}
        for gi, cell in enumerate(row):
            if gi >= len(grid_to_logical) or not cell or not cell.strip():
                continue
            logical.setdefault(grid_to_logical[gi], cell.strip())
        bbox = None
        if ri < len(table.rows):
            bbox = table.rows[ri].bbox
        rows.append((logical, bbox))
    return rows, table


def _collect_label_events(
    page: pdfplumber.Page,
    page_idx: int,
    table: object,
    label_col: tuple[str, float, float],
    label_re: re.Pattern[str],
) -> list[dict]:
    """Collect subtype/section label events on one page.

    Returns events sorted by y: {label, page, top, bottom, kind} where kind is
    'cell' (merged grid cell with precise span) or 'free' (text line).
    """
    _, lo, hi = label_col
    _, table_top, _, table_bottom = table.bbox
    events: list[dict] = []

    # 1. Merged grid cells inside the label column.
    for c in table.cells:
        if c[0] < lo - 2 or c[2] > hi + 2:
            continue
        text = (page.crop(c).extract_text() or "").strip()
        if not text:
            continue
        matches = [m.group(0) for m in label_re.finditer(text)]
        if not matches:
            continue
        # A cell may contain a section header plus a subtype label; the
        # subtype label is the more specific (later) match.
        events.append(
            {"label": matches[-1], "page": page_idx,
             "top": c[1], "bottom": c[3], "kind": "cell"}
        )

    # 2. Free text lines in the label column.
    words = [
        w for w in page.extract_words()
        if lo <= w["x0"] < hi and table_top <= w["top"] < table_bottom
    ]
    for ln in _group_words_into_lines(words):
        text = _line_text(ln).strip()
        if not label_re.match(text):
            continue
        top = min(w["top"] for w in ln)
        bottom = max(w["bottom"] for w in ln)
        # Skip free-text lines already covered by a cell event with the same
        # label.
        duplicate = any(
            ev["label"] == text.strip() and ev["kind"] == "cell"
            and ev["top"] - 2 <= top and bottom <= ev["bottom"] + 2
            for ev in events
        )
        if not duplicate:
            events.append(
                {"label": text.strip(), "page": page_idx,
                 "top": top, "bottom": bottom, "kind": "free"}
            )

    events.sort(key=lambda ev: (ev["page"], ev["top"]))
    return events


def _group_words_into_lines(words: list[dict]) -> list[list[dict]]:
    """Group words into visual lines using a small y tolerance."""
    words = sorted(words, key=lambda w: (w["top"], w["x0"]))
    lines: list[list[dict]] = []
    cur: list[dict] = []
    last_top: float | None = None
    for w in words:
        if last_top is None or abs(w["top"] - last_top) <= 3:
            cur.append(w)
        else:
            lines.append(cur)
            cur = [w]
        last_top = w["top"]
    if cur:
        lines.append(cur)
    return lines


def _line_text(line: list[dict]) -> str:
    return " ".join(w["text"] for w in sorted(line, key=lambda w: w["x0"]))


def _event_blocks(
    events: list[dict],
) -> tuple[list[tuple[int, float]], list[tuple[int, float] | None]]:
    """Compute (start, end) positions for each label event.

    Positions are (page_idx, y) tuples compared in reading order. Cell events
    use their precise span; free events span from the previous event's block
    end to the next event's block start (open-ended for the last event).
    """
    n = len(events)
    starts: list[tuple[int, float]] = []
    ends: list[tuple[int, float] | None] = []
    for i, ev in enumerate(events):
        if ev["kind"] == "cell":
            starts.append((ev["page"], ev["top"]))
            ends.append((ev["page"], ev["bottom"]))
        else:
            start = ends[i - 1] if i > 0 else (ev["page"], ev["top"] - 1.0)
            starts.append(start)
            if i + 1 < n:
                nxt = events[i + 1]
                ends.append((nxt["page"], nxt["top"]))
            else:
                ends.append(None)  # open-ended: extends to the table end
    return starts, ends


def _assign_labels(
    data_rows: list[tuple[dict, tuple | None, int]],
    events: list[dict],
) -> list[str]:
    """Assign a subtype label to each data row.

    Rows outside any label block carry the previous label forward; rows before
    the first label block get an empty label (fail-closed downstream).
    """
    starts, _ends = _event_blocks(events)
    labels: list[str] = []
    active = ""
    for _logical, bbox, page_idx in data_rows:
        if bbox is None:
            labels.append(active)
            continue
        center = (bbox[1] + bbox[3]) / 2
        # Blocks are ordered and non-overlapping, so the last event whose
        # block starts at or before this row is the containing block (or the
        # nearest preceding block, whose label carries forward through gaps).
        chosen: int | None = None
        for i in range(len(events)):
            if starts[i] <= (page_idx, center):
                chosen = i
        if chosen is not None:
            active = events[chosen]["label"]
        labels.append(active)
    return labels


def _leftmost_comment_lines(
    page: pdfplumber.Page,
    table: object,
    label_col: tuple[str, float, float],
    label_re: re.Pattern[str],
) -> list[dict]:
    """Extract non-label text lines from the human table's leftmost region.

    Returns {text, top, bottom} dicts in reading order. These are the row
    comments that share the leftmost column with the subtype labels.
    """
    _, lo, hi = label_col
    _, table_top, _, table_bottom = table.bbox
    words = [
        w for w in page.extract_words()
        if lo <= w["x0"] < hi and table_top <= w["top"] < table_bottom
    ]
    out: list[dict] = []
    for ln in _group_words_into_lines(words):
        text = _line_text(ln).strip()
        if label_re.match(text):
            continue
        out.append({
            "text": text,
            "top": min(w["top"] for w in ln),
            "bottom": max(w["bottom"] for w in ln),
        })
    return out


def _assign_comments(
    comment_lines: list[dict],
    data_bboxes: list[tuple | None],
) -> dict[int, list[str]]:
    """Assign comment lines to the data row containing each line's center."""
    assigned: dict[int, list[str]] = {}
    for line in comment_lines:
        center = (line["top"] + line["bottom"]) / 2
        best_i, best_d = None, 1e9
        for i, bbox in enumerate(data_bboxes):
            if bbox is None:
                continue
            if bbox[1] <= center <= bbox[3]:
                best_i, best_d = i, 0.0
                break
            d = min(abs(center - bbox[1]), abs(center - bbox[3]))
            if d < best_d:
                best_d, best_i = d, i
        if best_i is not None and best_d <= 30.0:
            assigned.setdefault(best_i, []).append(line["text"])
    return assigned


def _parse_nai_like(
    path: str,
    bounds: list[tuple[str, float, float]],
    label_col: tuple[str, float, float],
    label_re: re.Pattern[str],
    with_comments: bool,
) -> ParsedTable:
    """Shared parser for the human and avian NAI marker tables."""
    with _pdfplumber().open(path) as pdf:
        _check_text_layer(pdf, path)
        date = footer_date(pdf)

        page_segments: list[tuple[int, list[tuple[dict, tuple | None, int]], dict[int, list[str]]]] = []
        all_events: list[dict] = []
        for page_idx, page in enumerate(pdf.pages):
            grid_rows, table = _grid_rows(page, bounds)
            if table is None:
                continue
            events = _collect_label_events(page, page_idx, table, label_col, label_re)
            all_events.extend(events)

            # Data rows: logical cells with a substitution column value that
            # is not the header marker. Whitespace is normalized first so
            # header cells with line wraps or footnote suffixes still match.
            data_rows: list[tuple[dict, tuple | None, int]] = []
            for logical, bbox in grid_rows:
                substitution = " ".join(
                    (logical.get("substitution") or "").split()
                )
                if substitution and HEADER_SUBSTITUTION not in substitution:
                    data_rows.append((logical, bbox, page_idx))

            comments_by_row: dict[int, list[str]] = {}
            if with_comments:
                comment_lines = _leftmost_comment_lines(page, table, label_col, label_re)
                comments_by_row = _assign_comments(
                    comment_lines, [bbox for _, bbox, _ in data_rows]
                )
            page_segments.append((page_idx, data_rows, comments_by_row))

        all_events.sort(key=lambda ev: (ev["page"], ev["top"]))
        flat_rows = [r for _, seg_rows, _ in page_segments for r in seg_rows]
        labels = _assign_labels(flat_rows, all_events)

        rows: list[MarkerRow] = []
        offset = 0
        for page_idx, seg_rows, comments_by_row in page_segments:
            for i, (logical, _bbox, _page) in enumerate(seg_rows):
                rows.append(MarkerRow(
                    page=page_idx + 1,
                    subtype=labels[offset + i],
                    substitution=logical.get("substitution", ""),
                    n2_numbering=logical.get("n2_numbering", ""),
                    oseltamivir=logical.get("oseltamivir", ""),
                    zanamivir=logical.get("zanamivir", ""),
                    peramivir=logical.get("peramivir", ""),
                    laninamivir=logical.get("laninamivir", ""),
                    origin=logical.get("origin", "").replace("\n", " "),
                    references=logical.get("references", "").replace("\n", " "),
                    comments=" ".join(comments_by_row.get(i, [])),
                ))
            offset += len(seg_rows)
    return ParsedTable(footer_date=date, rows=rows)


def parse_human_nai(path: str) -> ParsedTable:
    """Parse the WHO human NAI marker table."""
    return _parse_nai_like(path, HUMAN_X_BOUNDS, HUMAN_LABEL_COL,
                           HUMAN_LABEL_RE, with_comments=True)


def parse_avian_nai(path: str) -> ParsedTable:
    """Parse the WHO avian NAI marker table."""
    return _parse_nai_like(path, AVIAN_X_BOUNDS, AVIAN_LABEL_COL,
                           AVIAN_LABEL_RE, with_comments=False)


def parse_pa(path: str) -> ParsedTable:
    """Parse the WHO PA (baloxavir) marker table."""
    with _pdfplumber().open(path) as pdf:
        _check_text_layer(pdf, path)
        date = footer_date(pdf)

        page_segments: list[tuple[int, list[tuple[dict, tuple | None, int]]]] = []
        all_events: list[dict] = []
        for page_idx, page in enumerate(pdf.pages):
            grid_rows, table = _grid_rows(page, PA_X_BOUNDS)
            if table is None:
                continue
            events = _collect_label_events(page, page_idx, table, PA_LABEL_COL,
                                           PA_LABEL_RE)
            all_events.extend(events)
            data_rows: list[tuple[dict, tuple | None, int]] = []
            for logical, bbox in grid_rows:
                substitution = " ".join(
                    (logical.get("substitution") or "").split()
                )
                if substitution and HEADER_SUBSTITUTION not in substitution:
                    data_rows.append((logical, bbox, page_idx))
            page_segments.append((page_idx, data_rows))

        all_events.sort(key=lambda ev: (ev["page"], ev["top"]))
        flat_rows = [r for _, seg_rows in page_segments for r in seg_rows]
        labels = _assign_labels(flat_rows, all_events)

        rows: list[MarkerRow] = []
        offset = 0
        for page_idx, seg_rows in page_segments:
            for i, (logical, _bbox, _page) in enumerate(seg_rows):
                rows.append(MarkerRow(
                    page=page_idx + 1,
                    subtype=labels[offset + i],
                    substitution=logical.get("substitution", ""),
                    fold_change=logical.get("fold_change", ""),
                    classification=logical.get("classification", ""),
                    origin=logical.get("origin", "").replace("\n", " "),
                    references=logical.get("references", "").replace("\n", " "),
                ))
            offset += len(seg_rows)
    return ParsedTable(footer_date=date, rows=rows)
