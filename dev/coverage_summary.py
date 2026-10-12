#!/usr/bin/env python3

"""
Combine Python (coverage.py) and C++ (gcovr) Cobertura XML reports
into a single compact text summary similar to pytest-cov output.

Usage:
  python dev/coverage_summary.py [PY_XML] [CPP_XML] [--html OUT_HTML]

Defaults:
  PY_XML  = build/coverage/coverage-python.xml
  CPP_XML = build/coverage/coverage-cpp.xml

A missing, unreadable or malformed report is an error (exit 1), never zero coverage, and so is a
report that counts no line, which would otherwise show as 100%: gcovr writes one, and exits 0,
when it reads no profile data (the profiles are not where it looks, or the filters match no
file), and coverage.py when the files it measured hold no statement. A report carries both
branch totals or neither, and covers at most what it counts. Each total is a count written in
decimal digits, optionally with a zero fraction ("12" or "12.0"), and at most 2**63 - 1 (far
beyond any real report); anything else (a sign, a nonzero fraction, an exponent, NaN or
infinity) is malformed. The counts are read and compared exactly.
"""

from __future__ import annotations

import argparse
import html
import os
import re
import sys
import xml.etree.ElementTree as ET
from typing import Optional, Tuple

# A count: decimal digits, optionally with a zero fraction (some tools write "12.0").
_COUNT = re.compile(r"([0-9]+)(?:\.0*)?")
_MAX_COUNT = 2**63 - 1


class ReportError(Exception):
    """A coverage report is missing, unreadable or malformed."""


def parse_cobertura_totals(
    xml_path: str,
) -> Tuple[int, int, Optional[int], Optional[int]]:
    """Return (lines_valid, lines_covered, branches_valid, branches_covered).

    The line totals are required; the branch totals are both None when the report has none.
    Raises ReportError when the file is missing or unreadable, is not a Cobertura report,
    holds a total that is not a count (the module docstring's spelling and bound), has one
    branch total without the other, covers more lines or branches than it counts, or counts
    no line (it measured nothing).
    """
    try:
        root = ET.parse(xml_path).getroot()
    except (OSError, ET.ParseError) as exc:
        raise ReportError(f"{xml_path}: cannot read the report: {exc}") from exc
    if root.tag != "coverage":
        raise ReportError(f"{xml_path}: not a Cobertura report (root <{root.tag}>)")

    def _to_int(name: str, required: bool) -> Optional[int]:
        value = root.attrib.get(name)
        if value is None:
            if required:
                raise ReportError(f"{xml_path}: no {name} total")
            return None
        match = _COUNT.fullmatch(value)
        if match is None:
            raise ReportError(
                f"{xml_path}: {name}={value!r} is not a count (a nonnegative integer)"
            )
        digits = match.group(1).lstrip("0") or "0"
        if len(digits) > len(str(_MAX_COUNT)) or int(digits) > _MAX_COUNT:
            raise ReportError(
                f"{xml_path}: {name}={value!r} is above the largest count, 2**63 - 1"
            )
        return int(digits)

    lines_valid = _to_int("lines-valid", True)
    lines_covered = _to_int("lines-covered", True)
    branches_valid = _to_int("branches-valid", False)
    branches_covered = _to_int("branches-covered", False)
    assert lines_valid is not None and lines_covered is not None
    if lines_covered > lines_valid:
        raise ReportError(f"{xml_path}: {lines_covered} lines covered of {lines_valid}")
    if lines_valid == 0:
        raise ReportError(
            f"{xml_path}: the report counts no line (no coverage data was read, or the"
            " measured files hold no statement)"
        )
    if branches_valid is None or branches_covered is None:
        if branches_valid is not None or branches_covered is not None:
            raise ReportError(
                f"{xml_path}: branches-valid and branches-covered must both be present or both absent"
            )
    elif branches_covered > branches_valid:
        raise ReportError(
            f"{xml_path}: {branches_covered} branches covered of {branches_valid}"
        )
    return lines_valid, lines_covered, branches_valid, branches_covered


def format_percent(numerator: int, denominator: int) -> str:
    if denominator <= 0:
        return "100%" if numerator == 0 else "0%"
    pct = 100.0 * float(numerator) / float(denominator)
    return f"{pct:.0f}%"


def format_text_table(
    py_totals: Tuple[int, int, Optional[int], Optional[int]],
    cpp_totals: Tuple[int, int, Optional[int], Optional[int]],
) -> str:
    py_lines_valid, py_lines_cov, py_br_valid, py_br_cov = py_totals
    cpp_lines_valid, cpp_lines_cov, cpp_br_valid, cpp_br_cov = cpp_totals

    total_lines_valid = py_lines_valid + cpp_lines_valid
    total_lines_cov = py_lines_cov + cpp_lines_cov

    # Branch columns when either report has branch totals; a report without them counts 0 of 0.
    have_branches = py_br_valid is not None or cpp_br_valid is not None
    total_br_valid = (py_br_valid or 0) + (cpp_br_valid or 0) if have_branches else 0
    total_br_cov = (py_br_cov or 0) + (cpp_br_cov or 0) if have_branches else 0

    header = "================== combined coverage =================="
    if have_branches:
        lines = [
            header,
            f"{'Name':<10}{'Stmts':>8}{'Miss':>8}{'Cover':>8}{'Br':>8}{'BrMiss':>10}{'BrCover':>10}",
            "-" * 62,
            f"{'python':<10}{py_lines_valid:>8}{(py_lines_valid - py_lines_cov):>8}{format_percent(py_lines_cov, py_lines_valid):>8}{(py_br_valid or 0):>8}{((py_br_valid or 0) - (py_br_cov or 0)):>10}{format_percent((py_br_cov or 0), (py_br_valid or 0)):>10}",
            f"{'c++':<10}{cpp_lines_valid:>8}{(cpp_lines_valid - cpp_lines_cov):>8}{format_percent(cpp_lines_cov, cpp_lines_valid):>8}{(cpp_br_valid or 0):>8}{((cpp_br_valid or 0) - (cpp_br_cov or 0)):>10}{format_percent((cpp_br_cov or 0), (cpp_br_valid or 0)):>10}",
            "-" * 62,
            f"{'TOTAL':<10}{total_lines_valid:>8}{(total_lines_valid - total_lines_cov):>8}{format_percent(total_lines_cov, total_lines_valid):>8}{total_br_valid:>8}{(total_br_valid - total_br_cov):>10}{format_percent(total_br_cov, total_br_valid):>10}",
        ]
    else:
        lines = [
            header,
            f"{'Name':<10}{'Stmts':>8}{'Miss':>8}{'Cover':>8}",
            "-" * 42,
            f"{'python':<10}{py_lines_valid:>8}{(py_lines_valid - py_lines_cov):>8}{format_percent(py_lines_cov, py_lines_valid):>8}",
            f"{'c++':<10}{cpp_lines_valid:>8}{(cpp_lines_valid - cpp_lines_cov):>8}{format_percent(cpp_lines_cov, cpp_lines_valid):>8}",
            "-" * 42,
            f"{'TOTAL':<10}{total_lines_valid:>8}{(total_lines_valid - total_lines_cov):>8}{format_percent(total_lines_cov, total_lines_valid):>8}",
        ]
    return "\n".join(lines)


def format_html_table(text_table: str, title: str = "Combined Coverage") -> str:
    pre = html.escape(text_table)
    return f"""<!doctype html><html><head><meta charset='utf-8'><title>{html.escape(title)}</title>
<style>body{{font-family: ui-monospace, SFMono-Regular, Menlo, monospace; padding: 16px; white-space: pre;}}</style>
</head><body><h2>{html.escape(title)}</h2><pre>{pre}</pre></body></html>"""


def main(argv: Optional[list[str]] = None) -> int:
    parser = argparse.ArgumentParser(
        description="Combine the Python and C++ Cobertura reports into one summary."
    )
    parser.add_argument(
        "py_xml", nargs="?", default="build/coverage/coverage-python.xml"
    )
    parser.add_argument("cpp_xml", nargs="?", default="build/coverage/coverage-cpp.xml")
    parser.add_argument(
        "--html", metavar="OUT_HTML", help="also write the table as HTML"
    )
    args = parser.parse_args(argv)

    try:
        py_totals = parse_cobertura_totals(args.py_xml)
        cpp_totals = parse_cobertura_totals(args.cpp_xml)
    except ReportError as exc:
        print(f"coverage_summary.py: {exc}", file=sys.stderr)
        return 1
    text = format_text_table(py_totals, cpp_totals)
    print(text)
    html_out: Optional[str] = args.html
    if html_out:
        out_dir = os.path.dirname(html_out)
        try:
            if out_dir:
                os.makedirs(out_dir, exist_ok=True)
            with open(html_out, "w", encoding="utf-8") as f:
                f.write(format_html_table(text))
        except OSError as exc:
            print(
                f"coverage_summary.py: cannot write {html_out}: {exc}", file=sys.stderr
            )
            return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
