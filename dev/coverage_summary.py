#!/usr/bin/env python3

"""
Combine Python (coverage.py) and C++ (gcovr) Cobertura XML reports
into a single compact text summary similar to pytest-cov output.

Usage:
  python dev/coverage_summary.py [PY_XML] [CPP_XML] [--html OUT_HTML]

Defaults:
  PY_XML  = build/coverage/coverage-python.xml
  CPP_XML = build/coverage/coverage-cpp.xml
"""

from __future__ import annotations

import html
import os
import sys
import xml.etree.ElementTree as ET
from typing import Optional, Tuple


def parse_cobertura_totals(
    xml_path: str,
) -> Tuple[int, int, Optional[int], Optional[int]]:
    """Return (lines_valid, lines_covered, branches_valid, branches_covered).

    Missing values are returned as None. If the file is missing, returns zeros/None.
    """
    if not os.path.exists(xml_path):
        return 0, 0, None, None

    try:
        root = ET.parse(xml_path).getroot()
    except Exception:
        return 0, 0, None, None

    def _to_int(value: Optional[str]) -> Optional[int]:
        if value is None:
            return None
        try:
            # Some tools may write floats; coerce safely
            return int(float(value))
        except Exception:
            return None

    lines_valid = _to_int(root.attrib.get("lines-valid")) or 0
    lines_covered = _to_int(root.attrib.get("lines-covered")) or 0
    branches_valid = _to_int(root.attrib.get("branches-valid"))
    branches_covered = _to_int(root.attrib.get("branches-covered"))

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

    # Branch totals only if present in both
    have_branches = (py_br_valid is not None and py_br_cov is not None) or (
        cpp_br_valid is not None and cpp_br_cov is not None
    )
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


def main() -> int:
    py_xml = (
        sys.argv[1]
        if len(sys.argv) > 1 and not sys.argv[1].startswith("--")
        else "build/coverage/coverage-python.xml"
    )
    cpp_xml = (
        sys.argv[2]
        if len(sys.argv) > 2 and not sys.argv[2].startswith("--")
        else "build/coverage/coverage-cpp.xml"
    )

    html_out: Optional[str] = None
    for arg in sys.argv[1:]:
        if arg.startswith("--html="):
            html_out = arg.split("=", 1)[1]

    py_totals = parse_cobertura_totals(py_xml)
    cpp_totals = parse_cobertura_totals(cpp_xml)
    text = format_text_table(py_totals, cpp_totals)
    print(text)
    if html_out:
        os.makedirs(os.path.dirname(html_out), exist_ok=True)
        with open(html_out, "w", encoding="utf-8") as f:
            f.write(format_html_table(text))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
