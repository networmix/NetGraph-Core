"""The rules of ``dev/coverage_summary.py``: what it accepts, every report it refuses, and its
exit status, on small synthetic Cobertura reports.

Nothing is built. Run first by ``make cov``, or directly:

    python -I -B dev/test_coverage_summary.py
"""

from __future__ import annotations

import contextlib
import io
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import coverage_summary as cs


def report(**totals: str) -> str:
    attrs = "".join(f' {k.replace("_", "-")}="{v}"' for k, v in totals.items())
    return (
        f'<?xml version="1.0" ?>\n<coverage{attrs} version="7"><packages/></coverage>\n'
    )


class Summary(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.dir = Path(self.tmp.name)

    def tearDown(self):
        self.tmp.cleanup()

    def write(self, name: str, text: str) -> str:
        path = self.dir / name
        path.write_text(text)
        return str(path)

    def main(self, *argv: str) -> tuple[int, str, str]:
        out, err = io.StringIO(), io.StringIO()
        with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
            status = cs.main(list(argv))
        return status, out.getvalue(), err.getvalue()

    def test_totals_read_exactly(self):
        for totals, want in (
            (dict(lines_valid="10", lines_covered="7"), (10, 7, None, None)),
            (
                dict(
                    lines_valid="12.0",
                    lines_covered="0012",
                    branches_valid="4.00",
                    branches_covered="0",
                ),
                (12, 12, 4, 0),
            ),
            # the largest count, and one below: compared exactly, not as floats
            (
                dict(lines_valid=str(2**63 - 1), lines_covered=str(2**63 - 2)),
                (2**63 - 1, 2**63 - 2, None, None),
            ),
        ):
            path = self.write("r.xml", report(**totals))
            self.assertEqual(cs.parse_cobertura_totals(path), want, totals)

    def test_refusals(self):
        good = dict(lines_valid="10", lines_covered="7")
        for name, text, why in (
            ("missing", None, "cannot read the report"),
            ("not xml", "<coverage", "cannot read the report"),
            ("other root", "<report/>", "not a Cobertura report (root <report>)"),
            ("no lines-valid", report(lines_covered="1"), "no lines-valid total"),
            ("no lines-covered", report(lines_valid="1"), "no lines-covered total"),
            (
                "zero lines",
                report(lines_valid="0", lines_covered="0"),
                "counts no line",
            ),
            (
                "lines above",
                report(lines_valid="7", lines_covered="8"),
                "8 lines covered of 7",
            ),
            (
                "branches above",
                report(**good, branches_valid="2", branches_covered="3"),
                "3 branches covered of 2",
            ),
            (
                "one branch total",
                report(**good, branches_valid="2"),
                "must both be present or both absent",
            ),
            (
                "the other branch total",
                report(**good, branches_covered="0"),
                "must both be present or both absent",
            ),
            (
                "above the largest count",
                report(lines_valid=str(2**63), lines_covered="1"),
                "above the largest count",
            ),
            *(
                (
                    f"spelling {value}",
                    report(lines_valid=value, lines_covered="1"),
                    f"lines-valid={value!r} is not a count",
                )
                for value in ("-1", "+1", "1.5", "1e3", "nan", "inf", "", " 1", "0x10")
            ),
        ):
            path = (
                str(self.dir / "absent.xml")
                if text is None
                else self.write(f"{len(name)}.xml", text)
            )
            with self.assertRaises(cs.ReportError, msg=name) as ctx:
                cs.parse_cobertura_totals(path)
            self.assertIn(why, str(ctx.exception), name)
            self.assertTrue(str(ctx.exception).startswith(path), name)

    def test_exit_status_and_outputs(self):
        py = self.write("py.xml", report(lines_valid="10", lines_covered="5"))
        cpp = self.write("cpp.xml", report(lines_valid="30", lines_covered="29"))
        bad = self.write("bad.xml", report(lines_valid="0", lines_covered="0"))
        html = self.dir / "sub" / "c.html"
        status, out, err = self.main(py, cpp, f"--html={html}")
        self.assertEqual((status, err), (0, ""))
        rows = out.splitlines()
        self.assertEqual(rows[3].split(), ["python", "10", "5", "50%"])
        self.assertEqual(rows[4].split(), ["c++", "30", "1", "97%"])
        self.assertEqual(rows[6].split(), ["TOTAL", "40", "6", "85%"])
        self.assertIn("TOTAL", html.read_text())
        # branch columns when one report has branch totals; the other counts 0 of 0
        branched = self.write(
            "br.xml",
            report(
                lines_valid="10",
                lines_covered="5",
                branches_valid="4",
                branches_covered="1",
            ),
        )
        status, out, err = self.main(branched, cpp)
        self.assertEqual(status, 0, err)
        self.assertEqual(out.splitlines()[6].split()[-3:], ["4", "3", "25%"])
        # any refused report: exit 1, the reason on stderr, no table
        for argv in ((bad, cpp), (py, bad), (py, str(self.dir / "absent.xml"))):
            status, out, err = self.main(*argv)
            self.assertEqual((status, out), (1, ""), argv)
            self.assertTrue(err.startswith("coverage_summary.py: "), argv)
        # an HTML output that cannot be written: exit 1
        blocker = self.write("file", "")
        status, out, err = self.main(py, cpp, f"--html={blocker}/c.html")
        self.assertEqual(status, 1)
        self.assertIn(f"cannot write {blocker}/c.html", err)


if __name__ == "__main__":
    unittest.main()
