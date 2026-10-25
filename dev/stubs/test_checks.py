"""Negative controls of the two package checks of ``bash dev/stubs/run.sh check``.

* check_contents.py: its sdist rule equals pyproject.toml's file list; the files of a checkout
  it takes as intended are the ones git tracks or would track; and on synthetic member lists,
  every kind of file it must refuse (a stray or ignored file, a left-out developer file, a
  missing file, a second extension module, ...) is refused, while the intended lists pass.
* check_clients.py: Pyright is replaced by a small program that prints a real run's output for
  the clients and then ends abnormally in one of its runs (killed by a signal, a fatal exit
  status, output that is not its report, a report that disagrees with the status); each must
  stop the check, and the same program ending normally must pass it. The same program also
  breaks each comparison rule in turn (search paths without the installation or with the source
  tree, a diagnostic on the valid client, a required case without its mark, a diagnostic count
  differing from the marks, a control line with no diagnostic, a valid client that does not
  run): each must fail the check with that rule's message and no other.

Nothing is built or installed. Run by ``bash dev/stubs/run.sh check``, or directly:

    python -I -B dev/stubs/test_checks.py
"""

from __future__ import annotations

import contextlib
import io
import os
import re
import subprocess
import sys
import tempfile
import tomllib
import unittest
from pathlib import Path
from unittest import mock

sys.path.insert(0, str(Path(__file__).resolve().parent))

import check_clients
import check_contents

ROOT = Path(__file__).resolve().parents[2]

# A checkout's files as git lists them: the intended sdist files and every kind it leaves out.
SHIPPED = [
    ".coveragerc",
    ".gitignore",
    "CHANGELOG.md",
    "CMakeLists.txt",
    "CONTRIBUTING.md",
    "LICENSE",
    "README.md",
    "pyproject.toml",
    "bindings/python/module.cpp",
    "docs/semantics.md",
    "docs/formal/theorems/round-robin.md",
    "include/netgraph/core/graph.hpp",
    "proofs/lean-toolchain",
    "proofs/lake-manifest.json",
    "proofs/NetGraphCore.lean",
    "proofs/Checks/audit.py",
    "python/netgraph_core/__init__.py",
    "python/netgraph_core/__init__.pyi",
    "python/netgraph_core/py.typed",
    "src/graph.cpp",
    "tests/py/test_graph.py",
    "tests/data/policy_fixtures.json",
    "third_party/pybind11/LICENSE",
]
LEFT_OUT = [
    "AGENTS.md",
    "CMakeLists.txt.user",
    "CLAUDE.md",
    "Makefile",
    ".pre-commit-config.yaml",
    ".github/workflows/tests.yml",
    ".superset/config.json",
    "dev/run-checks.sh",
    "dev/proofs/README.md",
    "tests/reference/policy/README.md",
    "proofs/.lake/build/lib/lean/NetGraphCore.olean",
    "proofs/.lake/packages/mathlib/Mathlib.lean",
]
INFO = "netgraph_core-0.12.0.dist-info/"
WHEEL = [
    "_netgraph_core.cpython-314-darwin.so",
    "netgraph_core/__init__.py",
    "netgraph_core/__init__.pyi",
    "netgraph_core/py.typed",
] + [INFO + n for n in sorted(check_contents.DIST_INFO_FILES)]
SDIST = SHIPPED + ["PKG-INFO"]


class Contents(unittest.TestCase):
    def problems(self, wheel=WHEEL, sdist=SDIST, intended=SHIPPED + LEFT_OUT):
        return check_contents.problems(wheel, sdist, intended)

    def refused(self, want: str, **lists) -> None:
        got = self.problems(**lists)
        self.assertIn(want, got)

    def test_the_rule_is_pyproject_toml_s(self):
        with open(ROOT / "pyproject.toml", "rb") as f:
            sdist = tomllib.load(f)["tool"]["scikit-build"]["sdist"]
        self.assertEqual(sdist["inclusion-mode"], "explicit")
        self.assertEqual(
            sorted(p.lstrip("/") for p in sdist["include"]),
            sorted(check_contents.SDIST_TOP),
        )
        left_out = [p.lstrip("/") for p in sdist["exclude"] if p.startswith("/")]
        self.assertEqual(sorted(left_out), sorted(check_contents.SDIST_LEFT_OUT))

    def test_the_intended_lists_pass(self):
        self.assertEqual(self.problems(), [])
        self.assertEqual(
            [n for n in SHIPPED + LEFT_OUT if check_contents.shipped(n)], SHIPPED
        )

    def test_the_sdist_refuses_every_other_file(self):
        for stray in (
            ".playwright-mcp/console.log",
            "notes.txt",
            "tests/py/__pycache__/test_graph.cpython-314.pyc",
            "python/netgraph_core/_netgraph_core.cpython-314-darwin.so",
            *LEFT_OUT,
        ):
            self.refused(
                f"the sdist carries {stray}, not an intended file",
                sdist=SDIST + [stray],
            )
        # A stray file the checkout holds: only the files the rule ships are intended.
        self.refused(
            "the sdist carries notes.txt, not an intended file",
            sdist=SDIST + ["notes.txt"],
            intended=SHIPPED + LEFT_OUT + ["notes.txt"],
        )

    def test_the_sdist_refuses_a_missing_file(self):
        for name in SDIST:
            rest = [n for n in SDIST if n != name]
            self.refused(f"the sdist lacks {name}", sdist=rest)
        # Required even when the checkout lost it.
        self.refused(
            "the sdist lacks python/netgraph_core/py.typed",
            sdist=[n for n in SDIST if not n.endswith("py.typed")],
            intended=[n for n in SHIPPED if not n.endswith("py.typed")],
        )

    def test_the_wheel_refuses_every_other_file(self):
        for stray in (
            "netgraph_core/notes.txt",
            "netgraph_core/__pycache__/__init__.cpython-314.pyc",
            INFO + "entry_points.txt",
            "tests/py/test_graph.py",
        ):
            self.refused(
                f"the wheel carries {stray}, not an intended file",
                wheel=WHEEL + [stray],
            )
        for name in WHEEL:
            if name.startswith("netgraph_core/") or name.startswith(INFO):
                rest = [n for n in WHEEL if n != name]
                self.refused(f"the wheel lacks {name}", wheel=rest)
        two = WHEEL + ["_netgraph_core.cpython-313-darwin.so"]
        self.assertTrue(
            any("2 extension modules" in p for p in self.problems(wheel=two))
        )
        none = [n for n in WHEEL if not n.startswith("_netgraph_core.")]
        self.assertTrue(
            any("0 extension modules" in p for p in self.problems(wheel=none))
        )
        other = WHEEL + ["netgraph_core-0.11.0.dist-info/METADATA"]
        self.assertTrue(any("2 .dist-info" in p for p in self.problems(wheel=other)))

    def test_no_docs_module(self):
        self.refused(
            "the wheel carries netgraph_core/_docs.py",
            wheel=WHEEL + ["netgraph_core/_docs.py"],
        )
        docs = "python/netgraph_core/_docs.py"
        self.refused(
            f"the sdist carries {docs}",
            sdist=SDIST + [docs],
            wheel=WHEEL + ["netgraph_core/_docs.py"],
            intended=SHIPPED + [docs],
        )

    def test_checkout_files(self):
        """Tracked files and new files git does not ignore; not an ignored file, nor a tracked
        file deleted from the work tree."""
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)

            def git(*args: str) -> None:
                subprocess.run(
                    ["git", "-C", tmp, *args], check=True, capture_output=True
                )

            git("init", "-q")
            files = {
                ".gitignore": "__pycache__/\n",
                "src/a.cpp": "",
                "src/gone.cpp": "",
                "a dir/with space.txt": "",
            }
            for name, text in files.items():
                (root / name).parent.mkdir(parents=True, exist_ok=True)
                (root / name).write_text(text)
            git("add", "-A")
            (root / "src/gone.cpp").unlink()
            (root / "src/new.cpp").write_text("")
            (root / "src/__pycache__").mkdir()
            (root / "src/__pycache__/m.pyc").write_text("")
            self.assertEqual(
                check_contents.checkout_files(root),
                [".gitignore", "a dir/with space.txt", "src/a.cpp", "src/new.cpp"],
            )


# Stands in for Pyright: in the project directory it is given, it reports one diagnostic per
# rule of each "# expect:" marker and one for the control's planted line, as a real run does,
# and prints the search paths of FAKE_PATHS (os.pathsep-separated). FAKE_BREAK names the run
# that ends abnormally or reports wrongly (verbose, report or control) and FAKE_HOW how: "drop"
# loses the run's first diagnostic, "extra-valid" and "extra-invalid" add one on line 1 of that
# client.
FAKE_PYRIGHT = r"""
import json, os, re, signal, sys
from pathlib import Path

args = sys.argv[1:]
work = Path(args[args.index("--project") + 1])
run = "verbose" if "--verbose" in args else ("control" if work.name == "control" else "report")
how = os.environ["FAKE_HOW"] if run == os.environ["FAKE_BREAK"] else "normal"
def diag(f, n, rule):
    span = {"line": n, "character": 0}
    return {"file": str(f), "severity": "error", "rule": rule,
            "message": "a diagnostic", "range": {"start": span, "end": span}}
diags = []
for f in sorted(work.glob("*client.py")):
    for n, line in enumerate(f.read_text().splitlines()):
        m = re.search(r"#\s*expect:\s*(\w+(?:\s*,\s*\w+)*)\s*$", line)
        rules = re.split(r"\s*,\s*", m.group(1)) if m else []
        if line == "oops: str = g.num_edges()":
            rules = ["reportAssignmentType"]
        diags += [diag(f, n, rule) for rule in rules]
if how == "drop":
    del diags[0]
if how.startswith("extra-"):
    diags.append(diag(work / (how[len("extra-"):] + "_client.py"), 0, "reportExtra"))
errors = len(diags)
status = 1 if errors else 0
if run == "verbose":
    paths = os.environ["FAKE_PATHS"].split(os.pathsep)
    print("Search paths for file:///\n" + "\n".join("  " + p for p in paths))
else:
    report = {"generalDiagnostics": diags, "summary": {"errorCount": errors}}
    if how == "no-diagnostics":
        del report["generalDiagnostics"]
    if how == "no-count":
        del report["summary"]["errorCount"]
    if how == "text-count":
        report["summary"]["errorCount"] = str(errors)
    if how == "disagree":
        status = 1 - status
    print("[1, 2]" if how == "list" else
          "Pyright crashed" if how == "not-json" else json.dumps(report))
sys.stdout.flush()
if how == "signal":
    os.kill(os.getpid(), signal.SIGTERM)
sys.exit({"exit-2": 2, "exit-3": 3, "exit-4": 4}.get(how, status))
"""


class Pyright(unittest.TestCase):
    """``check_clients.main`` with Pyright replaced by FAKE_PYRIGHT (a real child process) and
    the installation by a directory of its own; the valid client's run ends as `client_run`
    says (status, stderr), and the clients come from `clients` (tests/typing by default)."""

    def setUp(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.tmp = Path(os.path.realpath(tmp.name))
        self.site = self.tmp / "site"
        (self.site / "netgraph_core").mkdir(parents=True)
        self.fake = self.tmp / "fake_pyright.py"
        self.fake.write_text(FAKE_PYRIGHT)

    def main(
        self,
        brk: str = "none",
        how: str = "normal",
        paths: list[Path] | None = None,
        clients: Path = check_clients.CLIENTS,
        client_run: tuple[int, str] = (0, ""),
    ) -> tuple[int | str | None, str]:
        real_run = subprocess.run

        def run(cmd, **kwargs):
            if cmd[1:3] == ["-m", "pyright"]:
                return real_run(
                    [sys.executable, "-I", str(self.fake), *cmd[3:]], **kwargs
                )
            return subprocess.CompletedProcess(cmd, client_run[0], "", client_run[1])

        searched = os.pathsep.join(str(p) for p in (paths or [self.site]))
        env = {"FAKE_BREAK": brk, "FAKE_HOW": how, "FAKE_PATHS": searched}
        out, err = io.StringIO(), io.StringIO()
        with (
            mock.patch.object(
                check_clients, "installation", lambda p: self.site / "netgraph_core"
            ),
            mock.patch.object(check_clients, "CLIENTS", clients),
            mock.patch.object(check_clients.subprocess, "run", run),
            mock.patch.dict(os.environ, env),
            mock.patch.object(
                sys, "argv", ["check_clients.py", "t", str(self.tmp / "w")]
            ),
            contextlib.redirect_stdout(out),
            contextlib.redirect_stderr(err),
        ):
            try:
                status = check_clients.main()
            except SystemExit as exc:
                status = exc.code
        return status, out.getvalue() + err.getvalue()

    def test_a_normal_end_passes(self):
        status, out = self.main("none", "normal")
        self.assertEqual(status, 0, out)
        self.assertIn("invalid client: ", out)

    def test_every_abnormal_end_stops_the_check(self):
        cases = {
            "signal": "pyright was killed by signal 15",
            "exit-2": "pyright failed (exit 2)",
            "exit-3": "pyright failed (exit 3)",
            "exit-4": "pyright failed (exit 4)",
            "not-json": "pyright gave no JSON (exit ",
            "list": "pyright's report is not one of a complete run",
            "no-diagnostics": "pyright's report is not one of a complete run",
            "no-count": "pyright's report is not one of a complete run",
            "text-count": "pyright's report is not one of a complete run",
            "disagree": "pyright's report is not one of a complete run",
        }
        for brk in ("verbose", "report", "control"):
            for how, message in cases.items():
                if brk == "verbose" and how not in (
                    "signal",
                    "exit-2",
                    "exit-3",
                    "exit-4",
                ):
                    continue  # the verbose run prints a log, not a report
                with self.subTest(run=brk, how=how):
                    status, out = self.main(brk, how)
                    self.assertIsInstance(status, str, out)
                    self.assertIn(message, status)

    def refused(self, problem: str, **kwargs) -> None:
        """The check fails with exactly this one problem."""
        status, out = self.main(**kwargs)
        self.assertEqual(status, 1, out)
        self.assertEqual(out, f"dev/stubs/check_clients.py [t]: {problem}\n")

    def test_the_search_paths_must_hold_the_installation(self):
        elsewhere = self.tmp / "elsewhere"
        elsewhere.mkdir()
        self.refused(
            f"Pyright's search paths do not hold {self.site} (pyright-verbose.log)",
            paths=[elsewhere],
        )

    def test_the_search_paths_must_not_hold_the_source_tree(self):
        self.refused(
            "Pyright searches the source tree's python/ directory",
            paths=[self.site, check_clients.ROOT / "python"],
        )

    def test_a_diagnostic_on_the_valid_client_fails(self):
        self.refused(
            "valid_client.py:1: 1 unexpected diagnostic(s) reportExtra",
            brk="report",
            how="extra-valid",
        )

    def test_a_required_case_without_its_mark_fails(self):
        # a copy of the clients in which the first required case's line lost its mark (and
        # so, from the fake as from Pyright, its diagnostic)
        case = check_clients.REQUIRED[0]
        clients = self.tmp / "clients"
        clients.mkdir()
        for name in ("valid_client.py", "invalid_client.py"):
            lines = (check_clients.CLIENTS / name).read_text().splitlines()
            if name == "invalid_client.py":
                lines = [
                    check_clients.EXPECT.sub("", t) if case in t else t for t in lines
                ]
            (clients / name).write_text("\n".join(lines) + "\n")
        self.refused(
            f"invalid_client.py lost the marked case {case!r}", clients=clients
        )

    def test_a_count_differing_from_the_marks_fails(self):
        # the fake's first diagnostic is the first rule of the invalid client's first mark
        lines = (check_clients.CLIENTS / "invalid_client.py").read_text().splitlines()
        line, mark = next(
            (i, m)
            for i, t in enumerate(lines, 1)
            if (m := check_clients.EXPECT.search(t))
        )
        rules = re.split(r"\s*,\s*", mark.group(1))
        n = rules.count(rules[0])
        with self.subTest("a marked diagnostic missing"):
            self.refused(
                f"invalid_client.py:{line}: {rules[0]} expected {n}, got {n - 1}",
                brk="report",
                how="drop",
            )
        with self.subTest("an unmarked diagnostic"):
            self.refused(
                "invalid_client.py:1: reportExtra expected 0, got 1",
                brk="report",
                how="extra-invalid",
            )

    def test_the_control_line_must_give_its_one_diagnostic(self):
        planted = (check_clients.CLIENTS / "valid_client.py").read_text()
        last = len(planted.splitlines()) + 1
        self.refused(
            f"control: the planted line {last} gave {{}}", brk="control", how="drop"
        )

    def test_the_valid_client_must_run(self):
        self.refused(
            "valid_client.py does not run:\nTraceback: boom\n",
            client_run=(1, "Traceback: boom\n"),
        )


if __name__ == "__main__":
    unittest.main()
