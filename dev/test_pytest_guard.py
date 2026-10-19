"""The rules of ``dev/pytest_guard.py``: every check it makes on a run (the hook setting, the
imported modules, the skipped tests and their reasons) refuses its defect on a synthetic record;
real runs pass when clean (an unrelated plugin's warning on stderr included) and fail on one
unexpected skip, of a test or of a whole module while collected.

Run by ``make check-ci``, or directly (the end-to-end controls need the development build):

    python -I -B dev/test_pytest_guard.py
"""

from __future__ import annotations

import contextlib
import io
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import pytest_guard as guard

OUTSIDE = ["/elsewhere/_netgraph_core.so", "/elsewhere/netgraph_core/__init__.py"]
REPLAY = guard.EXPECTED_OFF_SKIPPED_PREFIX + "x]"
REASON = "Skipped: needs the private hooks (NETGRAPH_CORE_PRIVATE_HOOKS)"


def off_run(**change):
    """A correct hooks-off record, with `change` applied."""
    skipped = {name: REASON for name in guard.EXPECTED_OFF_SKIPS}
    skipped[REPLAY] = REASON
    record = dict(
        has_hooks=False,
        modules=list(OUTSIDE),
        seen=set(skipped) | {"tests/py/test_units.py::test_ok"},
        skipped=skipped,
    )
    record.update(change)
    return record


class Checks(unittest.TestCase):
    def test_hooks_on(self):
        ok = dict(has_hooks=True, modules=OUTSIDE, seen={"a"}, skipped={})
        self.assertEqual(guard.problems("on", **ok), [])
        self.assertEqual(
            guard.problems("on", **{**ok, "skipped": {"a": REASON}}),
            ["1 tests skipped in a hooks-on run: ['a']"],
        )
        self.assertEqual(
            guard.problems("on", **{**ok, "has_hooks": False}),
            ["the extension lacks the private hooks (expected the editable build)"],
        )
        self.assertEqual(
            guard.problems("on", **{**ok, "modules": [None, OUTSIDE[1]]}),
            [f"netgraph_core or _netgraph_core was not imported: {[None, OUTSIDE[1]]}"],
        )

    def test_hooks_off(self):
        self.assertEqual(guard.problems("off", **off_run()), [])
        unexpected = off_run()
        unexpected["skipped"] = {**unexpected["skipped"], "tests/py/t.py::u": REASON}
        self.assertEqual(
            guard.problems("off", **unexpected),
            ["skipped, not a known hook-only test: tests/py/t.py::u"],
        )
        missing = off_run()
        first = guard.EXPECTED_OFF_SKIPS[0]
        del missing["skipped"][first]
        self.assertEqual(
            guard.problems("off", **missing), [f"expected to skip, did not: {first}"]
        )
        other = off_run()
        other["skipped"][first] = "Skipped: something else"
        self.assertEqual(
            guard.problems("off", **other),
            [f"skipped for another reason: {first}: Skipped: something else"],
        )
        no_replay = off_run()
        del no_replay["skipped"][REPLAY]
        no_replay["seen"].discard(REPLAY)
        self.assertEqual(
            guard.problems("off", **no_replay),
            [f"no test {guard.EXPECTED_OFF_SKIPPED_PREFIX}...] ran"],
        )
        self.assertEqual(
            guard.problems("off", **off_run(has_hooks=True)),
            [
                "the extension carries the private hooks (expected a wheel or sdist build)"
            ],
        )
        inside = [str(guard.ROOT / "python" / "_netgraph_core.so"), OUTSIDE[1]]
        self.assertEqual(
            guard.problems("off", **off_run(modules=inside)),
            [f"imported from the source tree, not an installation: {inside[:1]}"],
        )

    def test_usage(self):
        with contextlib.redirect_stderr(io.StringIO()):
            self.assertEqual(guard.main([]), 2)
            self.assertEqual(guard.main(["both"]), 2)


def development_build():
    try:
        import _netgraph_core
    except ImportError:
        return False
    return hasattr(_netgraph_core, "_fail_native_allocation")


HEADER = "import _netgraph_core\nimport netgraph_core\nimport pytest\n\n"
PASSING = HEADER + "def test_ok():\n    assert True\n"
# A plugin of the environment that writes a warning when loaded (as an installed pytest plugin
# may): stderr of a clean run is not the guard's.
UNRELATED_PLUGIN = "import sys\n\nsys.stderr.write('an unrelated plugin warns\\n')\n"


@unittest.skipUnless(development_build(), "needs the development build (private hooks)")
class Runs(unittest.TestCase):
    """Real runs of the guard in hooks-on mode on small test files: clean, then one test that
    skips (as a test that silently needs a hook would on Windows), then one module skipped while
    collected (as a module-level platform skip would)."""

    def run_guard(self, files):
        """The guard's status and stderr on a directory of `files` (name: text), with the plugin
        guard_unrelated_plugin loaded. The directory is the run's rootdir, so a node id starts at the
        module's file name wherever the temporary directory is. Each test module's name is used
        once in this file: the runs share one interpreter, which keeps an imported test module."""
        with tempfile.TemporaryDirectory() as tmp:
            for name, text in {
                **files,
                "guard_unrelated_plugin.py": UNRELATED_PLUGIN,
            }.items():
                (Path(tmp) / name).write_text(text)
            err = io.StringIO()
            sys.path.insert(0, tmp)
            try:
                with (
                    contextlib.redirect_stdout(io.StringIO()),
                    contextlib.redirect_stderr(err),
                ):
                    status = guard.main(
                        ["on", tmp, "--rootdir", tmp, "-q", "-p", "no:cacheprovider"]
                        + ["-p", "guard_unrelated_plugin"]
                    )
            finally:
                sys.path.remove(tmp)
                sys.modules.pop("guard_unrelated_plugin", None)
            return status, err.getvalue()

    def guard_lines(self, err):
        return [ln for ln in err.splitlines() if ln.startswith("dev/pytest_guard.py")]

    def test_a_clean_run_passes(self):
        status, err = self.run_guard(
            {"test_guard_clean.py": PASSING, "test_guard_other.py": PASSING}
        )
        self.assertEqual((status, self.guard_lines(err)), (0, []))
        self.assertIn("an unrelated plugin warns", err)

    def test_an_unexpected_skip_fails(self):
        body = (
            HEADER
            + "def test_ok():\n    pass\n\n\ndef test_skips():\n    pytest.skip('planted')\n"
        )
        status, err = self.run_guard({"test_guard_skip.py": body})
        self.assertEqual(status, 1)
        self.assertEqual(len(self.guard_lines(err)), 1)
        self.assertIn("1 tests skipped in a hooks-on run", err)
        self.assertIn("test_guard_skip.py::test_skips", err)

    def test_a_module_skipped_while_collected_fails(self):
        module = "import pytest\n\npytest.skip('planted', allow_module_level=True)\n"
        status, err = self.run_guard(
            {"test_guard_kept.py": PASSING, "test_guard_platform.py": module}
        )
        self.assertEqual(status, 1)
        self.assertEqual(len(self.guard_lines(err)), 1)
        self.assertIn("1 tests skipped in a hooks-on run", err)
        self.assertIn("['test_guard_platform.py']", err)


if __name__ == "__main__":
    unittest.main()
