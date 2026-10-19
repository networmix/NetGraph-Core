"""Negative controls of generate.py's rules (its module docstring, steps 2 to 6) and of its
check of the installation.

A small synthetic module stands in for the compiled one: a raw stub as pybind11-stubgen writes
it, the module objects it describes and reviewed overrides. The stub built from them is checked
once, then each refusal of steps 2 to 5 is provoked once and ``build_stub`` must refuse with its
message. ``--write`` and ``--check`` run on that module (step 6: a changed, shortened, longer,
empty or missing stub fails the check with the difference), and the installation check refuses a
package or an extension outside site-packages and an extension with the private hooks. Nothing
is built or installed. Run by ``bash dev/stubs/run.sh check``, or directly:

    python -I -B dev/stubs/test_generate.py
"""

from __future__ import annotations

import contextlib
import io
import os
import sys
import sysconfig
import tempfile
import types
import unittest
from pathlib import Path
from unittest import mock

sys.path.insert(0, str(Path(__file__).resolve().parent))

import generate

RAW = '''"""raw"""
from __future__ import annotations
import typing
__all__: list[str] = ["Graph", "Result", "f"]
class Graph:
    def __init__(self, n: typing.Any) -> None:
        """Graph doc"""
    def __eq__(self, other: typing.Any) -> bool:
        ...
    def _hook(self) -> None:
        ...
class Result:
    def value(self) -> int:
        ...
def f(x: typing.Any, y: int = 1) -> typing.Any:
    ...
def _private() -> None:
    ...
'''

OVERRIDES = """import typing

class Graph:
    # reason: an integer
    # generated: def __init__(self, n: typing.Any) -> None:
    def __init__(self, n: typing.SupportsIndex) -> None: ...

class Result:
    # reason: built by the module only
    # generated: (none)
    def __init__(self, _no_public_constructor: typing.Never, /) -> None: ...

# reason: integers in and out
# generated: def f(x: typing.Any, y: int = 1) -> typing.Any:
def f(x: int, y: int = 1) -> int: ...
"""


EQ_DOC = "__eq__(self: Graph, other: object, /) -> bool"


def module(
    result_constructible: bool = False, eq_doc: str = EQ_DOC
) -> types.ModuleType:
    """The runtime side: __eq__ takes its operand positional-only (``eq_doc``, pybind11's
    signature line), Result has no constructor (pybind11's TypeError) unless
    ``result_constructible``."""

    class Graph:
        def __init__(self, n):
            self.n = n

        def __eq__(self, other):
            return self is other

        __eq__.__doc__ = eq_doc

    class Result:
        def __init__(self):
            if not result_constructible:
                raise TypeError("Result: No constructor defined!")

        def value(self):
            return 0

    def f(x, y=1):
        return x

    pkg = types.ModuleType("fake", "The package docstring.")
    pkg.__dict__.update(
        Graph=Graph, Result=Result, f=f, __all__=["Graph", "Result", "f"]
    )
    return pkg


class Rules(unittest.TestCase):
    def setUp(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        root = Path(tmp.name)
        saved = generate.ROOT, generate.OVERRIDES
        self.addCleanup(lambda: setattr(generate, "ROOT", saved[0]))
        self.addCleanup(lambda: setattr(generate, "OVERRIDES", saved[1]))
        generate.ROOT, generate.OVERRIDES = root, root / "overrides.pyi"

    def build(self, raw=RAW, overrides=OVERRIDES, pkg=None) -> str:
        generate.OVERRIDES.write_text(overrides)
        return generate.build_stub(raw, pkg or module())

    def refused(self, why: str, **kwargs) -> None:
        with self.assertRaises(generate.Failure) as ctx:
            self.build(**kwargs)
        self.assertIn(why, str(ctx.exception))

    def test_the_valid_stub(self):
        stub = self.build()
        for line in (
            "def __init__(self, n: typing.SupportsIndex) -> None:",  # an override
            "def __init__(self, _no_public_constructor: typing.Never, /) -> None:",  # marker
            "def __eq__(self, other: typing.Any, /) -> bool:",  # positional-only restored
            "def f(x: int, y: int=1) -> int:",  # before ruff formats it
        ):
            self.assertIn(line, stub)
        self.assertIn('"""Graph doc"""', stub)  # the generated body kept
        self.assertNotIn("_hook", stub)  # private members and names left out
        self.assertNotIn("_private", stub)

    def test_an_override_changes_annotations_only(self):
        for changed in ("def f(z: int, y: int = 1)", "def f(x: int, y: int = 2)"):
            overrides = OVERRIDES.replace("def f(x: int, y: int = 1)", changed)
            self.refused("an override changes only annotations", overrides=overrides)

    def test_the_quoted_generated_line(self):
        raw = RAW.replace("y: int = 1) -> typing.Any", "y: float = 1) -> typing.Any")
        self.refused("the generated declaration changed", raw=raw)
        no_reason = OVERRIDES.replace("# reason: integers in and out\n", "")
        self.refused("an override needs a '# reason:' line", overrides=no_reason)
        stale = (
            OVERRIDES
            + "\n# reason: gone\n# generated: def g() -> None:\ndef g() -> None: ...\n"
        )
        self.refused(
            "overrides without a generated declaration: ['g']", overrides=stale
        )

    def test_the_marker_matches_the_runtime_both_ways(self):
        pkg = module(result_constructible=True)
        self.refused(
            "Result: marked not constructible, but the module constructs it", pkg=pkg
        )
        start = OVERRIDES.index("class Result:")
        unmarked = (
            OVERRIDES[:start] + OVERRIDES[OVERRIDES.index("# reason: integers") :]
        )
        self.refused("Result: the module refuses to construct it", overrides=unmarked)

    def test_no_any_outside_its_domain(self):
        overrides = OVERRIDES.replace("y: int = 1) -> int", "y: int = 1) -> typing.Any")
        self.refused(
            "typing.Any left in the stub (add a reviewed override): f -> ...",
            overrides=overrides,
        )

    def test_positional_only_names_match(self):
        raw = RAW.replace(
            "__eq__(self, other: typing.Any)", "__eq__(self, rhs: typing.Any)"
        )
        self.refused("the module makes ['self', 'other'] positional-only", raw=raw)

    def test_the_overrides_file(self):
        f_line = "# generated: def f(x: typing.Any, y: int = 1) -> typing.Any:\n"
        f_def = "def f(x: int, y: int = 1) -> int: ..."
        overload = f"@typing.overload\n{f_def}\n"
        for overrides, why in (
            (OVERRIDES.replace(f_line, ""), "needs exactly one '# generated:' line"),
            (OVERRIDES.replace(f_line, f_line * 2), "needs exactly one '# generated:'"),
            (
                OVERRIDES.replace(
                    f"{f_def}\n", overload + "# reason: again\n" + overload
                ),
                "'# reason:'/'# generated:' belong above the first overload only",
            ),
            (
                OVERRIDES.replace(f"{f_def}\n", overload + f"{f_def}\n"),
                "repeated declarations must all be @typing.overload",
            ),
            (
                OVERRIDES.replace("class Graph:", f"{overload}\nclass Graph:"),
                "f: overloads must be adjacent",
            ),
            (
                OVERRIDES.replace("class Result:\n", "class Result:\n    x: int\n"),
                "only method overrides inside a class",
            ),
            (OVERRIDES + "x = 1\n", "unexpected statement"),
        ):
            self.refused(why, overrides=overrides)

    def test_what_an_override_replaces(self):
        overloads = RAW.replace(
            "def f(x: typing.Any, y: int = 1) -> typing.Any:\n    ...\n",
            "@typing.overload\ndef f(x: typing.Any, y: int = 1) -> typing.Any:\n    ...\n"
            "@typing.overload\ndef f(x: str) -> str:\n    ...\n",
        )
        self.refused(
            "f: generated as overloads; an override replaces one declaration",
            raw=overloads,
        )
        none = OVERRIDES.replace(
            "# generated: def __init__(self, n: typing.Any) -> None:",
            "# generated: (none)",
        )
        self.refused(
            "Graph.__init__: the override says '(none)', but Graph.__init__ is generated",
            overrides=none,
        )
        decorated = OVERRIDES.replace("def f(x: int", "@staticmethod\ndef f(x: int")
        self.refused(
            "f: decorators ['staticmethod'] differ from the generated []",
            overrides=decorated,
        )
        marker = OVERRIDES.replace(
            "_no_public_constructor: typing.Never", "_no_public_constructor: int"
        )
        self.refused(
            "Result.__init__: not generated; the only addition allowed is the constructor"
            " marker",
            overrides=marker,
        )

    def test_positional_only_signature_count(self):
        doc = (
            "__eq__(*args, **kwargs)\nOverloaded function.\n\n"
            "1. __eq__(self: Graph, other: object, /) -> bool\n\n"
            "2. __eq__(self: Graph, other: int, /) -> bool\n"
        )
        self.refused(
            "Graph.__eq__: 2 module signatures, 1 stub declarations",
            pkg=module(eq_doc=doc),
        )

    def test_no_any_in_parameters_or_attributes(self):
        raw = RAW.replace(
            "def value(self) -> int:", "def value(self, x: typing.Any) -> int:"
        )
        self.refused(
            "typing.Any left in the stub (add a reviewed override): Result.value(x)",
            raw=raw,
        )
        raw = RAW.replace("class Result:\n", "class Result:\n    y: typing.Any\n")
        self.refused(
            "typing.Any left in the stub (add a reviewed override): Result.y", raw=raw
        )
        # a module attribute, exported
        pkg = module()
        pkg.__all__ = [*pkg.__all__, "z"]
        pkg.z = 0
        self.refused(
            "typing.Any left in the stub (add a reviewed override): z",
            raw=RAW + "z: typing.Any\n",
            pkg=pkg,
        )

    def test_no_any_in_any_parameter_position(self):
        pkg = module()
        pkg.__all__ = [*pkg.__all__, "g"]
        pkg.g = lambda *args, **kwargs: 0
        for signature, where in (
            ("x: typing.Any, /", "g(x)"),  # positional-only
            ("x: typing.Any", "g(x)"),  # positional
            ("*args: typing.Any", "g(args)"),
            ("*, x: typing.Any", "g(x)"),  # keyword-only
            ("**kwargs: typing.Any", "g(kwargs)"),
        ):
            with self.subTest(signature):
                self.refused(
                    f"typing.Any left in the stub (add a reviewed override): {where}",
                    raw=RAW + f"def g({signature}) -> int:\n    ...\n",
                    pkg=pkg,
                )
        self.refused(
            "typing.Any left in the stub (add a reviewed override): g -> ...",
            raw=RAW + "def g(*args: int) -> typing.Any:\n    ...\n",
            pkg=pkg,
        )
        # typed variadic parameters are kept as generated
        stub = self.build(
            raw=RAW + "def g(*args: int, **kwargs: str) -> int:\n    ...\n", pkg=pkg
        )
        self.assertIn("def g(*args: int, **kwargs: str) -> int:", stub)

    def test_no_any_in_a_nested_class(self):
        raw = RAW.replace(
            "class Result:\n",
            "class Result:\n    class Part:\n"
            "        def h(self, x: typing.Any) -> int:\n            ...\n",
        )
        self.refused(
            "typing.Any left in the stub (add a reviewed override): Result.Part.h(x)",
            raw=raw,
        )

    def test_exports_match(self):
        self.refused(
            "_netgraph_core.g is public but not exported",
            raw=RAW + "def g() -> None:\n    ...\n",
        )
        raw = RAW.replace("def f(", "def _f(")
        self.refused("exported by netgraph_core, not generated: ['f']", raw=raw)


class CommandLine(unittest.TestCase):
    """``main`` on the synthetic module (the installation, pybind11-stubgen and ruff replaced):
    ``--write`` writes the stub, ``--check`` passes on an equal stub and fails on any other
    with the difference."""

    def setUp(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        root = Path(tmp.name)
        names = ("ROOT", "STUB", "OVERRIDES", "installed_package", "run_stubgen")
        saved = {name: getattr(generate, name) for name in names + ("ruff_format",)}
        self.addCleanup(lambda: [setattr(generate, n, v) for n, v in saved.items()])
        generate.ROOT, generate.STUB = root, root / "stub.pyi"
        generate.OVERRIDES = root / "overrides.pyi"
        generate.OVERRIDES.write_text(OVERRIDES)
        ext = types.SimpleNamespace(__file__="the compiled module")
        generate.installed_package = lambda: (module(), ext)
        generate.run_stubgen = lambda raw_dir: RAW
        generate.ruff_format = lambda text: text
        self.want = generate.build_stub(RAW, module())

    def main(self, *argv: str) -> tuple[int, str, str]:
        out, err = io.StringIO(), io.StringIO()
        with (
            mock.patch.object(sys, "argv", ["generate.py", *argv]),
            contextlib.redirect_stdout(out),
            contextlib.redirect_stderr(err),
        ):
            status = generate.main()
        return status, out.getvalue(), err.getvalue()

    def drift(self, committed: str | None) -> str:
        if committed is None:
            generate.STUB.unlink(missing_ok=True)
        else:
            generate.STUB.write_text(committed)
        status, out, err = self.main("--check")
        self.assertEqual(status, 1, out)
        self.assertIn("stub.pyi differs from a fresh generation", err)
        return out

    def test_write_then_check(self):
        status, out, _ = self.main("--write")
        self.assertEqual(status, 0)
        self.assertEqual(generate.STUB.read_text(), self.want)
        status, out, _ = self.main("--check")
        self.assertEqual(status, 0)
        self.assertIn("stub.pyi is up to date", out)

    def test_check_fails_on_any_difference(self):
        changed = self.want.replace("-> int:", "-> str:", 1)
        self.assertNotEqual(changed, self.want)
        out = self.drift(changed)
        self.assertIn("\n-", out)
        self.assertIn("\n+", out)
        self.assertIn("-> str:", out)
        self.drift(self.want + "\n")
        self.drift(self.want[:-1])
        self.drift("")
        self.assertIn("(generated)", self.drift(None))


class InstalledPackage(unittest.TestCase):
    """``installed_package``: both modules imported from this interpreter's site-packages and
    the extension without the private hooks, else a Failure."""

    def check(self, pkg_file: str, ext_file: str, hooks: bool = False):
        pkg = types.ModuleType("netgraph_core")
        ext = types.ModuleType("_netgraph_core")
        pkg.__file__, ext.__file__ = pkg_file, ext_file
        if hooks:
            ext._fail_native_allocation = lambda: 0
        with mock.patch.dict(
            sys.modules, {"netgraph_core": pkg, "_netgraph_core": ext}
        ):
            return generate.installed_package()

    def test_refusals(self):
        site = os.path.realpath(sysconfig.get_paths()["platlib"])
        inside = os.path.join(site, "netgraph_core", "__init__.py")
        ext = os.path.join(site, "_netgraph_core.so")
        with tempfile.TemporaryDirectory() as tmp:
            outside = os.path.join(os.path.realpath(tmp), "python", "netgraph_core.py")
            for pkg_file, ext_file, name in (
                (outside, ext, "netgraph_core"),
                (inside, outside, "_netgraph_core"),
            ):
                with self.assertRaises(generate.Failure) as ctx:
                    self.check(pkg_file, ext_file)
                self.assertIn(
                    f"{name} is imported from {outside}, not from an installation",
                    str(ctx.exception),
                )
        with self.assertRaises(generate.Failure) as ctx:
            self.check(inside, ext, hooks=True)
        self.assertIn(
            "the installed extension carries the private hooks", str(ctx.exception)
        )
        pkg, got = self.check(inside, ext)
        self.assertEqual((pkg.__file__, got.__file__), (inside, ext))


if __name__ == "__main__":
    unittest.main()
