"""Generate the package's type stub, python/netgraph_core/__init__.pyi, or check it for drift.

Run by ``bash dev/stubs/run.sh`` with the interpreter of an environment holding an installed,
non-editable netgraph_core (a wheel or an sdist build: the shipped configuration, without the
private hooks) and the pinned pybind11-stubgen and ruff of pyproject.toml's ``dev`` extra.

The steps, each failing loudly:

1. pybind11-stubgen writes the stub of the compiled module ``_netgraph_core``.
2. The package stub keeps the declarations of ``netgraph_core.__all__`` (the private names of the
   compiled module and private members are not part of the package), under the package's own
   docstring, ``__all__`` and ``__version__``.
3. The reviewed overrides of dev/stubs/overrides.pyi replace generated declarations. An override
   may change annotations only: its parameter names, kinds and defaults must equal the generated
   ones, the exact generated line it replaces is quoted above it (``# generated:``) with its
   reason (``# reason:``), and a changed generated line fails until the override is reviewed
   again. Overloads replace one generated declaration. The one addition allowed is the marker
   of a class the module refuses to construct (``def __init__(self, _no_public_constructor:
   typing.Never, /)``): a type checker then rejects ``Cls()``; the marker must match the runtime
   both ways.
4. No ``typing.Any`` remains, except where the module accepts any object: the operand of
   ``__eq__``/``__ne__`` and the memo of ``__deepcopy__``.
5. The module's positional-only parameters (``/`` in its signatures, which pybind11-stubgen
   3.0.0 drops) are positional-only in the stub too, wherever a parameter other than the
   receiver is affected (the enums' comparison operand and ``__setstate__``'s state).
6. ruff (pinned) formats the result; ``--write`` writes it, ``--check`` compares it with the
   committed stub and prints the difference.

Usage: python -I dev/stubs/generate.py (--write | --check) [--raw DIR]
"""

from __future__ import annotations

import argparse
import ast
import copy
import difflib
import importlib
import inspect
import os
import subprocess
import sys
import sysconfig
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
STUB = ROOT / "python/netgraph_core/__init__.pyi"
OVERRIDES = ROOT / "dev/stubs/overrides.pyi"
HEADER = "# Generated from the compiled module; do not edit.\n"
# The imports, grouped as ruff's isort rules want them (ast.unparse keeps no blank lines).
IMPORTS = """from __future__ import annotations

import collections.abc
import fractions
import typing

import numpy
import numpy.typing
"""
MARKER_PARAM = "_no_public_constructor"
# Where the module accepts any object, typing.Any is the exact domain: the comparison operand
# (pybind11 names it "other" on enums, "arg0" on the value classes, which define __eq__ only)
# and copy.deepcopy's memo.
ANY_ALLOWED = {
    ("__eq__", "other"),
    ("__eq__", "arg0"),
    ("__ne__", "other"),
    ("__deepcopy__", "memo"),
}


class Failure(Exception):
    pass


def installed_package():
    """netgraph_core and _netgraph_core, imported from this interpreter's site-packages (an
    installation, not the source tree or an editable build's redirect)."""
    pkg = importlib.import_module("netgraph_core")
    ext = importlib.import_module("_netgraph_core")
    roots = {os.path.realpath(sysconfig.get_paths()[k]) for k in ("purelib", "platlib")}
    for mod in (pkg, ext):
        where = os.path.realpath(mod.__file__ or "")
        if not any(where.startswith(r + os.sep) for r in roots):
            raise Failure(
                f"{mod.__name__} is imported from {where}, not from an installation in "
                "this interpreter's site-packages (build a wheel: bash dev/stubs/run.sh)"
            )
    if hasattr(ext, "_fail_native_allocation"):
        raise Failure(
            "the installed extension carries the private hooks; generate from the "
            "shipped configuration (a wheel or an sdist build)"
        )
    return pkg, ext


def run_stubgen(raw_dir: Path) -> str:
    cmd = [
        sys.executable,
        "-I",
        "-m",
        "pybind11_stubgen",
        "--exit-code",
        "-o",
        str(raw_dir),
        "_netgraph_core",
    ]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0 or proc.stderr.strip():
        raise Failure(
            f"pybind11-stubgen failed ({proc.returncode}):\n{proc.stderr}{proc.stdout}"
        )
    return (raw_dir / "_netgraph_core.pyi").read_text(encoding="utf-8")


def is_private(name: str) -> bool:
    return name.startswith("_") and not (name.startswith("__") and name.endswith("__"))


def params(fn: ast.FunctionDef) -> list[tuple[str, str, str | None]]:
    """(name, kind, default) of every parameter, annotations left out."""
    a = fn.args
    pos = [*a.posonlyargs, *a.args]
    defaults: list[ast.expr | None] = [None] * (len(pos) - len(a.defaults))
    defaults += list(a.defaults)
    out = []
    for i, (arg, d) in enumerate(zip(pos, defaults, strict=True)):
        kind = "positional-only" if i < len(a.posonlyargs) else "positional"
        out.append((arg.arg, kind, None if d is None else ast.unparse(d)))
    if a.vararg:
        out.append((a.vararg.arg, "*args", None))
    for arg, d in zip(a.kwonlyargs, a.kw_defaults, strict=True):
        out.append((arg.arg, "keyword-only", None if d is None else ast.unparse(d)))
    if a.kwarg:
        out.append((a.kwarg.arg, "**kwargs", None))
    return out


def is_marker(fn: ast.FunctionDef) -> bool:
    a = fn.args
    return (
        fn.name == "__init__"
        and not fn.decorator_list
        and [x.arg for x in a.posonlyargs] == ["self", MARKER_PARAM]
        and not (a.args or a.vararg or a.kwonlyargs or a.kwarg or a.defaults)
        and ast.unparse(a.posonlyargs[1].annotation or ast.Constant(None))
        == "typing.Never"
        and ast.unparse(fn.returns or ast.Constant(None)) == "None"
    )


def comment_block(lines: list[str], first_line: int) -> list[str]:
    """The comment lines directly above line ``first_line`` (1-based)."""
    out: list[str] = []
    i = first_line - 2
    while i >= 0 and lines[i].strip().startswith("#"):
        out.insert(0, lines[i].strip())
        i -= 1
    return out


class Override:
    def __init__(self, qual: str, defs: list[ast.FunctionDef], lines: list[str]):
        self.qual = qual
        self.defs = defs
        first = defs[0]
        start = min([first.lineno, *(d.lineno for d in first.decorator_list)])
        comments = comment_block(lines, start)
        self.reasons = [
            c.removeprefix("# reason:").strip()
            for c in comments
            if c.startswith("# reason:")
        ]
        gen = [
            c.removeprefix("# generated:").strip()
            for c in comments
            if c.startswith("# generated:")
        ]
        where = f"{OVERRIDES.relative_to(ROOT)}:{first.lineno} {qual}"
        if not self.reasons or not all(self.reasons):
            raise Failure(f"{where}: an override needs a '# reason:' line above it")
        if len(gen) != 1:
            raise Failure(
                f"{where}: an override needs exactly one '# generated:' line above it "
                "(the generated line it replaces, or '(none)' for the constructor marker)"
            )
        self.generated = gen[0]
        for d in defs[1:]:
            stray = comment_block(
                lines, min([d.lineno, *(x.lineno for x in d.decorator_list)])
            )
            if any(c.startswith(("# reason:", "# generated:")) for c in stray):
                raise Failure(
                    f"{where}: '# reason:'/'# generated:' belong above the first overload only"
                )
        if len(defs) > 1 and not all(
            any(ast.unparse(x) == "typing.overload" for x in d.decorator_list)
            for d in defs
        ):
            raise Failure(
                f"{where}: repeated declarations must all be @typing.overload"
            )


def load_overrides() -> dict[str, Override]:
    text = OVERRIDES.read_text(encoding="utf-8")
    lines = text.splitlines()
    tree = ast.parse(text)
    groups: dict[str, list[ast.FunctionDef]] = {}
    order: list[str] = []

    def add(qual: str, fn: ast.FunctionDef) -> None:
        if qual in groups and order[-1] != qual:
            raise Failure(
                f"{OVERRIDES.name}:{fn.lineno} {qual}: overloads must be adjacent"
            )
        if qual not in groups:
            order.append(qual)
        groups.setdefault(qual, []).append(fn)

    for node in tree.body:
        if isinstance(node, ast.FunctionDef):
            add(node.name, node)
        elif isinstance(node, ast.ClassDef):
            for item in node.body:
                if isinstance(item, ast.FunctionDef):
                    add(f"{node.name}.{item.name}", item)
                elif not (
                    isinstance(item, ast.Expr) and isinstance(item.value, ast.Constant)
                ):
                    raise Failure(
                        f"{OVERRIDES.name}:{item.lineno}: only method overrides inside a class"
                    )
        elif not isinstance(node, (ast.Import, ast.ImportFrom)) and not (
            isinstance(node, ast.Expr) and isinstance(node.value, ast.Constant)
        ):
            raise Failure(f"{OVERRIDES.name}:{node.lineno}: unexpected statement")
    return {q: Override(q, groups[q], lines) for q in order}


def generated_line(raw_lines: list[str], fn: ast.FunctionDef) -> str:
    return raw_lines[fn.lineno - 1].strip()


def body_of(fn: ast.FunctionDef) -> list[ast.stmt]:
    """The generated body (the docstring or ``...``)."""
    return copy.deepcopy(fn.body)


def with_body(fn: ast.FunctionDef, body: list[ast.stmt]) -> ast.FunctionDef:
    new = copy.deepcopy(fn)
    new.body = body
    return new


def apply_overrides(
    container: list[ast.stmt],
    owner: str,
    raw_lines: list[str],
    overrides: dict[str, Override],
    used: set[str],
) -> list[ast.stmt]:
    out: list[ast.stmt] = []
    funcs = [n for n in container if isinstance(n, ast.FunctionDef)]
    for node in container:
        if not isinstance(node, ast.FunctionDef):
            out.append(node)
            continue
        qual = f"{owner}.{node.name}" if owner else node.name
        ov = overrides.get(qual)
        if ov is None:
            out.append(node)
            continue
        if sum(n.name == node.name for n in funcs) > 1:
            raise Failure(
                f"{qual}: generated as overloads; an override replaces one declaration"
            )
        used.add(qual)
        if ov.generated == "(none)":
            raise Failure(
                f"{qual}: the override says '(none)', but {qual} is generated"
            )
        got = generated_line(raw_lines, node)
        if got != ov.generated:
            raise Failure(
                f"{qual}: the generated declaration changed; review the override again\n"
                f"  override quotes: {ov.generated}\n  generated now:   {got}"
            )
        want = params(node)
        for d in ov.defs:
            if params(d) != want:
                raise Failure(
                    f"{qual}: an override changes only annotations; parameters "
                    f"{params(d)} differ from the generated {want}"
                )
        decos = {ast.unparse(x) for x in node.decorator_list}
        for d in ov.defs:
            have = {ast.unparse(x) for x in d.decorator_list} - {"typing.overload"}
            if have != decos:
                raise Failure(
                    f"{qual}: decorators {sorted(have)} differ from the generated {sorted(decos)}"
                )
        for i, d in enumerate(ov.defs):
            body = body_of(node) if i == 0 else [ast.Expr(ast.Constant(Ellipsis))]
            out.append(with_body(d, body))
    return out


def add_markers(
    cls: ast.ClassDef, overrides: dict[str, Override], used: set[str]
) -> None:
    qual = f"{cls.name}.__init__"
    ov = overrides.get(qual)
    if ov is None or any(
        isinstance(n, ast.FunctionDef) and n.name == "__init__" for n in cls.body
    ):
        return
    used.add(qual)
    if ov.generated != "(none)" or len(ov.defs) != 1 or not is_marker(ov.defs[0]):
        raise Failure(
            f"{qual}: not generated; the only addition allowed is the constructor marker "
            f"'def __init__(self, {MARKER_PARAM}: typing.Never, /) -> None' "
            "with '# generated: (none)'"
        )
    marker = with_body(ov.defs[0], [ast.Expr(ast.Constant(Ellipsis))])
    # After the class docstring, before the first member (the generated order is kept).
    at = (
        1
        if cls.body
        and isinstance(cls.body[0], ast.Expr)
        and isinstance(cls.body[0].value, ast.Constant)
        else 0
    )
    cls.body.insert(at, marker)


def refuses_construction(cls: type) -> bool:
    try:
        cls()
    except TypeError as exc:
        return "No constructor defined" in str(exc)
    return False


def positional_only(doc: str | None) -> list[list[str]]:
    """Per pybind11 signature of a docstring, the names its ``/`` makes positional-only."""
    lines = (doc or "").strip("\n").split("\n")
    if len(lines) > 1 and lines[1].strip() == "Overloaded function.":
        sigs = [ln.split(". ", 1)[1] for ln in lines if ln[:1].isdigit() and ". " in ln]
    else:
        sigs = lines[:1]
    out = []
    for sig in sigs:
        if "(" not in sig or ")" not in sig:
            out.append([])
            continue
        inner = sig[sig.index("(") + 1 : sig.rindex(")")]
        parts, depth, cur = [], 0, ""
        for ch in inner:
            depth += ch in "([{"
            depth -= ch in ")]}"
            if ch == "," and depth == 0:
                parts.append(cur.strip())
                cur = ""
            else:
                cur += ch
        parts.append(cur.strip())
        names = [p.split(":")[0].split("=")[0].strip() for p in parts]
        out.append(names[: names.index("/")] if "/" in names else [])
    return out


def restore_positional_only(defs: list[ast.FunctionDef], obj, qual: str) -> None:
    """pybind11-stubgen 3.0.0 drops the ``/`` of the module's signatures. Where the module makes
    a parameter other than the receiver positional-only, the stub declares the same leading
    parameters positional-only (the enums' ``__eq__``/``__ne__`` operand, ``__setstate__``'s
    state); a lone positional-only receiver is left as generated."""
    sigs = positional_only(getattr(obj, "__doc__", None))
    if not any(set(names) - {"self"} for names in sigs):
        return
    if len(sigs) == 1:
        sigs = sigs * len(defs)
    if len(sigs) != len(defs):
        raise Failure(
            f"{qual}: {len(sigs)} module signatures, {len(defs)} stub declarations"
        )
    for names, fn in zip(sigs, defs, strict=True):
        a = fn.args
        pos = [*a.posonlyargs, *a.args]
        if [x.arg for x in pos[: len(names)]] != names:
            raise Failure(
                f"{qual}: the module makes {names} positional-only; the stub starts with "
                f"{[x.arg for x in pos[: len(names)]]}"
            )
        a.posonlyargs, a.args = pos[: len(names)], pos[len(names) :]


def apply_positional_only(decls: list[ast.stmt], pkg) -> None:
    def groups(body: list[ast.stmt]) -> dict[str, list[ast.FunctionDef]]:
        out: dict[str, list[ast.FunctionDef]] = {}
        for n in body:
            if isinstance(n, ast.FunctionDef) and "property" not in {
                ast.unparse(d) for d in n.decorator_list
            }:
                out.setdefault(n.name, []).append(n)
        return out

    for name, defs in groups(decls).items():
        restore_positional_only(defs, getattr(pkg, name), name)
    for node in decls:
        if isinstance(node, ast.ClassDef):
            cls = getattr(pkg, node.name)
            for name, defs in groups(node.body).items():
                if name in vars(cls) and not any(is_marker(d) for d in defs):
                    restore_positional_only(
                        defs, getattr(cls, name), f"{node.name}.{name}"
                    )


def check_markers(classes: dict[str, ast.ClassDef], pkg) -> None:
    for name, node in classes.items():
        runtime = getattr(pkg, name)
        if issubclass(runtime, BaseException):
            continue
        marked = any(isinstance(n, ast.FunctionDef) and is_marker(n) for n in node.body)
        refuses = refuses_construction(runtime)
        if marked and not refuses:
            raise Failure(
                f"{name}: marked not constructible, but the module constructs it"
            )
        if refuses and not marked:
            raise Failure(
                f"{name}: the module refuses to construct it (No constructor defined); "
                "add the constructor marker to dev/stubs/overrides.pyi"
            )


def check_no_any(tree: ast.Module) -> None:
    problems = []

    def visit(fn: ast.FunctionDef, owner: str) -> None:
        a = fn.args
        # Every place a parameter can be annotated: positional-only, positional, *args,
        # keyword-only and **kwargs.
        every = [*a.posonlyargs, *a.args, a.vararg, *a.kwonlyargs, a.kwarg]
        for arg in (x for x in every if x is not None):
            if arg.annotation is not None and "typing.Any" in ast.unparse(
                arg.annotation
            ):
                if (fn.name, arg.arg) not in ANY_ALLOWED:
                    problems.append(f"{owner}{fn.name}({arg.arg})")
        if fn.returns is not None and "typing.Any" in ast.unparse(fn.returns):
            problems.append(f"{owner}{fn.name} -> ...")

    def scan(body: list[ast.stmt], owner: str) -> None:
        for node in body:
            if isinstance(node, ast.FunctionDef):
                visit(node, owner)
            elif isinstance(node, ast.ClassDef):
                scan(node.body, f"{owner}{node.name}.")
            elif isinstance(node, ast.AnnAssign) and "typing.Any" in ast.unparse(
                node.annotation
            ):
                problems.append(f"{owner}{ast.unparse(node.target)}")

    scan(tree.body, "")
    if problems:
        raise Failure(
            "typing.Any left in the stub (add a reviewed override): "
            + ", ".join(problems)
        )


def build_stub(raw: str, pkg) -> str:
    raw_lines = raw.splitlines()
    tree = ast.parse(raw)
    exported = list(pkg.__all__)
    overrides = load_overrides()
    used: set[str] = set()
    decls: list[ast.stmt] = []
    seen: set[str] = set()
    classes: dict[str, ast.ClassDef] = {}
    for node in tree.body:
        if isinstance(node, (ast.ClassDef, ast.FunctionDef)):
            name = node.name
        elif isinstance(node, ast.AnnAssign) and isinstance(node.target, ast.Name):
            name = node.target.id
        else:
            continue  # the module docstring, imports and the module's own __all__
        if name == "__all__":
            continue
        if name not in exported:
            if not is_private(name):
                raise Failure(
                    f"_netgraph_core.{name} is public but not exported by netgraph_core"
                )
            continue
        seen.add(name)
        if isinstance(node, ast.ClassDef):
            node.body = [
                n
                for n in node.body
                if not (isinstance(n, (ast.FunctionDef,)) and is_private(n.name))
                and not (
                    isinstance(n, ast.AnnAssign)
                    and isinstance(n.target, ast.Name)
                    and is_private(n.target.id)
                )
            ]
            node.body = apply_overrides(
                node.body, node.name, raw_lines, overrides, used
            )
            add_markers(node, overrides, used)
            classes[name] = node
        decls.append(node)
    decls = apply_overrides(decls, "", raw_lines, overrides, used)
    missing = [n for n in exported if n != "__version__" and n not in seen]
    if missing:
        raise Failure(f"exported by netgraph_core, not generated: {missing}")
    unused = [q for q in overrides if q not in used]
    if unused:
        raise Failure(f"overrides without a generated declaration: {unused}")
    check_markers(classes, pkg)
    apply_positional_only(decls, pkg)

    doc = ast.Expr(ast.Constant(inspect.cleandoc(pkg.__doc__ or "")))
    body: list[ast.stmt] = ast.parse(f"__all__: list[str] = {exported!r}").body
    body += ast.parse("__version__: str").body
    body += decls
    module = ast.Module(body=body, type_ignores=[])
    check_no_any(module)
    return f"{HEADER}{ast.unparse(ast.Module([doc], []))}\n\n{IMPORTS}\n{ast.unparse(module)}\n"


def ruff_format(text: str) -> str:
    proc = subprocess.run(
        [
            sys.executable,
            "-m",
            "ruff",
            "format",
            "--stdin-filename",
            str(STUB.relative_to(ROOT)),
            "-",
        ],
        input=text,
        capture_output=True,
        text=True,
        cwd=ROOT,
    )
    if proc.returncode != 0:
        raise Failure(f"ruff format failed:\n{proc.stderr}")
    return proc.stdout


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n", 1)[0])
    mode = ap.add_mutually_exclusive_group(required=True)
    mode.add_argument("--write", action="store_true", help="write the stub")
    mode.add_argument("--check", action="store_true", help="fail on drift")
    ap.add_argument(
        "--raw", type=Path, help="keep the raw pybind11-stubgen output here"
    )
    args = ap.parse_args()
    try:
        pkg, ext = installed_package()
        print(f"generating from {ext.__file__}")
        with tempfile.TemporaryDirectory() as tmp:
            raw_dir = args.raw or Path(tmp)
            raw = run_stubgen(raw_dir)
        stub = ruff_format(build_stub(raw, pkg))
    except Failure as exc:
        print(f"dev/stubs/generate.py: {exc}", file=sys.stderr)
        return 1
    if args.write:
        STUB.write_text(stub, encoding="utf-8")
        print(f"wrote {STUB.relative_to(ROOT)} ({len(stub.splitlines())} lines)")
        return 0
    committed = STUB.read_text(encoding="utf-8") if STUB.exists() else ""
    if committed == stub:
        print(
            f"{STUB.relative_to(ROOT)} is up to date ({len(stub.splitlines())} lines)"
        )
        return 0
    rel = str(STUB.relative_to(ROOT))
    sys.stdout.writelines(
        difflib.unified_diff(
            committed.splitlines(keepends=True),
            stub.splitlines(keepends=True),
            f"{rel} (committed)",
            f"{rel} (generated)",
        )
    )
    print(
        f"dev/stubs/generate.py: {rel} differs from a fresh generation "
        "(regenerate: bash dev/stubs/run.sh generate)",
        file=sys.stderr,
    )
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
