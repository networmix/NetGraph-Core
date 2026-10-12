"""Check the files of the wheel and the sdist against the files intended for them.

Run by ``bash dev/stubs/run.sh check`` with the wheel and the sdist it built from this checkout:

* the sdist holds exactly the intended source files, plus the PKG-INFO the backend writes. The
  intended files are the checkout's files that git tracks or would track (new files git does
  not ignore) under the top-level entries of SDIST_TOP, without tests/reference/. That rule is
  the one of pyproject.toml's sdist.include and sdist.exclude, repeated here so that a change
  to that list is caught. A file git ignores (a cache, a compiled module, a local tool's
  output) or a file outside those entries is refused, and so is an intended file left out;
* the wheel holds exactly the package files of the sdist's python/netgraph_core/ under
  netgraph_core/, one compiled extension module, and its .dist-info directory with the
  metadata files and both licence notices (the project's and pybind11's);
* both carry the stub, py.typed and the licence notices, and neither carries a _docs.py.

The one file git cannot tell apart from an intended one is a new file it does not ignore inside
a shipped directory: ``git status`` lists it, and CI builds from a clean checkout.

Usage: python -I dev/stubs/check_contents.py WHEEL SDIST WORK_DIR
"""

from __future__ import annotations

import re
import subprocess
import sys
import tarfile
import zipfile
from collections.abc import Iterable
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
# The top-level files and directories the sdist ships (a directory ends with "/") and the
# directories under them it leaves out.
SDIST_TOP = (
    ".coveragerc",
    ".gitignore",
    "CHANGELOG.md",
    "CMakeLists.txt",
    "CONTRIBUTING.md",
    "LICENSE",
    "README.md",
    "pyproject.toml",
    "bindings/",
    "docs/",
    "include/",
    "python/",
    "src/",
    "tests/",
    "third_party/",
)
SDIST_LEFT_OUT = ("tests/reference/",)
PACKAGE = "python/netgraph_core/"
NEED = {
    "wheel": [
        "netgraph_core/__init__.py",
        "netgraph_core/__init__.pyi",
        "netgraph_core/py.typed",
    ],
    "sdist": [
        "python/netgraph_core/__init__.py",
        "python/netgraph_core/__init__.pyi",
        "python/netgraph_core/py.typed",
        "LICENSE",
        "third_party/pybind11/LICENSE",
    ],
}
EXTENSION = re.compile(r"_netgraph_core\.[^/]+\.(so|pyd)")
DIST_INFO = re.compile(r"netgraph_core-[^/]+\.dist-info/")
DIST_INFO_FILES = {
    "METADATA",
    "WHEEL",
    "RECORD",
    "licenses/LICENSE",
    "licenses/third_party/pybind11/LICENSE",
}


def shipped(name: str) -> bool:
    """Whether the sdist rule selects the checkout path ``name``."""
    if name.startswith(SDIST_LEFT_OUT):
        return False
    return any(
        name.startswith(top) if top.endswith("/") else name == top for top in SDIST_TOP
    )


def checkout_files(root: Path) -> list[str]:
    """The checkout's files that git tracks or would track (a new file it does not ignore),
    as they are on disk now: a tracked file deleted from the work tree is not one."""
    cmd = ["git", "-C", str(root), "ls-files", "-z", "--cached", "--others"]
    try:
        out = subprocess.run(
            [*cmd, "--exclude-standard"], capture_output=True, check=True
        ).stdout
    except (OSError, subprocess.CalledProcessError) as exc:
        sys.exit(
            f"dev/stubs/check_contents.py: cannot list the checkout's files ({exc})"
        )
    names = {n for n in out.decode("utf-8").split("\0") if n}
    return sorted(n for n in names if (root / n).is_file())


def problems(
    wheel: Iterable[str], sdist: Iterable[str], intended: Iterable[str]
) -> list[str]:
    """What is wrong with the member names of the wheel and of the sdist (the sdist's without
    their top directory), given the checkout files ``intended`` for the sdist."""
    wheel, sdist = set(wheel), set(sdist)
    want_sdist = {n for n in intended if shipped(n)} | {"PKG-INFO"}
    package = {n[len(PACKAGE) :] for n in want_sdist if n.startswith(PACKAGE)}
    out = []
    for kind, names in (("wheel", wheel), ("sdist", sdist)):
        out += [f"the {kind} lacks {n}" for n in NEED[kind] if n not in names]
        out += [
            f"the {kind} carries {n}" for n in sorted(names) if n.endswith("_docs.py")
        ]
    out += [
        f"the sdist carries {n}, not an intended file"
        for n in sorted(sdist - want_sdist)
    ]
    out += [f"the sdist lacks {n}" for n in sorted(want_sdist - sdist)]

    extensions = sorted(n for n in wheel if EXTENSION.fullmatch(n))
    if len(extensions) != 1:
        out.append(f"the wheel holds {len(extensions)} extension modules: {extensions}")
    infos = {m.group(0) for n in wheel if (m := DIST_INFO.match(n))}
    if len(infos) != 1:
        out.append(
            f"the wheel holds {len(infos)} .dist-info directories: {sorted(infos)}"
        )
    want_wheel = {f"netgraph_core/{n}" for n in package} | set(extensions)
    want_wheel |= {i + n for i in infos for n in DIST_INFO_FILES}
    out += [
        f"the wheel carries {n}, not an intended file"
        for n in sorted(wheel - want_wheel)
    ]
    out += [f"the wheel lacks {n}" for n in sorted(want_wheel - wheel)]
    return out


def main() -> int:
    if len(sys.argv) != 4:
        sys.exit("usage: check_contents.py WHEEL SDIST WORK_DIR")
    wheel, sdist, work = sys.argv[1], sys.argv[2], Path(sys.argv[3])
    with zipfile.ZipFile(wheel) as z:
        wheel_names = sorted(n for n in z.namelist() if not n.endswith("/"))
    with tarfile.open(sdist) as t:
        members = [m for m in t if not m.isdir()]
    tops = {m.name.split("/", 1)[0] for m in members}
    sdist_names = sorted(m.name.split("/", 1)[1] for m in members if "/" in m.name)
    found = list(problems(wheel_names, sdist_names, checkout_files(ROOT)))
    if len(tops) != 1 or len(sdist_names) != len(members):
        found.append(
            f"the sdist's files are not all under one top directory: {sorted(tops)}"
        )
    for kind, names in (("wheel", wheel_names), ("sdist", sdist_names)):
        (work / f"contents-{kind}.txt").write_text("\n".join(names) + "\n")
        print(f"{kind}: {len(names)} files ({work}/contents-{kind}.txt)")
    if found:
        for p in found:
            print(f"dev/stubs/check_contents.py: {p}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
