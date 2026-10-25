"""Check exact source/wheel payloads and bind release evidence to published bytes."""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import tarfile
import zipfile
from pathlib import Path, PurePosixPath

ROOT = Path(__file__).resolve().parents[2]


def sha(path: Path) -> str:
    return hashlib.file_digest(path.open("rb"), "sha256").hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sdist", type=Path, required=True)
    parser.add_argument("--wheel", type=Path, action="append", default=[])
    parser.add_argument("--receipt", type=Path, required=True)
    parser.add_argument("--extract", type=Path)
    parser.add_argument("--verify", action="store_true")
    args = parser.parse_args()
    if args.verify:
        report = json.loads(args.receipt.read_text())
        if report["sdist"] != {"name": args.sdist.name, "sha256": sha(args.sdist)}:
            raise ValueError(
                "published sdist differs from the archive whose proofs were checked"
            )
        if {path.name: sha(path) for path in args.wheel} != report["wheels"]:
            raise ValueError("published wheels differ from the checked release wheels")
        print("Published artifact identities: PASS")
        return
    spec = importlib.util.spec_from_file_location(
        "contents", ROOT / "dev/stubs/check_contents.py"
    )
    contents = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(contents)
    with tarfile.open(args.sdist) as archive:
        members = [member for member in archive if not member.isdir()]
        tops = {PurePosixPath(member.name).parts[0] for member in members}
        if len(tops) != 1:
            raise ValueError("sdist has no single root")
        source_names = []
        for member in members:
            name = PurePosixPath(member.name)
            if (
                not member.isfile()
                or ".." in name.parts
                or name.is_absolute()
                or len(name.parts) < 2
            ):
                raise ValueError(f"unsafe/nonregular sdist entry: {member.name}")
            source_names.append(str(PurePosixPath(*name.parts[1:])))
        if len(set(source_names)) != len(source_names):
            raise ValueError("duplicate sdist entry")
        intended = contents.checkout_files(ROOT)
        expected = {name for name in intended if contents.shipped(name)} | {"PKG-INFO"}
        if set(source_names) != expected:
            raise ValueError(
                f"sdist payload differs: missing={sorted(expected - set(source_names))}, extra={sorted(set(source_names) - expected)}"
            )
        if args.extract:
            if args.extract.exists() and any(args.extract.iterdir()):
                raise ValueError("extract destination must be empty")
            args.extract.mkdir(parents=True, exist_ok=True)
            for member, relative in zip(members, source_names, strict=True):
                output = args.extract / relative
                output.parent.mkdir(parents=True, exist_ok=True)
                output.write_bytes(archive.extractfile(member).read())
                output.chmod(member.mode & 0o777)
    for wheel in args.wheel:
        with zipfile.ZipFile(wheel) as archive:
            names = [name for name in archive.namelist() if not name.endswith("/")]
        problems = contents.problems(names, source_names, intended)
        if problems:
            raise ValueError("; ".join(problems))
    report = {
        "sdist": {"name": args.sdist.name, "sha256": sha(args.sdist)},
        "wheels": {path.name: sha(path) for path in args.wheel},
        "sdist_members": sorted(source_names),
    }
    args.receipt.parent.mkdir(parents=True, exist_ok=True)
    args.receipt.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(
        f"Proof distribution payload: PASS ({len(source_names)} source files; {len(args.wheel)} wheels)"
    )


if __name__ == "__main__":
    main()
