"""Negative controls of source safety, actual Lean axiom checks, replay and the bridge.

Mutations live only in a disposable proof project. --preflight-only runs the
standard-library source, comparison and artifact controls without a Lean toolchain.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
PROOFS = ROOT / "proofs"
MODEL = "proofs/NetGraphCore/Apportionment/RoundRobin.lean"


def load(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--preflight-only", action="store_true")
    args = parser.parse_args()
    tracked = {}
    for directory, children, files in os.walk(PROOFS):
        children[:] = [name for name in children if name != ".lake"]
        for name in files:
            path = Path(directory) / name
            tracked[str(path.relative_to(ROOT))] = path.read_bytes()
    work_parent = Path(os.environ.get("PROOFS_OUTPUT", str(ROOT / "build/proofs")))
    work_parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="controls-", dir=work_parent) as temporary:
        root = Path(temporary)
        shutil.copytree(PROOFS, root / "proofs", ignore=shutil.ignore_patterns(".lake"))
        shutil.copytree(ROOT / "dev/proofs", root / "dev/proofs")
        (root / "proofs/.lake").mkdir()
        (root / "proofs/.lake/packages").symlink_to(
            PROOFS / ".lake/packages", target_is_directory=True
        )
        baseline = {
            str(path.relative_to(root)): path.read_bytes()
            for path in root.rglob("*")
            if path.is_file() and ".lake" not in path.parts
        }
        environment = {
            **os.environ,
            "PROOFS_PYTHON": sys.executable,
            "PROOFS_OUTPUT": str(root / "evidence"),
        }

        def run(label: str, expected: str | None = None, mode: str = "check") -> None:
            command = (
                [sys.executable, "-I", "-B", "proofs/Checks/audit.py", "--preflight"]
                if mode == "preflight"
                else ["bash", "dev/proofs/run.sh", mode]
            )
            result = subprocess.run(
                command,
                cwd=root,
                env=environment,
                capture_output=True,
                text=True,
                timeout=300,
            )
            output = result.stdout + result.stderr
            (work_parent / f"control-{label}.log").write_text(
                output + f"\nSUMMARY exit={result.returncode}\n"
            )
            if expected is None:
                if result.returncode:
                    raise AssertionError(
                        f"{label}: clean control failed: {output[-6000:]}"
                    )
            elif result.returncode == 0 or expected not in output:
                raise AssertionError(
                    f"{label}: wrong gate outcome; wanted {expected!r}: {output[-6000:]}"
                )
            print(
                f"{label}: {'PASS' if expected is None else 'rejected as required'}",
                flush=True,
            )

        def restore() -> None:
            for path in root.rglob("*"):
                if (
                    path.is_file()
                    and ".lake" not in path.parts
                    and not path.is_relative_to(root / "evidence")
                ):
                    if str(path.relative_to(root)) not in baseline:
                        path.unlink()
            for name, data in baseline.items():
                (root / name).write_bytes(data)

        def edit(name: str, old: str, new: str) -> None:
            path = root / name
            text = path.read_text()
            if old not in text:
                raise AssertionError(f"plant not applied: {name}: {old}")
            path.write_text(text.replace(old, new, 1))

        run("source-clean", mode="preflight")
        for label, command in [
            ("sorry", "example : True := by sorry"),
            ("admit", "example : True := by admit"),
            ("project-axiom", "axiom unfinished : False"),
            ("unsafe", "unsafe def unchecked : Nat := 0"),
            ("native-decide", "example : True := by native_decide"),
            ("kernel-skip", "set_option debug.skipKernelTC true"),
            ("unchecked-command", "run_cmd pure ()"),
            ("unchecked-elaborator", 'elab "unchecked" : command => pure ()'),
        ]:
            path = root / MODEL
            path.write_text(path.read_text() + "\n" + command + "\n")
            run(label, "forbidden proof command or kernel bypass", mode="preflight")
            restore()
        with (root / MODEL).open("a") as stream:
            stream.write(
                "\n/- sorry /- unsafe -/ axiom -/\n-- set_option debug.skipKernelTC true\n"
            )
            stream.write(
                '\n/-- A string is data, not a proof command. -/\ndef textExample : String := "sorry unsafe axiom"\n'
            )
        run("quoted-and-commented", mode="preflight")
        restore()
        (root / "proofs/Checks/Forgotten.lean").write_text("module\n")
        run(
            "unimported-module",
            "aggregation root does not cover local modules",
            mode="preflight",
        )
        restore()
        run("source-restored", mode="preflight")

        if not args.preflight_only:
            run("clean")
            # This postulate is outside the mathematical-source scanner. The built
            # environment must expose the project axiom or its transitive use.
            (root / "proofs/Checks/Untrusted.lean").write_text(
                "module\npublic section\nnamespace Untrusted\n"
                "axiom forge {p : Prop} : p\nend Untrusted\nend\n"
            )
            edit(MODEL, "module\n", "module\nimport Checks.Untrusted\n")
            path = root / MODEL
            text = path.read_text()
            start = text.index("theorem native_increment_bounds")
            body = text.index(" := by", start) + len(" := by")
            path.write_text(
                text[:body]
                + "\n  have _ := hm\n  have _ := hi\n  have _ := hk\n  have _ := hv\n  exact Untrusted.forge\n\nend\nend NetGraphCore.Apportionment\n"
            )
            edit(
                "proofs/Checks.lean",
                "public import Checks.Main",
                "public import Checks.Main\npublic import Checks.Untrusted",
            )
            run("imported-axiom", "axiom")
            restore()
            (root / "proofs/Checks/Untrusted.lean").write_text(
                "module\npublic theorem unfinished : True := by sorry\n"
            )
            edit(
                "proofs/Checks.lean",
                "public import Checks.Main",
                "public import Checks.Main\npublic import Checks.Untrusted",
            )
            run("imported-sorry", "declaration uses `sorry`")
            restore()
            # Corrupt the checker environment, where metaprogramming is allowed.
            # The axiom audit alone passes; fresh kernel replay must reject it.
            (root / "proofs/Checks/Untrusted.lean").write_text(
                "module\npublic import Lean\nopen Lean Elab Command\n"
                "set_option debug.skipKernelTC true in\n"
                "run_cmd liftCoreM <| addDecl (.thmDecl {\n"
                "  name := `Untrusted.bad, levelParams := [], type := .const ``False [],\n"
                "  value := .const ``True.intro [] })\n"
            )
            edit(
                "proofs/Checks.lean",
                "public import Checks.Main",
                "public import Checks.Main\npublic import Checks.Untrusted",
            )
            run("kernel-plant-audit")
            run("kernel-plant-replay", "declaration type mismatch", mode="replay")
            restore()
            main_path = root / "proofs/Checks/Main.lean"
            text = main_path.read_text()
            start = text.index("public unsafe def main")
            main_path.write_text(
                text[:start]
                + "public unsafe def main (_args : List String) : IO UInt32 := pure 0\n"
            )
            run("disabled-audit", "missing proof evidence file")
            restore()
            run("restored-clean")
    for name, data in tracked.items():
        if (ROOT / name).read_bytes() != data:
            raise AssertionError(f"control changed product input: {name}")
    bridge = load(ROOT / "dev/proofs/bridge.py", "bridge_controls")
    for wrong in ([3, 3, 4], [4.0, 3, 3], [10, 0, 0]):
        try:
            bridge.exact(wrong, [4, 3, 3], "negative control")
        except ValueError:
            continue
        raise AssertionError(f"bridge accepted wrong result {wrong}")
    # The publish-time verifier must reject changed bytes of the already checked archive.
    with tempfile.TemporaryDirectory(
        prefix="artifact-identity-", dir=work_parent
    ) as temporary:
        work = Path(temporary)
        archive = work / "candidate.tar.gz"
        archive.write_bytes(b"identity-control fixture; not a real source archive")
        receipt = work / "receipt.json"
        receipt.write_text(
            json.dumps(
                {
                    "sdist": {
                        "name": archive.name,
                        "sha256": hashlib.sha256(archive.read_bytes()).hexdigest(),
                    },
                    "wheels": {},
                }
            )
        )
        command = [
            sys.executable,
            "-I",
            "-B",
            str(ROOT / "dev/proofs/package.py"),
            "--verify",
            "--sdist",
            str(archive),
            "--receipt",
            str(receipt),
        ]
        subprocess.run(command, check=True, capture_output=True)
        archive.write_bytes(archive.read_bytes() + b"changed")
        wrong = subprocess.run(command, capture_output=True, text=True)
        if wrong.returncode == 0 or "published sdist differs" not in wrong.stderr:
            raise AssertionError(
                "publish identity check accepted changed archive bytes"
            )
    print(
        "Proof gate controls: PASS; source inputs unchanged; bridge and publication identity reject wrong results."
        + (
            " Lean build/replay controls not run (--preflight-only)."
            if args.preflight_only
            else ""
        )
    )


if __name__ == "__main__":
    main()
