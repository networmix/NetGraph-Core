"""Real-discovery negative controls for the proof acceptance entry point.

Each mutation is made in a disposable copy. Dependencies are shared read-only by convention;
their checked-out revisions and clean state are checked by the audit. Product inputs are not
modified. Compilable weakened statements/definitions must fail the reviewed-binding check.
"""

from __future__ import annotations

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
    inventory = json.loads((PROOFS / "theorems.json").read_text())
    tracked = {MODEL: (ROOT / MODEL).read_bytes()}
    for record in inventory["theorems"]:
        for ref in [
            *record["semantics"],
            *record["implementation"]["sources"],
            record["explanation"],
        ]:
            tracked[ref["file"]] = (ROOT / ref["file"]).read_bytes()
    work_parent = ROOT / "build/proofs"
    work_parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="controls-", dir=work_parent) as temporary:
        root = Path(temporary)
        shutil.copytree(PROOFS, root / "proofs", ignore=shutil.ignore_patterns(".lake"))
        shutil.copytree(ROOT / "dev/proofs", root / "dev/proofs")
        for name, data in tracked.items():
            path = root / name
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(data)
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
        environment.pop("PROOFS_RECORD", None)

        def run(label: str, expected: str | None = None, mode: str = "check") -> None:
            result = subprocess.run(
                ["bash", "dev/proofs/run.sh", mode],
                cwd=root,
                env=environment,
                capture_output=True,
                text=True,
                timeout=300,
            )
            output = result.stdout + result.stderr
            (work_parent / f"control-{label}.log").write_text(output)
            if expected is None:
                if result.returncode:
                    raise AssertionError(f"clean control failed: {output[-6000:]}")
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

        def save_inventory(data: dict) -> None:
            (root / "proofs/theorems.json").write_text(
                json.dumps(data, indent=2) + "\n"
            )

        run("clean")
        edit(
            MODEL,
            "@[expose] def amountLimit : Nat := 2 ^ 62",
            "@[expose] def amountLimit : Nat := 2 ^ 61",
        )
        run("weakened-definition", "reviewed statement/definition digest drift")
        restore()
        edit(
            MODEL,
            "previous % members + extra < 2 ^ 63 ∧\n      0 ≤ increments previous extra members recipient ∧\n      increments previous extra members recipient < (amountLimit : Int)",
            "previous % members + extra < 2 ^ 63 ∧\n      0 ≤ increments previous extra members recipient",
        )
        run("weakened-statement", "reviewed statement/definition digest drift")
        restore()
        # An equivalent synthetic candidate imports a mathematical constant from an
        # owned checker module, outside the mathematical namespace. Its definition
        # must enter the fingerprint too. Candidate fingerprints are test setup in
        # this disposable directory; the repository's reviewed inventory is untouched.
        (root / "proofs/Checks/Scalar.lean").write_text(
            "module\npublic section\nnamespace Support\n"
            "/-- The synthetic candidate uses the same amount limit as the real model. -/\n"
            "@[expose] def importedLimit : Nat := 2 ^ 62\nend Support\nend\n"
        )
        edit(MODEL, "module\n", "module\npublic import Checks.Scalar\n")
        edit(
            MODEL,
            "def amountLimit : Nat := 2 ^ 62",
            "def amountLimit : Nat := Support.importedLimit",
        )
        edit(
            MODEL,
            "norm_num [amountLimit]",
            "norm_num [amountLimit, Support.importedLimit]",
        )
        changed = json.loads((root / "proofs/theorems.json").read_text())
        changed["modules"].append("Checks/Scalar.lean")
        save_inventory(changed)
        edit(
            "proofs/Checks.lean",
            "public import Checks.Main",
            "public import Checks.Main\npublic import Checks.Scalar",
        )
        built = subprocess.run(
            ["lake", "--wfail", "build"],
            cwd=root / "proofs",
            env=environment,
            capture_output=True,
            text=True,
            timeout=300,
        )
        (work_parent / "control-imported-definition-setup.log").write_text(
            built.stdout + built.stderr
        )
        built.check_returncode()
        exported = subprocess.run(
            ["lake", "exe", "proof-audit", "--metadata"],
            cwd=root / "proofs",
            env=environment,
            capture_output=True,
            text=True,
            check=True,
            timeout=300,
        )
        metadata = {
            item["name"]: item for item in json.loads(exported.stdout)["declarations"]
        }
        audit = load(root / "proofs/Checks/audit.py", "imported_definition_controls")
        for record in changed["theorems"]:
            record["binding"] = audit.binding(record["declaration"], metadata)
        if not any(
            "Support.importedLimit" in r["binding"]["local_closure"]
            for r in changed["theorems"]
        ):
            raise AssertionError(
                "owned imported definition omitted from actual constant closure"
            )
        save_inventory(changed)
        run("imported-definition-clean")
        edit("proofs/Checks/Scalar.lean", "2 ^ 62", "2 ^ 61")
        run(
            "imported-definition-weakened", "reviewed statement/definition digest drift"
        )
        restore()
        (root / "proofs/Checks/Scalar.lean").write_text(
            "module\npublic def undocumentedImported : Nat := 0\n"
        )
        changed = json.loads((root / "proofs/theorems.json").read_text())
        changed["modules"].append("Checks/Scalar.lean")
        save_inventory(changed)
        edit(
            "proofs/Checks.lean",
            "public import Checks.Main",
            "public import Checks.Main\npublic import Checks.Scalar",
        )
        run("checker-helper-doc", "missing declaration doc comment")
        restore()
        edit(
            MODEL,
            "namespace NetGraphCore.Apportionment",
            "namespace NetGraphCore.Apportionment\nset_option debug.skipKernelTC true",
        )
        run("kernel-skip", "forbidden proof command or kernel bypass")
        restore()
        # An imported postulate is deliberately outside the mathematical source scanner.
        # Transitive axiom inspection must reject its use even though Lean accepts it.
        (root / "proofs/Checks/Untrusted.lean").write_text(
            "module\npublic section\nnamespace Untrusted\n"
            "/-- A deliberately untrusted postulate for the transitive-axiom control. -/\n"
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
        changed = json.loads((root / "proofs/theorems.json").read_text())
        changed["modules"].append("Checks/Untrusted.lean")
        save_inventory(changed)
        edit(
            "proofs/Checks.lean",
            "public import Checks.Main",
            "public import Checks.Main\npublic import Checks.Untrusted",
        )
        run("imported-axiom", "unaccepted transitive axioms")
        restore()
        (root / "proofs/Checks/Untrusted.lean").write_text(
            "module\npublic theorem unfinished : True := by sorry\n"
        )
        changed = json.loads((root / "proofs/theorems.json").read_text())
        changed["modules"].append("Checks/Untrusted.lean")
        save_inventory(changed)
        edit(
            "proofs/Checks.lean",
            "public import Checks.Main",
            "public import Checks.Main\npublic import Checks.Untrusted",
        )
        run("imported-sorry", "declaration uses `sorry`")
        restore()
        # A trusted-looking axiom list cannot establish kernel validity. Deliberately
        # corrupt an imported environment outside the mathematical-source scanner.
        # The ordinary audit succeeds, but fresh replay must reject the false proof.
        (root / "proofs/Checks/Untrusted.lean").write_text(
            "module\npublic import Lean\nopen Lean Elab Command\n"
            "set_option debug.skipKernelTC true in\n"
            "run_cmd liftCoreM <| addDecl (.thmDecl {\n"
            "  name := `Untrusted.bad, levelParams := [], type := .const ``False [],\n"
            "  value := .const ``True.intro [] })\n"
        )
        changed = json.loads((root / "proofs/theorems.json").read_text())
        changed["modules"].append("Checks/Untrusted.lean")
        save_inventory(changed)
        edit(
            "proofs/Checks.lean",
            "public import Checks.Main",
            "public import Checks.Main\npublic import Checks.Untrusted",
        )
        run("kernel-plant-audit")
        run("kernel-plant-replay", "declaration type mismatch", mode="replay")
        restore()
        (root / "proofs/Checks/Forgotten.lean").write_text("module\n")
        run("unimported-module", "module inventory differs")
        restore()
        changed = json.loads((root / "proofs/theorems.json").read_text())
        changed["theorems"] = []
        save_inventory(changed)
        run("empty-inventory", "empty or unsupported theorem inventory")
        restore()
        changed = json.loads((root / "proofs/theorems.json").read_text())
        changed["theorems"].pop()
        save_inventory(changed)
        run("missing-theorem", "inventory differs from theorem ID markers")
        restore()
        edit("docs/semantics.md", "Both are monotone:", "Both are not monotone:")
        run("semantic-drift", "applicability digest drift")
        restore()
        with (root / "src/route_kernels.cpp").open("a") as stream:
            stream.write("\n// changed mapped input\n")
        run("native-drift", "applicability digest drift")
        restore()
        explanation = "docs/formal/theorems/round-robin.md"
        edit(explanation, '<a id="NGC-AP-004"></a>', "")
        changed = json.loads((root / "proofs/theorems.json").read_text())
        new_hash = hashlib.sha256((root / explanation).read_bytes()).hexdigest()
        for record in changed["theorems"]:
            record["mapped_inputs"]["explanation"]["sha256"] = new_hash
        save_inventory(changed)
        run("missing-explanation-anchor", "missing explanation anchor")
        restore()
        text = (root / MODEL).read_text()
        start = text.index("/-- The share computed")
        end = text.index("-/", start) + 2
        (root / MODEL).write_text(text[:start] + text[end:])
        run("missing-doc-comment", "missing declaration doc comment")
        restore()
        edit(
            MODEL,
            "end\nend NetGraphCore.Apportionment",
            "private def forgotten : Nat := 0\n\nend\nend NetGraphCore.Apportionment",
        )
        run("private-helper-doc", "missing declaration doc comment")
        restore()
        for label, declaration in [
            (
                "protected-helper-doc",
                "protected theorem undocumented : (1 : Nat) = 1 := rfl",
            ),
            ("indented-helper-doc", "  theorem undocumented : (1 : Nat) = 1 := rfl"),
            ("instance-doc", "instance undocumented : Inhabited (Fin 1) := ⟨0⟩"),
            ("structure-doc", "structure Undocumented where\n  value : Nat"),
            ("opaque-doc", "opaque undocumented : Nat := 0"),
            (
                "where-helper-doc",
                "/-- The outer value has a doc, but its written helper does not. -/\ndef undocumentedWhere : Nat := aux 1\nwhere\n  aux (n : Nat) : Nat := n",
            ),
        ]:
            edit(
                MODEL,
                "end\nend NetGraphCore.Apportionment",
                declaration + "\n\nend\nend NetGraphCore.Apportionment",
            )
            run(label, "missing declaration doc comment")
            restore()
        with (root / "proofs/NetGraphCore.lean").open("a") as stream:
            stream.write(
                "\nnamespace NetGraphCore\npublic def undocumented : Nat := 0\nend NetGraphCore\n"
            )
        run("root-helper-doc", "missing declaration doc comment")
        restore()
        with (root / MODEL).open("a") as stream:
            stream.write(
                "\nnamespace Elsewhere\n/-- An owned helper in the wrong namespace. -/\npublic def misplaced : Nat := 0\nend Elsewhere\n"
            )
        run("misplaced-namespace", "mathematical declaration outside project namespace")
        restore()
        edit(
            MODEL,
            "end\nend NetGraphCore.Apportionment",
            "/-- A documented protected helper, referencing NGC-AP-001 in prose. -/\n"
            "protected theorem documented : (1 : Nat) = 1 := rfl\n"
            "  /-- An indented helper with its own explanation. -/\n"
            "  theorem documentedIndented : (1 : Nat) = 1 := rfl\n\n"
            "/-- A documented record; its constructor is generated by Lean. -/\n"
            "structure Documented where\n  /-- The stored whole-unit count. -/\n  value : Nat\n\n"
            "/-- A documented value with a written local helper. -/\n"
            "def documentedWhere : Nat := aux 1\nwhere\n"
            "  /-- Return the supplied whole-unit count unchanged. -/\n"
            "  aux (n : Nat) : Nat := n\n\n"
            "end\nend NetGraphCore.Apportionment",
        )
        run("documented-variants")
        restore()
        main_path = root / "proofs/Checks/Main.lean"
        text = main_path.read_text()
        start = text.index("public unsafe def main")
        main_path.write_text(
            text[:start]
            + "public unsafe def main (_args : List String) : IO UInt32 := pure 0\n"
        )
        # Even an explicitly rebound, disabled audit must not turn a missing report into green.
        changed = json.loads((root / "proofs/theorems.json").read_text())
        changed["audit_sources"]["Checks/Main.lean"] = hashlib.sha256(
            main_path.read_bytes()
        ).hexdigest()
        save_inventory(changed)
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
        "Proof gate controls: PASS; original mapped inputs unchanged; bridge rejects wrong order, type and distribution."
    )


if __name__ == "__main__":
    main()
