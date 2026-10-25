"""Packaged proof inventory and evidence checks. Python never proves a theorem.

Lean supplies structural declarations and transitive axioms from its built environment.
This file binds them to reviewed inputs, explanations and independently recorded examples.
It uses only Python's standard library and works without Git or repository developer tools.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import subprocess
import sys
from pathlib import Path

PROOFS = Path(__file__).resolve().parents[1]
ROOT = PROOFS.parent
ALLOWED_AXIOMS = {"propext", "Classical.choice", "Quot.sound"}
AUDIT_INPUTS = ["Checks/Main.lean", "Checks/Examples.lean", "Checks/audit.py"]


def fail(message: str) -> None:
    raise ValueError(message)


def digest(data: bytes | str) -> str:
    return hashlib.sha256(data.encode() if isinstance(data, str) else data).hexdigest()


def canonical(value: object) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True)


def unique_object(pairs: list[tuple[str, object]]) -> dict:
    result = {}
    for key, value in pairs:
        if key in result:
            fail(f"duplicate JSON key: {key}")
        result[key] = value
    return result


def read_json(path: Path) -> dict:
    return json.loads(path.read_text(), object_pairs_hook=unique_object)


def source_file(name: str) -> Path:
    path = ROOT / name
    if (
        not path.is_file()
        or path.is_symlink()
        or not path.resolve().is_relative_to(ROOT.resolve())
    ):
        fail(f"missing or external mapped input: {name}")
    return path


def section_text(name: str, heading: str) -> str:
    lines = source_file(name).read_text().splitlines(keepends=True)
    starts = [i for i, line in enumerate(lines) if line.rstrip() == heading]
    if len(starts) != 1:
        fail(f"semantic heading missing or ambiguous: {name}: {heading}")
    start = starts[0]
    level = len(heading) - len(heading.lstrip("#"))
    end = len(lines)
    fence = None
    for i in range(start + 1, len(lines)):
        marker = re.match(r"^ {0,3}(`{3,}|~{3,})", lines[i])
        if fence:
            if re.fullmatch(
                r" {0,3}" + re.escape(fence[0]) + "{" + str(fence[1]) + r",}\s*",
                lines[i],
            ):
                fence = None
            continue
        if marker:
            fence = (marker[1][0], len(marker[1]))
            continue
        match = re.match(r"^ {0,3}(#+) ", lines[i])
        if match and len(match[1]) <= level:
            end = i
            break
    return "".join(lines[start:end])


def lean_sources() -> list[str]:
    result = []
    for directory, children, files in os.walk(PROOFS):
        children[:] = [name for name in children if name != ".lake"]
        for name in files:
            if name.endswith(".lean"):
                result.append(str((Path(directory) / name).relative_to(PROOFS)))
    return sorted(result)


def stripped_lean(text: str) -> str:
    """Remove strings and nested comments before checking source commands."""
    result, i, depth = [], 0, 0
    while i < len(text):
        if depth:
            if text.startswith("/-", i):
                depth += 1
                i += 2
            elif text.startswith("-/", i):
                depth -= 1
                i += 2
            else:
                result.append("\n" if text[i] == "\n" else " ")
                i += 1
        elif text.startswith("/-", i):
            depth, i = 1, i + 2
        elif text.startswith("--", i):
            end = text.find("\n", i)
            i = len(text) if end < 0 else end
        elif text[i] == '"':
            i += 1
            while i < len(text):
                if text[i] == "\\":
                    i += 2
                elif text[i] == '"':
                    i += 1
                    break
                else:
                    i += 1
            result.append(" ")
        else:
            result.append(text[i])
            i += 1
    if depth:
        fail("unclosed Lean comment")
    return "".join(result)


def input_state(inventory: dict) -> dict:
    if inventory.get("schema") != 1 or not inventory.get("theorems"):
        fail("empty or unsupported theorem inventory")
    if lean_sources() != sorted(inventory["modules"]):
        fail("maintained Lean module inventory differs (including unimported sources)")
    # Every maintained library module must be explicitly imported by its aggregation root.
    for root, prefix in [
        ("NetGraphCore.lean", "NetGraphCore/"),
        ("Checks.lean", "Checks/"),
    ]:
        expected = sorted(
            name[:-5].replace("/", ".")
            for name in inventory["modules"]
            if name.startswith(prefix)
        )
        imports = re.findall(
            r"^(?:public )?import ([A-Za-z0-9_.]+)\s*$",
            (PROOFS / root).read_text(),
            re.M,
        )
        if sorted(imports) != expected:
            fail(f"aggregation root does not cover maintained modules: {root}")
    inputs = {}
    for name in [
        "lean-toolchain",
        "lakefile.toml",
        "lake-manifest.json",
        *AUDIT_INPUTS,
    ]:
        inputs[f"proofs/{name}"] = digest((PROOFS / name).read_bytes())
    for name in inventory["modules"]:
        path = PROOFS / name
        inputs[f"proofs/{name}"] = digest(path.read_bytes())
        if name.startswith("NetGraphCore") or name == "Checks/Examples.lean":
            code = stripped_lean(path.read_text())
            if re.search(
                r"\b(sorry|admit|axiom|unsafe|native_decide)\b|\bdebug\.|#eval\b|^\s*(?:set_option|initialize|run_cmd|elab|macro|attribute)\b",
                code,
                re.M,
            ):
                fail(f"forbidden proof command or kernel bypass: {name}")
            if any(
                value.strip() != "expose"
                for value in re.findall(r"@\[([^\]]*)\]", code)
            ):
                fail(f"unreviewed mathematical declaration attribute: {name}")
    return inputs


def record_inputs(record: dict, inputs: dict) -> dict:
    result = {"semantics": [], "implementation": [], "explanation": {}}
    for ref in record["semantics"]:
        value = digest(section_text(ref["file"], ref["heading"]))
        result["semantics"].append({**ref, "sha256": value})
        inputs[ref["file"]] = digest(source_file(ref["file"]).read_bytes())
    for ref in record["implementation"]["sources"]:
        value = digest(source_file(ref["file"]).read_bytes())
        result["implementation"].append({**ref, "sha256": value})
        inputs[ref["file"]] = value
    ref = record["explanation"]
    text = source_file(ref["file"]).read_text()
    if f'<a id="{ref["anchor"]}"></a>' not in text:
        fail(f"missing explanation anchor: {ref['anchor']}")
    result["explanation"] = {**ref, "sha256": digest(text)}
    inputs[ref["file"]] = digest(text)
    return result


def preflight(inventory: dict, recording: bool = False) -> dict:
    inputs = input_state(inventory)
    pins = {
        name: inputs[f"proofs/{name}"]
        for name in ["lean-toolchain", "lakefile.toml", "lake-manifest.json"]
    }
    tools = {name: inputs[f"proofs/{name}"] for name in AUDIT_INPUTS}
    if not recording and (
        inventory.get("pins") != pins or inventory.get("audit_sources") != tools
    ):
        fail("toolchain/dependency/audit-source digest drift")
    ids, names = set(), set()
    for record in inventory["theorems"]:
        if not record["semantics"]:
            fail(f"accepted theorem has no semantic contract reference: {record['id']}")
        if record["id"] in ids or record["declaration"] in names:
            fail("duplicate accepted ID or declaration")
        ids.add(record["id"])
        names.add(record["declaration"])
        for field in ("module", "claim", "domain"):
            if not isinstance(record.get(field), str) or not record[field].strip():
                fail(f"missing theorem {field}: {record['id']}")
        mapping = record_inputs(record, inputs)
        if not recording and mapping != record.get("mapped_inputs"):
            fail(f"explanation or applicability digest drift: {record['id']}")
        if recording:
            record["mapped_inputs"] = mapping
        if record["implementation"]["status"] not in {"unconnected", "tested"}:
            fail("unsupported native connection claim; no refinement is supplied")
        if (
            record["implementation"]["status"] == "tested"
            and record["implementation"].get("bridge") != "round-robin"
        ):
            fail("connected declaration has no supported bridge")
    if recording:
        inventory["pins"], inventory["audit_sources"] = pins, tools
    return inputs


def binding(declaration: str, metadata: dict) -> dict:
    if declaration not in metadata:
        fail(f"unresolved Lean declaration: {declaration}")
    root = metadata[declaration]
    if root["value"] is None:
        fail(f"proof body unavailable to audit: {declaration}")
    axioms = sorted(root["axioms"])
    if set(axioms) - ALLOWED_AXIOMS:
        fail(f"unaccepted transitive axioms: {declaration}: {axioms}")
    todo, seen, closure = [declaration], set(), {}
    while todo:
        name = todo.pop()
        if name in seen or name not in metadata:
            continue
        seen.add(name)
        item = metadata[name]
        fields = {"kind": item["kind"], "type": item["type"]}
        if item["kind"] != "theorem":
            if item["value"] is None:
                fail(f"definition body unavailable to audit: {name}")
            fields["value"] = item["value"]
        closure[name] = digest(canonical(fields))
        todo.extend(item["uses"])
    if root["kind"] != "theorem":
        fail(f"accepted declaration is not a theorem: {declaration}")
    return {
        "type_sha256": digest(root["type"]),
        "local_closure": closure,
        "axioms": axioms,
    }


def audit(data: dict, inventory: dict, recording: bool) -> dict:
    inputs = preflight(inventory, recording)
    toolchain = (PROOFS / "lean-toolchain").read_text().strip()
    if not re.fullmatch(r"leanprover/lean4:v\d+\.\d+\.\d+", toolchain) or data[
        "lean_version"
    ] != toolchain.removeprefix("leanprover/lean4:v"):
        fail("unexpected Lean runtime or toolchain")
    dependencies = read_json(PROOFS / "lake-manifest.json")["packages"]
    for package in dependencies:
        directory = PROOFS / ".lake/packages" / package["name"]
        actual = subprocess.run(
            ["git", "-C", str(directory), "rev-parse", "HEAD"],
            capture_output=True,
            text=True,
            check=True,
        ).stdout.strip()
        dirty = subprocess.run(
            [
                "git",
                "-C",
                str(directory),
                "status",
                "--porcelain",
                "--untracked-files=normal",
            ],
            capture_output=True,
            text=True,
            check=True,
        ).stdout
        if actual != package["rev"] or dirty:
            fail(f"dependency source drift: {package['name']}")
    metadata = {item["name"]: item for item in data["declarations"]}
    mathematical_modules = {
        name[:-5].replace("/", ".")
        for name in inventory["modules"]
        if name == "NetGraphCore.lean" or name.startswith("NetGraphCore/")
    }
    tagged = []
    for item in data["declarations"]:
        if not item["authored"]:
            continue
        name = item["user_name"]
        if not item["doc"] or not item["doc"].strip():
            fail(f"missing declaration doc comment: {name}")
        if item["source_module"] not in mathematical_modules:
            continue
        if not name.startswith("NetGraphCore."):
            fail(f"mathematical declaration outside project namespace: {name}")
        marker = re.match(r"^(NGC-[A-Z]+-\d{3}):", item["doc"].strip())
        if marker:
            tagged.append((marker[1], item["name"]))
    expected = [(r["id"], r["declaration"]) for r in inventory["theorems"]]
    if sorted(tagged) != sorted(expected):
        fail("accepted inventory differs from theorem ID markers in built declarations")
    expected_examples = read_json(PROOFS / "cases.json")
    if canonical(data["examples"]) != canonical(expected_examples):
        fail("model example inputs or exact outputs differ from independent cases")
    inputs["proofs/cases.json"] = digest((PROOFS / "cases.json").read_bytes())
    case_ids = {item["id"] for item in expected_examples["cases"]}
    if len(case_ids) != len(expected_examples["cases"]):
        fail("duplicate executable example ID")
    for record in inventory["theorems"]:
        if not record["examples"] or set(record["examples"]) - case_ids:
            fail(f"missing executable explanation example: {record['id']}")
        if record["module"] != metadata[record["declaration"]]["source_module"]:
            fail(f"theorem module differs from built declaration: {record['id']}")
        actual = binding(record["declaration"], metadata)
        if recording:
            record["binding"] = actual
        elif actual != record.get("binding"):
            fail(f"reviewed statement/definition digest drift: {record['id']}")
    revision, dirty = None, None
    if (ROOT / ".git").exists():
        revision = subprocess.run(
            ["git", "-C", str(ROOT), "rev-parse", "HEAD"],
            capture_output=True,
            text=True,
            check=True,
        ).stdout.strip()
        dirty = bool(
            subprocess.run(
                [
                    "git",
                    "-C",
                    str(ROOT),
                    "status",
                    "--porcelain",
                    "--untracked-files=normal",
                ],
                capture_output=True,
                text=True,
                check=True,
            ).stdout
        )
    return {
        "status": "PASS",
        "scope": "model proofs; no native refinement",
        "lean_version": data["lean_version"],
        "inventory_sha256": digest(canonical(inventory)),
        "source_revision": revision,
        "worktree_dirty": dirty,
        "source_inputs": dict(sorted(inputs.items())),
        "modules": inventory["modules"],
        "accepted": [record["id"] for record in inventory["theorems"]],
        "bindings": {
            record["id"]: record["binding"] for record in inventory["theorems"]
        },
        "example_count": len(case_ids),
        "public_sequence_steps": len(expected_examples["public_sequence"]),
        "dependencies": {p["name"]: p["rev"] for p in dependencies},
        "dependency_origin": os.environ.get("PROOFS_DEPENDENCY_ORIGIN", "unrecorded"),
        "native_bridge": "not run here",
        "kernel_replay": "not run by this command",
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--preflight", action="store_true")
    parser.add_argument("--record", action="store_true")
    parser.add_argument("--verify-report", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    inventory_path = PROOFS / "theorems.json"
    inventory = read_json(inventory_path)
    if args.verify_report:
        if not args.verify_report.is_file():
            fail("missing proof evidence file")
        inputs = preflight(inventory)
        inputs["proofs/cases.json"] = digest((PROOFS / "cases.json").read_bytes())
        report = read_json(args.verify_report)
        if (
            report.get("status") != "PASS"
            or report.get("accepted") != [r["id"] for r in inventory["theorems"]]
            or report.get("inventory_sha256") != digest(canonical(inventory))
            or report.get("source_inputs") != dict(sorted(inputs.items()))
        ):
            fail("missing, incomplete or stale proof evidence")
        print("Current proof evidence: PASS")
        return
    if args.preflight:
        preflight(inventory)
        print("Proof input preflight: PASS")
        return
    if args.record and (os.environ.get("CI") or os.environ.get("PROOFS_RECORD") != "1"):
        fail("recording requires explicit PROOFS_RECORD=1 outside CI")
    data = json.load(sys.stdin, object_pairs_hook=unique_object)
    report = audit(data, inventory, args.record)
    if args.record:
        inventory_path.write_text(
            json.dumps(inventory, indent=2, sort_keys=True) + "\n"
        )
        print("Recorded candidate bindings; independent review is still required.")
        return
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps(report, sort_keys=True))


if __name__ == "__main__":
    try:
        main()
    except (ValueError, KeyError, TypeError, OSError) as error:
        print(f"proof-audit: {error}", file=sys.stderr)
        sys.exit(1)
