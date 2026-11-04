"""Check proof safety and exact model examples exported from Lean's built environment.

Python never proves a theorem. Lean builds the proofs and leanchecker replays them;
this checker rejects forbidden source commands and transitive axioms, and compares
executable model examples with independently derived cases.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import sys
from pathlib import Path
from typing import NoReturn

PROOFS = Path(__file__).resolve().parents[1]
ALLOWED_AXIOMS = {"propext", "Classical.choice", "Quot.sound"}


def fail(message: str) -> NoReturn:
    raise ValueError(message)


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


def preflight() -> list[str]:
    sources = lean_sources()
    # Every local proof/check module must be imported by its aggregation root.
    for root, prefix in [
        ("NetGraphCore.lean", "NetGraphCore/"),
        ("Checks.lean", "Checks/"),
    ]:
        expected = sorted(
            name[:-5].replace("/", ".") for name in sources if name.startswith(prefix)
        )
        imports = re.findall(
            r"^(?:public )?import ([A-Za-z0-9_.]+)\s*$",
            (PROOFS / root).read_text(),
            re.M,
        )
        if sorted(imports) != expected:
            fail(f"aggregation root does not cover local modules: {root}")
    for name in sources:
        if name.startswith("NetGraphCore") or name == "Checks/Examples.lean":
            code = stripped_lean((PROOFS / name).read_text())
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
    return sources


def checked_examples(examples: dict) -> dict:
    expected = read_json(PROOFS / "cases.json")
    if not expected.get("cases") or not expected.get("public_sequence"):
        fail("empty model examples")
    ids = [item["id"] for item in expected["cases"]]
    if len(set(ids)) != len(ids):
        fail("duplicate executable example ID")
    if canonical(examples) != canonical(expected):
        fail("model example inputs or exact outputs differ from independent cases")
    return expected


def audit(data: dict) -> dict:
    modules = {name[:-5].replace("/", ".") for name in preflight()}
    mathematical = {
        name
        for name in modules
        if name == "NetGraphCore" or name.startswith("NetGraphCore.")
    }
    theorems = {}
    seen = set()
    for item in data["declarations"]:
        name = item["name"]
        if name in seen or item["source_module"] not in modules:
            fail(f"duplicate or external project declaration: {name}")
        seen.add(name)
        if item["kind"] == "axiom":
            fail(f"project axiom: {name}")
        in_namespace = item["user_name"].startswith("NetGraphCore.")
        if item["source_module"] in mathematical and not in_namespace:
            fail(f"mathematical declaration outside project namespace: {name}")
        if not in_namespace:
            continue
        axioms = sorted(item["axioms"])
        if set(axioms) - ALLOWED_AXIOMS:
            fail(f"unaccepted transitive axioms: {name}: {axioms}")
        if item["kind"] == "theorem":
            theorems[name] = axioms
    if not theorems:
        fail("no project theorems in built environment")
    examples = checked_examples(data["examples"])
    return {
        "status": "PASS",
        "scope": "model proofs; no native refinement",
        "lean_version": data["lean_version"],
        "theorems": dict(sorted(theorems.items())),
        "examples": examples,
        "dependency_origin": os.environ.get("PROOFS_DEPENDENCY_ORIGIN", "unrecorded"),
        "native_bridge": "not run here",
        "kernel_replay": "not run by this command",
    }


def verify_report(path: Path) -> dict:
    if not path.is_file():
        fail("missing proof evidence file")
    report = read_json(path)
    theorems = report.get("theorems")
    if report.get("status") != "PASS" or not isinstance(theorems, dict) or not theorems:
        fail("missing or incomplete proof evidence")
    for name, axioms in theorems.items():
        if not name or not isinstance(axioms, list) or set(axioms) - ALLOWED_AXIOMS:
            fail("invalid theorem or transitive axioms in proof evidence")
    checked_examples(report["examples"])
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--preflight", action="store_true")
    parser.add_argument("--modules", action="store_true")
    parser.add_argument("--verify-report", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.modules:
        print(json.dumps([name[:-5].replace("/", ".") for name in preflight()]))
        return
    if args.preflight:
        preflight()
        print("Proof source preflight: PASS")
        return
    if args.verify_report:
        verify_report(args.verify_report)
        print("Proof evidence: PASS")
        return
    report = audit(json.load(sys.stdin, object_pairs_hook=unique_object))
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
