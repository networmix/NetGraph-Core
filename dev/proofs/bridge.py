"""Compare audited model cases with a specifically selected installed native wheel.

The public calls cover ordered source injection and cumulative offers beyond 2^64.
No claim of universal native refinement follows from these finite comparisons.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import importlib.util
import json
import sys
import tomllib
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def exact(actual: object, expected: object, label: str) -> None:
    if type(actual) is not type(expected):
        raise ValueError(f"{label}: scalar/container type differs")
    if isinstance(expected, list):
        if len(actual) != len(expected):
            raise ValueError(f"{label}: length differs")
        for i, (left, right) in enumerate(zip(actual, expected, strict=True)):
            exact(left, right, f"{label}[{i}]")
    elif actual != expected:
        raise ValueError(f"{label}: {actual!r} != {expected!r}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--proof-report", type=Path, required=True)
    parser.add_argument("--wheel", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    spec = importlib.util.spec_from_file_location(
        "proof_audit", ROOT / "proofs/Checks/audit.py"
    )
    audit = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(audit)
    inventory = audit.read_json(ROOT / "proofs/theorems.json")
    inputs = audit.preflight(inventory)
    inputs["proofs/cases.json"] = sha((ROOT / "proofs/cases.json").read_bytes())
    proof = audit.read_json(args.proof_report)
    if (
        proof.get("status") != "PASS"
        or proof.get("inventory_sha256") != audit.digest(audit.canonical(inventory))
        or proof.get("source_inputs") != dict(sorted(inputs.items()))
        or proof.get("accepted") != [r["id"] for r in inventory["theorems"]]
    ):
        raise ValueError("bridge requires complete current proof evidence")

    import _netgraph_core as native

    import netgraph_core as ngc

    expected_version = tomllib.loads((ROOT / "pyproject.toml").read_text())["project"][
        "version"
    ]
    if importlib.metadata.version("netgraph-core") != expected_version:
        raise ValueError("installed package version differs from this source tree")
    # Check bytes, not just a version string: the selected wheel is the one under test.
    with zipfile.ZipFile(args.wheel) as archive:
        for module in (ngc, native):
            path = Path(module.__file__).resolve()
            if not path.is_relative_to(Path(sys.prefix).resolve()):
                raise ValueError(
                    f"module is not installed in the selected environment: {path}"
                )
            member = "netgraph_core/__init__.py" if module is ngc else path.name
            if sha(path.read_bytes()) != sha(archive.read(member)):
                raise ValueError(
                    f"installed module differs from selected wheel: {member}"
                )
    cases = audit.read_json(ROOT / "proofs/cases.json")
    checked = []
    for case in cases["cases"]:
        if case["kind"] != "rr":
            continue
        k = case["members"]
        graph = ngc.Graph(k + 1, list(range(k)), [k] * k, [100] * k, [1] * k, 1.0)
        policy = ngc.FlowPolicy(
            ngc.FlowLedger(graph),
            list(reversed(range(k))),
            k,
            ngc.FlowPolicyConfig(
                routing=ngc.ShortestPaths(), admission=ngc.Admission.LOSSY
            ),
        )
        result = policy.place(case["total"])
        exact(
            [row.source for row in result.by_source],
            list(range(k)),
            case["id"] + " source order",
        )
        exact([row.offered for row in result.by_source], case["expected"], case["id"])
        checked.append(case["id"])

    graph = ngc.Graph(4, [0, 1, 2], [3, 3, 3], [1, 1, 1], [1, 2, 3], 1.0)
    sp = ngc.FlowPolicy(
        ngc.FlowLedger(graph),
        [2, 0, 1],
        3,
        ngc.FlowPolicyConfig(
            routing=ngc.ShortestPaths(), admission=ngc.Admission.LOSSY
        ),
    )
    program = sp.replace(0).program
    pinned = ngc.FlowPolicy(
        ngc.FlowLedger(graph),
        [2, 0, 1],
        3,
        ngc.FlowPolicyConfig(
            routing=ngc.Pinned(program), admission=ngc.Admission.LOSSY
        ),
    )
    policies = [("ShortestPaths", sp), ("Pinned", pinned)]
    for kind, policy in policies:
        previous = [0, 0, 0]
        for index, case in enumerate(cases["public_sequence"]):
            result = policy.place(case["offer"])
            exact(
                [row.source for row in result.by_source],
                [0, 1, 2],
                f"{kind} source order {index}",
            )
            actual = [row.offered for row in result.by_source]
            exact(result.offered, case["total"], f"{kind} total {index}")
            exact(actual, case["expected"], f"{kind} shares {index}")
            if index == 1:
                exact(
                    [a - b for a, b in zip(actual, previous, strict=True)],
                    [0, 0, 1],
                    f"{kind} continues the circle",
                )
            previous = actual
    if cases["public_sequence"][-1]["total"] <= 2**64:
        raise ValueError("public sequence no longer exercises the high word")
    report = {
        "status": "PASS",
        "scope": "finite native comparisons, no universal refinement",
        "proof_report_sha256": sha(args.proof_report.read_bytes()),
        "inventory_sha256": proof["inventory_sha256"],
        "source_inputs": proof["source_inputs"],
        "interpreter": sys.executable,
        "package": ngc.__file__,
        "extension": native.__file__,
        "version": importlib.metadata.version("netgraph-core"),
        "extension_sha256": sha(Path(native.__file__).read_bytes()),
        "wheel_sha256": sha(args.wheel.read_bytes()),
        "round_robin_cases": checked,
        "public_policies": [kind for kind, _ in policies],
        "sequence_steps_per_policy": len(cases["public_sequence"]),
        "largest_cumulative_offer": cases["public_sequence"][-1]["total"],
        "not_covered": [
            "all machine inputs",
            "mod_u128 refinement",
            "near-2^127 native execution",
        ],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(
        f"Native bridge: PASS ({len(checked)} RR cases, {len(policies)} policies, {len(cases['public_sequence'])} steps each)"
    )


if __name__ == "__main__":
    main()
