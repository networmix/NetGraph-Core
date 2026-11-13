"""Conservative self-weight partition for one requested one-worker capture.

Stripped solver subphases remain solver_subphase_unresolved until matching wheel
call sites are identified. Never assign that work to per-pair fixed costs.
"""

import argparse
import hashlib
import json
import re
from collections import Counter
from pathlib import Path

from evidence import equal, normalize
from original_analyze_native_sample import parse

BASE = Path(__file__).resolve().parents[1]
EFFORT = BASE.parents[1]


def classify(chain):
    frames = [node.frame for node in chain]
    native = any("(in _netgraph_core." in frame for frame in frames)
    if not native:
        if (
            any("textiowrapper_readline" in f.lower() for f in frames)
            and frames[-1].split("  (in ", 1)[0] == "read"
        ):
            return "readiness_wait", False
        return "python_runtime_outside_core", False
    # Deep named branches take priority over their outer caller.
    for category, tokens in [
        (
            "solver_residual_search",
            ("MinCostFlowSolver::dijkstra", "MinCostFlowSolver::close_sources"),
        ),
        (
            "solver_blocking_flow_dinic",
            (
                "MinCostFlowSolver::drain",
                "MinCostFlowSolver::build_levels",
                "MinCostFlowSolver::find_path",
            ),
        ),
        (
            "solver_potentials",
            ("MinCostFlowSolver::settle", "MinCostFlowSolver::update_all"),
        ),
        ("solver_subphase_unresolved", ("solve_min_cost_flow(", "MinCostFlowSolver::")),
        (
            "roster_support_search",
            (
                "detail::static_costs(",
                "detail::static_support(",
                "detail::tunnel_roster(",
                "detail::run_search(",
            ),
        ),
        (
            "result_record_construction",
            (
                "detail::derived_entry(",
                "LedgerAccess::commit_batch(",
                "PlacementResult::",
                "FlowRecord::",
            ),
        ),
        (
            "ledger_policy_construction_cleanup",
            (
                "FlowLedger::FlowLedger(",
                "FlowPolicy::FlowPolicy(",
                "PolicyAccess::prepare(",
                "PolicyAccess::construct(",
                "FlowPolicy::~FlowPolicy(",
                "FlowLedger::~FlowLedger(",
                "LedgerAccess::detach(",
                "FlowPolicy::remove(",
            ),
        ),
        (
            "binding_conversion",
            (
                "pybind11::detail::type_caster",
                "pybind11::cast",
                "pybind11::detail::list_caster",
            ),
        ),
    ]:
        if any(token in frame for frame in frames for token in tokens):
            return category, True
    return "unresolved_native_bindings", True


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--capture", type=Path, required=True)
    parser.add_argument("--workload", choices=("backbone", "nsfnet"), required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    text = (args.capture / "sample.txt").read_text()
    nodes = parse(text)
    roots = [n for n in nodes if n.parent is None]
    assert len(roots) == 1 and "com.apple.main-thread" in roots[0].frame
    grouped = Counter()
    leaves = []
    function_self = Counter()
    for i, node in enumerate(nodes):
        own = node.count - sum(nodes[c].count for c in node.children)
        if not own:
            continue
        chain = []
        cursor = i
        while cursor is not None:
            chain.append(nodes[cursor])
            cursor = nodes[cursor].parent
        chain.reverse()
        category, native = classify(chain)
        grouped[category] += own
        function_self[node.frame.split("  (in ", 1)[0]] += own
        leaves.append(
            {
                "sample_line": node.line,
                "self_samples": own,
                "category": category,
                "native_context": native,
                "stack_lines": [n.line for n in chain],
                "stack": [n.frame for n in chain],
            }
        )
    assert sum(grouped.values()) == roots[0].count
    # Independent top-of-stack check on every visible named native function.
    top = text.split("Sort by top of stack, same collapsed", 1)[1].split(
        "Binary Images:", 1
    )[0]
    crosschecks = {}
    for line in top.splitlines():
        match = re.fullmatch(r"\s*(netgraph::.+?)  \(in .+?\)\s+(\d+)\s*", line)
        if match:
            name, count = match[1], int(match[2])
            assert function_self[name] == count, (name, count, function_self[name])
            crosschecks[name] = count
    expected = (
        BASE / "backbone-wall-1-01/results.json"
        if args.workload == "backbone"
        else EFFORT / "scratch/orient-consumers/nsfnet-1-02/results.json"
    )
    equal(
        normalize(json.loads(expected.read_text())),
        normalize(json.loads((args.capture / "results.json").read_text())),
    )
    index = json.loads((args.capture / "index.json").read_text())
    assert index["sample_exit"] == index["worker_exit"] == 0
    native_total = sum(row["self_samples"] for row in leaves if row["native_context"])
    rows = [
        {
            "category": category,
            "samples": count,
            "share_all": count / roots[0].count,
            "share_native": count / native_total
            if category not in ("readiness_wait", "python_runtime_outside_core")
            else None,
        }
        for category, count in grouped.items()
    ]
    args.out.mkdir()
    summary = {
        "exploratory": True,
        "workload": args.workload,
        "samples": roots[0].count,
        "native_samples": native_total,
        "categories": rows,
        "top_of_stack_checks": crosschecks,
        "sample_sha256": hashlib.sha256(
            (args.capture / "sample.txt").read_bytes()
        ).hexdigest(),
        "summary": "self weights conserved; named top-of-stack crosschecks and complete export comparison passed; stripped solver subphases require call-site attribution",
    }
    (args.out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    with (args.out / "self-weights.jsonl").open("x") as stream:
        for row in leaves:
            stream.write(json.dumps(row) + "\n")
        stream.write(
            json.dumps(
                {
                    "summary": "every self weight assigned once",
                    "samples": roots[0].count,
                }
            )
            + "\n"
        )
    print(json.dumps(summary), flush=True)


if __name__ == "__main__":
    main()
