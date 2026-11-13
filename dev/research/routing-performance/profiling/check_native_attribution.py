"""Separating controls for same-wheel range attribution and batch denominator."""

import json
from pathlib import Path

from native_attribution import batch_context, classify, verify_witnesses
from original_analyze_native_sample import Node

BASE = Path(__file__).resolve().parents[1]


def frame(text):
    return Node(1, 0, 1, text, None)


def native_pc(pc):
    return frame(
        f"???  (in _netgraph_core.test)  load address 0x100000000 + {hex(pc)}  [0x{0x100000000 + pc:x}]"
    )


solver = frame(
    "netgraph::core::detail::solve_min_cost_flow(  (in _netgraph_core.test) + 2480"
)
cases = []
for pc, expected in [
    (0xAB34, "solver_residual_dijkstra"),
    (0xA75C, "solver_potentials"),
    (0xAFA0, "solver_blocking_flow_dinic"),
    (0xB840, "solver_residual_closure_cut"),
    (0xA31C, "solver_setup_output_other"),
]:
    assert classify([solver, native_pc(pc)]) == (expected, True)
    cases.append({"pc": hex(pc), "category": expected})
assert classify([solver, native_pc(0xA75C), native_pc(0xAB34)]) == (
    "solver_residual_dijkstra",
    True,
)
assert classify([native_pc(0x13325C)]) == ("binding_conversion_publication", True)
evaluator = frame(
    "netgraph::core::detail::evaluate_pairs_in_call(  (in _netgraph_core.test) + 10"
)
assert classify([native_pc(0x13325C), evaluator, native_pc(0x52E84)]) == (
    "unresolved_native_bindings",
    True,
)
assert batch_context([native_pc(0x13325C)])
assert not batch_context([native_pc(0x12CB08)])
assert not batch_context(
    [frame("???  (in other.so) load address 0x100000000 + 0x13325c")]
)
assert classify([native_pc(0x12CB08)]) == ("binding_conversion_publication", True)
# Complete real partitions and scaled denominator must conserve their counts.
for workload in ("nsfnet", "backbone"):
    x = json.loads((BASE / f"native-{workload}-analysis-03/summary.json").read_text())
    assert sum(row["samples"] for row in x["categories"]) == x["samples"]
    assert sum(x["evaluate_pairs_categories"].values()) == x["evaluate_pairs_samples"]
    if workload == "nsfnet":
        value = x["nsfnet_derived_per_evaluation"]
        assert (
            abs(sum(value["categories_us"].values()) - value["total_us_per_evaluation"])
            < 1e-12
        )
        assert (
            abs(
                value["solver_total_us"]
                + value["non_solver_us"]
                - value["total_us_per_evaluation"]
            )
            < 1e-12
        )
report = {
    "phase_controls": cases,
    "instruction_witnesses": len(verify_witnesses()["instruction_witnesses"]),
    "summary": "phase distinctions, closest-ancestor priority, binding/native separation, exact batch denominator and real partition/time conservation passed",
}
with (BASE / "native-attribution-controls-02.json").open("x") as stream:
    json.dump(report, stream, indent=2)
print(json.dumps(report), flush=True)
