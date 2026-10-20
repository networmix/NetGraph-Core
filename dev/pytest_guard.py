"""The Python suite with the check of what it imported and skipped: dev/run-checks.sh runs it for
`pytest on` and `pytest off`, and the Windows CI jobs run it directly (no shell there):

    python dev/pytest_guard.py on|off [PYTEST_ARGS...]

on: the extension carries the private hooks (the editable development build) and no test skips.
off: an installation without them (a wheel, or a build from the sdist), imported from outside
the source tree, and exactly the tests that need a hook skip (EXPECTED_OFF_SKIPS, and every
replay of the reference fixtures). A small plugin records the skipped node ids (tests, and modules
or classes skipped while collected) and the imported modules; the arguments are pytest's. Exit status: pytest's when it fails, else 1 when a check
fails (each problem printed), else 0."""

from __future__ import annotations

import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

# The Python tests that need the private hooks and skip without them, by pytest node id; every
# replay of the fixtures produced by the reference under tests/reference/policy/ (the
# parametrizations of test_fixture) skips as well.
EXPECTED_OFF_SKIPS = [
    "tests/py/test_apportionment.py::test_apportionment_entries",
    (
        "tests/py/test_shortest_path_cache.py::test_frozen_statistics_survive_python_allocation_failures"
    ),
    "tests/py/test_shortest_path_cache.py::test_the_frozen_sweep_rejects_a_wrong_call",
    "tests/py/test_flow_ledger.py::test_a_flow_id_read_passes_an_allocation_failure_through",
    "tests/py/test_flow_policy.py::test_evaluator_results_survive_python_allocation_failures",
    "tests/py/test_flow_policy.py::test_owner_hook_is_the_owner_of_the_records",
    "tests/py/test_flow_policy.py::test_python_publication_is_all_or_nothing[False-costonly]",
    "tests/py/test_flow_policy.py::test_python_publication_is_all_or_nothing[False-global]",
    "tests/py/test_flow_policy.py::test_python_publication_is_all_or_nothing[True-costonly]",
    "tests/py/test_flow_policy.py::test_python_publication_is_all_or_nothing[True-global]",
    "tests/py/test_flow_policy_fixtures.py::test_a_corrupted_expectation_fails[ledger]",
    "tests/py/test_flow_policy_fixtures.py::test_a_corrupted_expectation_fails[policy]",
    "tests/py/test_flow_policy_fixtures.py::test_a_corrupted_expectation_fails[result]",
    "tests/py/test_flow_policy_fixtures.py::test_a_corrupted_expectation_fails[trace]",
    (
        "tests/py/test_flow_policy_fixtures.py::test_an_expected_value_of_another_type_fails[bool as integer]"
    ),
    (
        "tests/py/test_flow_policy_fixtures.py::test_an_expected_value_of_another_type_fails[integer as float]"
    ),
    "tests/py/test_flow_policy_fixtures.py::test_every_other_replay_clause_reports_its_operation",
    "tests/py/test_getter_allocation.py::test_every_read_survives_python_allocation_failures",
    "tests/py/test_getter_allocation.py::test_every_number_read_raises_memory_error",
    "tests/py/test_getter_allocation.py::test_a_failed_integer_read_raises_memory_error",
    "tests/py/test_ksp.py::test_python_result_lists_are_all_or_nothing",
    (
        "tests/py/test_number_types.py::test_the_converter_raises_memory_error_on_every_python_allocation_failure"
    ),
    "tests/py/test_number_types.py::test_a_cold_conversion_raises_memory_error",
    "tests/py/test_native_allocation.py::test_allocation_hooks_preserve_exception_identity",
    (
        "tests/py/test_native_allocation.py::test_native_allocation_failure_raises_and_changes_nothing[Eligibility.intersect]"
    ),
    (
        "tests/py/test_native_allocation.py::test_native_allocation_failure_raises_and_changes_nothing[FlowLedger.reset]"
    ),
    (
        "tests/py/test_native_allocation.py::test_native_allocation_failure_raises_and_changes_nothing[FlowRecord.route]"
    ),
    (
        "tests/py/test_native_allocation.py::test_native_allocation_failure_raises_and_changes_nothing[PlacementResult.program]"
    ),
    (
        "tests/py/test_native_allocation.py::test_native_allocation_failure_raises_and_changes_nothing[evaluate_pairs.global]"
    ),
    (
        "tests/py/test_native_allocation.py::test_native_allocation_failure_raises_and_changes_nothing[evaluate_pairs.global_equal]"
    ),
    (
        "tests/py/test_native_allocation.py::test_native_allocation_failure_raises_and_changes_nothing[capacity.frozen_cache]"
    ),
    (
        "tests/py/test_native_allocation.py::test_native_allocation_failure_raises_and_changes_nothing[capacity.global]"
    ),
    (
        "tests/py/test_native_allocation.py::test_native_allocation_failure_raises_and_changes_nothing[capacity.pinned_lossy]"
    ),
    (
        "tests/py/test_native_allocation.py::test_native_allocation_failure_raises_and_changes_nothing[evaluate.global]"
    ),
    (
        "tests/py/test_native_allocation.py::test_native_allocation_failure_raises_and_changes_nothing[forward.costonly_pooled_lossy_black_holes]"
    ),
    (
        "tests/py/test_native_allocation.py::test_native_allocation_failure_raises_and_changes_nothing[from_dag]"
    ),
    (
        "tests/py/test_native_allocation.py::test_native_allocation_failure_raises_and_changes_nothing[from_routes]"
    ),
    (
        "tests/py/test_native_allocation.py::test_native_allocation_failure_raises_and_changes_nothing[ksp]"
    ),
    (
        "tests/py/test_native_allocation.py::test_native_allocation_failure_raises_and_changes_nothing[place.costonly]"
    ),
    (
        "tests/py/test_native_allocation.py::test_native_allocation_failure_raises_and_changes_nothing[place.costonly_pooled_lossy]"
    ),
    (
        "tests/py/test_native_allocation.py::test_native_allocation_failure_raises_and_changes_nothing[place.global]"
    ),
    (
        "tests/py/test_native_allocation.py::test_native_allocation_failure_raises_and_changes_nothing[place.global_equal]"
    ),
    (
        "tests/py/test_native_allocation.py::test_native_allocation_failure_raises_and_changes_nothing[place.pinned_split_route]"
    ),
    (
        "tests/py/test_native_allocation.py::test_native_allocation_failure_raises_and_changes_nothing[place_traced.global]"
    ),
    (
        "tests/py/test_native_allocation.py::test_native_allocation_failure_raises_and_changes_nothing[replace.global_support]"
    ),
    (
        "tests/py/test_native_allocation.py::test_native_allocation_failure_raises_and_changes_nothing[replace.pinned_tunnels_lossy]"
    ),
    (
        "tests/py/test_native_allocation.py::test_native_allocation_failure_raises_and_changes_nothing[resolve_to_paths]"
    ),
    (
        "tests/py/test_native_allocation.py::test_native_allocation_failure_raises_and_changes_nothing[to_units.list]"
    ),
    (
        "tests/py/test_native_allocation.py::test_native_allocation_failure_raises_and_changes_nothing[sensitivity.global]"
    ),
    (
        "tests/py/test_native_allocation.py::test_native_allocation_failure_raises_and_changes_nothing[sensitivity.pinned]"
    ),
    (
        "tests/py/test_native_allocation.py::test_native_allocation_failure_raises_and_changes_nothing[spf.full]"
    ),
    (
        "tests/py/test_native_allocation.py::test_native_allocation_failure_raises_and_changes_nothing[spf.single_path]"
    ),
    (
        "tests/py/test_native_allocation.py::test_native_allocation_failure_raises_and_changes_nothing[spf.targets_residual]"
    ),
    "tests/py/test_native_allocation.py::test_the_row_checks_reject_a_wrong_call",
    "tests/py/test_programs.py::test_from_routes_and_getters_survive_python_allocation_failures",
    "tests/py/test_programs.py::test_operations_are_all_or_nothing[costonly forward]",
    "tests/py/test_programs.py::test_operations_are_all_or_nothing[costonly replace]",
    "tests/py/test_programs.py::test_operations_are_all_or_nothing[global evaluate]",
    "tests/py/test_programs.py::test_operations_are_all_or_nothing[pinned replace]",
    "tests/py/test_search.py::test_python_result_construction_is_all_or_nothing",
    "tests/py/test_search.py::test_the_python_allocation_sweep_rejects_a_wrong_call",
    "tests/py/test_distributed_te.py::test_traced_place_checks_the_operation_first",
    (
        "tests/py/test_distributed_te.py::test_t1_getters_survive_python_allocation_failures[ca.tunnel_count]"
    ),
    (
        "tests/py/test_distributed_te.py::test_t1_getters_survive_python_allocation_failures[fg.record(fid)]"
    ),
    (
        "tests/py/test_distributed_te.py::test_t1_getters_survive_python_allocation_failures[fg.records()]"
    ),
    (
        "tests/py/test_distributed_te.py::test_t1_getters_survive_python_allocation_failures[fg.reserved()]"
    ),
    (
        "tests/py/test_distributed_te.py::test_t1_getters_survive_python_allocation_failures[fg.unreserved()]"
    ),
    (
        "tests/py/test_distributed_te.py::test_t1_getters_survive_python_allocation_failures[p.record(fid)]"
    ),
    (
        "tests/py/test_distributed_te.py::test_t1_getters_survive_python_allocation_failures[p.records()]"
    ),
    (
        "tests/py/test_distributed_te.py::test_t1_getters_survive_python_allocation_failures[r.plan.design_by_head]"
    ),
    "tests/py/test_distributed_te.py::test_t1_getters_survive_python_allocation_failures[r.plan]",
    (
        "tests/py/test_distributed_te.py::test_t1_getters_survive_python_allocation_failures[r.program.tunnels]"
    ),
    (
        "tests/py/test_distributed_te.py::test_t1_getters_survive_python_allocation_failures[repr(abs3)]"
    ),
    (
        "tests/py/test_distributed_te.py::test_t1_getters_survive_python_allocation_failures[repr(st2)]"
    ),
    (
        "tests/py/test_distributed_te.py::test_t1_getters_survive_python_allocation_failures[st2.cost_ceiling.factor]"
    ),
    (
        "tests/py/test_distributed_te.py::test_t1_getters_survive_python_allocation_failures[st2.cost_ceiling]"
    ),
    (
        "tests/py/test_distributed_te.py::test_t1_getters_survive_python_allocation_failures[cfg.routing]"
    ),
    "tests/py/test_distributed_te.py::test_t1_getters_survive_python_allocation_failures[repr(r)]",
    "tests/py/test_distributed_te.py::test_t1_getters_survive_python_allocation_failures[repr(sp)]",
    (
        "tests/py/test_distributed_te.py::test_t1_getters_survive_python_allocation_failures[t0.status]"
    ),
    "tests/py/test_distributed_te.py::test_t1_results_survive_python_allocation_failures[evaluate]",
    "tests/py/test_distributed_te.py::test_t1_results_survive_python_allocation_failures[forward]",
    "tests/py/test_distributed_te.py::test_t1_results_survive_python_allocation_failures[replace]",
    "tests/py/test_centralized_te_tunnels.py::test_traced_place_checks_the_operation_first",
    (
        "tests/py/test_centralized_te_tunnels.py::test_t2_getters_survive_python_allocation_failures[gl.tunnel_count]"
    ),
    (
        "tests/py/test_centralized_te_tunnels.py::test_t2_getters_survive_python_allocation_failures[repr(r)]"
    ),
    (
        "tests/py/test_centralized_te_tunnels.py::test_t2_getters_survive_python_allocation_failures[r.plan.max_fractional_bandwidth_per_tunnel]"
    ),
    (
        "tests/py/test_centralized_te_tunnels.py::test_t2_getters_survive_python_allocation_failures[r.plan.headroom_penalty]"
    ),
    (
        "tests/py/test_centralized_te_tunnels.py::test_t2_getters_survive_python_allocation_failures[r.plan.heads_with_capacity]"
    ),
    (
        "tests/py/test_centralized_te_tunnels.py::test_t2_getters_survive_python_allocation_failures[r.plan.total_tunnel_cost]"
    ),
    (
        "tests/py/test_centralized_te_tunnels.py::test_t2_getters_survive_python_allocation_failures[r.plan]"
    ),
    (
        "tests/py/test_centralized_te_tunnels.py::test_t2_getters_survive_python_allocation_failures[r.program.tunnels]"
    ),
    (
        "tests/py/test_centralized_te_tunnels.py::test_t2_results_survive_python_allocation_failures[evaluate]"
    ),
    (
        "tests/py/test_centralized_te_tunnels.py::test_t2_results_survive_python_allocation_failures[forward]"
    ),
    (
        "tests/py/test_centralized_te_tunnels.py::test_t2_results_survive_python_allocation_failures[replace]"
    ),
]
EXPECTED_OFF_SKIPPED_PREFIX = "tests/py/test_flow_policy_fixtures.py::test_fixture["


class Record:
    """The run as the checks read it: the node ids seen, the skipped ones with their reasons."""

    def __init__(self):
        self.seen, self.skipped = set(), {}

    def skip(self, report):
        r = report.longrepr
        self.skipped[report.nodeid] = str(r[-1] if isinstance(r, tuple) else r)

    def pytest_runtest_logreport(self, report):
        self.seen.add(report.nodeid)
        if report.skipped:
            self.skip(report)

    def pytest_collectreport(self, report):
        # A module or a class skipped while collected (pytest.skip(allow_module_level=True),
        # pytest.importorskip) runs no test: its collection report is its only record.
        if report.skipped:
            self.skip(report)


def problems(hooks, has_hooks, modules, seen, skipped, root=ROOT):
    """The checks of one run (hooks: "on" or "off"; modules: the files of _netgraph_core and
    netgraph_core as imported, None when not imported); returns the problems found."""
    out = []
    if None in modules:
        out.append(f"netgraph_core or _netgraph_core was not imported: {modules}")
    if hooks == "on":
        if not has_hooks:
            out.append(
                "the extension lacks the private hooks (expected the editable build)"
            )
        if skipped:
            out.append(
                f"{len(skipped)} tests skipped in a hooks-on run: {sorted(skipped)}"
            )
        return out
    if has_hooks:
        out.append(
            "the extension carries the private hooks (expected a wheel or sdist build)"
        )
    inside = [
        m
        for m in modules
        if m and os.path.realpath(m).startswith(os.path.realpath(root) + os.sep)
    ]
    if inside:
        out.append(f"imported from the source tree, not an installation: {inside}")
    replays = {n for n in seen if n.startswith(EXPECTED_OFF_SKIPPED_PREFIX)}
    if not replays:
        out.append(f"no test {EXPECTED_OFF_SKIPPED_PREFIX}...] ran")
    want = set(EXPECTED_OFF_SKIPS) | replays
    got = set(skipped)
    out += [f"expected to skip, did not: {n}" for n in sorted(want - got)]
    out += [f"skipped, not a known hook-only test: {n}" for n in sorted(got - want)]
    out += [
        f"skipped for another reason: {n}: {r}"
        for n, r in sorted(skipped.items())
        if "needs the private hooks" not in r
    ]
    return out


def main(argv):
    if not argv or argv[0] not in ("on", "off"):
        print(
            "usage: python dev/pytest_guard.py on|off [PYTEST_ARGS...]", file=sys.stderr
        )
        return 2
    import pytest

    hooks, record = argv[0], Record()
    status = pytest.main(argv[1:], plugins=[record])
    if status != 0:
        return int(status)
    ext = sys.modules.get("_netgraph_core")
    pkg = sys.modules.get("netgraph_core")
    modules = [getattr(ext, "__file__", None), getattr(pkg, "__file__", None)]
    found = problems(
        hooks,
        hasattr(ext, "_fail_native_allocation"),
        modules,
        record.seen,
        record.skipped,
    )
    for p in found:
        print(f"dev/pytest_guard.py: {p}", file=sys.stderr)
    if found:
        return 1
    if hooks == "off":
        print(
            f"hooks off: exactly the {len(record.skipped)} expected hook-only tests skipped;"
            f" imported {modules[0]}"
        )
    else:
        print(f"hooks on: no test skipped; imported {modules[0]}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
