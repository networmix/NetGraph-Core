"""Small real Core streams that separate fields, ownership and failure state."""

import argparse
import gc
from pathlib import Path

import netgraph_core as core

from events import Recorder, provenance


def small_calls():
    graph = core.Graph(4, [0, 1, 0, 2], [1, 3, 2, 3], [10, 10, 2, 2], [1, 1, 1, 1], 1.0)
    edge_ids = graph.input_positions()
    assert edge_ids.tolist() == [0, 2, 1, 3]
    mask = graph.eligibility(excluded_edges=[1])
    assert core.Eligibility.intersect(mask, mask).uid() == mask.uid()
    ledger = core.FlowLedger(graph)
    cfg = core.FlowPolicyConfig(routing=core.ShortestPaths())
    policy = core.FlowPolicy(ledger, 0, 3, cfg)
    selections = [
        {},
        *[
            {flag: True}
            for flag in (
                "edge_flows",
                "residual",
                "saturated_boundary_edges",
                "edge_reservations",
                "min_cut",
                "program",
                "by_source",
            )
        ],
        {"program": True, "by_source": True},
    ]
    for selection in selections:
        result = core.apply_actions(
            ledger,
            [core.Replace(policy, float("inf"))],
            outputs=core.OutputOptions(**selection),
        )[0]
        assert result.delivered == 5
        assert ledger.carried().tolist() == [3, 2, 3, 2]
        policy.remove()
    try:
        ledger.reset()
    except ValueError as exc:
        assert (
            str(exc) == "FlowLedger.reset is rejected while a policy is bound (1 bound)"
        )
    else:
        raise AssertionError("reset accepted while bound")
    del policy
    ledger.reset()
    wanted = core.OutputOptions(
        **dict.fromkeys(
            (
                "edge_flows",
                "residual",
                "saturated_boundary_edges",
                "edge_reservations",
                "min_cut",
                "program",
                "by_source",
            ),
            True,
        )
    )
    controller = core.FlowPolicyConfig(
        routing=core.CentralizedTE(source_allocation=core.SourceAllocation.OPTIMIZED)
    )
    policy = core.FlowPolicy(ledger, 0, 3, controller)
    result = policy.place(float("inf"), wanted)
    assert (result.delivered, result.total_cost, result.min_cut) == (12, 24, [0, 1])
    try:
        policy.place(1)
    except ValueError as exc:
        assert (
            str(exc) == "finite placement after infinite placement requires a new epoch"
        )
    else:
        raise AssertionError("finite placement accepted after infinity")
    policy.remove()
    assert policy.place(1).delivered == 1
    del policy
    # Orphans remain after owner destruction, then reset retires them.
    assert len(ledger.records()) == 1
    entry = ledger.records()[0]
    assert ledger.record(entry.flow_id).flow_id == entry.flow_id
    del entry
    ledger.reset()
    cache = core.ShortestPathCache(graph)
    opts = core.EvaluationOptions(cfg, cache=cache, outputs=wanted)
    assert len(core.evaluate_pairs(graph, [(0, 3), (1, 3)], opts, workers=1)) == 2
    assert (cache.stats().hits, cache.stats().misses, cache.stats().entries) == (
        1,
        1,
        1,
    )
    cache.freeze()
    core.evaluate_pairs(graph, [(0, 3), (1, 3)], opts, workers=2)
    try:
        core.evaluate_demand(graph, 0, 2, opts)
    except core.ShortestPathCacheMiss:
        pass
    else:
        raise AssertionError("frozen cache miss accepted")
    cache.thaw()
    cache.clear()
    for routing in (
        core.DistributedTE(core.TunnelCount(2)),
        core.CentralizedTE(
            source_allocation=core.SourceAllocation.EQUAL,
            tunnel_count=core.TunnelCount(2),
        ),
    ):
        policy = core.FlowPolicy(ledger, 0, 3, core.FlowPolicyConfig(routing=routing))
        result = policy.replace(4, design_demand=4, outputs=wanted)
        assert result.plan is not None and result.program is not None
        policy.remove()
        del policy
    core.evaluate_pairs(graph, [], core.EvaluationOptions(controller), workers=2)
    excluded = graph.eligibility(excluded_nodes=[0])
    core.evaluate_demand(
        graph,
        0,
        3,
        core.EvaluationOptions(controller, eligibility=excluded, outputs=wanted),
    )
    search = core.spf(
        graph,
        range(1),
        core.SpfOptions(core.Direction.FORWARD, residual=ledger.residual()),
    )
    core.resolve_to_paths(graph, search.dag, 0, 3)
    core.ksp(graph, 0, 3, core.KspOptions(3))
    core.evaluate_edge_failures(
        graph, 0, 3, core.EvaluationOptions(controller, outputs=wanted), workers=2
    )

    # Cumulative Lossy counters cross 64 bits without narrowing per-call amounts.
    wide = core.Graph(2, [0], [1], [2**61], [1], 2**-12)
    wide_ledger = core.FlowLedger(wide)
    lossy = core.FlowPolicy(
        wide_ledger,
        0,
        1,
        core.FlowPolicyConfig(
            routing=core.ShortestPaths(), admission=core.Admission.LOSSY
        ),
    )
    for _ in range(10):
        wide_result = lossy.place(2**61)
    assert wide_result.offered == 10 * 2**61
    assert wide_result.delivered == 2**61
    assert wide_result.edge_drops == [(0, 9 * 2**61)]
    lossy.remove()
    priced = core.Graph(2, [0], [1], [2**60], [17], 1.0)
    priced_result = core.evaluate_demand(
        priced, 0, 1, core.EvaluationOptions(controller)
    )
    assert priced_result.total_cost == 17 * 2**60
    fractional = core.Graph(2, [0], [1], [3], [1], 1.0)
    fractional_policy = core.FlowPolicy(
        core.FlowLedger(fractional),
        0,
        1,
        core.FlowPolicyConfig(
            routing=core.CentralizedTE(
                source_allocation=core.SourceAllocation.EQUAL,
                tunnel_count=core.TunnelCount(2),
            )
        ),
    )
    fraction_result = fractional_policy.replace(float("inf"), outputs=wanted)
    from fractions import Fraction

    assert fraction_result.plan.max_fractional_bandwidth_per_tunnel == Fraction(3, 2)
    fractional_policy.remove()


def lookup_preflight():
    """Check every record field on hand-derived one- and two-source placements."""
    checked = 0
    for heads, cost in ((1, 2), (2, 1)):
        graph = core.Graph(
            heads + 1,
            list(range(heads)),
            [heads] * heads,
            [5] * heads,
            [cost] * heads,
            1.0,
        )
        ledger = core.FlowLedger(graph)
        policy = core.FlowPolicy(
            ledger,
            list(range(heads)),
            heads,
            core.FlowPolicyConfig(routing=core.ShortestPaths()),
        )
        assert policy.place(3 * heads).delivered == 3 * heads
        for owner in (ledger, policy):
            records = owner.records()
            assert len(records) == heads
            for source, record in enumerate(records):
                value = owner.record(record.flow_id)
                assert type(value.flow_id) is tuple and value.flow_id == record.flow_id
                assert all(type(x) is int for x in value.flow_id)
                assert value.flow_id[0] == ledger.uid()
                assert (
                    type(value.policy_id) is int and value.policy_id == record.policy_id
                )
                assert type(value.sources) is list and type(value.destinations) is list
                assert value.sources == [source] and value.destinations == [heads]
                assert all(type(x) is int for x in value.sources + value.destinations)
                assert type(value.admission) is core.Admission
                assert (
                    value.admission == core.Admission.LOSSLESS
                    and value.entry_id is None
                )
                assert (
                    type(value.min_route_cost) is int and type(value.delivered) is int
                )
                assert (
                    value.min_route_cost == cost
                    and value.offered is None
                    and value.delivered == 3
                )
                assert (
                    value.edge_flows == [(source, 3)]
                    and value.edge_drops == []
                    and value.edge_reservations == []
                )
                assert (
                    type(value.edge_flows) is list
                    and type(value.edge_drops) is list
                    and type(value.edge_reservations) is list
                )
                assert type(value.edge_flows[0]) is tuple
                assert all(type(x) is int for x in value.edge_flows[0])
                assert (
                    type(value.route.graph_uid()) is int
                    and value.route.graph_uid() == graph.uid()
                )
                assert value.route.parent_offsets().tolist() == [0] * (heads + 1) + [1]
                assert value.route.parents().tolist() == [source]
                assert value.route.via_edges().tolist() == [source]
                for array in (
                    value.route.parent_offsets(),
                    value.route.parents(),
                    value.route.via_edges(),
                ):
                    assert array.dtype.name == "int32" and not array.flags.writeable
                checked += 1
        assert ledger.carried().tolist() == [3] * heads
        assert ledger.residual().tolist() == [2] * heads
        assert ledger.reserved().tolist() == [0] * heads
        assert ledger.unreserved().tolist() == [5] * heads
        policy.remove()
        assert ledger.records() == [] and policy.records() == []
    return {
        "lookup_records": checked,
        "summary": "PASS single-source and two-source lookup preflight, all fields",
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("out", type=Path)
    args = parser.parse_args()
    provenance(core)
    recorder = Recorder(core, args.out)
    recorder.install()
    recorder.phase = "run"
    small_calls()
    lookup_preflight()
    gc.collect()
    recorder.finish({"workload": "controls", "workers": 2})


if __name__ == "__main__":
    main()
