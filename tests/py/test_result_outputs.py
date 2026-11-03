"""Completed-result retention is independent of live policy state and numeric answers."""

import itertools

import numpy as np
import pytest
from helpers import build, plain

import netgraph_core as ngc

FLAGS = tuple(itertools.product((False, True), repeat=2))


def fixture(lossy=False):
    graph = build(2, [(0, 1, 10, 3)])
    ledger = ngc.FlowLedger(graph)
    config = ngc.FlowPolicyConfig(
        routing=ngc.ShortestPaths(),
        admission=ngc.Admission.LOSSY if lossy else ngc.Admission.LOSSLESS,
    )
    return graph, ledger, ngc.FlowPolicy(ledger, 0, 1, config), config


def source_values(result):
    rows = result.by_source
    assert rows is not None
    return [
        (
            r.source,
            r.offered,
            r.admitted,
            r.delivered,
            r.source_drops,
            r.total_edge_drops,
        )
        for r in rows
    ]


def aggregates(result):
    value = plain(result)
    del value["program"], value["by_source"]
    return value


def record_values(records):
    fields = (
        "flow_id",
        "policy_id",
        "sources",
        "destinations",
        "admission",
        "entry_id",
        "route",
        "min_route_cost",
        "offered",
        "delivered",
        "edge_flows",
        "edge_drops",
        "edge_reservations",
    )
    return [
        tuple(plain(getattr(record, name)) for name in fields) for record in records
    ]


def invoke(kind, program, by_source):
    lossy = kind == "simultaneous"
    graph, ledger, policy, config = fixture(lossy)
    outputs = ngc.OutputOptions(program=program, by_source=by_source)
    if kind in ("place", "replace", "evaluate", "forward"):
        if kind == "forward":
            policy.replace(2)
        result = getattr(policy, kind)(4, outputs=outputs)
    elif kind in ("ordered", "simultaneous"):
        mode = getattr(ngc.BatchMode, kind.upper())
        result = ngc.apply_actions(
            ledger, [ngc.Replace(policy, 4)], mode=mode, outputs=outputs
        )[0]
    else:
        options = ngc.EvaluationOptions(config, offer=4, outputs=outputs)
        if kind == "demand":
            result = ngc.evaluate_demand(graph, 0, 1, options)
        elif kind.startswith("pairs"):
            result = ngc.evaluate_pairs(
                graph, [(0, 1)], options, workers=int(kind[-1])
            )[0]
        else:
            sensitivity = ngc.evaluate_edge_failures(
                graph, 0, 1, options, workers=int(kind[-1])
            )
            assert [(r.edge, r.delivery_reduction) for r in sensitivity.entries] == [
                (0, 4)
            ]
            result = sensitivity.baseline
    return result, policy, lossy


@pytest.mark.parametrize(
    "kind",
    (
        "place",
        "replace",
        "evaluate",
        "forward",
        "ordered",
        "simultaneous",
        "demand",
        "pairs1",
        "pairs2",
        "failures1",
        "failures2",
    ),
)
def test_every_producer_honors_each_flag_without_changing_aggregates(kind):
    # On cap10/cost3, offer4 delivers4 without a drop. Flags only choose retention.
    full, _, _ = invoke(kind, True, True)
    for program, by_source in FLAGS:
        result, policy, lossy = invoke(kind, program, by_source)
        assert aggregates(result) == aggregates(full)
        assert result.delivered == result.admitted == result.requested == 4
        assert result.call_delivered == 4
        assert (
            result.call_undelivered
            == result.call_source_drops
            == result.source_drops
            == 0
        )
        assert result.edge_drops == result.call_edge_drops == []
        assert (result.program is not None) is program
        if by_source:
            assert source_values(result) == [(0, 4 if lossy else None, 4, 4, 0, 0)]
        else:
            assert result.by_source is None
        for name in (
            "edge_flows",
            "residual",
            "saturated_boundary_edges",
            "min_cut",
            "edge_reservations",
        ):
            assert getattr(result, name) is None
        if kind in ("place", "replace", "forward", "ordered", "simultaneous"):
            assert policy.program is not None
            assert policy.records()[0].delivered == 4


def test_omitted_options_zero_calls_and_later_opt_in_keep_live_state():
    _, ledger, policy, _ = fixture()
    for program, by_source in FLAGS:
        outputs = ngc.OutputOptions(program=program, by_source=by_source)
        result = policy.place(0, outputs=outputs)
        assert result.program is None
        assert (result.by_source is not None) is by_source
        if by_source:
            assert source_values(result) == [(0, None, 0, 0, 0, 0)]
    result = policy.replace(4)
    assert result.program is result.by_source is None
    installed = policy.program
    records = record_values(policy.records())
    for program, by_source in FLAGS:
        result = policy.place(
            0, outputs=ngc.OutputOptions(program=program, by_source=by_source)
        )
        assert result.delivered == 4 and result.call_delivered == 0
        assert result.delivery_delta == 0
        assert not (
            result.routes_changed
            or result.weights_changed
            or result.reservations_changed
        )
        assert (result.program is installed) if program else result.program is None
        assert (result.by_source is not None) is by_source
        assert (
            record_values(policy.records())
            == records
            == record_values(ledger.records())
        )
    snapshot = policy.replace(
        4, outputs=ngc.OutputOptions(program=True, by_source=True)
    )
    policy.remove()
    del policy, ledger, installed
    assert snapshot.program is not None
    assert snapshot.program.sources == (0,)
    assert source_values(snapshot) == [(0, None, 4, 4, 0, 0)]


@pytest.mark.parametrize("sources", (1, 2))
def test_lossy_rows_are_epoch_totals_even_after_calls_without_rows(sources):
    if sources == 1:
        _, _, policy, _ = fixture(True)
        expected_rows = [(0, 12, 12, 10, 0, 2)]
        expected_delivery, expected_drops = 10, [(0, 2)]
    else:
        graph = build(4, [(0, 2, 10, 1), (1, 2, 10, 1), (2, 3, 6, 1)])
        policy = ngc.FlowPolicy(
            ngc.FlowLedger(graph),
            [0, 1],
            3,
            ngc.FlowPolicyConfig(
                routing=ngc.ShortestPaths(), admission=ngc.Admission.LOSSY
            ),
        )
        # First RR(5)=(3,2) delivers all five. The next increments are (3,4),
        # and the remaining shared unit goes to source1: final delivery (3,3).
        expected_rows = [(0, 6, 6, 3, 0, 3), (1, 6, 6, 3, 0, 3)]
        expected_delivery, expected_drops = 6, [(2, 6)]
    first = policy.place(5)
    assert first.program is first.by_source is None
    result = policy.place(7, outputs=ngc.OutputOptions(by_source=True))
    assert (result.offered, result.delivered, result.edge_drops) == (
        12,
        expected_delivery,
        expected_drops,
    )
    rows = source_values(result)
    assert rows == expected_rows
    assert sum(row[1] for row in rows) == result.offered
    assert sum(row[2] for row in rows) == result.admitted
    assert sum(row[3] for row in rows) == result.delivered
    assert sum(row[4] for row in rows) == result.source_drops
    assert sum(row[5] for row in rows) == sum(value for _, value in result.edge_drops)
    assert result.program is None
    after = policy.place(0)
    assert after.by_source is after.program is None
    assert aggregates(after)["offered"] == 12
    assert sum(record.delivered for record in policy.records()) == expected_delivery


@pytest.mark.parametrize("name", ("program", "by_source"))
def test_new_flags_use_strict_booleans_and_are_read_only(name):
    for value in (True, False, np.bool_(True), np.bool_(False)):
        outputs = ngc.OutputOptions(**{name: value})
        assert getattr(outputs, name) is bool(value)
        with pytest.raises(AttributeError):
            setattr(outputs, name, False)
    for value in (0, 1, 0.0, 1.0, None, "true"):
        with pytest.raises(TypeError):
            ngc.OutputOptions(**{name: value})


def test_new_flags_do_not_change_refusal_order_or_state():
    _, _, policy, _ = fixture()
    for program, by_source in FLAGS:
        outputs = ngc.OutputOptions(program=program, by_source=by_source)
        with pytest.raises(ValueError) as exc:
            policy.forward(-1, outputs=outputs)
        assert (
            str(exc.value) == "offer must be an integer in [0, 2^62) or permitted +inf"
        )
        assert policy.program is None and policy.records() == []
        with pytest.raises(ValueError) as exc:
            policy.forward(4, outputs=outputs)
        assert str(exc.value) == "forward requires an installed program"
    assert policy.replace(4).delivered == 4


MODELS = (
    "equal",
    "capacity",
    "explicit",
    "snapshot",
    "optimized",
    "optimized_equal",
    "global",
    "global_equal",
    "global_weighted",
    "supported",
    "t1",
    "t2",
    "routes",
    "dag",
)


def model_config(graph, name, lossy):
    dag = ngc.PredecessorDAG.from_csr(graph, [0, 0, 1], [0], [0])
    if name in ("equal", "capacity", "explicit", "snapshot"):
        kw = {}
        if name != "equal":
            kw["allocation"] = ngc.RouteAllocation.WEIGHTED
        if name == "explicit":
            kw.update(weight_basis=ngc.WeightBasis.EXPLICIT, edge_weights=[1])
        if name == "snapshot":
            kw["weight_basis"] = ngc.WeightBasis.RESIDUAL_SNAPSHOT
        routing = ngc.ShortestPaths(**kw)
    elif name.startswith("optimized"):
        kw = dict(allocation=ngc.RouteAllocation.OPTIMIZED)
        if name.endswith("equal"):
            kw["source_allocation"] = ngc.SourceAllocation.EQUAL
        routing = ngc.ShortestPaths(**kw)
    elif name in ("global", "global_equal", "global_weighted", "supported", "t2"):
        kw = dict(source_allocation=ngc.SourceAllocation.OPTIMIZED)
        if name in ("global_equal", "t2"):
            kw["source_allocation"] = ngc.SourceAllocation.EQUAL
        if name == "global_weighted":
            kw.update(
                source_allocation=ngc.SourceAllocation.WEIGHTED, source_weights=[1, 0]
            )
        if name == "supported":
            kw["support"] = dag
        if name == "t2":
            kw["tunnel_count"] = ngc.TunnelCount(1)
        routing = ngc.CentralizedTE(**kw)
    elif name == "t1":
        routing = ngc.DistributedTE(ngc.TunnelCount(1))
    else:
        if name == "dag":
            program = ngc.ForwardingProgram.from_dag(graph, 0, 1, dag)
        else:
            program = ngc.ForwardingProgram.from_routes(
                graph, 0, 1, [ngc.RouteDefinition(0, ngc.Path(graph, [[0]]))]
            )
        routing = ngc.Pinned(program)
    return ngc.FlowPolicyConfig(
        routing=routing,
        admission=ngc.Admission.LOSSY if lossy else ngc.Admission.LOSSLESS,
    )


@pytest.mark.parametrize("name", MODELS)
@pytest.mark.parametrize("operation", ("replace", "evaluate"))
def test_every_model_retains_aggregate_values_and_requested_rows(name, operation):
    graph = build(2, [(0, 1, 10, 3)])
    for lossy in (False, True):
        if lossy and name in (
            "optimized",
            "optimized_equal",
            "global",
            "global_equal",
            "global_weighted",
            "supported",
        ):
            continue
        config = model_config(graph, name, lossy)
        expected = None
        for program, by_source in FLAGS:
            ledger = ngc.FlowLedger(graph)
            policy = ngc.FlowPolicy(ledger, 0, 1, config)
            extra = {"design_demand": 4} if name in ("t1", "t2") else {}
            result = getattr(policy, operation)(
                4,
                outputs=ngc.OutputOptions(program=program, by_source=by_source),
                **extra,
            )
            assert result.delivered == result.admitted == result.requested == 4
            assert result.source_drops == 0 and result.edge_drops == []
            if expected is None:
                expected = aggregates(result)
            assert aggregates(result) == expected
            assert (result.program is not None) is program
            assert (result.by_source is not None) is by_source
            if by_source:
                assert source_values(result) == [(0, 4 if lossy else None, 4, 4, 0, 0)]
            if name in ("t1", "t2"):
                assert result.plan is not None
                assert len(result.by_tunnel) == 1
                assert result.by_tunnel[0].delivered == 4
            if operation == "evaluate":
                assert policy.program is None and policy.records() == []
            else:
                assert policy.program is not None
                assert sum(r.delivered for r in policy.records()) == 4


def test_masked_roster_keeps_zero_rows_only_when_requested():
    graph = build(2, [(0, 1, 10, 3)])
    mask = graph.eligibility(excluded_nodes=[0])
    config = ngc.FlowPolicyConfig(
        routing=ngc.ShortestPaths(), admission=ngc.Admission.LOSSY
    )
    for program, by_source in FLAGS:
        policy = ngc.FlowPolicy(ngc.FlowLedger(graph), 0, 1, config, mask)
        result = policy.replace(
            4, outputs=ngc.OutputOptions(program=program, by_source=by_source)
        )
        assert result.offered == result.admitted == result.source_drops == 4
        assert result.delivered == 0 and result.edge_drops == []
        assert (result.program is not None) is program
        assert (result.by_source is not None) is by_source
        if by_source:
            assert source_values(result) == [(0, 0, 0, 0, 0, 0)]


def test_two_cyclic_joint_participants_honor_each_flag():
    # The two shared links carry H(5,[21,5])=(4,1), H(6,[4,15])=(1,5).
    # A's source deliveries are (1,0), B's (1,0); each source's rest drops.
    graph = build(
        10,
        [
            (0, 6, 100, 1),
            (1, 6, 100, 1),
            (2, 4, 100, 1),
            (3, 4, 100, 1),
            (6, 7, 5, 1),
            (7, 4, 100, 1),
            (4, 5, 6, 1),
            (5, 6, 100, 1),
            (5, 8, 100, 1),
            (7, 9, 100, 1),
        ],
    )
    config = ngc.FlowPolicyConfig(
        routing=ngc.ShortestPaths(), admission=ngc.Admission.LOSSY
    )
    expected_rows = [
        [(0, 11, 11, 1, 0, 10), (1, 10, 10, 0, 0, 10)],
        [(2, 8, 8, 1, 0, 7), (3, 7, 7, 0, 0, 7)],
    ]
    expected_aggregates = None
    for program, by_source in FLAGS:
        ledger = ngc.FlowLedger(graph)
        policies = [
            ngc.FlowPolicy(ledger, [0, 1], 8, config),
            ngc.FlowPolicy(ledger, [2, 3], 9, config),
        ]
        results = ngc.apply_actions(
            ledger,
            [ngc.Replace(policies[0], 21), ngc.Replace(policies[1], 15)],
            mode=ngc.BatchMode.SIMULTANEOUS,
            outputs=ngc.OutputOptions(program=program, by_source=by_source),
        )
        assert len(results) == 2
        values = [aggregates(result) for result in results]
        if expected_aggregates is None:
            expected_aggregates = values
        assert values == expected_aggregates
        for result, rows, offer in zip(results, expected_rows, (21, 15), strict=True):
            assert result.call_delivered == result.delivered == 1
            assert result.offered == result.admitted == offer
            assert result.source_drops == result.call_source_drops == 0
            assert sum(amount for _, amount in result.edge_drops) == offer - 1
            assert (result.program is not None) is program
            assert (result.by_source is not None) is by_source
            if by_source:
                assert source_values(result) == rows
        assert [record.delivered for record in ledger.records()] == [1, 0, 1, 0]
        assert all(policy.program is not None for policy in policies)


@pytest.mark.parametrize("lossy", (False, True))
def test_wide_epoch_totals_do_not_depend_on_output_flags(lossy):
    maximum = (1 << 62) - 1
    requested = 5 * maximum  # 23058430092136939515, above 64 bits.
    for program, by_source in FLAGS:
        _, _, policy, _ = fixture(lossy)
        for _ in range(4):
            policy.place(maximum)
        result = policy.place(
            maximum,
            outputs=ngc.OutputOptions(program=program, by_source=by_source),
        )
        assert result.requested == requested
        assert result.delivered == 10
        assert result.admitted == (requested if lossy else 10)
        assert result.offered == (requested if lossy else None)
        assert result.unadmitted == (0 if lossy else requested - 10)
        assert result.source_drops == 0
        assert result.edge_drops == ([(0, requested - 10)] if lossy else [])
        assert (result.program is not None) is program
        assert (result.by_source is not None) is by_source
        if by_source:
            assert source_values(result) == [
                (
                    0,
                    requested if lossy else None,
                    requested if lossy else 10,
                    10,
                    0,
                    requested - 10 if lossy else 0,
                )
            ]
        assert policy.records()[0].delivered == 10
