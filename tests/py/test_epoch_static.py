"""Static epoch preparation is reused as a value, without retaining mutable epochs."""

import itertools

import numpy as np
import pytest
from helpers import ALL_OUT, build, dag_of, plain
from private_hooks import needs_hooks

import netgraph_core as ngc


def cases():
    for pooled, allocation, tie, multipath, parallel, lossy in itertools.product(
        (False, True),
        ("equal", "capacity", "explicit"),
        ("FEWER_HOPS", "HIGHER_CAPACITY"),
        (False, True),
        (False, True),
        (False, True),
    ):
        if tie == "HIGHER_CAPACITY" and (
            (multipath and parallel) or (pooled and not multipath)
        ):
            continue
        yield pooled, allocation, tie, multipath, parallel, lossy


def fixture(pooled, allocation, tie, multipath, parallel, lossy):
    g = build(
        6,
        [
            (0, 2, 20, 1),
            (1, 2, 20, 1),
            (2, 3, 4, 1),
            (2, 3, 10, 1),
            (2, 4, 20, 1),
            (3, 5, 20, 1),
            (4, 5, 20, 1),
        ],
    )
    kw = dict(
        tie_break=getattr(ngc.RouteTieBreak, tie),
        multipath=multipath,
        parallel_edges=parallel,
    )
    if allocation != "equal":
        kw.update(allocation=ngc.RouteAllocation.WEIGHTED)
        if allocation == "capacity":
            kw.update(weight_basis=ngc.WeightBasis.CAPACITY)
        else:
            kw.update(
                weight_basis=ngc.WeightBasis.EXPLICIT,
                edge_weights=np.arange(7, dtype=np.int64) % 3,
            )
    cfg = ngc.FlowPolicyConfig(
        routing=ngc.ShortestPaths(**kw),
        admission=ngc.Admission.LOSSY if lossy else ngc.Admission.LOSSLESS,
    )
    ledger = ngc.FlowLedger(g)
    p = ngc.FlowPolicy(ledger, [0, 1] if pooled else [0], 5, cfg)
    return g, ledger, p


@pytest.mark.parametrize("case", list(cases()))
def test_static_identity_all_preparation_points_and_unchanged_values(case):
    g, ledger, p = fixture(*case)
    assert p.evaluate(0).program is not None
    assert p.program is None and not p.records()
    assert p.place(0).program is None
    first = p.place(6, ALL_OUT)
    before_ids = [r.flow_id for r in p.records()]
    assert p.evaluate(6, outputs=ALL_OUT).program is first.program
    assert [r.flow_id for r in p.records()] == before_ids
    p.remove()
    assert p.program is None
    assert p.place(0).program is None
    with pytest.raises(ValueError, match="forward requires an installed program"):
        p.forward(1)
    second = p.place(6, ALL_OUT)
    assert second.program is first.program
    assert plain(second) == plain(first)
    assert all(r.flow_id[1] > before_ids[-1][1] for r in p.records())
    replaced = p.replace(6, outputs=ALL_OUT)
    assert replaced.program is first.program
    assert not replaced.routes_changed and not replaced.weights_changed
    p.remove()
    ordered = ngc.apply_actions(ledger, [ngc.Place(p, 6)], outputs=ALL_OUT)[0]
    assert ordered.program is first.program and plain(ordered) == plain(first)
    ordered = ngc.apply_actions(ledger, [ngc.Replace(p, 6)], outputs=ALL_OUT)[0]
    assert ordered.program is first.program and plain(ordered) == plain(replaced)
    if case[-1]:
        simultaneous = ngc.apply_actions(
            ledger,
            [ngc.Replace(p, 6)],
            mode=ngc.BatchMode.SIMULTANEOUS,
            outputs=ALL_OUT,
        )[0]
        assert simultaneous.program is first.program
        assert plain(simultaneous) == plain(replaced)
    del g  # policy and ledger keep the graph owner
    p.remove()
    assert ledger.records() == []


@pytest.mark.parametrize("mode", ("available", "headroom", "snapshot", "optimized"))
def test_excluded_preparation_is_new_and_reads_changed_residual(mode):
    g = build(2, [(0, 1, 10, 1), (0, 1, 20 if mode == "headroom" else 10, 1)])
    kw = {}
    if mode in ("available", "headroom"):
        kw.update(
            parallel_edges=False,
            tie_break=getattr(
                ngc.RouteTieBreak,
                "HIGHER_AVAILABLE_CAPACITY"
                if mode == "available"
                else "HIGHER_HEADROOM",
            ),
        )
    elif mode == "snapshot":
        kw.update(
            allocation=ngc.RouteAllocation.WEIGHTED,
            weight_basis=ngc.WeightBasis.RESIDUAL_SNAPSHOT,
        )
    else:
        kw.update(allocation=ngc.RouteAllocation.OPTIMIZED)
    ledger = ngc.FlowLedger(g)
    p = ngc.FlowPolicy(
        ledger, 0, 1, ngc.FlowPolicyConfig(routing=ngc.ShortestPaths(**kw))
    )
    first = p.place(2)
    program = ngc.ForwardingProgram.from_routes(
        g, [0], [1], [ngc.RouteDefinition(0, dag_of(g, [0]))]
    )
    other = ngc.FlowPolicy(
        ledger, 0, 1, ngc.FlowPolicyConfig(routing=ngc.Pinned(program))
    )
    other.place(1 if mode == "headroom" else 8)
    second = p.replace(1 if mode in ("available", "headroom") else 12)
    assert second.program is not first.program
    if mode in ("available", "headroom"):
        assert list(p.records()[0].route.via_edges()) == [1]
    else:
        assert second.delivered == 12
    p.remove()
    assert p.place(1).program is not first.program


@needs_hooks
@pytest.mark.parametrize("pooled", (False, True))
def test_reused_prepare_trace_is_non_null_and_equal_to_fresh(pooled):
    _, _, p = fixture(pooled, "equal", "FEWER_HOPS", True, True, False)
    first, trace = p._place_traced(3, ALL_OUT)
    p.remove()
    second, repeated = p._place_traced(3, ALL_OUT)
    assert second.program is first.program
    assert plain(trace)["prepare"] is not None
    assert plain(trace)["prepare"] == plain(repeated)["prepare"]
    _, increment = p._place_traced(1, ALL_OUT)
    assert plain(increment)["prepare"] is None
