"""The policy engine, the DagCache and the evaluators through the bindings (APPROACH 8 to 11):
configuration objects and their structural validation, PlacementResult fields and outputs,
lifetimes, the evaluators with workers, and the cache statistics. Prototype 2's fixtures
replay in ``test_flow_policy_fixtures.py``."""

from __future__ import annotations

import gc
import math

import _netgraph_core as ext
import numpy as np
import pytest
from private_hooks import HOOKS, needs_hooks

import netgraph_core as ngc

INF = math.inf
L62 = 1 << 62
T = ngc.Termination


def build(n, edges, q=1.0):
    return ngc.Graph(
        n,
        [e[0] for e in edges],
        [e[1] for e in edges],
        [e[2] for e in edges],
        [e[3] for e in edges],
        q,
    )


def ids(g, positions):
    inv = [0] * g.num_edges()
    for e, pos in enumerate(g.input_index()):
        inv[int(pos)] = e
    return sorted(inv[p] for p in positions)


def dag(g, edges):
    per = [[] for _ in range(g.num_nodes())]
    for e in sorted(edges):
        per[int(g.edge_dst()[e])].append(e)
    offsets, parents, via = [0], [], []
    for v in range(g.num_nodes()):
        for e in per[v]:
            parents.append(int(g.edge_src()[e]))
            via.append(e)
        offsets.append(len(via))
    return ngc.PredDAG.from_csr(g, offsets, parents, via)


def sel(
    tie=ngc.EdgeTieBreak.DETERMINISTIC, prio=ngc.EdgePriority.CAPACITY, mp=True, me=True
):
    return ngc.PathSelection(tie, prio, mp, me)


def config(
    routing, placement=ngc.FlowPlacement.EQUAL_BALANCED, admission=None, selection=None
):
    return ngc.FlowPolicyConfig(
        routing=routing,
        placement=placement,
        admission=admission or ngc.Admission.LOSSLESS,
        selection=selection or sel(),
    )


GLOBAL = config(ngc.Global(), ngc.FlowPlacement.MAX_FLOW)


def diamond(a=5, b=3, c=1):
    # input edges: 0: 0->1, 1: 1->3 (cap a, cost 1); 2: 0->2, 3: 2->3 (cap b, cost c)
    return build(4, [(0, 1, a, 1), (1, 3, a, 1), (0, 2, b, c), (2, 3, b, c)])


def test_configurations_are_structural_and_explicit():
    with pytest.raises(TypeError):
        ngc.FlowPolicyConfig(
            ngc.CostOnly(), ngc.FlowPlacement.MAX_FLOW, ngc.Admission.LOSSLESS, sel()
        )
    with pytest.raises(TypeError):
        ngc.FlowPolicyConfig(
            routing=ngc.CostOnly(), placement=ngc.FlowPlacement.MAX_FLOW
        )
    with pytest.raises(TypeError, match="incompatible"):
        ngc.CostOnly(flow_count=2)  # a field outside its routing variant (rule 13)
    with pytest.raises(TypeError, match="incompatible"):
        ngc.Global(max_flow_count=4)
    with pytest.raises(TypeError, match="routing"):
        config("CostOnly")
    ca = ngc.CapacityAware(
        flow_count=2,
        max_flow_count=4,
        max_path_cost=ngc.Stretch(1.5),
        min_flow_bandwidth=3,
        reoptimize=True,
        diminishing=ngc.DiminishingReturns(2, 0.25),
    )
    assert (ca.flow_count, ca.max_flow_count, ca.min_flow_bandwidth, ca.reoptimize) == (
        2,
        4,
        3,
        True,
    )
    assert ca.max_path_cost == ngc.Stretch(1.5)
    assert (ca.diminishing.window, ca.diminishing.fraction) == (2, 0.25)
    assert ngc.CapacityAware().max_path_cost is None
    assert ngc.CapacityAware(
        max_path_cost=ngc.Absolute(7)
    ).max_path_cost == ngc.Absolute(7)
    with pytest.raises(TypeError):
        ngc.CapacityAware(2)  # keyword-only
    for bad in (2**40, True, 1.5):
        with pytest.raises(ValueError, match="flow_count"):
            ngc.CapacityAware(flow_count=bad)
    with pytest.raises(ValueError, match="min_flow_bandwidth"):
        ngc.CapacityAware(min_flow_bandwidth=2**80)
    with pytest.raises(TypeError, match="PredDAG"):
        ngc.Pinned([[0, 1]])
    cfg = config(ca)
    assert isinstance(cfg.routing, ngc.CapacityAware)
    assert cfg.routing.flow_count == 2


def test_validation_messages_before_binding():
    g = diamond()
    fg = ngc.FlowLedger(g)

    def rejects(cfg, msg, src=0, dst=3):
        with pytest.raises(ValueError, match=msg):
            ngc.FlowPolicy(fg, src, dst, cfg)

    rejects(
        config(ngc.CostOnly(), ngc.FlowPlacement.MAX_FLOW, ngc.Admission.LOSSY), "Lossy"
    )
    rejects(config(ngc.CapacityAware(), admission=ngc.Admission.LOSSY), "Lossy")
    rejects(config(ngc.Global()), "MaxFlow")
    rejects(config(ngc.CapacityAware(flow_count=0)), "flow_count")
    rejects(config(ngc.CapacityAware(flow_count=2, max_flow_count=1)), "max_flow_count")
    rejects(config(ngc.CapacityAware(max_path_cost=ngc.Stretch(0.5))), "Stretch")
    rejects(config(ngc.CapacityAware(max_path_cost=ngc.Stretch(math.nan))), "Stretch")
    rejects(
        config(ngc.CapacityAware(diminishing=ngc.DiminishingReturns(1, 1.0))),
        "fraction",
    )
    rejects(config(ngc.Pinned([])), "empty bundle list")
    b = dag(g, ids(g, [0, 1]))
    rejects(config(ngc.Pinned([b, b])), "duplicate bundle")
    rejects(config(ngc.Pinned([dag(g, ids(g, [0, 3]))])), "malformed bundle")
    rejects(config(ngc.Pinned([b])), "single endpoints", src=[0, 1])
    rejects(config(ngc.CostOnly()), "both endpoint sets", src=[0, 3])
    rejects(config(ngc.CostOnly()), "empty", src=[])
    rejects(config(ngc.CostOnly()), "out of node range", src=[L62])
    with pytest.raises(ValueError, match="foreign"):
        ngc.FlowPolicy(
            fg, 0, 3, config(ngc.CostOnly()), diamond().eligibility(excluded_edges=[0])
        )
    with pytest.raises(ValueError, match="foreign"):
        ngc.FlowPolicy(fg, 0, 3, config(ngc.CostOnly()), None, ngc.DagCache(diamond()))
    fg.reset()  # nothing bound


def test_enum_values_outside_their_enumeration_are_rejected():
    # Rule 13: an enum built from an integer outside the enumeration is a domain error named
    # by its field, at every constructor that takes it, so it never reaches a policy or record.
    with pytest.raises(
        ValueError, match="tie_break: 99 is not a value of EdgeTieBreak"
    ):
        ngc.PathSelection(ngc.EdgeTieBreak(99), ngc.EdgePriority.CAPACITY)
    with pytest.raises(ValueError, match="priority: 99 is not a value of EdgePriority"):
        ngc.PathSelection(ngc.EdgeTieBreak.DETERMINISTIC, ngc.EdgePriority(99))
    for routing in (ngc.CostOnly(), ngc.CapacityAware(), ngc.Global()):
        with pytest.raises(
            ValueError, match="placement: 99 is not a value of FlowPlacement"
        ):
            config(routing, ngc.FlowPlacement(99))
        with pytest.raises(
            ValueError, match="admission: 99 is not a value of Admission"
        ):
            config(routing, ngc.FlowPlacement.MAX_FLOW, ngc.Admission(99))
    for valid in ngc.FlowPlacement.__members__.values():
        assert config(ngc.CostOnly(), valid).placement == valid


def test_placement_result_fields_and_outputs():
    g = diamond()
    fg = ngc.FlowLedger(g)
    pol = ngc.FlowPolicy(fg, [0], [3], config(ngc.CostOnly()))
    r = pol.place(4)
    assert (r.placed, r.remaining, r.termination) == (4, 0, T.COMPLETE)
    assert r.local_outcomes == [(0, T.COMPLETE)]
    assert (r.routes_changed, r.weights_changed, r.drops) == (True, False, [])
    assert r.cost_distribution is None
    assert (
        r.edge_flows is None
        and r.residual is None
        and r.bottleneck is None
        and r.min_cut is None
    )
    r = pol.place(
        INF, ngc.OutputOptions(edge_flows=True, residual=True, bottleneck=True)
    )
    assert (r.placed, r.remaining, r.termination) == (3, INF, T.EXHAUSTED)
    assert r.edge_flows.dtype == np.int64 and r.residual.dtype == np.int64
    assert list(r.residual) == list(fg.residual())
    own = [0] * 4
    for rec in pol.records():
        for e, x in rec.edges:
            own[e] += x
    assert list(r.edge_flows) == own
    assert r.bottleneck == fg.bottleneck(0)
    assert r.min_cut is None  # Global only
    assert pol.placed_demand() == 7
    with pytest.raises(ValueError, match="integer units or"):
        pol.place(1.5)
    with pytest.raises(ValueError, match=r"after \+inf"):
        pol.place(np.int64(1))
    pol.remove()
    assert pol.records() == [] and list(fg.allocated()) == [0, 0, 0, 0]
    assert pol.place(np.int64(2)).placed == 2
    # Global: the marginal curve; min_cut with bottleneck.
    fg2 = ngc.FlowLedger(diamond(5, 3, 2))
    gl = ngc.FlowPolicy(fg2, 0, 3, GLOBAL)
    r = gl.place(INF, ngc.OutputOptions(bottleneck=True))
    assert r.cost_distribution == [(2, 5), (4, 3)]
    assert sorted(r.min_cut) == ids(diamond(5, 3, 2), [0, 2])
    assert [rec.delivered for rec in gl.records()] == [8]


def test_lifetimes_keep_alive_and_orphans():
    g = diamond()
    fg = ngc.FlowLedger(g)
    el = g.eligibility(excluded_edges=ids(g, [2]))
    cache = ngc.DagCache(g)
    pol = ngc.FlowPolicy(fg, 0, 3, config(ngc.CostOnly()), el, cache)
    del fg, el, cache
    gc.collect()
    assert (
        pol.place(10).placed == 5
    )  # the ledger, the Eligibility and the cache are alive
    fg2 = ngc.FlowLedger(g)
    p = ngc.FlowPolicy(fg2, 0, 3, config(ngc.CostOnly()))
    p.place(2)
    del p
    gc.collect()
    if HOOKS:
        assert fg2._state()["orphans"] == [0]
    assert fg2.records()[0].delivered == 2
    with pytest.raises(ValueError, match="bound"):
        q = ngc.FlowPolicy(fg2, 0, 3, config(ngc.CostOnly()))
        fg2.reset()
    del q
    gc.collect()
    fg2.reset()
    assert fg2.records() == []


def test_pinned_usable_bundles_and_records():
    g = diamond()
    fg = ngc.FlowLedger(g)
    r1, r2 = dag(g, ids(g, [0, 1])), dag(g, ids(g, [2, 3]))
    el = g.eligibility(excluded_edges=ids(g, [3]))
    pol = ngc.FlowPolicy(fg, 0, 3, config(ngc.Pinned([r1, r2])), el)
    assert pol.usable_bundles() == 1
    r = pol.place(9)
    assert (r.placed, r.termination) == (5, T.EXHAUSTED)
    rec = pol.records()[0]
    assert rec.bundle == r1 and pol.record(rec.flow_id).delivered == 5
    # record() reads the policy's own records only (section 8: a policy's records).
    other = ngc.FlowPolicy(fg, 0, 3, config(ngc.Pinned([r1])))
    other.place(1)
    with pytest.raises(ValueError, match="not a record of this policy"):
        pol.record(other.records()[0].flow_id)
    assert fg.record(other.records()[0].flow_id).delivered == 0
    with pytest.raises(ValueError, match="another ledger"):
        pol.record((ngc.FlowLedger(g).uid(), 0))
    pol.remove()
    with pytest.raises(ValueError, match="retired or unknown"):
        pol.record(rec.flow_id)
    all_down = ngc.FlowPolicy(
        fg, 0, 3, config(ngc.Pinned([r2])), g.eligibility(excluded_edges=ids(g, [3]))
    )
    assert all_down.usable_bundles() == 0
    r = all_down.place(1)
    assert (r.placed, r.termination, r.local_outcomes) == (0, T.NO_ROUTE, [])


def test_dag_cache_statistics_and_projection():
    g = diamond()
    cache = ngc.DagCache(g)
    fg = ngc.FlowLedger(g)
    a = ngc.FlowPolicy(fg, 0, 3, config(ngc.CostOnly()), None, cache)
    b = ngc.FlowPolicy(
        fg,
        0,
        3,
        config(
            ngc.CostOnly(),
            selection=sel(
                ngc.EdgeTieBreak.PREFER_HIGHER_PRIORITY,
                ngc.EdgePriority.RESIDUAL,
                False,
                False,
            ),
        ),
        None,
        cache,
    )
    s = cache.stats()
    assert (s.hits, s.misses, s.entries, s.live_values) == (1, 1, 1, 1)
    assert s.bytes > 0
    assert a.place(2).placed == 2 and b.place(2).placed == 2
    cache.clear()
    s = cache.stats()
    assert (s.hits, s.misses, s.entries, s.bytes) == (1, 1, 0, 0)
    assert s.live_values == 1  # held by the policies
    assert cache.graph_uid() == g.uid()


def test_capacity_batch_and_workers(monkeypatch):
    g = diamond(5, 3)
    opts = ngc.CapacityOptions(GLOBAL, outputs=ngc.OutputOptions(edge_flows=True))
    r = ngc.capacity(g, 0, 3, opts)
    assert (r.placed, r.termination, r.cost_distribution) == (8, T.EXHAUSTED, [(2, 8)])
    assert int(np.sum(r.edge_flows)) == 16
    common = g.eligibility(excluded_edges=ids(g, [3]))
    pairs = [(0, 3), ([0], [3]), (1, 3), (0, [1, 2])]
    per = [None, g.eligibility(excluded_edges=ids(g, [0])), None, None]
    want = [5, 0, 5, 8]
    for workers in (1, 3, None):
        got = ngc.batch_capacity(
            g, pairs, ngc.CapacityOptions(GLOBAL, elig=common), per, workers
        )
        assert [x.placed for x in got] == want
    monkeypatch.setenv("NGRAPH_CORE_BATCH_THREADS", "4")
    got = ngc.batch_capacity(g, pairs, ngc.CapacityOptions(GLOBAL, elig=common), per)
    assert [x.placed for x in got] == want
    assert ngc.batch_capacity(g, [], ngc.CapacityOptions(GLOBAL)) == []
    with pytest.raises(ValueError, match="node 7"):
        ngc.batch_capacity(
            g, [(0, 3), (0, 7), (0, 8)], ngc.CapacityOptions(GLOBAL), workers=4
        )
    with pytest.raises(ValueError, match="per_pair_elig has 1 entries"):
        ngc.batch_capacity(g, pairs, ngc.CapacityOptions(GLOBAL), [None])
    with pytest.raises(ValueError, match="Pinned"):
        ngc.batch_capacity(
            g, pairs, ngc.CapacityOptions(config(ngc.Pinned([dag(g, ids(g, [0, 1]))])))
        )
    cache = ngc.DagCache(g)
    eb = ngc.CapacityOptions(config(ngc.CostOnly()), cache=cache)
    with pytest.raises(ValueError, match="workers == 1"):
        ngc.batch_capacity(g, pairs, eb, workers=2)
    got = ngc.batch_capacity(
        g, pairs, eb
    )  # None with a cache: 1 worker, even with the env set
    assert [x.placed for x in got] == [7, 7, 5, 7]  # EB: [4, 3] over links of 5 and 3
    with pytest.raises(ValueError, match="workers"):
        ngc.batch_capacity(g, pairs, ngc.CapacityOptions(GLOBAL), workers=0)


def test_sensitivity_candidates_losses_and_unevaluable():
    # Global: an unsaturated critical edge appears (rule 1).
    g = build(3, [(0, 1, 9, 1), (1, 2, 4, 1), (0, 2, 2, 5)])
    s = ngc.sensitivity(g, 0, 2, ngc.SensitivityOptions(GLOBAL))
    assert s.baseline.placed == 6
    pos = g.input_index()
    assert {int(pos[e.edge]): e.loss for e in s.entries} == {0: 4, 1: 4, 2: 2}
    assert all(
        e.status == ngc.SensitivityStatus.EVALUATED and e.reason is None
        for e in s.entries
    )
    # A policy model: signed losses; zero-capacity edges absent.
    p = build(2, [(0, 1, 1, 1), (0, 1, 10, 1), (0, 1, 0, 1)])
    opts = ngc.SensitivityOptions(config(ngc.CostOnly()))
    s = ngc.sensitivity(p, 0, 1, opts, workers=2)
    assert s.baseline.placed == 2
    pos = p.input_index()
    assert {int(pos[e.edge]): e.loss for e in s.entries} == {0: -8, 1: 1}
    with pytest.raises(ValueError, match="both endpoint sets"):
        ngc.sensitivity(p, 0, 0, opts)


@needs_hooks
def test_sensitivity_reports_a_model_exception_as_unevaluable():
    # No valid configuration raises under a deletion; the private hook injects one (rule 1).
    p = build(2, [(0, 1, 1, 1), (0, 1, 10, 1), (0, 1, 0, 1)])
    opts = ngc.SensitivityOptions(config(ngc.CostOnly()))
    ext._sensitivity_fail_edge(ids(p, [1])[0])
    try:
        s = ngc.sensitivity(p, 0, 1, opts)
    finally:
        ext._sensitivity_fail_edge(None)
    bad = [e for e in s.entries if e.status == ngc.SensitivityStatus.UNEVALUABLE]
    assert len(bad) == 1 and bad[0].loss is None and "injected" in bad[0].reason


def _state_of(policy, ledger):
    return (
        policy.placed_demand(),
        len(ledger.records()),
        policy._state()["requested"],
        [int(x) for x in ledger.allocated()],
    )


@needs_hooks
@pytest.mark.parametrize("routing", ["costonly", "global"])
@pytest.mark.parametrize("traced", [False, True])
def test_python_publication_is_all_or_nothing(routing, traced):
    # review-1c F2: every Python object of a placement exists before anything is published, so
    # the k-th Python allocation failing raises with the ledger and the policy unchanged, never
    # crashes, and a retry places once.
    cfg = GLOBAL if routing == "global" else config(ngc.CostOnly())
    g = build(2, [(0, 1, 10, 1)])
    ledger = ngc.FlowLedger(g)
    want = ngc.FlowPolicy(ngc.FlowLedger(g), 0, 1, cfg).place(3).cost_distribution
    reached = 0
    for k in range(400):
        ledger = ngc.FlowLedger(g)
        policy = ngc.FlowPolicy(ledger, 0, 1, cfg)
        fn = policy._place_traced if traced else policy.place
        before = _state_of(policy, ledger)
        hit, _, value, error = ext._fail_python_allocation(fn, (3,), k)
        if not hit:
            assert error is None
            r = value[0] if traced else value
            assert (r.placed, r.cost_distribution) == (3, want)
            break
        reached += 1
        if (
            error is None
        ):  # a failure inside Python that the call absorbed (none expected)
            continue
        assert isinstance(error, (MemoryError, RuntimeError)), repr(error)
        assert _state_of(policy, ledger) == before
        r = fn(3)
        r = r[0] if traced else r
        assert r.placed == 3 and policy.placed_demand() == 3
    else:
        pytest.fail("the call never completed")
    assert reached >= 3  # not vacuous


@needs_hooks
def test_evaluator_results_survive_python_allocation_failures():
    # The evaluators publish nothing; their result objects are made through the checked path
    # too, so a failed Python allocation raises instead of crashing.
    g = diamond(5, 3)
    calls = [
        (ngc.capacity, (g, 0, 3, ngc.CapacityOptions(GLOBAL))),
        (ngc.batch_capacity, (g, [(0, 3), (1, 3)], ngc.CapacityOptions(GLOBAL))),
        (ngc.sensitivity, (g, 0, 3, ngc.SensitivityOptions(GLOBAL))),
    ]
    for fn, args in calls:
        for k in range(400):
            hit, _, _, error = ext._fail_python_allocation(fn, args, k)
            if not hit:
                assert error is None
                break
            assert error is None or isinstance(error, (MemoryError, RuntimeError))
        else:
            pytest.fail(f"{fn.__name__} never completed")


def test_batch_reports_the_lowest_failing_pair():
    # review-1c F3: the arguments first, then each pair completely in pair order, across the
    # Python/C++ boundary.
    g = build(2, [(0, 1, 3, 7)])
    opts = ngc.CapacityOptions(GLOBAL)
    with pytest.raises(ValueError, match="node 0 is in both endpoint sets"):
        ngc.batch_capacity(g, [(0, 0), (2, 1)], opts, workers=2)
    with pytest.raises(ValueError, match="node 2 out of node range"):
        ngc.batch_capacity(g, [(0, 1), (2, 1), (0, 0)], opts)
    other = build(2, [(0, 1, 3, 7)]).eligibility()
    with pytest.raises(ValueError, match=r"per_pair_elig\[0\]"):
        ngc.batch_capacity(g, [(0, 1), (2, 1)], opts, per_pair_elig=[other, None])
    with pytest.raises(TypeError, match="per_pair_elig entries"):
        ngc.batch_capacity(g, [(0, 1), (2, 1)], opts, per_pair_elig=[3, None])
    with pytest.raises(ValueError, match="node 2 out of node range"):
        ngc.batch_capacity(g, [(0, 1), (2, 1)], opts, per_pair_elig=[None, 3])
    # An argument error comes before every pair error.
    with pytest.raises(ValueError, match="workers"):
        ngc.batch_capacity(g, [(0, 0)], opts, workers=0)


def test_batch_checks_an_explicit_empty_eligibility_list():
    # review-1c F4: a given per_pair_elig has one entry per pair, also when it is empty.
    g = build(2, [(0, 1, 3, 7)])
    opts = ngc.CapacityOptions(GLOBAL)
    with pytest.raises(ValueError, match="per_pair_elig has 0 entries for 1 pairs"):
        ngc.batch_capacity(g, [(0, 1)], opts, per_pair_elig=[])
    assert ngc.batch_capacity(g, [], opts, per_pair_elig=[]) == []
    assert [
        r.placed for r in ngc.batch_capacity(g, [(0, 1)], opts, per_pair_elig=[None])
    ] == [3]


def test_global_zero_call_reports_the_cumulative_curve():
    # review-1c F5: 9.1's curve is cumulative (its amounts sum to placed_demand()); a zero call
    # reports it ([] before the first placement and after remove()), changes nothing and claims
    # no cut.
    g = build(2, [(0, 1, 3, 7)])
    ledger = ngc.FlowLedger(g)
    p = ngc.FlowPolicy(ledger, 0, 1, GLOBAL)
    assert p.place(0).cost_distribution == []
    assert p.place(2).cost_distribution == [(7, 2)]
    z = p.place(0, ngc.OutputOptions(edge_flows=True, residual=True, bottleneck=True))
    assert (z.placed, z.cost_distribution, z.min_cut) == (0, [(7, 2)], None)
    assert p.placed_demand() == 2 and [int(x) for x in ledger.allocated()] == [2]
    p.remove()
    assert p.place(0).cost_distribution == []
