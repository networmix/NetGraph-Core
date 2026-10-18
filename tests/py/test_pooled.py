"""Pooled CostOnly sources (docs/semantics.md, section 7): one split per shared router,
attribution to sources by highest averages, the equal origin's entry that cannot pass anything,
Lossy forwarding through a source, weights read once per epoch, the first-hop experiment,
increments and epochs. The values are derived by hand (dev/design step-2 LEDGER rows P1 to P17
and I1 to I4); the random tests compare every field with the independent reference of
``origin_oracle.py`` and check that each source's record uses only its own route."""

from __future__ import annotations

import random

import pytest
from helpers import build, dag_of, exact, sel
from origin_oracle import INF, PooledReference

import netgraph_core as ngc

T = ngc.Termination
EB = ngc.FlowPlacement.EQUAL_BALANCED
W = ngc.FlowPlacement.WEIGHTED
LOSSY = ngc.Admission.LOSSY
LOSSLESS = ngc.Admission.LOSSLESS
EDGES = ngc.OutputOptions(edge_flows=True)

# The graphs of the ledger, every edge written in EdgeId order: (u, v, capacity, cost).
POOL = [(0, 2, 10, 1), (1, 2, 10, 1), (2, 3, 1, 1), (2, 3, 1, 1)]  # {0, 1} -> 3
B1 = [(0, 3, 2, 1), (1, 2, 10, 1), (2, 3, 10, 2)]  # {0, 1} -> 3
B3Z = [(0, 2, 0, 1), (1, 2, 10, 1)]  # {0, 1} -> 2
B3P = [(0, 2, 10, 1), (1, 2, 0, 1)]  # {0, 1} -> 2
CHAIN = [(0, 2, 3, 1), (1, 0, 10, 1)]  # {0, 1} -> 2: source 1's route passes source 0
SNAP = [(0, 2, 10, 1), (1, 2, 10, 1), (2, 3, 10, 1), (2, 3, 10, 1)]  # {0, 1} -> 3
PRIV = [(0, 2, 4, 1), (1, 2, 6, 1)]  # {0, 1} -> 2
INC = [(0, 2, 1, 1), (1, 2, 3, 1)]  # {0, 1} -> 2
ONE = [(0, 1, 10, 1)]
RC = [(0, 1, 1, 1), (0, 1, 1, 1)]


def cfg(placement=EB, admission=LOSSLESS, selection=None):
    return ngc.FlowPolicyConfig(
        routing=ngc.CostOnly(),
        placement=placement,
        admission=admission,
        selection=selection,
    )


def policy(edges, sources, dest, config, n=None, residual=None):
    g = build(n or 1 + max(max(u, v) for u, v, *_ in edges), edges)
    ledger = ngc.FlowLedger(g, residual)
    return g, ledger, ngc.FlowPolicy(ledger, sources, dest, config)


def rows(r):
    """by_source as (source, offered, admitted, delivered, source_loss, dropped) tuples."""
    return [
        (x.source, x.offered, x.admitted, x.delivered, x.source_loss, x.dropped)
        for x in r.by_source
    ]


def records(p):
    """Per record: (its one source, its edges, its delivery)."""
    return [(rec.sources, rec.edges, rec.delivered) for rec in p.records()]


def test_b_pooled():
    # P1: RR(2, 2) = (1, 1) injected; router 2 holds 2 and splits (1, 1) over e2, e3; e2's unit
    # goes to source 0 by H(1, [1, 1]) = [1, 0], e3's to source 1. Split separately, both units
    # would take e2 (capacity 1) and only 1 would fit.
    _, _, p = policy(POOL, [0, 1], 3, cfg())
    r = p.place(2, EDGES)
    assert (r.placed, r.termination, r.edge_flows.tolist()) == (
        2,
        T.COMPLETE,
        [1, 1, 1, 1],
    )
    assert rows(r) == [(0, None, 1, 1, 0, 0), (1, None, 1, 1, 0, 0)]
    assert records(p) == [([0], [(0, 1), (2, 1)], 1), ([1], [(1, 1), (3, 1)], 1)]
    # P2: 3 would give e2 two units: the unlimited request admits 2.
    _, _, p = policy(POOL, [0, 1], 3, cfg())
    r = p.place(INF)
    assert (r.placed, r.remaining, r.termination) == (2, INF, T.EXHAUSTED)
    assert (r.requested, r.unadmitted, r.offered) == (INF, None, None)
    # P3: capacity weights (1, 1) are equal: highest averages is the round robin.
    _, _, p = policy(POOL, [0, 1], 3, cfg(W))
    assert p.place(2).placed == 2


def test_lossy_attribution_by_arrivals():
    # P4: injections (2, 2); router 2 holds 4: (2, 2) over e2, e3, each attributed
    # H(2, [2, 2]) = [1, 1] (then [1, 1] remain for e3); each edge carries 1 of its 2 arrivals,
    # H(1, [1, 1]) = [1, 0]: source 0 delivers 2, source 1 drops one unit on each.
    _, _, p = policy(POOL, [0, 1], 3, cfg(admission=LOSSY))
    r = p.place(4)
    assert (r.placed, r.termination, r.call_drops) == (2, T.COMPLETE, [(2, 1), (3, 1)])
    assert rows(r) == [(0, 2, 2, 2, 0, 0), (1, 2, 2, 0, 0, 2)]
    assert (r.offered, r.admitted, r.unadmitted, r.drops) == (4, 4, 0, [(2, 1), (3, 1)])


@pytest.mark.parametrize(
    ("edges", "lossless"),
    [(B3Z, 0), (B3P, 1)],
    ids=["B3-zero-first", "B3-positive-first"],
)
def test_b3(edges, lossless):
    # P5, P6: Lossless admits a total only when it gives the source behind the zero edge
    # nothing: RR(1) gives source 0 the unit, RR(2) gives source 1 one.
    _, _, p = policy(edges, [0, 1], 2, cfg())
    r = p.place(10)
    assert (r.placed, r.termination) == (lossless, T.EXHAUSTED)
    assert [rec.delivered for rec in p.records()] == [lossless, 0]
    # P9: the same with weights: a share without a usable route is not admitted.
    _, _, p = policy(edges, [0, 1], 2, cfg(W))
    assert p.place(10).placed == lossless


def test_lossy_entry_that_cannot_pass():
    # P7: EB keeps the zero edge and drops source 0's 5 there.
    _, _, p = policy(B3Z, [0, 1], 2, cfg(admission=LOSSY))
    r = p.place(10)
    assert (r.placed, r.termination, r.call_drops, r.call_source_loss) == (
        5,
        T.COMPLETE,
        [(0, 5)],
        0,
    )
    # P8: W gives the zero edge weight 0: source 0 has no member, its 5 is lost at the source.
    _, _, p = policy(B3Z, [0, 1], 2, cfg(W, LOSSY))
    r = p.place(10)
    assert (r.placed, r.termination, r.call_drops, r.call_source_loss) == (
        5,
        T.COMPLETE,
        [],
        5,
    )
    assert rows(r) == [(0, 5, 5, 0, 5, 0), (1, 5, 5, 5, 0, 0)]
    assert (r.offered, r.source_loss, r.delivered) == (10, 5, 5)


def test_equal_origin_limits_the_total():
    # P10: RR(5, 2) = (3, 2) would put 3 on e0 (capacity 2): 4, rows (2, 2).
    _, _, p = policy(B1, [0, 1], 3, cfg())
    r = p.place(12)
    assert (r.placed, r.termination) == (4, T.EXHAUSTED)
    assert [x.delivered for x in r.by_source] == [2, 2]
    # P11: Lossy offers (6, 6): source 0 drops 4 on e0.
    _, _, p = policy(B1, [0, 1], 3, cfg(admission=LOSSY))
    r = p.place(12)
    assert (r.placed, r.call_drops) == (8, [(0, 4)])


def test_transit_through_a_source():
    # P12: injections (5, 5); e1 carries source 1's 5 to node 0, which holds 10 for e0 (3):
    # carried H(3, [5, 5]) = [2, 1].
    _, _, p = policy(CHAIN, [0, 1], 2, cfg(admission=LOSSY))
    r = p.place(10)
    assert (r.placed, r.call_drops) == (3, [(0, 7)])
    assert rows(r) == [(0, 5, 5, 2, 0, 3), (1, 5, 5, 1, 0, 4)]
    assert records(p) == [([0], [(0, 2)], 2), ([1], [(0, 1), (1, 5)], 1)]
    # P13: Lossless: e0 carries the whole total, so 3, as RR(3) = (2, 1).
    _, _, p = policy(CHAIN, [0, 1], 2, cfg())
    r = p.place(INF)
    assert (r.placed, [x.delivered for x in r.by_source]) == (3, [2, 1])
    # P17: the first-hop experiment: source 0 injects e0's 3, source 1 e1's 10; node 0 adds
    # source 1's 10 of transit on e0: arrivals (3, 10), carried H(3, [3, 10]) = [0, 3].
    _, _, p = policy(CHAIN, [0, 1], 2, cfg(admission=LOSSY))
    r = p.place(INF)
    assert (r.placed, r.termination, r.requested, r.offered) == (
        3,
        T.EXHAUSTED,
        INF,
        13,
    )
    assert r.drops == [(0, 10)]
    assert rows(r) == [(0, 3, 3, 0, 0, 3), (1, 10, 10, 3, 0, 7)]


def test_residual_weights_on_the_common_view():
    # P14: weights read once, on the residual at the first positive call's entry: e2 2, e3 10;
    # router 2: H(12, [2, 10]) = (2, 10); RR(12) = (6, 6); each member attributed by highest
    # averages over what the sources still hold: e2 H(2, [6, 6]) = [1, 1], e3 H(10, [5, 5]).
    pick = sel(prio=ngc.EdgePriority.RESIDUAL)
    _, _, p = policy(SNAP, [0, 1], 3, cfg(W, selection=pick), residual=[10, 10, 2, 10])
    r = p.place(INF)
    assert (r.placed, [x.delivered for x in r.by_source]) == (12, [6, 6])
    assert records(p) == [
        ([0], [(0, 6), (2, 1), (3, 5)], 6),
        ([1], [(1, 6), (2, 1), (3, 5)], 6),
    ]


def test_weights_frozen_for_the_epoch():
    # P15: another owner holds 8 of e2 at the first call: weights (2, 10) for the epoch, also
    # after the owner leaves (weights read again would give [6, 6, 6, 6]).
    g = build(4, SNAP)
    ledger = ngc.FlowLedger(g)
    pinned = ngc.FlowPolicyConfig(routing=ngc.Pinned([dag_of(g, [2])]), placement=EB)
    other = ngc.FlowPolicy(ledger, 2, 3, pinned)
    assert other.place(8).placed == 8
    pick = sel(prio=ngc.EdgePriority.RESIDUAL)
    p = ngc.FlowPolicy(ledger, [0, 1], 3, cfg(W, selection=pick))
    assert p.place(2).placed == 2
    other.remove()
    r = p.place(10, EDGES)
    assert (r.placed, r.edge_flows.tolist()) == (10, [6, 6, 2, 10])


def test_first_hop_experiment():
    # P16: every source injects its first edges' residual: 4 + 6.
    _, _, p = policy(PRIV, [0, 1], 2, cfg(admission=LOSSY))
    r = p.place(INF)
    assert (r.placed, r.termination, r.requested, r.offered) == (
        10,
        T.EXHAUSTED,
        INF,
        10,
    )
    for finite in (1, 0):
        with pytest.raises(
            ValueError,
            match=exact(
                "finite placement after infinite placement requires a new epoch"
            ),
        ):
            p.place(finite)


def test_inc():
    # I1: RR(3) = (2, 1) does not fit e0 (1): call 1 admits 2; call 2 targets 3 again: 0.
    _, _, p = policy(INC, [0, 1], 2, cfg())
    a, b = p.place(3), p.place(1)
    assert (a.placed, a.termination, b.placed, b.termination) == (
        2,
        T.EXHAUSTED,
        0,
        T.EXHAUSTED,
    )
    assert (b.requested, b.admitted, b.unadmitted, b.offered) == (4, 2, 2, None)
    # I2: Lossy: call 1 injects (2, 1) and drops 1 on e0; call 2 injects RR(4) - RR(3) = (0, 1).
    _, _, p = policy(INC, [0, 1], 2, cfg(admission=LOSSY))
    a, b = p.place(3), p.place(1)
    assert (a.placed, a.call_drops, b.placed, b.call_drops) == (2, [(0, 1)], 1, [])
    assert (b.requested, b.offered, b.delivered, b.drops, b.unadmitted) == (
        4,
        4,
        3,
        [(0, 1)],
        0,
    )
    assert rows(b) == [(0, 2, 2, 1, 0, 1), (1, 2, 2, 2, 0, 0)]
    # I3: another owner holds 8 of 10: call 1 admits 2; it leaves and places 5; call 2 targets
    # 3 on the free 3 plus the policy's own 2.
    g = build(2, ONE)
    ledger = ngc.FlowLedger(g)
    other = ngc.FlowPolicy(ledger, 0, 1, cfg())
    other.place(8)
    p = ngc.FlowPolicy(ledger, 0, 1, cfg())
    assert p.place(5).placed == 2
    other.remove()
    other.place(5)
    r = p.place(1)
    assert (r.placed, r.requested, r.admitted, r.unadmitted) == (1, 6, 3, 3)
    # I4: two parallel unit edges: Lossy 1 then 1 drops the second unit on e0 (RR(1) again);
    # one call of 2 carries both; Lossless 1 then 1 fills both.
    _, _, p = policy(RC, 0, 1, cfg(admission=LOSSY))
    assert [p.place(1).placed, p.place(1).placed] == [1, 0]
    _, _, p = policy(RC, 0, 1, cfg(admission=LOSSY))
    assert p.place(2).placed == 2
    _, _, p = policy(RC, 0, 1, cfg())
    assert [p.place(1).placed, p.place(1, EDGES).edge_flows.tolist()] == [1, [1, 1]]


def random_case(rng):
    """A small graph (cycles and parallel edges allowed), usually several sources, and edges
    biased toward the destinations so that most sources have a route."""
    n = rng.randint(3, 7)
    nodes = list(range(n))
    rng.shuffle(nodes)
    k = rng.randint(2, min(3, n - 1)) if rng.random() < 0.8 else 1
    sources = nodes[:k]
    dests = nodes[k : k + rng.randint(1, min(2, n - k))]
    edges = []
    for _ in range(rng.randint(n, 2 * n + 3)):
        u = rng.randrange(n)
        v = rng.choice(dests) if rng.random() < 0.3 else rng.randrange(n)
        if u != v:
            edges.append((u, v, rng.choice([0, 1, 2, 3, 5, 8]), rng.randint(1, 3)))
    if not edges:
        edges.append((sources[0], dests[0], 1, 1))
    return n, edges, sources, dests


def compare(where, r, want, p):
    got = {
        "placed": r.placed,
        "termination": r.termination.name,
        "call_drops": r.call_drops,
        "call_source_loss": r.call_source_loss,
        "requested": r.requested,
        "offered": r.offered,
        "admitted": r.admitted,
        "delivered": r.delivered,
        "unadmitted": r.unadmitted,
        "source_loss": r.source_loss,
        "drops": r.drops,
        "by_source": rows(r),
        "records": {rec.sources[0]: (rec.edges, rec.delivered) for rec in p.records()},
    }
    for key, value in want.items():
        assert got[key] == value, (where, key, got[key], value)


@pytest.mark.parametrize("seed", range(400))
def test_random_against_reference(seed):
    # P19: random graphs with masks, zero-capacity edges, nonuniform costs and a partial
    # residual; EB and W (capacity or residual weights), both admissions, multipath or the
    # canonical path; a series of calls. Every field equals the reference, the accounting
    # identities hold, and (P18) each record's edges lie in its own route.
    rng = random.Random(seed)
    n, edges, sources, dests = random_case(rng)
    g = build(n, edges)
    m = g.num_edges()
    cap = [int(x) for x in g.capacity_units()]
    residual = [rng.randint(0, c) if rng.random() < 0.3 else c for c in cap]
    xn = [
        v
        for v in range(n)
        if v not in sources and v not in dests and rng.random() < 0.15
    ]
    xe = [e for e in range(m) if rng.random() < 0.1]
    kernel = rng.choice(["EB", "W"])
    lossy = rng.random() < 0.5
    multipath = rng.random() < 0.7
    priority = rng.choice(["CAPACITY", "RESIDUAL"]) if kernel == "W" else "CAPACITY"
    selection = sel(
        prio=getattr(ngc.EdgePriority, priority),
        mp=multipath,
    )
    config = cfg(EB if kernel == "EB" else W, LOSSY if lossy else LOSSLESS, selection)
    ledger = ngc.FlowLedger(g, residual)
    elig = g.eligibility(excluded_nodes=xn, excluded_edges=xe)
    p = ngc.FlowPolicy(ledger, sources, dests, config, elig)
    ref = PooledReference(
        g, sources, dests, kernel, lossy, residual, multipath, priority, xn, xe
    )
    calls = [rng.choice([0, 1, 2, 3, 5, 9]) for _ in range(rng.randint(1, 4))]
    if rng.random() < 0.3:
        calls.append(INF)
    for i, volume in enumerate(calls):
        r = p.place(volume)
        want = ref.place(volume)
        compare((seed, i, volume), r, want, p)
        if lossy:
            assert r.offered == r.source_loss + sum(x for _, x in r.drops) + r.delivered
        elif volume != INF and not ref.infinite:
            assert r.requested == r.unadmitted + r.delivered
        for field in ("admitted", "delivered", "source_loss", "dropped"):
            total = sum(getattr(x, field) for x in r.by_source)
            state = r.drops if field == "dropped" else getattr(r, field)
            if field == "dropped":
                state = sum(x for _, x in state)
            if ref.roster or not lossy:
                assert total == state, (seed, field)
        for rec in p.records():
            route = {int(e) for e in rec.bundle.via_edges()}
            assert {e for e, _ in rec.edges} <= route, (seed, rec.sources)
            assert {e for e, _ in rec.drops} <= route, (seed, rec.sources)
        allocated = ledger.allocated().tolist()
        assert all(0 <= a <= c for a, c in zip(allocated, residual, strict=True)), seed
