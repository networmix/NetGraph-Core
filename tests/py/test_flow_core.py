"""The private flow-core entry `_netgraph_core._place_global_units` (APPROACH 2.2): the payload
schema, the review witnesses, input validation and randomized certificates."""

from __future__ import annotations

import random

import _netgraph_core as core
import numpy as np
import pytest
from flow_oracle import NEVER_REACHED, certify

import netgraph_core as ngc

COUNTERS = {
    "dijkstra_runs",
    "dijkstra_pops",
    "settled",
    "arc_scans",
    "dinic_phases",
    "augmentations",
}
CHAIN = [(0, 1, 1, 1), (1, 2, 1, 1), (2, 3, 1, 1), (0, 2, 1, 3), (1, 3, 1, 3)]


def build(n, edges, q=1.0):
    return ngc.Graph(
        n,
        [e[0] for e in edges],
        [e[1] for e in edges],
        [e[2] for e in edges],
        [e[3] for e in edges],
        q,
    )


def by_input(g, values):
    out = [0] * len(values)
    for e, x in enumerate(values):
        out[int(g.input_index()[e])] = x
    return out


def to_ids(g, values_in):
    return [values_in[int(i)] for i in g.input_index()]


def solve(g, sources, dests, cap, view=None, elig=None):
    view = g.capacity_units().tolist() if view is None else view
    return core._place_global_units(g, view, sources, dests, elig, cap)


def test_payload_schema():
    g = build(4, CHAIN)
    p = solve(g, [0], [3], None)
    assert set(p) == {
        "placed",
        "cost_distribution",
        "edge_flows",
        "min_cut",
        "potentials",
        "counters",
    }
    assert type(p["placed"]) is int
    assert all(type(c) is int and type(u) is int for c, u in p["cost_distribution"])
    assert isinstance(p["edge_flows"], list) and all(
        type(x) is int for x in p["edge_flows"]
    )
    assert isinstance(p["min_cut"], list) and all(type(x) is int for x in p["min_cut"])
    assert all(type(x) is int for x in p["potentials"])
    assert set(p["counters"]) == COUNTERS and all(
        type(x) is int for x in p["counters"].values()
    )
    assert p == solve(
        g, [0], [3], None
    )  # identical calls, identical payloads (counters included)


def test_chain_plus_bypasses():
    g = build(4, CHAIN)
    one = solve(g, [0], [3], 1)
    assert one["placed"] == 1 and one["cost_distribution"] == [[3, 1]]
    assert by_input(g, one["edge_flows"]) == [1, 1, 1, 0, 0]
    assert one["min_cut"] is None
    two = solve(g, [0], [3], 2)
    assert two["cost_distribution"] == [[3, 1], [5, 1]]
    assert by_input(g, two["edge_flows"]) == [1, 0, 1, 1, 1]
    assert sorted(int(g.input_index()[e]) for e in two["min_cut"]) == [0, 3]
    assert two["potentials"] == [0, 2, 3, 5]
    inf = solve(g, [0], [3], None)
    assert (
        inf["edge_flows"] == two["edge_flows"] and inf["counters"]["dijkstra_runs"] == 3
    )
    for cap, p in ((1, one), (2, two), (None, inf)):
        assert certify(g, g.capacity_units().tolist(), [0], [3], cap, p) == []
    masked = solve(
        g,
        [0],
        [3],
        None,
        elig=g.eligibility(excluded_edges=[to_ids(g, list(range(5))).index(1)]),
    )
    assert masked["cost_distribution"] == [[4, 2]]


def test_zero_cap_and_positive_cap_on_zero_view():
    g = build(3, [(0, 1, 5, 1), (1, 2, 5, 1)])
    z = solve(g, [0], [2], 0)
    assert z == {
        "placed": 0,
        "cost_distribution": [],
        "edge_flows": [0, 0],
        "min_cut": None,
        "potentials": [NEVER_REACHED] * 3,
        "counters": dict.fromkeys(COUNTERS, 0),
    }
    p = solve(g, [0], [2], 4, view=to_ids(g, [0, 5]))
    assert p["placed"] == 0 and p["min_cut"] == [] and p["potentials"] == [0, -1, -1]
    assert p["counters"]["dijkstra_runs"] == 1


def test_legal_maximum_and_anycast():
    big = (1 << 62) - 1
    g = build(2, [(0, 1, big, 1)])
    p = solve(g, [0], [1], None)
    assert (
        p["placed"] == big
        and p["cost_distribution"] == [[1, big]]
        and p["min_cut"] == [0]
    )
    a = build(3, [(0, 1, 1, 1), (0, 2, 1, 2)])
    q = solve(a, [0], [1, 2], None)
    assert q["placed"] == 2 and q["cost_distribution"] == [[1, 1], [2, 1]]
    assert q["potentials"] == [0, 2, 2]
    assert solve(a, np.array([0, 0]), np.array([2, 1, 2]), np.int64(5)) == solve(
        a, [0], [1, 2], 5
    )


def test_release_witness_fresh_solve():
    g = build(3, [(0, 1, 1, 1), (0, 1, 2, 3), (1, 2, 2, 1)])
    p = solve(g, [0], [2], 1, view=to_ids(g, [1, 1, 1]))
    assert by_input(g, p["edge_flows"]) == [1, 0, 1]
    assert p["cost_distribution"] == [[2, 1]]


@pytest.mark.parametrize(
    "kwargs, pattern",
    [
        ({"view": [1]}, "view has length 1, expected 2"),
        ({"view": [1, 6]}, "capacity_units"),
        ({"view": [1.0, 1.0]}, "view"),
        ({"view": [-1, 0]}, "outside"),
        ({"sources": []}, "sources: the node set is empty"),
        ({"destinations": []}, "destinations: the node set is empty"),
        ({"sources": [3]}, "out of node range"),
        ({"sources": [0, 2]}, "both endpoint sets"),
        ({"sources": [0.0]}, "sources"),
        ({"cap": -1}, "cap"),
        ({"cap": 1 << 63}, "cap"),
        ({"cap": (1 << 63) - 1}, "cap"),
        ({"cap": 1.0}, "cap"),
        ({"cap": True}, "cap"),
    ],
)
def test_invalid_inputs(kwargs, pattern):
    g = build(3, [(0, 1, 5, 1), (1, 2, 5, 1)])
    call = {
        "graph": g,
        "view": [5, 5],
        "sources": [0],
        "destinations": [2],
        "elig": None,
        "cap": 3,
    }
    call.update(kwargs)
    with pytest.raises(ValueError, match=pattern):
        core._place_global_units(**call)


def test_foreign_eligibility_and_largest_cap():
    g = build(3, [(0, 1, 5, 1), (1, 2, 5, 1)])
    h = build(3, [(0, 1, 5, 1), (1, 2, 5, 1)])
    with pytest.raises(ValueError, match="foreign Eligibility"):
        solve(g, [0], [2], 1, elig=h.eligibility())
    assert solve(g, [0], [2], (1 << 63) - 2)["placed"] == 5


def test_random_graphs_certified():
    rng = random.Random(2026)
    positive = 0
    for _ in range(400):
        n = rng.randint(2, 5)
        edges = [
            (rng.randrange(n), rng.randrange(n), rng.randint(0, 2), rng.randint(1, 5))
            for _ in range(rng.randint(1, 6))
        ]
        g = build(n, edges)
        nodes = list(range(n))
        rng.shuffle(nodes)
        ns = rng.randint(1, min(2, n - 1))
        sources, dests = nodes[:ns], nodes[ns : ns + rng.randint(1, min(2, n - ns))]
        view = [
            c if rng.random() < 0.7 else rng.randint(0, c)
            for c in g.capacity_units().tolist()
        ]
        cap = rng.choice([None, 0, 1, rng.randint(0, 6), (1 << 63) - 2])
        elig = None
        if rng.random() < 0.4:
            elig = g.eligibility(
                excluded_nodes=[rng.randrange(n)] if rng.random() < 0.3 else [],
                excluded_edges=[rng.randrange(len(edges))]
                if rng.random() < 0.6
                else [],
            )
        p = solve(g, sources, dests, cap, view=view, elig=elig)
        fails = certify(
            g,
            view,
            sources,
            dests,
            cap,
            p,
            None if elig is None else elig.node_allowed(),
            None if elig is None else elig.edge_allowed(),
        )
        assert fails == [], (edges, sources, dests, view, cap, fails)
        positive += p["placed"] > 0
    assert positive > 50


def test_entry_validates_in_the_core_order_before_narrowing():
    # review-1a D1 for the private entry: view length, view entries, sources, destinations,
    # overlap, cap, each checked on the unnarrowed value.
    g = build(3, [(0, 1, 5, 1), (1, 2, 5, 1)])

    def call(**kw):
        args = {
            "view": [5, 5],
            "sources": [0],
            "destinations": [2],
            "elig": None,
            "cap": 3,
        }
        args.update(kw)
        return core._place_global_units(g, **args)

    with pytest.raises(ValueError, match="view has length 1, expected 2"):
        call(view=[2**70], sources=[2**40])
    with pytest.raises(
        ValueError, match=r"view\[0\] = 1180591620717411303424 is outside"
    ):
        call(view=[2**70, 5], sources=[2**40])
    with pytest.raises(
        ValueError, match=r"sources: node 1099511627776 out of node range"
    ):
        call(sources=[2**40], destinations=[-(2**70)])
    with pytest.raises(ValueError, match=r"destinations: node -1180591620717411303424"):
        call(destinations=[-(2**70)], cap=-1)
    with pytest.raises(ValueError, match="both endpoint sets"):
        call(sources=[0, 2], cap=2**70)
    with pytest.raises(ValueError, match="cap 1180591620717411303424 is outside"):
        call(cap=2**70)
    with pytest.raises(ValueError, match="cap -1180591620717411303424 is outside"):
        call(cap=-(2**70))
    with pytest.raises(
        ValueError, match=r"view\[1\] = 18446744073709551615 is outside"
    ):
        call(view=np.array([5, 2**64 - 1], dtype=np.uint64))


def test_entry_view_rank_is_checked_with_its_length():
    g = build(3, [(0, 1, 5, 1), (1, 2, 5, 1)])
    with pytest.raises(
        ValueError, match=r"view must be a 1-D array \(got 2 dimensions\)"
    ):
        core._place_global_units(g, np.array([[5, 5]]), [2**40], [2], None, 3)
    with pytest.raises(
        ValueError, match=r"view must be a 1-D array \(got 0 dimensions\)"
    ):
        core._place_global_units(g, np.array(5), [0], [2], None, 3)
    with pytest.raises(
        ValueError, match=r"view must be a 1-D sequence of integers \(entry 1"
    ):
        core._place_global_units(g, [5, [5]], [0], [2], None, 3)
    with pytest.raises(ValueError, match="sources node ids must be a 1-D array"):
        core._place_global_units(g, [5, 5], np.array([[0]]), [2], None, -1)
