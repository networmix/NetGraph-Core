"""The private flow-core entry `_netgraph_core._solve_min_cost_flow`: the payload
schema, the hand-derived witnesses, input validation and randomized certificates."""

from __future__ import annotations

import random

import _netgraph_core as core
import numpy as np
import pytest
from flow_oracle import (
    NEVER_REACHED,
    certify,
    exhaustive,
    independent_curve,
    points_value,
)
from helpers import CHAIN, build, by_input, exact, ids, to_ids

COUNTERS = {
    "dijkstra_runs",
    "dijkstra_pops",
    "settled",
    "arc_scans",
    "dinic_phases",
    "augmentations",
}


def allowed(eligibility):
    """The certificate's masks of an Eligibility (None: everything allowed)."""
    if eligibility is None:
        return None, None
    return eligibility.node_allowed(), eligibility.edge_allowed()


def solve(g, sources, dests, cap, view=None, eligibility=None):
    view = g.edge_capacity_units().tolist() if view is None else view
    return core._solve_min_cost_flow(g, view, sources, dests, eligibility, cap)


def test_payload_schema():
    g = build(4, CHAIN)
    p = solve(g, [0], [3], None)
    assert set(p) == {
        "flow_value",
        "marginal_cost_curve",
        "edge_flows",
        "min_cut",
        "potentials",
        "counters",
    }
    assert type(p["flow_value"]) is int
    assert all(type(c) is int and type(u) is int for c, u in p["marginal_cost_curve"])
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
    assert one["flow_value"] == 1 and one["marginal_cost_curve"] == [[3, 1]]
    assert by_input(g, one["edge_flows"]) == [1, 1, 1, 0, 0]
    assert one["min_cut"] is None
    two = solve(g, [0], [3], 2)
    assert two["marginal_cost_curve"] == [[3, 1], [5, 1]]
    assert by_input(g, two["edge_flows"]) == [1, 0, 1, 1, 1]
    assert sorted(int(g.input_positions()[e]) for e in two["min_cut"]) == [0, 3]
    assert two["potentials"] == [0, 2, 3, 5]
    inf = solve(g, [0], [3], None)
    assert (
        inf["edge_flows"] == two["edge_flows"] and inf["counters"]["dijkstra_runs"] == 3
    )
    for cap, p in ((1, one), (2, two), (None, inf)):
        assert certify(g, g.edge_capacity_units().tolist(), [0], [3], cap, p) == []
    mask = g.eligibility(excluded_edges=[ids(g, [1])[0]])
    masked = solve(g, [0], [3], None, eligibility=mask)
    assert masked["marginal_cost_curve"] == [[4, 2]]
    view = g.edge_capacity_units().tolist()
    assert certify(g, view, [0], [3], None, masked, *allowed(mask)) == []
    # Every cost times a factor: the curve's costs and the potentials scale, the volumes, flows,
    # cut and operation counters do not.
    for factor in (2, 7, 1000003):
        scaled_g = build(4, [(u, v, c, w * factor) for u, v, c, w in CHAIN])
        scaled = solve(scaled_g, [0], [3], 2)
        assert scaled["marginal_cost_curve"] == [[3 * factor, 1], [5 * factor, 1]]
        assert scaled["potentials"] == [0, 2 * factor, 3 * factor, 5 * factor]
        for key in ("flow_value", "edge_flows", "min_cut", "counters"):
            assert scaled[key] == two[key], (factor, key)
        assert (
            certify(
                scaled_g, scaled_g.edge_capacity_units().tolist(), [0], [3], 2, scaled
            )
            == []
        )


def test_zero_cap_and_positive_cap_on_zero_view():
    g = build(3, [(0, 1, 5, 1), (1, 2, 5, 1)])
    z = solve(g, [0], [2], 0)
    assert z == {
        "flow_value": 0,
        "marginal_cost_curve": [],
        "edge_flows": [0, 0],
        "min_cut": None,
        "potentials": [NEVER_REACHED] * 3,
        "counters": dict.fromkeys(COUNTERS, 0),
    }
    assert certify(g, [5, 5], [0], [2], 0, z) == []
    zero_view = to_ids(g, [0, 5])
    p = solve(g, [0], [2], 4, view=zero_view)
    # The closure {0} is separated from 2 only by the edge 0 -> 1 of view 0 (id 0): it is the
    # cut, of value 0.
    assert (
        p["flow_value"] == 0 and p["min_cut"] == [0] and p["potentials"] == [0, -1, -1]
    )
    assert p["counters"]["dijkstra_runs"] == 1
    assert certify(g, zero_view, [0], [2], 4, p) == []


def test_legal_maximum_and_anycast():
    big = (1 << 62) - 1
    g = build(2, [(0, 1, big, 1)])
    p = solve(g, [0], [1], None)
    assert (
        p["flow_value"] == big
        and p["marginal_cost_curve"] == [[1, big]]
        and p["min_cut"] == [0]
    )
    assert certify(g, [big], [0], [1], None, p) == []
    a = build(3, [(0, 1, 1, 1), (0, 2, 1, 2)])
    q = solve(a, [0], [1, 2], None)
    assert q["flow_value"] == 2 and q["marginal_cost_curve"] == [[1, 1], [2, 1]]
    assert q["potentials"] == [0, 2, 2]
    assert certify(a, [1, 1], [0], [1, 2], None, q) == []
    five = solve(a, np.array([0, 0]), np.array([2, 1, 2]), np.int64(5))
    assert five == solve(a, [0], [1, 2], 5)
    assert certify(a, [1, 1], [0], [1, 2], 5, five) == []


def test_fresh_solve_uses_the_released_cheap_edge():
    # Two parallel edges 0 -> 1 (cost 1 with 1 unit, cost 3 with 2 units), then 1 -> 2. A view
    # in which the cheap edge has its unit back (say, released by an earlier flow) is solved
    # from scratch: the unit goes over the cheap edge, at cost 2.
    g = build(3, [(0, 1, 1, 1), (0, 1, 2, 3), (1, 2, 2, 1)])
    view = to_ids(g, [1, 1, 1])
    p = solve(g, [0], [2], 1, view=view)
    assert by_input(g, p["edge_flows"]) == [1, 0, 1]
    assert p["marginal_cost_curve"] == [[2, 1]]
    assert certify(g, view, [0], [2], 1, p) == []


# A type the position never takes (a float or a bool for an integer, None for the graph) is a
# TypeError, a value outside the domain a ValueError.
TYPE_DEFECTS = {
    "sources node ids[0] must be an integer (got float)",
    "cap must be an integer number of units in [0, 2^63 - 2] or None",
    "graph must be a Graph",
}


@pytest.mark.parametrize(
    "kwargs, message",
    [
        ({"sources": []}, "sources: the node set is empty"),
        ({"destinations": []}, "destinations: the node set is empty"),
        ({"sources": [3]}, "sources: node 3 out of node range [0, 3)"),
        (
            {"sources": [0, 2]},
            "node 2 is in both endpoint sets: sources and destinations must not overlap",
        ),
        ({"sources": [0.0]}, "sources node ids[0] must be an integer (got float)"),
        ({"cap": -1}, "cap -1 is outside [0, 2^63 - 2]"),
        ({"cap": 1 << 63}, "cap 9223372036854775808 is outside [0, 2^63 - 2]"),
        ({"cap": (1 << 63) - 1}, "cap 9223372036854775807 is outside [0, 2^63 - 2]"),
        (
            {"cap": 1.0},
            "cap must be an integer number of units in [0, 2^63 - 2] or None",
        ),
        (
            {"cap": True},
            "cap must be an integer number of units in [0, 2^63 - 2] or None",
        ),
        ({"graph": None}, "graph must be a Graph"),
    ],
)
def test_invalid_inputs(kwargs, message):
    g = build(3, [(0, 1, 5, 1), (1, 2, 5, 1)])
    call = {
        "graph": g,
        "view": [5, 5],
        "sources": [0],
        "destinations": [2],
        "eligibility": None,
        "cap": 3,
    }
    call.update(kwargs)
    with pytest.raises((TypeError, ValueError), match=exact(message)) as error:
        core._solve_min_cost_flow(**call)
    assert error.type is (TypeError if message in TYPE_DEFECTS else ValueError)


def test_foreign_eligibility_and_largest_cap():
    g = build(3, [(0, 1, 5, 1), (1, 2, 5, 1)])
    h = build(3, [(0, 1, 5, 1), (1, 2, 5, 1)])
    with pytest.raises(
        ValueError, match=exact("foreign Eligibility: built for another graph")
    ):
        solve(g, [0], [2], 1, eligibility=h.eligibility())
    assert solve(g, [0], [2], (1 << 63) - 2)["flow_value"] == 5


def kinds_of(sources, dests, view, g, cap, eligibility, excluded_nodes, p):
    """The kinds of case a certified solve with a positive flow belongs to (none for no
    flow): several sources or destinations, a view below the capacities, a cap that stopped
    the solve (no cut) or a maximum flow proved by its cut, an Eligibility, an excluded node."""
    if p["flow_value"] == 0:
        return []
    return [
        kind
        for kind, present in (
            ("positive", True),
            ("several sources", len(sources) > 1),
            ("several destinations", len(dests) > 1),
            ("reduced view", view != g.edge_capacity_units().tolist()),
            ("stopped by the cap", cap is not None and p["min_cut"] is None),
            ("maximum with its cut", bool(p["min_cut"])),
            ("masked", eligibility is not None),
            ("node masked", bool(excluded_nodes)),
        )
        if present
    ]


def test_random_graphs_certified():
    rng = random.Random(2026)
    # The kinds of case the draws aim at, each to be certified at least this often (well below
    # of today's counts: a change of the draws that loses a kind fails here).
    floors = {
        "positive": 51,
        "several sources": 14,
        "several destinations": 11,
        "reduced view": 10,
        "stopped by the cap": 6,
        "maximum with its cut": 28,
        "masked": 16,
        "node masked": 1,
    }
    kinds = dict.fromkeys(floors, 0)
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
            for c in g.edge_capacity_units().tolist()
        ]
        cap = rng.choice([None, 0, 1, rng.randint(0, 6), (1 << 63) - 2])
        eligibility = None
        excluded_nodes = []
        if rng.random() < 0.4:
            excluded_nodes = [rng.randrange(n)] if rng.random() < 0.3 else []
            eligibility = g.eligibility(
                excluded_nodes=excluded_nodes,
                excluded_edges=[rng.randrange(len(edges))]
                if rng.random() < 0.6
                else [],
            )
        p = solve(g, sources, dests, cap, view=view, eligibility=eligibility)
        fails = certify(g, view, sources, dests, cap, p, *allowed(eligibility))
        assert fails == [], (edges, sources, dests, view, cap, fails)
        for kind in kinds_of(
            sources, dests, view, g, cap, eligibility, excluded_nodes, p
        ):
            kinds[kind] += 1
    assert all(kinds[k] >= floors[k] for k in floors), kinds


def test_random_larger_graphs_certified_by_the_independent_solver():
    # Past the exhaustive bound (more than 7 edges here), the certificate compares the curve
    # with the independent solver; that solver first matches the exhaustive optimum at every
    # value of a small graph, so the second tier's reference is itself checked.
    # The chain with bypasses at capacities 4: its second batch of 4 units cancels the first
    # batch's flow on 1 -> 2 (the backward residual arc), so the optimum is 4 units at cost 3,
    # then 4 at cost 5.
    wide_chain = build(4, [(u, v, 4, w) for u, v, _c, w in CHAIN])
    args = (
        4,
        wide_chain.edge_sources().tolist(),
        wide_chain.edge_destinations().tolist(),
        wide_chain.edge_costs().tolist(),
        wide_chain.edge_capacity_units().tolist(),
    )
    assert independent_curve(*args, [0], [3]) == [(0, 0), (4, 12), (8, 32)]
    p = solve(wide_chain, [0], [3], None)
    assert p["marginal_cost_curve"] == [[3, 4], [5, 4]]
    assert certify(wide_chain, args[4], [0], [3], None, p) == []
    # Small graphs of 4 to 7 edges of capacity 1 or 2, with one or two sources and destinations
    # and most with some edges excluded (the certificate passes both solvers its sets and its
    # edge mask): the independent solver's curve equals the exhaustive optimum at every value,
    # many of them after two augmentations or more.
    oracle_rng = random.Random(2028)
    # The kinds of case, each to be crossed at least this often (well below today's counts).
    floors = {
        "several augmentations": 45,
        "several sources": 60,
        "several destinations": 50,
        "edges excluded": 50,
    }
    crossed = dict.fromkeys(floors, 0)
    for _ in range(400):
        n = oracle_rng.randint(3, 4)
        small = [
            (
                oracle_rng.randrange(n),
                oracle_rng.randrange(n),
                oracle_rng.randint(1, 2),
                oracle_rng.randint(1, 5),
            )
            for _ in range(oracle_rng.randint(4, 7))
        ]
        args = (
            n,
            [u for u, _v, _c, _w in small],
            [v for _u, v, _c, _w in small],
            [w for _u, _v, _c, w in small],
            [c for _u, _v, c, _w in small],
        )
        nodes = list(range(n))
        oracle_rng.shuffle(nodes)
        ns = oracle_rng.randint(1, 2)
        nd = oracle_rng.randint(1, min(2, n - ns))
        sources, dests = nodes[:ns], nodes[ns : ns + nd]
        mask = None
        if oracle_rng.random() < 0.7:
            mask = [oracle_rng.random() < 0.8 for _ in small]
        case = (small, sources, dests, mask)
        best = exhaustive(*args, sources, dests, mask)
        ref = independent_curve(*args, sources, dests, mask)
        assert ref[-1][0] == max(best), case
        assert all(points_value(ref, x) == c for x, c in best.items()), case
        if max(best) > 0:
            crossed["several augmentations"] += len(ref) > 2
            crossed["several sources"] += len(sources) > 1
            crossed["several destinations"] += len(dests) > 1
            crossed["edges excluded"] += mask is not None and not all(mask)
    assert all(crossed[k] >= floors[k] for k in floors), crossed
    rng = random.Random(2027)
    # The kinds of case the draws aim at, each to be certified at least this often (well below
    # today's counts).
    floors = {
        "positive": 40,
        "several augmentations": 21,
        "cost factor above 1": 8,
        "several sources": 23,
        "several destinations": 22,
        "reduced view": 38,
        "stopped by the cap": 11,
        "maximum with its cut": 30,
        "masked": 21,
        "node masked": 3,
    }
    kinds = dict.fromkeys(floors, 0)
    for _ in range(120):
        n = rng.randint(3, 7)
        factor = rng.choice([1, 1, 3, 1000003])
        edges = [
            (
                rng.randrange(n),
                rng.randrange(n),
                rng.randint(0, 9),
                rng.randint(1, 9) * factor,
            )
            for _ in range(rng.randint(8, 14))
        ]
        g = build(n, edges)
        nodes = list(range(n))
        rng.shuffle(nodes)
        ns = rng.randint(1, 2)
        sources, dests = nodes[:ns], nodes[ns : ns + rng.randint(1, min(2, n - ns))]
        view = [
            c if rng.random() < 0.7 else rng.randint(0, c)
            for c in g.edge_capacity_units().tolist()
        ]
        cap = rng.choice([None, 1, rng.randint(0, 20), (1 << 63) - 2])
        eligibility = None
        excluded_nodes = []
        if rng.random() < 0.4:
            excluded_nodes = [rng.randrange(n)] if rng.random() < 0.3 else []
            eligibility = g.eligibility(
                excluded_nodes=excluded_nodes,
                excluded_edges=[rng.randrange(len(edges))],
            )
        p = solve(g, sources, dests, cap, view=view, eligibility=eligibility)
        fails = certify(g, view, sources, dests, cap, p, *allowed(eligibility))
        assert fails == [], (edges, sources, dests, view, cap, fails)
        for kind in kinds_of(
            sources, dests, view, g, cap, eligibility, excluded_nodes, p
        ):
            kinds[kind] += 1
        several = len(p["marginal_cost_curve"]) > 1
        kinds["several augmentations"] += several
        kinds["cost factor above 1"] += several and factor > 1
    assert all(kinds[k] >= floors[k] for k in floors), kinds


def test_entry_validates_in_the_core_order_before_narrowing():
    # The private entry validates in a fixed order: view length, view entries, sources,
    # destinations, overlap, cap, each checked on the unnarrowed value. The view's own checks,
    # entries before a later argument included, are the shared ones of test_edge_amounts.py.
    g = build(3, [(0, 1, 5, 1), (1, 2, 5, 1)])

    def call(**kw):
        args = {
            "view": [5, 5],
            "sources": [0],
            "destinations": [2],
            "eligibility": None,
            "cap": 3,
        }
        args.update(kw)
        return core._solve_min_cost_flow(g, **args)

    with pytest.raises(
        ValueError, match=exact("view has length 1, expected 2 (one entry per edge)")
    ):
        call(view=[2**70], sources=[2**40])
    with pytest.raises(
        ValueError, match=exact("sources: node 1099511627776 out of node range [0, 3)")
    ):
        call(sources=[2**40], destinations=[-(2**70)])
    with pytest.raises(
        ValueError,
        match=exact(
            "destinations: node -1180591620717411303424 out of node range [0, 3)"
        ),
    ):
        call(destinations=[-(2**70)], cap=-1)
    with pytest.raises(
        ValueError,
        match=exact("destinations: node 1099511627776 out of node range [0, 3)"),
    ):
        call(sources=[0, 2], destinations=[2, 2**40], cap=2**70)
    with pytest.raises(
        ValueError,
        match=exact(
            "node 2 is in both endpoint sets: sources and destinations must not overlap"
        ),
    ):
        call(sources=[0, 2], cap=2**70)
    with pytest.raises(
        ValueError, match=exact("cap 1180591620717411303424 is outside [0, 2^63 - 2]")
    ):
        call(cap=2**70)
    with pytest.raises(
        ValueError, match=exact("cap -1180591620717411303424 is outside [0, 2^63 - 2]")
    ):
        call(cap=-(2**70))


def test_entry_view_rank_is_checked_with_its_length():
    # A 0-D view (the 2-D and nested views are test_edge_amounts.py's), and a 2-D sources array.
    g = build(3, [(0, 1, 5, 1), (1, 2, 5, 1)])
    with pytest.raises(
        TypeError, match=exact("view must be a 1-D array (got 0 dimensions)")
    ):
        core._solve_min_cost_flow(g, np.array(5), [0], [2], None, 3)
    with pytest.raises(
        TypeError,
        match=exact("sources node ids must be a 1-D array (got 2 dimensions)"),
    ):
        core._solve_min_cost_flow(g, [5, 5], np.array([[0]]), [2], None, -1)


def test_certificate_rejects_mutated_payloads():
    # Each clean payload certifies first; changing one field must fail the clause it
    # breaks. These controls include payload shape and the destination potential
    # pi(t*) against the last curve cost, which the C++ certificate also checks.
    def rejected(g, cap, payload, mutate, dests=(1,), allowed=(None, None)):
        view = g.edge_capacity_units().tolist()
        assert certify(g, view, [0], list(dests), cap, payload, *allowed) == []
        mutate(payload)
        return certify(g, view, [0], list(dests), cap, payload, *allowed)

    one = build(2, [(0, 1, 3, 7)])
    # Underfill: placed 1 is certified for cap 1 and an early stop for cap 2.
    under = solve(one, [0], [1], 1)
    assert certify(one, [3], [0], [1], 1, under) == []
    assert "attained value: early stop" in certify(one, [3], [0], [1], 2, under)
    # The same underfill without a cap, in both branches (capacity 3 takes the exhaustive one,
    # capacity 4 the independent solver's): an uncapped solve must be maximal, so the two clauses
    # that prove it fail and nothing else does.
    for capacity, branch in ((3, "exhaustive"), (4, "independent")):
        single = build(2, [(0, 1, capacity, 7)])
        assert certify(
            single, [capacity], [0], [1], None, solve(single, [0], [1], 1)
        ) == [
            "uncapped solve left the destinations connected",
            f"placed 1 != min(cap, F_max) {capacity} ({branch})",
        ]
    assert "curve moment differs from the flow cost" in rejected(
        one,
        2,
        solve(one, [0], [1], 2),
        lambda p: p.update(marginal_cost_curve=[[6, 2]]),
    )
    assert rejected(
        one, None, solve(one, [0], [1], None), lambda p: p.update(min_cut=[])
    ) == ["min_cut [] != closure cut [0]"]
    # The curve-shape clauses, each the only failure: an entry of volume 0 before the real one,
    # and the one tier split in two at the same cost (total, moment and pi(t*) hold). Capacity 3
    # takes the exhaustive branch, capacity 4 the independent solver's.
    for capacity in (3, 4):
        single = build(2, [(0, 1, capacity, 7)])
        assert solve(single, [0], [1], None)["marginal_cost_curve"] == [[7, capacity]]
        assert rejected(
            single,
            None,
            solve(single, [0], [1], None),
            lambda p, c=capacity: p.update(marginal_cost_curve=[[6, 0], [7, c]]),
        ) == ["curve volume not positive"]
        assert rejected(
            single,
            None,
            solve(single, [0], [1], None),
            lambda p, c=capacity: p.update(marginal_cost_curve=[[7, 1], [7, c - 1]]),
        ) == ["curve costs not strictly increasing"]
    par = build(2, [(0, 1, 1, 1), (0, 1, 1, 3)])
    # Same moment as [[1, 1], [3, 1]], the breakpoint at 1 omitted.
    assert "curve differs from the exhaustive optimum at 1" in rejected(
        par,
        None,
        solve(par, [0], [1], None),
        lambda p: p.update(marginal_cost_curve=[[2, 2]]),
    )
    anycast = build(3, [(0, 1, 1, 1), (0, 2, 1, 2)])
    # Destination 1 lowered to 1 becomes pi(t*), so destination 2, which also receives, differs.
    assert "inflow destination 2 potential differs from pi(t*)" in rejected(
        anycast,
        None,
        solve(anycast, [0], [1, 2], None),
        lambda p: p["potentials"].__setitem__(
            1, 1
        ),  # an inflow destination below pi(t*)
        dests=(1, 2),
    )
    assert "closure vertex 1 holds NEVER_REACHED" in rejected(
        one,
        1,
        solve(one, [0], [1], 1),
        lambda p: p.update(potentials=[NEVER_REACHED, NEVER_REACHED]),
    )
    assert rejected(
        one, 0, solve(one, [0], [1], 0), lambda p: p["counters"].update(dijkstra_runs=1)
    ) == ["zero cap: the prescribed zero payload"]
    # The payload shape: one edge flow too many.
    assert rejected(
        one, None, solve(one, [0], [1], None), lambda p: p["edge_flows"].append(0)
    ) == ["shape"]
    # pi(t*) alone: the saturated edge leaves the destination outside the closure, so only the
    # last curve cost (7) pins its potential.
    assert rejected(
        one,
        None,
        solve(one, [0], [1], None),
        lambda p: p["potentials"].__setitem__(1, 8),
    ) == ["pi(t*) differs from the last curve cost"]
    # The flow clauses. Two parallel edges of capacity 1 and cost 1 carry one unit each; both
    # units on EdgeId 0 keep the value, conservation and the cost but exceed its view of 1.
    twin = build(2, [(0, 1, 1, 1), (0, 1, 1, 1)])
    over = solve(twin, [0], [1], None)
    assert over["edge_flows"] == [1, 1]
    assert "feasibility at 0" in rejected(
        twin, None, over, lambda p: p.update(edge_flows=[2, 0])
    )
    # EdgeId 1 masked: the solve uses EdgeId 0; moving its unit to EdgeId 1 puts flow on the
    # ineligible edge at the same cost.
    mask = twin.eligibility(excluded_edges=[1])
    masked = solve(twin, [0], [1], None, eligibility=mask)
    assert masked["edge_flows"] == [1, 0]
    assert "flow on ineligible edge 1" in rejected(
        twin,
        None,
        masked,
        lambda p: p.update(edge_flows=[0, 1]),
        allowed=(mask.node_allowed(), mask.edge_allowed()),
    )
    # The chain 0 -> 1 -> 2 carries 2 units; one unit less out of node 1 leaves it unbalanced.
    chain = build(3, [(0, 1, 2, 1), (1, 2, 2, 1)])
    out_of_1 = 0 if int(chain.edge_sources()[0]) == 1 else 1
    assert "conservation at 1" in rejected(
        chain,
        None,
        solve(chain, [0], [2], None),
        lambda p: p["edge_flows"].__setitem__(out_of_1, 1),
        dests=(2,),
    )
    # The same chain with capacities 5, capped at 1: every vertex is in the closure with
    # potentials (0, 1, 2). Raising pi(1) to 3 makes the reduced costs of 0 -> 1 (1 + 0 - 3)
    # and of the backward arc 2 -> 1 (-1 + 2 - 3) negative; pi(t*) = 2 is still the last curve
    # cost, so these two are the only failures.
    wide = build(3, [(0, 1, 5, 1), (1, 2, 5, 1)])
    capped = solve(wide, [0], [2], 1)
    assert capped["potentials"] == [0, 1, 2]
    assert rejected(
        wide, 1, capped, lambda p: p["potentials"].__setitem__(1, 3), dests=(2,)
    ) == ["negative reduced cost 0->1", "negative reduced cost 2->1"]
    # Past the exhaustive bound (capacities 4): the cut 0 -> 1 lets one unit through, and from
    # node 1 the direct edge to 3 costs 5, the detour through 2 costs 1 + 1. One unit sent the
    # direct way is feasible and conserved with the same cut and a matching pi(t*), but costs 6,
    # not the optimum 3; only the independent solver's curve catches it.
    detour = build(4, [(0, 1, 1, 1), (1, 3, 4, 5), (1, 2, 4, 1), (2, 3, 4, 1)])
    direct = [
        1 if (int(u), int(v)) in ((0, 1), (1, 3)) else 0
        for u, v in zip(detour.edge_sources(), detour.edge_destinations(), strict=True)
    ]

    def costlier(p):
        p.update(edge_flows=direct, marginal_cost_curve=[[6, 1]])
        p["potentials"][3] = 6

    assert rejected(
        detour, None, solve(detour, [0], [3], None), costlier, dests=(3,)
    ) == ["curve differs from the independent solver at 1"]
    # Two parallel edges of capacity 4 at costs 2 and 4: the curve [[3, 8]] keeps the moment 24,
    # the value and pi(t*) = 3, and equals the optimum at its own breakpoints 0 and 8; only the
    # independent solver's breakpoint at 4 (cost 8, not 12) shows it.
    par4 = build(2, [(0, 1, 4, 2), (0, 1, 4, 4)])

    def averaged(p):
        p.update(marginal_cost_curve=[[3, 8]])
        p["potentials"][1] = 3

    assert rejected(par4, None, solve(par4, [0], [1], None), averaged) == [
        "curve differs from the independent solver at 4"
    ]


def test_certificate_rejects_every_other_clause():
    # The clauses (and each condition of a clause with several) that the test above does not
    # reach, one control each, as in the C++ twin. A control expects the exact list where the
    # mutation can break one clause alone, else the clause it targets among those it also breaks.
    def mutated(g, sources, dests, cap, mutate, eligibility=None):
        view = g.edge_capacity_units().tolist()
        masks = allowed(eligibility)
        payload = solve(g, sources, dests, cap, eligibility=eligibility)
        assert certify(g, view, sources, dests, cap, payload, *masks) == []
        mutate(payload)
        return certify(g, view, sources, dests, cap, payload, *masks)

    def flows(g, by_position):
        """A mutation setting the edge flows of the given input positions."""

        def mutate(p):
            for pos, x in by_position.items():
                p["edge_flows"][ids(g, [pos])[0]] = x

        return mutate

    def potential(v, x):
        return lambda p: p["potentials"].__setitem__(v, x)

    # One edge 0 -> 1 of capacity 3 and cost 7. Uncapped, the edge is saturated, so the closure
    # is {0} and the destination's potential (7) is pinned only by the last curve cost; under
    # cap 1 the destination stays in the closure.
    one = build(2, [(0, 1, 3, 7)])
    assert mutated(one, [0], [1], None, lambda p: p["potentials"].append(0)) == [
        "shape"
    ]
    # The zero-cap payload, each of its fields alone (the counters are the test above's control).
    for field, value in (
        ("flow_value", 1),
        ("marginal_cost_curve", [[7, 1]]),
        ("min_cut", []),
        ("edge_flows", [1]),
        ("potentials", [0, NEVER_REACHED]),
    ):
        assert mutated(
            one, [0], [1], 0, lambda p, f=field, x=value: p.update({f: x})
        ) == ["zero cap: the prescribed zero payload"], field
    # A negative edge flow (the test above sends one above the view).
    assert "feasibility at 0" in mutated(one, [0], [1], None, flows(one, {0: -1}))
    # Sources {0, 1}, destination 2, edges 0 -> 1 and 1 -> 2: source 1 sends its unit directly.
    # Moving it to 0 -> 1 makes source 1 absorb a unit. With destinations {1, 2} and source 0
    # the unit stops at destination 1; moving it to 1 -> 2 makes destination 1 emit a unit.
    two = build(3, [(0, 1, 1, 1), (1, 2, 1, 1)])
    assert "source 1" in mutated(two, [0, 1], [2], None, flows(two, {0: 1, 1: 0}))
    assert "destination 1" in mutated(two, [0], [1, 2], None, flows(two, {0: 0, 1: 1}))
    # Excluded endpoints: 0 -> 2 and 1 -> 2 with sources {0, 1} and source 0 excluded (and with
    # it the edge 0 -> 2); a unit on 0 -> 2 is supplied by the excluded source. Mirrored, 0 -> 1
    # and 0 -> 2 with destinations {1, 2} and destination 2 excluded.
    into_2 = build(3, [(0, 2, 1, 1), (1, 2, 1, 1)])
    assert "source 0" in mutated(
        into_2,
        [0, 1],
        [2],
        None,
        flows(into_2, {0: 1}),
        into_2.eligibility(excluded_nodes=[0]),
    )
    out_of_0 = build(3, [(0, 1, 1, 1), (0, 2, 1, 1)])
    assert "destination 2" in mutated(
        out_of_0,
        [0],
        [1, 2],
        None,
        flows(out_of_0, {1: 1}),
        out_of_0.eligibility(excluded_nodes=[2]),
    )
    # The value against the supply alone and the inflow alone: on the chain 0 -> 1 -> 2 carrying
    # 2, one unit less into node 1 lowers the supply only, one unit less out of it the inflow only.
    chain = build(3, [(0, 1, 2, 1), (1, 2, 2, 1)])
    for pos in (0, 1):
        assert "value differs from supply or inflow" in mutated(
            chain, [0], [2], None, flows(chain, {pos: 1})
        ), pos
    # placed outside [0, cap]: above the cap, and negative without one.
    assert "placed outside [0, cap]" in mutated(
        one, [0], [1], 1, lambda p: p.update(flow_value=2)
    )
    assert "placed outside [0, cap]" in mutated(
        one, [0], [1], None, lambda p: p.update(flow_value=-1)
    )
    # The cut: present while cap 1 leaves the destination reachable, and its value (3) against
    # a lowered placed (which also breaks the value and the curve total).
    assert mutated(one, [0], [1], 1, lambda p: p.update(min_cut=[])) == [
        "min_cut on a connected solve"
    ]
    assert "cut value differs from placed" in mutated(
        one, [0], [1], None, lambda p: p.update(flow_value=2)
    )
    # The potentials: node 1's (outside the closure) below the sentinel (also below pi(t*)),
    # the root's (the closure), and no destination potential at all although 3 units arrive
    # (node 1 is reached at zero flow, so the sentinel also breaks that clause).
    assert "potential below the sentinel" in mutated(
        one, [0], [1], None, potential(1, NEVER_REACHED - 1)
    )
    assert mutated(one, [0], [1], None, potential(0, 1)) == ["root 0 potential 1"]
    assert mutated(one, [0], [1], None, potential(1, NEVER_REACHED)) == [
        "NEVER_REACHED differs from reachability at zero flow at node 1",
        "flow delivered but no destination reached",
    ]
    # The sentinel against reachability at zero flow, each direction alone. On 0 -> 1 and
    # 2 -> 1 nothing reaches node 2, which must hold NEVER_REACHED, not 0 or another value. On
    # 0 -> 1 (capacity 3, cost 7) and 1 -> 2, uncapped, the first search settles node 2 through
    # node 1 and the saturated edge then leaves both outside the closure: node 2 must not hold
    # NEVER_REACHED.
    unreached = build(3, [(0, 1, 1, 1), (2, 1, 1, 1)])
    for x in (0, 123456):
        assert mutated(unreached, [0], [1], None, potential(2, x)) == [
            "NEVER_REACHED differs from reachability at zero flow at node 2"
        ], x
    beyond = build(3, [(0, 1, 3, 7), (1, 2, 1, 1)])
    assert mutated(beyond, [0], [1], None, potential(2, NEVER_REACHED)) == [
        "NEVER_REACHED differs from reachability at zero flow at node 2"
    ]
    # The curve total: a free tier of 1 before the real one keeps the moment (21) and pi(t*),
    # and the curve then also differs from the optimum at 1.
    assert mutated(
        one, [0], [1], None, lambda p: p.update(marginal_cost_curve=[[0, 1], [7, 3]])
    ) == [
        "curve total differs from placed",
        "curve differs from the exhaustive optimum at 1",
    ]
