"""Graph on integer units (docs/semantics.md#2-the-graph): construction, validation order, EdgeIds,
read-only views, Eligibility."""

from __future__ import annotations

import gc
import weakref

import numpy as np
import pytest
from helpers import L62, build, by_input, exact

import netgraph_core as ngc

# (u, v, capacity in units, cost) in input order.
EDGES = [
    (2, 0, 5, 3),
    (0, 2, 4, 2),
    (0, 1, 3, 2),
    (1, 2, 2, 1),
    (0, 1, 1, 2),
    (0, 2, 6, 1),
    (2, 1, 7, 1),
]


def test_edge_id_order_and_input_positions():
    g = build(3, EDGES)
    order = sorted(
        range(len(EDGES)), key=lambda i: (EDGES[i][0], EDGES[i][3], EDGES[i][1], i)
    )
    assert g.input_positions().tolist() == order == [5, 2, 4, 1, 3, 6, 0]
    assert g.outgoing_offsets().tolist() == [0, 4, 5, 7]
    assert g.edge_sources().tolist() == [EDGES[i][0] for i in order]
    assert g.edge_destinations().tolist() == [EDGES[i][1] for i in order]
    assert g.edge_capacity_units().tolist() == [EDGES[i][2] for i in order]
    assert g.edge_costs().tolist() == [EDGES[i][3] for i in order]
    inv = np.empty(len(EDGES), dtype=np.int64)
    inv[g.input_positions()] = np.arange(len(EDGES))
    assert g.input_positions()[inv].tolist() == list(range(len(EDGES)))
    # Reverse CSR in (cost, source, input) order, as input positions.
    rev = []
    off, adj = g.incoming_offsets(), g.incoming_edges()
    for v in range(3):
        rev.append([int(g.input_positions()[e]) for e in adj[off[v] : off[v + 1]]])
    assert rev == [[0], [6, 2, 4], [5, 3, 1]]
    assert g.max_edge_cost() == 3
    assert g.edge_cost_gcd() == 1
    assert g.num_nodes() == 3 and g.num_edges() == 7
    assert g.resolution() == 1.0


def test_numpy_integer_inputs_and_exact_large_units():
    big = (1 << 53) + 1
    g = ngc.Graph(
        np.int64(2),
        np.array([0, 0], dtype=np.int32),
        np.array([1, 1], dtype=np.uint8),
        np.array([big, 7], dtype=np.int64),
        np.array([2, 1], dtype=np.uint64),
        2.0**-12,
    )
    assert g.edge_capacity_units().tolist() == [7, big]
    assert g.resolution() == 2.0**-12


def test_empty_and_isolated():
    g = ngc.Graph(0, [], [], [], [], 1)
    assert g.outgoing_offsets().tolist() == [0] and g.incoming_offsets().tolist() == [0]
    assert g.max_edge_cost() == 0 and g.edge_cost_gcd() == 0
    h = ngc.Graph(3, [], [], [], [], 1)
    assert h.outgoing_offsets().tolist() == [0, 0, 0, 0]


def test_uids_are_distinct():
    a, b = build(3, EDGES), build(3, EDGES)
    assert a.uid() != b.uid() and a.uid() > 0


def test_resolution_has_no_default():
    with pytest.raises(TypeError):
        ngc.Graph(2, [0], [1], [1], [1])


@pytest.mark.parametrize(
    "args, message",
    [
        # The resolution.
        (
            (2, [0], [1], [1], [1], 3.0),
            "resolution must be a power of two in [2^-32, 2^32]",
        ),
        (
            (2, [0], [1], [1], [1], 2.0**-33),
            "resolution must be a power of two in [2^-32, 2^32]",
        ),
        (
            (2, [0], [1], [1], [1], 0.0),
            "resolution must be a power of two in [2^-32, 2^32]",
        ),
        (
            (2, [0], [1], [1], [1], float("nan")),
            "resolution must be a power of two in [2^-32, 2^32]",
        ),
        # The unit values.
        (
            (2, [0], [1], [-1], [1], 1),
            "edge_capacity_units[0] = -1 is outside the units domain [0, 2^62)",
        ),
        (
            (2, [0], [1], [L62], [1], 1),
            "edge_capacity_units[0] = 4611686018427387904 is outside the units domain [0, 2^62)",
        ),
        (
            (2, [0], [1], [1.0], [1], 1),
            "edge_capacity_units[0] must be an integer (got float)",
        ),
        (
            (2, [0], [1], np.array([1.5]), [1], 1),
            "edge_capacity_units entries must be integers (got dtype float64)",
        ),
        (
            (2, [0], [1], [True], [1], 1),
            "edge_capacity_units[0] must be an integer (got bool)",
        ),
        # The shapes.
        ((-1, [], [], [], [], 1), "num_nodes = -1 is out of range [0, 2147483647]"),
        ((2.0, [0], [1], [1], [1], 1), "num_nodes must be an integer (got float)"),
        (
            (2, [0, 1], [1], [1], [1], 1),
            "edge_sources, edge_destinations, edge_capacity_units and edge_costs must have the same length (got 2, 1, 1, 1)",
        ),
        (
            (2, [0], [1], [1, 1], [1], 1),
            "edge_sources, edge_destinations, edge_capacity_units and edge_costs must have the same length (got 1, 1, 2, 1)",
        ),
        # Zero, negative and non-integer costs.
        (
            (2, [0], [1], [1], [0], 1),
            "edge_costs[0] = 0: costs must be positive (edge_costs >= 1)",
        ),
        (
            (2, [0], [1], [1], [-3], 1),
            "edge_costs[0] = -3: costs must be positive (edge_costs >= 1)",
        ),
        ((2, [0], [1], [1], [1.5], 1), "edge_costs[0] must be an integer (got float)"),
        # The remaining bounds.
        (
            (2, [0], [2], [1], [1], 1),
            "edge 0: endpoint (0, 2) out of node range [0, 2)",
        ),
        (
            (2, [0, 0], [1, 1], [L62 // 2, L62 // 2], [1, 1], 1),
            "the sum of all capacities reaches the bound of the units domain (2^62 units)",
        ),
        (
            (2, [0, 0], [1, 1], [1, 1], [L62 // 2, L62 // 2], 1),
            "the sum of all edge costs reaches the bound 2^62: path costs would overflow int64",
        ),
    ],
)
def test_validation_errors(args, message):
    # A value of a type its position never takes is a TypeError, every other defect a ValueError.
    cls = TypeError if message in TYPE_DEFECTS else ValueError
    with pytest.raises((TypeError, ValueError), match=exact(message)) as error:
        ngc.Graph(*args)
    assert error.type is cls


TYPE_DEFECTS = {
    "edge_capacity_units[0] must be an integer (got float)",
    "edge_capacity_units entries must be integers (got dtype float64)",
    "edge_capacity_units[0] must be an integer (got bool)",
    "num_nodes must be an integer (got float)",
    "edge_costs[0] must be an integer (got float)",
}


def test_validation_order():
    # The unit values before the shapes (the length message names edge_capacity_units too), the
    # shapes before the costs, the costs before the remaining bounds (endpoints and sums).
    with pytest.raises(
        ValueError,
        match=exact(
            "edge_capacity_units[0] = -1 is outside the units domain [0, 2^62)"
        ),
    ):
        ngc.Graph(2, [0, 1], [1], [-1], [1], 1)
    with pytest.raises(
        ValueError,
        match=exact(
            "edge_sources, edge_destinations, edge_capacity_units and edge_costs must have the same length (got 2, 1, 1, 1)"
        ),
    ):
        ngc.Graph(2, [0, 1], [1], [1], [0], 1)
    with pytest.raises(
        ValueError,
        match=exact("edge_costs[0] = 0: costs must be positive (edge_costs >= 1)"),
    ):
        ngc.Graph(2, [0], [9], [1], [0], 1)
    with pytest.raises(
        ValueError,
        match=exact("edge_costs[0] = 0: costs must be positive (edge_costs >= 1)"),
    ):
        ngc.Graph(2, [0, 0], [1, 1], [L62 // 2, L62 // 2], [0, 1], 1)
    with pytest.raises(
        ValueError, match=exact("edge 1: endpoint (5, 1) out of node range [0, 2)")
    ):
        ngc.Graph(2, [0, 5], [1, 1], [L62 // 2, L62 // 2], [1, 1], 1)
    # Within the shapes: num_nodes before the lengths.
    with pytest.raises(
        ValueError, match=exact("num_nodes = -1 is out of range [0, 2147483647]")
    ):
        ngc.Graph(-1, [0, 1], [1], [1], [0], 1)
    # The sums run edge by edge in input order: the first edge at which a running sum reaches
    # 2^62 reports that sum (here the cost sum at edge 1, the capacity sum only at edge 2), and
    # at the same edge the capacity sum comes first.
    with pytest.raises(
        ValueError,
        match=exact(
            "the sum of all edge costs reaches the bound 2^62: path costs would overflow int64"
        ),
    ):
        ngc.Graph(
            2, [0, 0, 0], [1, 1, 1], [1, L62 // 2, L62 // 2], [L62 // 2] * 2 + [1], 1
        )
    with pytest.raises(
        ValueError,
        match=exact(
            "the sum of all capacities reaches the bound of the units domain (2^62 units)"
        ),
    ):
        ngc.Graph(2, [0, 0], [1, 1], [L62 // 2] * 2, [L62 // 2] * 2, 1)


def test_bounds_just_below_are_accepted():
    g = ngc.Graph(
        2, [0, 0], [1, 1], [L62 // 2, L62 // 2 - 1], [L62 // 2, L62 // 2 - 1], 1
    )
    assert sum(int(x) for x in g.edge_capacity_units()) == L62 - 1
    assert sum(int(x) for x in g.edge_costs()) == L62 - 1


def test_eligibility():
    g = build(4, [(0, 1, 1, 1), (1, 2, 1, 1), (2, 3, 1, 1), (1, 1, 1, 2), (3, 0, 1, 1)])
    every = g.eligibility()
    assert every.uid() == 0 and every.graph_uid() == g.uid()
    assert every.node_allowed().tolist() == [True] * 4
    el = g.eligibility(excluded_nodes=[1, 1])
    assert el.uid() not in (0, every.uid())
    assert el.node_allowed().tolist() == [True, False, True, True]
    assert by_input(g, el.edge_allowed().tolist()) == [False, False, True, False, True]
    # Both flag arrays are read-only boolean views.
    for view in (el.node_allowed(), el.edge_allowed()):
        assert view.dtype == np.bool_
        assert not view.flags.writeable
        with pytest.raises(ValueError):
            view[0] = True
    # Built only by Graph.eligibility and Eligibility.intersect.
    with pytest.raises(TypeError):
        ngc.Eligibility()
    xe = g.eligibility(excluded_edges=np.array([4, 0, 4]))
    assert xe.edge_allowed().tolist() == [False, True, True, True, False]
    assert g.eligibility([], [0]).uid() != g.eligibility([], [0]).uid()
    with pytest.raises(
        ValueError, match=exact("excluded node 4 out of node range [0, 4)")
    ):
        g.eligibility(excluded_nodes=[4])
    with pytest.raises(
        ValueError, match=exact("excluded edge -1 out of edge range [0, 5)")
    ):
        g.eligibility(excluded_edges=[-1])
    with pytest.raises(
        TypeError, match=exact("excluded_edges[0] must be an integer (got float)")
    ):
        g.eligibility(excluded_edges=[0.0])


# The binding validates in the documented order on unnarrowed integers, so the
# reported defect does not depend on how the same defective input is represented.
@pytest.mark.parametrize(
    "args, message",
    [
        # A zero or negative cost comes before any remaining bound, even one that does not fit
        # the C++ type (endpoints beyond int32, a cost beyond int64).
        (
            (2, [2**31], [1], [1], [0], 1.0),
            "edge_costs[0] = 0: costs must be positive (edge_costs >= 1)",
        ),
        (
            (2, [-1], [1], [1], [0], 1.0),
            "edge_costs[0] = 0: costs must be positive (edge_costs >= 1)",
        ),
        (
            (2, [0], [2**31], [1], [-1], 1.0),
            "edge_costs[0] = -1: costs must be positive (edge_costs >= 1)",
        ),
        (
            (2, [0], [-(2**70)], [1], [-(2**70)], 1.0),
            "edge_costs[0] = -1180591620717411303424: costs must be positive (edge_costs >= 1)",
        ),
        (
            (2, [0, 0], [1, 1], [1, 1], [2**63, 0], 1.0),
            "edge_costs[1] = 0: costs must be positive (edge_costs >= 1)",
        ),
        # The unit values come before the costs, the resolution before everything.
        (
            (2, [0], [1], [2**62], [0], 1.0),
            "edge_capacity_units[0] = 4611686018427387904 is outside the units domain [0, 2^62)",
        ),
        (
            (2, [0], [1], [2**64], [0], 1.0),
            "edge_capacity_units[0] = 18446744073709551616 is outside the units domain [0, 2^62)",
        ),
        (
            (2, [0], [1], [-(2**70)], [0], 1.0),
            "edge_capacity_units[0] = -1180591620717411303424 is outside the units domain [0, 2^62)",
        ),
        (
            (2, [0], [1], [-1], [0], 3.0),
            "resolution must be a power of two in [2^-32, 2^32]",
        ),
        # The remaining bounds name the exact value: ids, then the sums.
        (
            (2, [2**40], [1], [1], [1], 1.0),
            "edge 0: endpoint (1099511627776, 1) out of node range [0, 2)",
        ),
        (
            (2, [0], [-(2**70)], [1], [1], 1.0),
            "edge 0: endpoint (0, -1180591620717411303424) out of node range [0, 2)",
        ),
        (
            (2, [0], [1], [1], [2**63], 1.0),
            "the sum of all edge costs reaches the bound 2^62: path costs would overflow int64",
        ),
        (
            (2, [0, 5], [1, 1], [1, 1], [2**64, 1], 1.0),
            "edge 1: endpoint (5, 1) out of node range [0, 2)",
        ),
        # The shapes.
        (
            (2**40, [], [], [], [], 1.0),
            "num_nodes = 1099511627776 is out of range [0, 2147483647]",
        ),
    ],
)
def test_validation_order_before_narrowing(args, message):
    with pytest.raises(ValueError, match=exact(message)):
        ngc.Graph(*args)


def test_validation_independent_of_representation():
    # The same defective graph as Python lists and as NumPy arrays (uint64 beyond int64).
    as_lists = (2, [2**63], [1], [1], [0], 1.0)
    as_arrays = (
        2,
        np.array([2**63], dtype=np.uint64),
        np.array([1], dtype=np.int32),
        np.array([1], dtype=np.int64),
        np.array([0], dtype=np.int64),
        1.0,
    )
    messages = []
    for args in (as_lists, as_arrays):
        with pytest.raises(ValueError) as info:
            ngc.Graph(*args)
        messages.append(str(info.value))
    assert (
        messages[0]
        == messages[1]
        == "edge_costs[0] = 0: costs must be positive (edge_costs >= 1)"
    )
    with pytest.raises(
        ValueError,
        match=exact(
            "edge 0: endpoint (9223372036854775808, 1) out of node range [0, 2)"
        ),
    ):
        ngc.Graph(2, np.array([2**63], dtype=np.uint64), [1], [1], [1], 1.0)


def test_eligibility_ids_checked_before_narrowing():
    g = build(3, EDGES)
    with pytest.raises(
        ValueError, match=exact("excluded node 1099511627776 out of node range [0, 3)")
    ):
        g.eligibility(excluded_nodes=[2**40])
    with pytest.raises(
        ValueError,
        match=exact("excluded edge -1180591620717411303424 out of edge range [0, 7)"),
    ):
        g.eligibility(excluded_edges=[-(2**70)])
    with pytest.raises(
        ValueError, match=exact("excluded node 5 out of node range [0, 3)")
    ):  # nodes before edges
        g.eligibility(excluded_nodes=[5], excluded_edges=[2**70])


ARG_INDEX = {
    "edge_sources": 1,
    "edge_destinations": 2,
    "edge_capacity_units": 3,
    "edge_costs": 4,
}
VALID = {
    "edge_sources": 0,
    "edge_destinations": 1,
    "edge_capacity_units": 1,
    "edge_costs": 0,
}  # the cost stays zero


@pytest.mark.parametrize("name", list(ARG_INDEX))
@pytest.mark.parametrize("ndim", [0, 2])
@pytest.mark.parametrize("cost", [0, -1])
def test_array_rank_is_a_shape_defect_before_costs(name, ndim, cost):
    # The shape step checks rank 1 for every array argument before the
    # cost step, so a 2-D or 0-D array is reported before a zero or negative cost (a rank is a
    # type: TypeError).
    args = [2, [0], [1], [1], [cost], 1.0]
    value = cost if name == "edge_costs" else VALID[name]
    bad = np.array(value, dtype=np.int64)
    args[ARG_INDEX[name]] = bad.reshape(1, 1) if ndim == 2 else bad
    pattern = exact(f"{name} must be a 1-D array (got {ndim} dimensions)")
    with pytest.raises(TypeError, match=pattern):
        ngc.Graph(*args)


@pytest.mark.parametrize("name", list(ARG_INDEX))
def test_nested_and_scalar_arguments_are_shape_defects(name):
    value = VALID[name]
    nested = [2, [0], [1], [1], [0], 1.0]
    nested[ARG_INDEX[name]] = [[value]]
    with pytest.raises(
        TypeError,
        match=exact(
            f"{name} must be a 1-D sequence of integers (entry 0 is a sequence)"
        ),
    ):
        ngc.Graph(*nested)
    scalar = [2, [0], [1], [1], [0], 1.0]
    scalar[ARG_INDEX[name]] = value
    with pytest.raises(
        TypeError, match=exact(f"{name} must be a 1-D sequence of integers (got int)")
    ):
        ngc.Graph(*scalar)


def test_shape_step_keeps_its_place_in_the_order():
    # The resolution and the unit values still come before the shapes, whatever the shape.
    with pytest.raises(
        ValueError, match=exact("resolution must be a power of two in [2^-32, 2^32]")
    ):
        ngc.Graph(2, np.array([[0]]), [1], [1], [0], 3.0)
    with pytest.raises(
        ValueError,
        match=exact(
            "edge_capacity_units[0] = 4611686018427387904 is outside the units domain [0, 2^62)"
        ),
    ):
        ngc.Graph(2, [0], [1], np.array([[2**62]], dtype=np.int64), [0], 1.0)
    with pytest.raises(
        ValueError,
        match=exact(
            "edge_capacity_units[0] = -1 is outside the units domain [0, 2^62)"
        ),
    ):
        ngc.Graph(2, [0], [1], np.array(-1, dtype=np.int64), [0], 1.0)
    with pytest.raises(
        ValueError,
        match=exact(
            "edge_capacity_units[1] = -5 is outside the units domain [0, 2^62)"
        ),
    ):
        ngc.Graph(2, [0], [1], [[1], -5], [0], 1.0)
    with pytest.raises(
        ValueError,
        match=exact("edge_capacity_units = -5 is outside the units domain [0, 2^62)"),
    ):
        ngc.Graph(2, [0], [1], -5, [0], 1.0)
    # With valid values the misshapen capacity is a shape defect, before the costs.
    with pytest.raises(
        TypeError,
        match=exact("edge_capacity_units must be a 1-D array (got 2 dimensions)"),
    ):
        ngc.Graph(2, [0], [1], np.array([[1]]), [0], 1.0)
    # Rank (a type) before lengths (a value).
    with pytest.raises(
        TypeError,
        match=exact("edge_destinations must be a 1-D array (got 2 dimensions)"),
    ):
        ngc.Graph(2, [0, 1], np.array([[1, 0]]), [1, 1], [1, 1], 1.0)
    # A valid graph from 1-D arrays of every integer kind still builds.
    g = ngc.Graph(
        2, np.array([0]), np.array([1], dtype=np.uint8), np.array([7]), (3,), 1.0
    )
    assert g.edge_capacity_units().tolist() == [7]


def test_repeated_views_share_storage_and_refuse_write_access():
    g = build(3, EDGES)
    for name in (
        "edge_capacity_units",
        "edge_costs",
        "edge_sources",
        "edge_destinations",
        "input_positions",
        "outgoing_offsets",
        "incoming_offsets",
        "incoming_edges",
    ):
        a, b = getattr(g, name)(), getattr(g, name)()
        assert np.shares_memory(a, b), name
        with pytest.raises(ValueError):
            a.setflags(write=True)
        assert not a.flags.writeable


VIEW_DTYPES = {
    "edge_capacity_units": np.int64,
    "edge_costs": np.int64,
    "edge_sources": np.int32,
    "edge_destinations": np.int32,
    "input_positions": np.int32,
    "outgoing_offsets": np.int32,
    "incoming_offsets": np.int32,
    "incoming_edges": np.int32,
}


@pytest.mark.parametrize(
    ("n", "edges"),
    [(0, []), (100, []), (3, EDGES)],
    ids=["empty", "isolated", "nonempty"],
)
def test_every_view_keeps_its_graph_alive_and_refuses_write_access(n, edges):
    # An empty view has no bytes to share, so the owner is checked by identity and a weak
    # reference rather than by np.shares_memory. The values stay readable after the caller
    # drops the graph.
    for name, dtype in VIEW_DTYPES.items():
        g = build(n, edges)
        ref = weakref.ref(g)
        view = getattr(g, name)()
        values = view.tolist()
        assert view.base is g, name
        assert view.dtype == dtype, name
        assert not view.flags.writeable, name
        with pytest.raises(ValueError):
            view[...] = 0
        with pytest.raises(ValueError):
            view.setflags(write=True)
        del g
        gc.collect()
        assert ref() is not None, name
        assert view.tolist() == values, name
        del view
        gc.collect()
        assert ref() is None, name


def _reference_layout(n, src, dst, cost):
    # EdgeIds: (source, cost, destination, input); the reverse slots: (destination, cost, source, input).
    i = np.arange(len(src))
    order = np.lexsort((i, dst, cost, src))
    pos = np.empty_like(order)
    pos[order] = np.arange(len(order))
    rev = pos[np.lexsort((i, src, cost, dst))]
    return order, np.searchsorted(src[order], np.arange(n + 1)), rev


@pytest.mark.parametrize(
    "shape",
    ["hash", "combined_keys", "all_cost_sort", "reversed", "edgeless"],
)
def test_construction_orders_match_a_reference(shape):
    # One input per construction path of the graph: distinct-cost hashing, the combined-key
    # passes (N * K < E), the all-cost sort (K > 65,536), input in reversed EdgeId order, no edges.
    rng = np.random.default_rng(7)
    n, e = 300, 1000  # N * K = 1500 > E: the separate passes
    if shape == "combined_keys":
        n, e = 20, 4000  # N * K = 100 < E
    if shape == "edgeless":
        e = 0
    src = rng.integers(0, n, e)
    dst = rng.integers(0, n, e)
    cost = rng.integers(1, 6, e)
    if shape == "all_cost_sort":
        e = 70000
        src, dst = rng.integers(0, n, e), rng.integers(0, n, e)
        cost = rng.permutation(np.arange(1, e + 1))
    cap = rng.integers(0, 1 << 40, e)
    if shape == "reversed":
        order, _, _ = _reference_layout(n, src, dst, cost)
        src, dst, cap, cost = (a[order[::-1]] for a in (src, dst, cap, cost))
    g = ngc.Graph(n, src, dst, cap, cost, 1)
    order, row, rev = _reference_layout(n, src, dst, cost)
    assert g.input_positions().tolist() == order.tolist()
    assert g.outgoing_offsets().tolist() == row.tolist()
    assert g.incoming_edges().tolist() == rev.tolist()
    assert (
        g.incoming_offsets().tolist()
        == np.searchsorted(np.sort(dst), np.arange(n + 1)).tolist()
    )
    assert g.edge_capacity_units().tolist() == cap[order].tolist()
    assert g.max_edge_cost() == (int(cost.max()) if e else 0)
