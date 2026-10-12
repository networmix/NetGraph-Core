"""Graph on integer units (APPROACH 3): construction, validation order, B' ids,
read-only views, Eligibility."""

from __future__ import annotations

import gc
import weakref

import numpy as np
import pytest

import netgraph_core as ngc

L62 = 1 << 62

# (u, v, capacity_units, cost) in input order.
EDGES = [
    (2, 0, 5, 3),
    (0, 2, 4, 2),
    (0, 1, 3, 2),
    (1, 2, 2, 1),
    (0, 1, 1, 2),
    (0, 2, 6, 1),
    (2, 1, 7, 1),
]


def build(n, edges, q=1.0):
    return ngc.Graph(
        n,
        [e[0] for e in edges],
        [e[1] for e in edges],
        [e[2] for e in edges],
        [e[3] for e in edges],
        q,
    )


def test_bprime_ids_and_input_index():
    g = build(3, EDGES)
    order = sorted(
        range(len(EDGES)), key=lambda i: (EDGES[i][0], EDGES[i][3], EDGES[i][1], i)
    )
    assert g.input_index().tolist() == order == [5, 2, 4, 1, 3, 6, 0]
    assert g.row_offsets().tolist() == [0, 4, 5, 7]
    assert g.edge_src().tolist() == [EDGES[i][0] for i in order]
    assert g.edge_dst().tolist() == [EDGES[i][1] for i in order]
    assert g.capacity_units().tolist() == [EDGES[i][2] for i in order]
    assert g.cost().tolist() == [EDGES[i][3] for i in order]
    inv = np.empty(len(EDGES), dtype=np.int64)
    inv[g.input_index()] = np.arange(len(EDGES))
    assert g.input_index()[inv].tolist() == list(range(len(EDGES)))
    # Reverse CSR in (cost, src, input) order, as input positions.
    rev = []
    off, adj = g.in_row_offsets(), g.in_adj_edge_index()
    for v in range(3):
        rev.append([int(g.input_index()[e]) for e in adj[off[v] : off[v + 1]]])
    assert rev == [[0], [6, 2, 4], [5, 3, 1]]
    assert g.max_cost() == 3
    assert g.cost_gcd() == 1
    assert g.num_nodes() == 3 and g.num_edges() == 7
    assert g.resolution() == 1.0


def test_views_are_read_only_and_keep_the_graph_alive():
    g = build(3, EDGES)
    views = [
        g.capacity_units(),
        g.cost(),
        g.edge_src(),
        g.edge_dst(),
        g.input_index(),
        g.row_offsets(),
        g.in_row_offsets(),
        g.in_adj_edge_index(),
    ]
    dtypes = [
        np.int64,
        np.int64,
        np.int32,
        np.int32,
        np.int32,
        np.int32,
        np.int32,
        np.int32,
    ]
    for view, dtype in zip(views, dtypes, strict=True):
        assert view.dtype == dtype
        assert not view.flags.writeable
        with pytest.raises(ValueError):
            view[0] = 0
    expect = [v.tolist() for v in views]
    del g
    gc.collect()
    assert [v.tolist() for v in views] == expect


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
    assert g.capacity_units().tolist() == [7, big]
    assert g.resolution() == 2.0**-12


def test_empty_and_isolated():
    g = ngc.Graph(0, [], [], [], [], 1)
    assert g.row_offsets().tolist() == [0] and g.in_row_offsets().tolist() == [0]
    assert g.max_cost() == 0 and g.cost_gcd() == 0
    h = ngc.Graph(3, [], [], [], [], 1)
    assert h.row_offsets().tolist() == [0, 0, 0, 0]


def test_uids_are_distinct():
    a, b = build(3, EDGES), build(3, EDGES)
    assert a.uid() != b.uid() and a.uid() > 0


@pytest.mark.parametrize(
    "args, pattern",
    [
        # The resolution.
        ((2, [0], [1], [1], [1], 3.0), "resolution"),
        ((2, [0], [1], [1], [1], 2.0**-33), "resolution"),
        ((2, [0], [1], [1], [1], 0.0), "resolution"),
        ((2, [0], [1], [1], [1], float("nan")), "resolution"),
        # The unit values.
        ((2, [0], [1], [-1], [1], 1), "capacity_units"),
        ((2, [0], [1], [L62], [1], 1), "capacity_units"),
        ((2, [0], [1], [1.0], [1], 1), "capacity_units"),
        ((2, [0], [1], np.array([1.5]), [1], 1), "capacity_units"),
        ((2, [0], [1], [True], [1], 1), "capacity_units"),
        # The shapes.
        ((-1, [], [], [], [], 1), "num_nodes"),
        ((2.0, [0], [1], [1], [1], 1), "num_nodes"),
        ((2, [0, 1], [1], [1], [1], 1), "same length"),
        ((2, [0], [1], [1, 1], [1], 1), "same length"),
        # Zero and negative costs (rule 28).
        ((2, [0], [1], [1], [0], 1), "cost"),
        ((2, [0], [1], [1], [-3], 1), "cost"),
        ((2, [0], [1], [1], [1.5], 1), "cost"),
        # The remaining bounds.
        ((2, [0], [2], [1], [1], 1), "out of node range"),
        ((2, [0, 0], [1, 1], [L62 // 2, L62 // 2], [1, 1], 1), "sum of all capacities"),
        ((2, [0, 0], [1, 1], [1, 1], [L62 // 2, L62 // 2], 1), "sum of all edge costs"),
    ],
)
def test_validation_errors(args, pattern):
    with pytest.raises(ValueError, match=pattern):
        ngc.Graph(*args)


def test_validation_order():
    # The unit values before the shapes, the shapes before the costs, the costs before the
    # remaining bounds (endpoints and sums).
    with pytest.raises(ValueError, match="capacity_units"):
        ngc.Graph(2, [0, 1], [1], [-1], [1], 1)
    with pytest.raises(ValueError, match="same length"):
        ngc.Graph(2, [0, 1], [1], [1], [0], 1)
    with pytest.raises(ValueError, match="cost"):
        ngc.Graph(2, [0], [9], [1], [0], 1)
    with pytest.raises(ValueError, match="cost"):
        ngc.Graph(2, [0, 0], [1, 1], [L62 // 2, L62 // 2], [0, 1], 1)
    with pytest.raises(ValueError, match="node range"):
        ngc.Graph(2, [0, 5], [1, 1], [L62 // 2, L62 // 2], [1, 1], 1)


def test_bounds_just_below_are_accepted():
    g = ngc.Graph(
        2, [0, 0], [1, 1], [L62 // 2, L62 // 2 - 1], [L62 // 2, L62 // 2 - 1], 1
    )
    assert sum(int(x) for x in g.capacity_units()) == L62 - 1
    assert sum(int(x) for x in g.cost()) == L62 - 1


def test_eligibility():
    g = build(4, [(0, 1, 1, 1), (1, 2, 1, 1), (2, 3, 1, 1), (1, 1, 1, 2), (3, 0, 1, 1)])
    every = g.eligibility()
    assert every.uid() == 0 and every.graph_uid() == g.uid()
    assert every.node_allowed().tolist() == [True] * 4
    el = g.eligibility(excluded_nodes=[1, 1])
    assert el.uid() not in (0, every.uid())
    assert el.node_allowed().tolist() == [True, False, True, True]
    allowed_in = [False] * 5
    for e, ok in enumerate(el.edge_allowed().tolist()):
        allowed_in[int(g.input_index()[e])] = ok
    assert allowed_in == [False, False, True, False, True]
    assert not el.edge_allowed().flags.writeable
    xe = g.eligibility(excluded_edges=np.array([4, 0, 4]))
    assert xe.edge_allowed().tolist() == [False, True, True, True, False]
    assert g.eligibility([], [0]).uid() != g.eligibility([], [0]).uid()
    with pytest.raises(ValueError, match="range"):
        g.eligibility(excluded_nodes=[4])
    with pytest.raises(ValueError, match="range"):
        g.eligibility(excluded_edges=[-1])
    with pytest.raises(ValueError):
        g.eligibility(excluded_edges=[0.0])


# review-1a D1: the binding validates in the contract's order on unnarrowed integers, so the
# reported defect does not depend on how the same defective input is represented.
@pytest.mark.parametrize(
    "args, pattern",
    [
        # A zero or negative cost comes before any remaining bound, even one that does not fit
        # the C++ type (endpoints beyond int32, a cost beyond int64).
        ((2, [2**31], [1], [1], [0], 1.0), r"cost\[0\] = 0: costs must be positive"),
        ((2, [-1], [1], [1], [0], 1.0), r"cost\[0\] = 0: costs must be positive"),
        ((2, [0], [2**31], [1], [-1], 1.0), r"cost\[0\] = -1: costs must be positive"),
        (
            (2, [0], [-(2**70)], [1], [-(2**70)], 1.0),
            r"cost\[0\] = -1180591620717411303424",
        ),
        (
            (2, [0, 0], [1, 1], [1, 1], [2**63, 0], 1.0),
            r"cost\[1\] = 0: costs must be positive",
        ),
        # The unit values come before the costs, the resolution before everything.
        (
            (2, [0], [1], [2**62], [0], 1.0),
            r"capacity_units\[0\] = 4611686018427387904",
        ),
        (
            (2, [0], [1], [2**64], [0], 1.0),
            r"capacity_units\[0\] = 18446744073709551616",
        ),
        (
            (2, [0], [1], [-(2**70)], [0], 1.0),
            r"capacity_units\[0\] = -1180591620717411303424",
        ),
        ((2, [0], [1], [-1], [0], 3.0), "resolution"),
        # The remaining bounds name the exact value: ids, then the sums.
        (
            (2, [2**40], [1], [1], [1], 1.0),
            r"edge 0: endpoint \(1099511627776, 1\) out of node range",
        ),
        (
            (2, [0], [-(2**70)], [1], [1], 1.0),
            r"endpoint \(0, -1180591620717411303424\)",
        ),
        ((2, [0], [1], [1], [2**63], 1.0), "sum of all edge costs"),
        ((2, [0, 5], [1, 1], [1, 1], [2**64, 1], 1.0), "out of node range"),
        # The shapes.
        ((2**40, [], [], [], [], 1.0), "num_nodes"),
    ],
)
def test_validation_order_before_narrowing(args, pattern):
    with pytest.raises(ValueError, match=pattern):
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
        messages[0] == messages[1] == "cost[0] = 0: costs must be positive (cost >= 1)"
    )
    with pytest.raises(ValueError, match=r"endpoint \(9223372036854775808, 1\)"):
        ngc.Graph(2, np.array([2**63], dtype=np.uint64), [1], [1], [1], 1.0)


def test_eligibility_ids_checked_before_narrowing():
    g = build(3, EDGES)
    with pytest.raises(
        ValueError, match=r"excluded node 1099511627776 out of node range \[0, 3\)"
    ):
        g.eligibility(excluded_nodes=[2**40])
    with pytest.raises(
        ValueError, match=r"excluded edge -1180591620717411303424 out of edge range"
    ):
        g.eligibility(excluded_edges=[-(2**70)])
    with pytest.raises(ValueError, match="excluded node"):  # nodes before edges
        g.eligibility(excluded_nodes=[5], excluded_edges=[2**70])


ARG_INDEX = {"src": 1, "dst": 2, "capacity_units": 3, "cost": 4}
VALID = {"src": 0, "dst": 1, "capacity_units": 1, "cost": 0}  # the cost stays zero


@pytest.mark.parametrize("name", list(ARG_INDEX))
@pytest.mark.parametrize("ndim", [0, 2])
@pytest.mark.parametrize("cost", [0, -1])
def test_array_rank_is_a_shape_defect_before_costs(name, ndim, cost):
    # review-1a-fixes D1: the shape step checks rank 1 for every array argument before the
    # cost step, so a 2-D or 0-D array is reported before a zero or negative cost.
    args = [2, [0], [1], [1], [cost], 1.0]
    value = cost if name == "cost" else VALID[name]
    bad = np.array(value, dtype=np.int64)
    args[ARG_INDEX[name]] = bad.reshape(1, 1) if ndim == 2 else bad
    pattern = rf"^{name} must be a 1-D array \(got {ndim} dimensions\)$"
    with pytest.raises(ValueError, match=pattern):
        ngc.Graph(*args)


@pytest.mark.parametrize("name", list(ARG_INDEX))
def test_nested_and_scalar_arguments_are_shape_defects(name):
    value = VALID[name]
    nested = [2, [0], [1], [1], [0], 1.0]
    nested[ARG_INDEX[name]] = [[value]]
    with pytest.raises(
        ValueError, match=rf"^{name} must be a 1-D sequence of integers \(entry 0"
    ):
        ngc.Graph(*nested)
    scalar = [2, [0], [1], [1], [0], 1.0]
    scalar[ARG_INDEX[name]] = value
    with pytest.raises(
        ValueError, match=rf"^{name} must be a 1-D sequence of integers \(got int\)$"
    ):
        ngc.Graph(*scalar)


def test_shape_step_keeps_its_place_in_the_order():
    # The resolution and the unit values still come before the shapes, whatever the shape.
    with pytest.raises(ValueError, match="resolution"):
        ngc.Graph(2, np.array([[0]]), [1], [1], [0], 3.0)
    with pytest.raises(
        ValueError, match=r"capacity_units\[0\] = 4611686018427387904 is outside"
    ):
        ngc.Graph(2, [0], [1], np.array([[2**62]], dtype=np.int64), [0], 1.0)
    with pytest.raises(ValueError, match=r"capacity_units\[0\] = -1 is outside"):
        ngc.Graph(2, [0], [1], np.array(-1, dtype=np.int64), [0], 1.0)
    with pytest.raises(ValueError, match=r"capacity_units\[1\] = -5 is outside"):
        ngc.Graph(2, [0], [1], [[1], -5], [0], 1.0)
    with pytest.raises(ValueError, match=r"capacity_units = -5 is outside"):
        ngc.Graph(2, [0], [1], -5, [0], 1.0)
    # With valid values the misshapen capacity is a shape defect, before the costs.
    with pytest.raises(ValueError, match="capacity_units must be a 1-D array"):
        ngc.Graph(2, [0], [1], np.array([[1]]), [0], 1.0)
    # Rank before lengths.
    with pytest.raises(ValueError, match="dst must be a 1-D array"):
        ngc.Graph(2, [0, 1], np.array([[1, 0]]), [1, 1], [1, 1], 1.0)
    # A valid graph from 1-D arrays of every integer kind still builds.
    g = ngc.Graph(
        2, np.array([0]), np.array([1], dtype=np.uint8), np.array([7]), (3,), 1.0
    )
    assert g.capacity_units().tolist() == [7]


def test_repeated_views_share_storage_and_refuse_write_access():
    g = build(3, EDGES)
    for name in (
        "capacity_units",
        "cost",
        "edge_src",
        "edge_dst",
        "input_index",
        "row_offsets",
        "in_row_offsets",
        "in_adj_edge_index",
    ):
        a, b = getattr(g, name)(), getattr(g, name)()
        assert np.shares_memory(a, b), name
        with pytest.raises(ValueError):
            a.setflags(write=True)
        assert not a.flags.writeable


VIEW_DTYPES = {
    "capacity_units": np.int64,
    "cost": np.int64,
    "edge_src": np.int32,
    "edge_dst": np.int32,
    "input_index": np.int32,
    "row_offsets": np.int32,
    "in_row_offsets": np.int32,
    "in_adj_edge_index": np.int32,
}


@pytest.mark.parametrize(
    ("n", "edges"),
    [(0, []), (100, []), (3, EDGES)],
    ids=["empty", "isolated", "nonempty"],
)
def test_every_view_keeps_its_graph_alive_and_refuses_write_access(n, edges):
    # An empty view has no bytes to share, so the owner is checked by identity and a weak
    # reference rather than by np.shares_memory.
    for name, dtype in VIEW_DTYPES.items():
        g = build(n, edges)
        ref = weakref.ref(g)
        view = getattr(g, name)()
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
        del view
        gc.collect()
        assert ref() is None, name


def _reference_layout(n, src, dst, cap, cost):
    # B' ids: (src, cost, dst, input); the reverse slots: (dst, cost, src, input).
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
    # One input per construction path of APPROACH 3: distinct-cost hashing, the combined-key
    # passes (N * K < E), the all-cost sort (K > 65,536), a reversed B' order, no edges.
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
        order, _, _ = _reference_layout(n, src, dst, cap, cost)
        src, dst, cap, cost = (a[order[::-1]] for a in (src, dst, cap, cost))
    g = ngc.Graph(n, src, dst, cap, cost, 1)
    order, row, rev = _reference_layout(n, src, dst, cap, cost)
    assert g.input_index().tolist() == order.tolist()
    assert g.row_offsets().tolist() == row.tolist()
    assert g.in_adj_edge_index().tolist() == rev.tolist()
    assert (
        g.in_row_offsets().tolist()
        == np.searchsorted(np.sort(dst), np.arange(n + 1)).tolist()
    )
    assert g.capacity_units().tolist() == cap[order].tolist()
    assert g.max_cost() == (int(cost.max()) if e else 0)
