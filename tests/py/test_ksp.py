"""KSP and resolve_to_paths through the binding (docs/semantics.md#k-shortest-paths): the
order, Expand and Collapse, the ceilings, k = 1 against the single path of spf, the validation
order on unnarrowed values, and the all-or-nothing construction of the result lists."""

from __future__ import annotations

from fractions import Fraction

import numpy as np
import pytest
from helpers import (
    CHAIN,
    FEWER,
    HIGHER_CAP,
    assert_read_only,
    build,
    exact,
    sweep_python_allocation,
)
from private_hooks import needs_hooks

import netgraph_core as ngc


def bgraph(n, rows):
    """A graph from rows (u, v, capacity, cost), sorted by (source, cost, destination) first: that is the
    EdgeId order, so each EdgeId is the row's position in the sorted list."""
    return build(n, sorted(rows, key=lambda r: (r[0], r[3], r[1])))


def ksp(g, s, t, k, *, parallel_edges=True, tie=FEWER, **kw):
    return ngc.ksp(
        g,
        s,
        t,
        ngc.KspOptions(k, tie_break=tie, parallel_edges=parallel_edges, **kw),
    )


def test_order_and_ceilings():
    g = bgraph(4, CHAIN)
    got = ksp(g, 0, 3, 8)
    assert [c for c, _ in got] == [3, 4, 4]
    assert all(isinstance(p, ngc.Path) and p.cost() == c for c, p in got)
    assert got[0][1].hops() == [[0], [2], [4]]
    assert [c for c, _ in ksp(g, 0, 3, 8, cost_ceiling=ngc.CostStretch(1.0))] == [3]
    assert [c for c, _ in ksp(g, 0, 3, 8, cost_ceiling=ngc.CostStretch(1.5))] == [
        3,
        4,
        4,
    ]
    assert [c for c, _ in ksp(g, 0, 3, 8, cost_ceiling=ngc.AbsoluteCost(3))] == [3]
    assert ksp(g, 0, 3, 8, cost_ceiling=ngc.AbsoluteCost(0)) == []
    assert len(ksp(g, 0, 3, 2)) == 2
    # The CostStretch baseline is static: gating the middle edge drops the cost-3 path only.
    residual = [1, 1, 0, 1, 1]
    assert ksp(g, 0, 3, 8, cost_ceiling=ngc.CostStretch(1.0), residual=residual) == []
    assert [c for c, _ in ksp(g, 0, 3, 8, residual=residual)] == [4, 4]
    excluded = g.eligibility([0], [])
    assert ksp(g, 0, 3, 8, eligibility=excluded) == []


def test_expand_collapse_and_the_single_path():
    pr = bgraph(3, [(0, 1, 1, 1), (0, 1, 1, 1), (1, 2, 1, 1), (1, 2, 1, 1)])
    assert len(ksp(pr, 0, 2, 8)) == 4
    assert len(ksp(pr, 0, 2, 8, parallel_edges=False)) == 1
    assert [p.hops() for _c, p in ksp(pr, 0, 2, 3)] == [
        [[0], [2]],
        [[0], [3]],
        [[1], [2]],
    ]
    three = bgraph(
        6,
        [
            (0, 2, 1, 1),
            (2, 3, 1, 1),
            (3, 5, 1, 1),
            (0, 1, 1, 1),
            (1, 4, 1, 1),
            (4, 5, 1, 1),
        ],
    )
    for tie in (FEWER, HIGHER_CAP):
        first = ksp(three, 0, 5, 1, tie=tie)[0][1]
        for direction, roots, targets in (
            (ngc.Direction.FORWARD, 0, 5),
            (ngc.Direction.REVERSE, 5, 0),
        ):
            options = ngc.SpfOptions(
                direction, tie_break=tie, multipath=False, targets=targets
            )
            r = ngc.spf(three, roots, options)
            assert r.dag == first.to_dag()
    # ksp has no multipath field: k sets the number of paths. Every tie_break is valid with
    # every parallel_edges (ksp always orders paths, so no tie rule is ineffective).
    with pytest.raises(TypeError):
        ngc.KspOptions(8, multipath=False)
    for g, t in ((pr, 2), (three, 5)):
        for tie, pe in ((FEWER, True), (HIGHER_CAP, True), (FEWER, False)):
            assert ksp(g, 0, t, 8, tie=tie, parallel_edges=pe)


def test_ksp_validation_order():
    # The residual's own checks are the shared ones of test_edge_amounts.py.
    g = bgraph(3, [(0, 1, 1, 1), (1, 2, 1, 1)])
    other = bgraph(3, [(0, 1, 1, 1)])
    with pytest.raises(
        ValueError, match=exact("source: node 3 out of node range [0, 3)")
    ):
        ksp(g, 3, 1, 1)
    with pytest.raises(
        ValueError, match=exact("destination: node -1 out of node range [0, 3)")
    ):
        ksp(g, 0, -1, 1)
    with pytest.raises(
        ValueError,
        match=exact(
            "destination: node 1180591620717411303424 out of node range [0, 3)"
        ),
    ):
        ksp(g, 0, 2**70, 1)
    with pytest.raises(
        ValueError,
        match=exact(
            "node 1 is in both endpoint sets: source and destination must not overlap"
        ),
    ):
        ksp(g, 1, 1, 1)
    with pytest.raises(ValueError, match=exact("k must be >= 1 (got 0)")):
        ksp(g, 0, 2, 0)
    with pytest.raises(ValueError, match=exact("k must be >= 1 (got -1)")):
        ksp(g, 0, 2, -1)
    with pytest.raises(
        ValueError, match=exact("k must be >= 1 (got -1180591620717411303424)")
    ):
        ksp(g, 0, 2, -(2**70))
    with pytest.raises(ValueError, match=exact("k must be an integer (got float)")):
        ksp(g, 0, 2, 1.5)
    absolute = "cost_ceiling AbsoluteCost cost must be >= 0 (got -1)"
    stretch = "cost_ceiling CostStretch factor must be >= 1 (got {})"
    with pytest.raises(ValueError, match=exact(absolute)):
        ksp(g, 0, 2, 1, cost_ceiling=ngc.AbsoluteCost(-1))
    with pytest.raises(ValueError, match=exact(stretch.format("1/2"))):
        ksp(g, 0, 2, 1, cost_ceiling=ngc.CostStretch(0.5))
    with pytest.raises(
        TypeError,
        match=exact("cost_ceiling must be AbsoluteCost, CostStretch or None"),
    ):
        ksp(g, 0, 2, 1, cost_ceiling=3)
    with pytest.raises(
        ValueError, match=exact("foreign Eligibility: built for another graph")
    ):
        ksp(g, 0, 2, 1, eligibility=other.eligibility([], [0]))
    # The order: the endpoints (source, destination, source == destination), k, the ceiling, the residual.
    with pytest.raises(
        ValueError, match=exact("source: node 3 out of node range [0, 3)")
    ):
        ksp(g, 3, -1, 0, cost_ceiling=ngc.AbsoluteCost(-1), residual=[9])
    with pytest.raises(
        ValueError, match=exact("destination: node -1 out of node range [0, 3)")
    ):
        ksp(g, 0, -1, 0, cost_ceiling=ngc.AbsoluteCost(-1), residual=[9])
    with pytest.raises(
        ValueError,
        match=exact(
            "node 2 is in both endpoint sets: source and destination must not overlap"
        ),
    ):
        ksp(g, 2, 2, 0, cost_ceiling=ngc.AbsoluteCost(-1), residual=[9])
    with pytest.raises(ValueError, match=exact("k must be >= 1 (got 0)")):
        ksp(g, 0, 2, 0, cost_ceiling=ngc.AbsoluteCost(-1), residual=[9])
    with pytest.raises(ValueError, match=exact(absolute)):
        ksp(g, 0, 2, 1, cost_ceiling=ngc.AbsoluteCost(-1), residual=[9])
    # A huge k lists every path (no error).
    assert len(ksp(g, 0, 2, 2**70)) == 1
    with pytest.raises(TypeError, match=exact("KspOptions expected (got SpfOptions)")):
        ngc.ksp(g, 0, 2, ngc.SpfOptions(ngc.Direction.FORWARD))
    with pytest.raises(TypeError):
        ngc.KspOptions()  # k is explicit
    with pytest.raises(TypeError):
        ngc.KspOptions(1, ngc.RouteTieBreak.FEWER_HOPS)  # the route choice is keywords
    for bad in (None, 0, "no"):
        with pytest.raises(TypeError):
            ngc.KspOptions(1, parallel_edges=bad)
    # KspOptions keeps its arguments as given (the same objects) and validates none of them:
    # ksp does, at the call, in the order above.
    k, ceiling, residual = 2**70, ngc.AbsoluteCost(-1), [9]
    o = ngc.KspOptions(k, cost_ceiling=ceiling, residual=residual)
    assert o.k is k and o.cost_ceiling is ceiling and o.residual is residual
    assert (o.tie_break, o.parallel_edges, o.eligibility) == (FEWER, True, None)
    assert_read_only(o, "k", "tie_break", "parallel_edges", "cost_ceiling")
    with pytest.raises(
        ValueError, match=exact("cost_ceiling AbsoluteCost cost must be >= 0 (got -1)")
    ):
        ngc.ksp(g, 0, 2, o)
    with pytest.raises(ValueError, match=exact("k must be >= 1 (got 0)")):
        ngc.ksp(g, 0, 2, ngc.KspOptions(0, cost_ceiling=ceiling))


def test_resolve_to_paths():
    # 0->1 twice (cost 1), 0->2 (1), 0->3 (5), 1->3 (1), 2->3 (1): ids in this order.
    g = bgraph(
        4,
        [
            (0, 1, 1, 1),
            (0, 1, 1, 1),
            (0, 2, 1, 1),
            (0, 3, 1, 5),
            (1, 3, 1, 1),
            (2, 3, 1, 1),
        ],
    )
    dag = ngc.PredecessorDAG.from_csr(
        g, [0, 0, 2, 3, 6], [0, 0, 0, 0, 1, 2], [0, 1, 2, 3, 4, 5]
    )
    grouped = ngc.resolve_to_paths(g, dag, 0, 3)
    assert [p.hops() for p in grouped] == [[[0, 1], [4]], [[2], [5]], [[3]]]
    assert [p.cost() for p in grouped] == [2, 2, 5]
    split = ngc.resolve_to_paths(g, dag, 0, 3, split_parallel_edges=True)
    assert [p.hops() for p in split] == [[[0], [4]], [[1], [4]], [[2], [5]], [[3]]]
    assert len(ngc.resolve_to_paths(g, dag, 0, 3, True, 2)) == 2
    assert ngc.resolve_to_paths(g, dag, 3, 0) == []
    other = bgraph(4, [(0, 1, 1, 1)])
    foreign = "foreign PredecessorDAG: built for another graph"
    overlap = "node 2 is in both endpoint sets: source and destination must not overlap"
    with pytest.raises(ValueError, match=exact(foreign)):
        ngc.resolve_to_paths(other, dag, 0, 3)
    with pytest.raises(
        ValueError, match=exact("source: node 4 out of node range [0, 4)")
    ):
        ngc.resolve_to_paths(g, dag, 4, 3)
    with pytest.raises(ValueError, match=exact(overlap)):
        ngc.resolve_to_paths(g, dag, 2, 2)
    with pytest.raises(ValueError, match=exact("max_paths must be >= 1 (got 0)")):
        ngc.resolve_to_paths(g, dag, 0, 3, False, 0)
    with pytest.raises(
        ValueError, match=exact("max_paths must be >= 1 (got -1180591620717411303424)")
    ):
        ngc.resolve_to_paths(g, dag, 0, 3, False, -(2**70))
    # A count above INT64_MAX lists every walk (no error).
    assert len(ngc.resolve_to_paths(g, dag, 0, 3, True, 2**70)) == 4
    # Several defects: the DAG's graph, source, destination, source == destination, then max_paths.
    with pytest.raises(ValueError, match=exact(foreign)):
        ngc.resolve_to_paths(other, dag, 9, -9, False, 0)
    with pytest.raises(
        ValueError, match=exact("source: node 9 out of node range [0, 4)")
    ):
        ngc.resolve_to_paths(g, dag, 9, -9, False, 0)
    with pytest.raises(
        ValueError, match=exact("destination: node -9 out of node range [0, 4)")
    ):
        ngc.resolve_to_paths(g, dag, 0, -9, False, 0)
    with pytest.raises(ValueError, match=exact(overlap)):
        ngc.resolve_to_paths(g, dag, 2, 2, False, 0)
    # A search's DAG resolves the same way from a reverse search: the two cost-2 paths (the
    # cost-5 edge is on no shortest path).
    fw = ngc.spf(g, 0, ngc.SpfOptions(ngc.Direction.FORWARD)).dag
    rv = ngc.spf(g, 3, ngc.SpfOptions(ngc.Direction.REVERSE)).dag
    forward = ngc.resolve_to_paths(g, fw, 0, 3)
    assert [(p.cost(), p.hops()) for p in forward] == [
        (2, [[0, 1], [4]]),
        (2, [[2], [5]]),
    ]
    assert ngc.resolve_to_paths(g, rv, 0, 3) == forward


@needs_hooks
def test_python_result_lists_are_all_or_nothing():
    # Every Python object allocation of the call fails once: the call raises (never crashes)
    # and a clean retry gives the same list (resolve_to_paths: the four split paths of
    # test_resolve_to_paths).
    g = bgraph(4, CHAIN)
    opts = ngc.KspOptions(8)
    split = bgraph(
        4,
        [
            (0, 1, 1, 1),
            (0, 1, 1, 1),
            (0, 2, 1, 1),
            (0, 3, 1, 5),
            (1, 3, 1, 1),
            (2, 3, 1, 1),
        ],
    )
    dag = ngc.PredecessorDAG.from_csr(
        split, [0, 0, 2, 3, 6], [0, 0, 0, 0, 1, 2], [0, 1, 2, 3, 4, 5]
    )
    for fn, args in (
        (ngc.ksp, (g, 0, 3, opts)),
        (ngc.resolve_to_paths, (split, dag, 0, 3, True)),
    ):
        failures = sweep_python_allocation(
            lambda fn=fn, args=args: (fn, args, list), 200
        )
        assert failures >= 3, fn.__name__
    # Costs beyond the cached small ints, the free lists emptied before each attempt: the
    # result's cost ints are new objects, and their failures raise MemoryError too.
    wide = bgraph(
        4, [(0, 1, 1, 2**30), (1, 2, 1, 2**30), (2, 3, 1, 2**30), (0, 3, 1, 2**32)]
    )
    assert [c for c, _ in ngc.ksp(wide, 0, 3, opts)] == [3 * 2**30, 2**32]
    failures = sweep_python_allocation(lambda: (ngc.ksp, (wide, 0, 3, opts), list), 200)
    assert failures >= 5


def test_tie_order():
    # TRI (three parallel cost-1 links of capacities 4, 10, 100 and a cost-2 link of 50) with
    # availability 4, 8, 5, 50: the cost-1 links in each rule's order, then the cost-2 link.
    g = build(2, [(0, 1, 4, 1), (0, 1, 10, 1), (0, 1, 100, 1), (0, 1, 50, 2)])
    tie = ngc.RouteTieBreak
    for t, want in (
        (tie.FEWER_HOPS, [0, 1, 2, 3]),
        (tie.HIGHER_CAPACITY, [2, 1, 0, 3]),
        (tie.HIGHER_AVAILABLE_CAPACITY, [1, 2, 0, 3]),
        (tie.HIGHER_HEADROOM, [0, 1, 2, 3]),  # 1, 4/5, 1/20
    ):
        paths = ngc.ksp(g, 0, 1, ngc.KspOptions(4, tie_break=t, residual=[4, 8, 5, 50]))
        assert [p.hops()[0][0] for _, p in paths] == want, t
        assert [c for c, _ in paths] == [1, 1, 1, 2]


def test_stretch_examples():
    # KSP-stretch-boundary / -below: links of costs 5 (residual 0, so gated) and 6; the static
    # baseline is 5. 1.2 is exactly 6/5 and admits cost 6 at equality; 11/10 does not.
    # KSP-stretch-numpy / -numpy-equal: base 10^8 and np.float32(1.2), also exactly 6/5, reject
    # 120000001 and admit 120000000 (the float's binary value would not separate them).
    def paths(cost0, cost1, factor):
        g = build(2, [(0, 1, 0, cost0), (0, 1, 10, cost1)])
        options = ngc.KspOptions(
            2, cost_ceiling=ngc.CostStretch(factor), residual=[0, 10]
        )
        return [(c, [h[0] for h in p.hops()]) for c, p in ngc.ksp(g, 0, 1, options)]

    assert paths(5, 6, 1.2) == [(6, [1])]
    assert paths(5, 6, Fraction(11, 10)) == []
    assert paths(10**8, 120000001, np.float32(1.2)) == []
    assert paths(10**8, 120000000, np.float32(1.2)) == [(120000000, [1])]


def test_stretch_below_one():
    # KSP-stretch-below-one and KSP-ceiling-before-residual: a factor below 1 is refused where the
    # ceiling is checked, after k and before the residual; the lowest numerator too.
    g = build(2, [(0, 1, 10, 5)])
    for factor, shown in ((Fraction(1, 2), "1/2"), (-(2**63), f"{-(2**63)}/1")):
        text = f"cost_ceiling CostStretch factor must be >= 1 (got {shown})"
        for residual in (None, [10], [99, 99]):
            options = ngc.KspOptions(
                2, cost_ceiling=ngc.CostStretch(factor), residual=residual
            )
            with pytest.raises(ValueError, match=exact(text)):
                ngc.ksp(g, 0, 1, options)
    with pytest.raises(ValueError, match=exact("k must be >= 1 (got 0)")):
        ngc.ksp(g, 0, 1, ngc.KspOptions(0, cost_ceiling=ngc.CostStretch(0.5)))


def test_stretch_wide_products():
    # BIG(2^60): links of costs 2^60 and 2^60 + 1. With (2^60 + 1) / 2^60 both qualify (equality);
    # with (2^61 + 1) / 2^61 only the first: (2^60 + 1) * 2^61 = 2^121 + 2^61 > (2^61 + 1) * 2^60.
    c0 = 2**60
    g = build(2, [(0, 1, 1, c0), (0, 1, 1, c0 + 1)])
    for factor, want in (
        (Fraction(2**60 + 1, 2**60), [c0, c0 + 1]),
        (Fraction(2**61 + 1, 2**61), [c0]),
    ):
        options = ngc.KspOptions(2, cost_ceiling=ngc.CostStretch(factor))
        assert [c for c, _ in ngc.ksp(g, 0, 1, options)] == want, factor


def test_stretch_max_components():
    # One link of cost 2^62 - 1 (c = c0) and the factor (2^63 - 1) / (2^63 - 2): the products stay
    # below 2^125 and the path qualifies.
    g = build(2, [(0, 1, 1, 2**62 - 1)])
    options = ngc.KspOptions(
        1, cost_ceiling=ngc.CostStretch(Fraction(2**63 - 1, 2**63 - 2))
    )
    assert [c for c, _ in ngc.ksp(g, 0, 1, options)] == [2**62 - 1]
