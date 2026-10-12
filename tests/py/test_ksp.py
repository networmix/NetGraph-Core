"""KSP and resolve_to_paths through the binding (docs/semantics.md#k-shortest-paths): the
order, Expand and Collapse, the ceilings, k = 1 against the single path of spf, the validation
order on unnarrowed values, and the all-or-nothing construction of the result lists."""

from __future__ import annotations

import pytest
from helpers import CHAIN, DET, PREFER, build, exact, sel, sweep_python_allocation
from private_hooks import needs_hooks

import netgraph_core as ngc


def bgraph(n, rows):
    """A graph from rows (u, v, capacity, cost), sorted by (src, cost, dst) first: that is the
    EdgeId order, so each EdgeId is the row's position in the sorted list."""
    return build(n, sorted(rows, key=lambda r: (r[0], r[3], r[1])))


def ksp(g, s, t, k, *, multi_edge=True, tie=DET, **kw):
    return ngc.ksp(g, s, t, ngc.KspOptions(k, sel(tie, me=multi_edge), **kw))


def test_order_and_ceilings():
    g = bgraph(4, CHAIN)
    got = ksp(g, 0, 3, 8)
    assert [c for c, _ in got] == [3, 4, 4]
    assert all(isinstance(p, ngc.Path) and p.cost() == c for c, p in got)
    assert got[0][1].hops() == [[0], [2], [4]]
    assert [c for c, _ in ksp(g, 0, 3, 8, max_path_cost=ngc.Stretch(1.0))] == [3]
    assert [c for c, _ in ksp(g, 0, 3, 8, max_path_cost=ngc.Stretch(1.5))] == [3, 4, 4]
    assert [c for c, _ in ksp(g, 0, 3, 8, max_path_cost=ngc.Absolute(3))] == [3]
    assert ksp(g, 0, 3, 8, max_path_cost=ngc.Absolute(0)) == []
    assert len(ksp(g, 0, 3, 2)) == 2
    # The Stretch baseline is static: gating the middle edge drops the cost-3 path only.
    residual = [1, 1, 0, 1, 1]
    assert ksp(g, 0, 3, 8, max_path_cost=ngc.Stretch(1.0), residual=residual) == []
    assert [c for c, _ in ksp(g, 0, 3, 8, residual=residual)] == [4, 4]
    excluded = g.eligibility([0], [])
    assert ksp(g, 0, 3, 8, elig=excluded) == []


def test_expand_collapse_and_the_single_path():
    pr = bgraph(3, [(0, 1, 1, 1), (0, 1, 1, 1), (1, 2, 1, 1), (1, 2, 1, 1)])
    assert len(ksp(pr, 0, 2, 8)) == 4
    assert len(ksp(pr, 0, 2, 8, multi_edge=False)) == 1
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
    for tie in (DET, PREFER):
        first = ksp(three, 0, 5, 1, tie=tie)[0][1]
        single = sel(tie, mp=False)
        for direction, roots, targets in (
            (ngc.Direction.FORWARD, 0, 5),
            (ngc.Direction.REVERSE, 5, 0),
        ):
            r = ngc.spf(
                three, roots, ngc.SpfOptions(direction, single, targets=targets)
            )
            assert r.dag == first.to_dag()
    # multipath is ignored: the same selection with multipath = False lists the same paths.
    for g, t in ((pr, 2), (three, 5)):
        for tie, me in ((DET, True), (PREFER, True), (DET, False)):
            want = ngc.ksp(g, 0, t, ngc.KspOptions(8, sel(tie, me=me)))
            assert want
            assert (
                ngc.ksp(g, 0, t, ngc.KspOptions(8, sel(tie, mp=False, me=me))) == want
            )


def test_ksp_validation_order():
    # The residual's own checks are the shared ones of test_edge_amounts.py.
    g = bgraph(3, [(0, 1, 1, 1), (1, 2, 1, 1)])
    other = bgraph(3, [(0, 1, 1, 1)])
    with pytest.raises(ValueError, match=exact("src: node 3 out of node range [0, 3)")):
        ksp(g, 3, 1, 1)
    with pytest.raises(
        ValueError, match=exact("dst: node -1 out of node range [0, 3)")
    ):
        ksp(g, 0, -1, 1)
    with pytest.raises(
        ValueError,
        match=exact("dst: node 1180591620717411303424 out of node range [0, 3)"),
    ):
        ksp(g, 0, 2**70, 1)
    with pytest.raises(
        ValueError,
        match=exact("node 1 is in both endpoint sets: src and dst must not overlap"),
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
    absolute = "max_path_cost Absolute cost must be >= 0 (got -1)"
    stretch = "max_path_cost Stretch factor must be finite and in [1, inf) (got {})"
    with pytest.raises(ValueError, match=exact(absolute)):
        ksp(g, 0, 2, 1, max_path_cost=ngc.Absolute(-1))
    with pytest.raises(ValueError, match=exact(stretch.format("0.500000"))):
        ksp(g, 0, 2, 1, max_path_cost=ngc.Stretch(0.5))
    with pytest.raises(ValueError, match=exact(stretch.format("nan"))):
        ksp(g, 0, 2, 1, max_path_cost=ngc.Stretch(float("nan")))
    with pytest.raises(
        TypeError,
        match=exact("max_path_cost must be Absolute, Stretch or None (got int)"),
    ):
        ksp(g, 0, 2, 1, max_path_cost=3)
    with pytest.raises(
        ValueError, match=exact("foreign Eligibility: built for another graph")
    ):
        ksp(g, 0, 2, 1, elig=other.eligibility([], [0]))
    # The order: the endpoints (src, dst, src == dst), k, the ceiling, the residual.
    with pytest.raises(ValueError, match=exact("src: node 3 out of node range [0, 3)")):
        ksp(g, 3, -1, 0, max_path_cost=ngc.Absolute(-1), residual=[9])
    with pytest.raises(
        ValueError, match=exact("dst: node -1 out of node range [0, 3)")
    ):
        ksp(g, 0, -1, 0, max_path_cost=ngc.Absolute(-1), residual=[9])
    with pytest.raises(
        ValueError,
        match=exact("node 2 is in both endpoint sets: src and dst must not overlap"),
    ):
        ksp(g, 2, 2, 0, max_path_cost=ngc.Absolute(-1), residual=[9])
    with pytest.raises(ValueError, match=exact("k must be >= 1 (got 0)")):
        ksp(g, 0, 2, 0, max_path_cost=ngc.Absolute(-1), residual=[9])
    with pytest.raises(ValueError, match=exact(absolute)):
        ksp(g, 0, 2, 1, max_path_cost=ngc.Absolute(-1), residual=[9])
    # A huge k lists every path (no error).
    assert len(ksp(g, 0, 2, 2**70)) == 1
    with pytest.raises(
        TypeError, match=exact("KspOptions expected (got PathSelection)")
    ):
        ngc.ksp(g, 0, 2, sel())
    with pytest.raises(TypeError):
        ngc.KspOptions(1)  # k and selection are explicit
    # KspOptions keeps its arguments as given (the same objects) and validates none of them:
    # ksp does, at the call, in the order above.
    k, ceiling, residual = 2**70, ngc.Absolute(-1), [9]
    o = ngc.KspOptions(k, sel(mp=False), max_path_cost=ceiling, residual=residual)
    assert o.k is k and o.max_path_cost is ceiling and o.residual is residual
    assert o.selection == sel(mp=False) and o.elig is None
    with pytest.raises(
        ValueError, match=exact("max_path_cost Absolute cost must be >= 0 (got -1)")
    ):
        ngc.ksp(g, 0, 2, o)
    with pytest.raises(ValueError, match=exact("k must be >= 1 (got 0)")):
        ngc.ksp(g, 0, 2, ngc.KspOptions(0, sel(), max_path_cost=ceiling))


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
    dag = ngc.PredDAG.from_csr(
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
    foreign = "foreign PredDAG: built for another graph"
    overlap = "node 2 is in both endpoint sets: src and dst must not overlap"
    with pytest.raises(ValueError, match=exact(foreign)):
        ngc.resolve_to_paths(other, dag, 0, 3)
    with pytest.raises(ValueError, match=exact("src: node 4 out of node range [0, 4)")):
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
    # Several defects: the DAG's graph, src, dst, src == dst, then max_paths.
    with pytest.raises(ValueError, match=exact(foreign)):
        ngc.resolve_to_paths(other, dag, 9, -9, False, 0)
    with pytest.raises(ValueError, match=exact("src: node 9 out of node range [0, 4)")):
        ngc.resolve_to_paths(g, dag, 9, -9, False, 0)
    with pytest.raises(
        ValueError, match=exact("dst: node -9 out of node range [0, 4)")
    ):
        ngc.resolve_to_paths(g, dag, 0, -9, False, 0)
    with pytest.raises(ValueError, match=exact(overlap)):
        ngc.resolve_to_paths(g, dag, 2, 2, False, 0)
    # A search's DAG resolves the same way from a reverse search: the two cost-2 paths (the
    # cost-5 edge is on no shortest path).
    fw = ngc.spf(g, 0, ngc.SpfOptions(ngc.Direction.FORWARD, sel())).dag
    rv = ngc.spf(g, 3, ngc.SpfOptions(ngc.Direction.REVERSE, sel())).dag
    forward = ngc.resolve_to_paths(g, fw, 0, 3)
    assert [(p.cost(), p.hops()) for p in forward] == [
        (2, [[0, 1], [4]]),
        (2, [[2], [5]]),
    ]
    assert ngc.resolve_to_paths(g, rv, 0, 3) == forward


@needs_hooks
def test_python_result_lists_are_all_or_nothing():
    # Every Python object allocation of the call fails once: the call raises (never crashes)
    # and a clean retry gives the same list.
    g = bgraph(4, CHAIN)
    opts = ngc.KspOptions(8, sel())
    dag = ngc.spf(g, 0, ngc.SpfOptions(ngc.Direction.FORWARD, sel())).dag
    for fn, args in (
        (ngc.ksp, (g, 0, 3, opts)),
        (ngc.resolve_to_paths, (g, dag, 0, 2, True)),
    ):
        failures = sweep_python_allocation(
            lambda fn=fn, args=args: (fn, args, list), 200
        )
        assert failures >= 3, fn.__name__
