"""KSP and resolve_to_paths through the binding (APPROACH 5, 6, rule 10): the order, Expand and
Collapse, the ceilings, k = 1 against the single path of spf, the validation order on
unnarrowed values, and the all-or-nothing construction of the result lists."""

from __future__ import annotations

import _netgraph_core as ext
import pytest
from private_hooks import needs_hooks

import netgraph_core as ngc

DET = ngc.EdgeTieBreak.DETERMINISTIC
PREFER = ngc.EdgeTieBreak.PREFER_HIGHER_PRIORITY
CAP = ngc.EdgePriority.CAPACITY


def bgraph(n, rows):
    """Rows (u, v, capacity, cost) in B' order (src, cost, dst): ids are input positions."""
    rows = sorted(rows, key=lambda r: (r[0], r[3], r[1]))
    return ngc.Graph(
        n,
        [r[0] for r in rows],
        [r[1] for r in rows],
        [r[2] for r in rows],
        [r[3] for r in rows],
        1.0,
    )


def ksp(g, s, t, k, *, multi_edge=True, tie=DET, **kw):
    sel = ngc.PathSelection(tie, CAP, multi_edge=multi_edge)
    return ngc.ksp(g, s, t, ngc.KspOptions(k, sel, **kw))


CHAIN = [(0, 1, 1, 1), (1, 2, 1, 1), (2, 3, 1, 1), (0, 2, 1, 3), (1, 3, 1, 3)]


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
        single = ngc.PathSelection(tie, CAP, multipath=False)
        for direction, roots, targets in (
            (ngc.Direction.FORWARD, 0, 5),
            (ngc.Direction.REVERSE, 5, 0),
        ):
            r = ngc.spf(
                three, roots, ngc.SpfOptions(direction, single, targets=targets)
            )
            assert r.dag == first.to_dag()


def test_ksp_validation_order():
    g = bgraph(3, [(0, 1, 1, 1), (1, 2, 1, 1)])
    other = bgraph(3, [(0, 1, 1, 1)])
    with pytest.raises(ValueError, match=r"^src: node 3 out of node range \[0, 3\)$"):
        ksp(g, 3, 1, 1)
    with pytest.raises(ValueError, match=r"^dst: node -1 out of node range"):
        ksp(g, 0, -1, 1)
    with pytest.raises(ValueError, match=r"^dst: node 1180591620717411303424 out of"):
        ksp(g, 0, 2**70, 1)
    with pytest.raises(
        ValueError, match="^node 1 is in both endpoint sets: src and dst"
    ):
        ksp(g, 1, 1, 1)
    with pytest.raises(ValueError, match=r"^k must be >= 1 \(got 0\)$"):
        ksp(g, 0, 2, 0)
    with pytest.raises(ValueError, match="k must be an integer"):
        ksp(g, 0, 2, 1.5)
    with pytest.raises(ValueError, match="^max_path_cost Absolute cost must be >= 0"):
        ksp(g, 0, 2, 1, max_path_cost=ngc.Absolute(-1))
    with pytest.raises(ValueError, match="^max_path_cost Stretch factor"):
        ksp(g, 0, 2, 1, max_path_cost=ngc.Stretch(0.5))
    with pytest.raises(ValueError, match="^max_path_cost Stretch factor"):
        ksp(g, 0, 2, 1, max_path_cost=ngc.Stretch(float("nan")))
    with pytest.raises(TypeError, match="max_path_cost must be Absolute, Stretch"):
        ksp(g, 0, 2, 1, max_path_cost=3)
    with pytest.raises(ValueError, match=r"^residual\[1\] = 2 is outside"):
        ksp(g, 0, 2, 1, residual=[1, 2])
    with pytest.raises(ValueError, match="^foreign Eligibility"):
        ksp(g, 0, 2, 1, elig=other.eligibility([], [0]))
    # The order: the endpoints, k, the ceiling, the residual.
    with pytest.raises(ValueError, match="both endpoint sets"):
        ksp(g, 2, 2, 0, max_path_cost=ngc.Absolute(-1), residual=[9])
    with pytest.raises(ValueError, match="^k must be"):
        ksp(g, 0, 2, 0, max_path_cost=ngc.Absolute(-1), residual=[9])
    with pytest.raises(ValueError, match="^max_path_cost"):
        ksp(g, 0, 2, 1, max_path_cost=ngc.Absolute(-1), residual=[9])
    # A huge k lists every path (no error).
    assert len(ksp(g, 0, 2, 2**70)) == 1
    with pytest.raises(TypeError, match="KspOptions expected"):
        ngc.ksp(g, 0, 2, ngc.PathSelection(DET, CAP))


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
    with pytest.raises(ValueError, match="^foreign PredDAG"):
        ngc.resolve_to_paths(other, dag, 0, 3)
    with pytest.raises(ValueError, match=r"^src: node 4 out of node range \[0, 4\)$"):
        ngc.resolve_to_paths(g, dag, 4, 3)
    with pytest.raises(
        ValueError, match="^node 2 is in both endpoint sets: src and dst"
    ):
        ngc.resolve_to_paths(g, dag, 2, 2)
    with pytest.raises(ValueError, match=r"^max_paths must be >= 1 \(got 0\)$"):
        ngc.resolve_to_paths(g, dag, 0, 3, False, 0)
    # A search's DAG resolves the same way from a reverse search.
    sel = ngc.PathSelection(DET, CAP)
    fw = ngc.spf(g, 0, ngc.SpfOptions(ngc.Direction.FORWARD, sel)).dag
    rv = ngc.spf(g, 3, ngc.SpfOptions(ngc.Direction.REVERSE, sel)).dag
    assert ngc.resolve_to_paths(g, fw, 0, 3) == ngc.resolve_to_paths(g, rv, 0, 3)


@needs_hooks
def test_python_result_lists_are_all_or_nothing():
    g = bgraph(4, CHAIN)
    opts = ngc.KspOptions(8, ngc.PathSelection(DET, CAP))
    want = ngc.ksp(g, 0, 3, opts)
    dag = ngc.spf(
        g, 0, ngc.SpfOptions(ngc.Direction.FORWARD, ngc.PathSelection(DET, CAP))
    ).dag
    want_paths = ngc.resolve_to_paths(g, dag, 0, 2, True)
    for fn, args, check in (
        (ngc.ksp, (g, 0, 3, opts), lambda r: r == want),
        (ngc.resolve_to_paths, (g, dag, 0, 2, True), lambda r: r == want_paths),
    ):
        failures = 0
        for k in range(200):
            hit, _seen, result, error = ext._fail_python_allocation(fn, args, k)
            if not hit:
                assert check(result)
                break
            failures += 1
            assert result is None and isinstance(error, (MemoryError, RuntimeError)), (
                error
            )
            assert check(fn(*args))
        assert failures >= 3
