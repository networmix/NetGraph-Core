"""The public search through the binding (``spf``, APPROACH 6): results, the DAG's canonical
form, the validation order on unnarrowed values, the enumerations, and the all-or-nothing
construction of the result objects."""

from __future__ import annotations

import _netgraph_core as ext
import numpy as np
import pytest
from private_hooks import needs_hooks

import netgraph_core as ngc

DET = ngc.EdgeTieBreak.DETERMINISTIC
PREFER = ngc.EdgeTieBreak.PREFER_HIGHER_PRIORITY
CAP = ngc.EdgePriority.CAPACITY
RES = ngc.EdgePriority.RESIDUAL
FWD = ngc.Direction.FORWARD
REV = ngc.Direction.REVERSE


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
    return sorted(int(e) for e, p in enumerate(g.input_index()) if p in positions)


def via(dag):
    return sorted(int(e) for e in dag.via_edges())


def opts(direction=FWD, selection=None, **kw):
    return ngc.SpfOptions(direction, selection or ngc.PathSelection(DET, CAP), **kw)


ANYCAST = [(0, 1, 1, 1), (1, 2, 1, 1), (0, 3, 1, 2), (3, 2, 1, 5), (1, 4, 1, 5)]


def test_full_search_and_anycast():
    g = build(5, ANYCAST)
    full = ngc.spf(g, 0, opts())
    assert full.costs.dtype == np.int64
    assert list(full.costs) == [0, 1, 2, 2, 6]
    assert full.target_cost is None
    assert not full.costs.flags.writeable
    any_ = ngc.spf(g, [0], opts(targets=[2, 4]))
    assert any_.costs is None
    assert any_.target_cost == 2
    assert via(any_.dag) == ids(g, {0, 1})
    rev = ngc.spf(g, [2, 4], opts(REV, targets=0))
    assert rev.dag == any_.dag
    assert rev.target_cost == 2
    none = ngc.spf(g, 2, opts(targets=0))
    assert none.target_cost is None
    assert list(none.dag.via_edges()) == []
    unreachable = ngc.spf(g, 2, opts())
    assert list(unreachable.costs) == [
        ngc.UNREACHABLE,
        ngc.UNREACHABLE,
        0,
        ngc.UNREACHABLE,
        ngc.UNREACHABLE,
    ]
    # Reading a field twice returns the same object: the result allocates nothing on access.
    assert full.dag is full.dag and full.costs is full.costs
    assert "target_cost=2" in repr(any_)


def test_reverse_full_search_is_forward_oriented():
    g = build(4, [(0, 1, 1, 1), (0, 1, 1, 1), (0, 2, 1, 2), (1, 3, 1, 1), (2, 3, 1, 1)])
    r = ngc.spf(g, 3, opts(REV))
    assert list(r.costs) == [2, 1, 1, 0]
    assert via(r.dag) == ids(g, {0, 1, 3, 4})
    off = list(r.dag.parent_offsets())
    assert off[0] == 0 and off[-1] == len(r.dag.via_edges())


def test_staged_single_path_review13_d3_both_directions():
    # s -> v (1, prio 1), v -> t (3, prio 1), v -> a -> b -> t (1 + 1 + 1, prio 2): the
    # bottleneck of every walk is 1, so the two-hop walk wins under both tie modes.
    g = build(5, [(0, 1, 1, 1), (1, 4, 1, 3), (1, 2, 2, 1), (2, 3, 2, 1), (3, 4, 2, 1)])
    for tie in (DET, PREFER):
        s = ngc.PathSelection(tie, CAP, multipath=False)
        f = ngc.spf(g, 0, opts(selection=s, targets=4))
        assert via(f.dag) == ids(g, {0, 1})
        assert f.target_cost == 4
        r = ngc.spf(g, 4, opts(REV, selection=s, targets=0))
        assert r.dag == f.dag


def test_residual_gating_and_priority():
    # Capacities 5; residuals by input [5, 4, 1, 3, 2]: 0->1->3 has bottleneck 4.
    g = build(4, [(0, 1, 5, 1), (1, 3, 5, 1), (0, 2, 5, 1), (2, 3, 5, 1), (0, 3, 5, 2)])
    by_input = [5, 4, 1, 3, 2]
    residual = np.array([by_input[int(p)] for p in g.input_index()], dtype=np.int64)
    s = ngc.PathSelection(PREFER, RES, multipath=False)
    r = ngc.spf(g, 0, opts(selection=s, targets=3, residual=residual))
    assert via(r.dag) == ids(g, {0, 1})
    # Gated at 4 units: only 0->1 (5) and 1->3 (4) are routable on that walk; the direct edge
    # (residual 2) is gated out, so the cost-2 walk is the two-hop one through node 1.
    gated = ngc.spf(g, 0, opts(targets=3, residual=residual, min_residual=4))
    assert gated.target_cost == 2
    assert via(gated.dag) == ids(g, {0, 1})
    # min_residual without a residual gates nothing.
    assert ngc.spf(g, 0, opts(min_residual=7)).costs[3] == 2


def test_validation_order_on_unnarrowed_values():
    g = build(3, [(0, 1, 2, 1), (1, 2, 2, 1)])
    other = build(3, [(0, 1, 2, 1)])
    with pytest.raises(ValueError, match=r"^roots: the node set is empty$"):
        ngc.spf(g, [], opts())
    with pytest.raises(ValueError, match=r"^roots: node 3 out of node range \[0, 3\)$"):
        ngc.spf(g, 3, opts())
    big = 2**70
    with pytest.raises(ValueError, match=rf"^roots: node {big} out of node range"):
        ngc.spf(g, [0, big], opts())
    with pytest.raises(ValueError, match=r"^targets: the node set is empty$"):
        ngc.spf(g, 0, opts(targets=[]))
    with pytest.raises(ValueError, match=r"^targets: node -1 out of node range"):
        ngc.spf(g, 0, opts(targets=[-1]))
    with pytest.raises(
        ValueError, match=r"^node 1 is in both endpoint sets: roots and targets"
    ):
        ngc.spf(g, [0, 1], opts(targets=[1]))
    with pytest.raises(ValueError, match="^multipath = false requires targets"):
        ngc.spf(g, 0, opts(selection=ngc.PathSelection(DET, CAP, multipath=False)))
    with pytest.raises(ValueError, match="^residual has length 1, expected 2"):
        ngc.spf(g, 0, opts(residual=[1]))
    with pytest.raises(ValueError, match=r"^residual\[1\] = 3 is outside"):
        ngc.spf(g, 0, opts(residual=[1, 3]))
    with pytest.raises(ValueError, match=r"^residual\[0\] = -" + str(big)):
        ngc.spf(g, 0, opts(residual=[-big, 1]))
    with pytest.raises(
        ValueError, match=r"^min_residual must be at least 1 unit \(got 0\)$"
    ):
        ngc.spf(g, 0, opts(min_residual=0))
    with pytest.raises(ValueError, match=r"min_residual must be an integer"):
        ngc.spf(g, 0, opts(min_residual=1.5))
    with pytest.raises(ValueError, match="^foreign Eligibility"):
        ngc.spf(g, 0, opts(elig=other.eligibility([], [0])))
    # Several defects: the contract's order (roots, targets, overlap, single path, residual,
    # min_residual, Eligibility).
    with pytest.raises(ValueError, match="^roots: the node set is empty$"):
        ngc.spf(g, [], opts(targets=[], residual=[1], min_residual=0))
    with pytest.raises(ValueError, match="both endpoint sets"):
        ngc.spf(g, 0, opts(targets=0, residual=[]))
    with pytest.raises(ValueError, match="requires targets"):
        ngc.spf(
            g,
            0,
            opts(selection=ngc.PathSelection(DET, CAP, multipath=False), residual=[]),
        )
    with pytest.raises(ValueError, match="^residual"):
        ngc.spf(g, 0, opts(residual=[1, 9], min_residual=0))
    with pytest.raises(ValueError, match="^min_residual"):
        ngc.spf(g, 0, opts(min_residual=0, elig=other.eligibility([], [0])))
    # A huge min_residual gates every residual edge (no error).
    assert list(ngc.spf(g, 0, opts(residual=[2, 2], min_residual=big)).costs)[1:] == [
        ngc.UNREACHABLE,
        ngc.UNREACHABLE,
    ]


def test_options_types_and_enumerations():
    g = build(2, [(0, 1, 1, 1)])
    o = opts(targets=[1], min_residual=3)
    assert o.direction == FWD and o.targets == [1] and o.min_residual == 3
    assert o.residual is None and o.elig is None
    with pytest.raises(ValueError, match="^direction: 7 is not a value of Direction$"):
        ngc.SpfOptions(ngc.Direction(7), ngc.PathSelection(DET, CAP))
    with pytest.raises(TypeError):
        ngc.SpfOptions(FWD, ngc.PathSelection(DET, CAP), fanout_edges=[0])
    with pytest.raises(TypeError, match=r"^spf\(\): incompatible function arguments"):
        ngc.spf(g, 0, opts(), fanout_edges=[0])
    with pytest.raises(TypeError, match="SpfOptions expected"):
        ngc.spf(g, 0, ngc.PathSelection(DET, CAP))
    with pytest.raises(TypeError):
        ngc.SpfOptions(FWD)
    # Roots as a NumPy array, with duplicates.
    r = ngc.spf(g, np.array([0, 0], dtype=np.int32), opts())
    assert list(r.costs) == [0, 1]


@needs_hooks
def test_python_result_construction_is_all_or_nothing():
    # Every Python object allocation of the call fails once: the call raises (never crashes)
    # and a clean retry gives the same result.
    g = build(5, ANYCAST)
    o = opts(targets=[2, 4])
    want = ngc.spf(g, 0, o)
    seen = None
    k = 0
    failures = 0
    while True:
        hit, seen, result, error = ext._fail_python_allocation(ngc.spf, (g, 0, o), k)
        if not hit:
            assert result is not None and result.dag == want.dag
            break
        failures += 1
        assert result is None and isinstance(error, (MemoryError, RuntimeError)), error
        again = ngc.spf(g, 0, o)
        assert again.dag == want.dag and again.target_cost == want.target_cost
        k += 1
    assert failures >= 3 and seen > 0
    full = opts()
    want_full = ngc.spf(g, 0, full)
    for k in range(40):
        hit, _seen, result, error = ext._fail_python_allocation(
            ngc.spf, (g, 0, full), k
        )
        if not hit:
            assert list(result.costs) == list(want_full.costs)
            break
        assert result is None and error is not None
