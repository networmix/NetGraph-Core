"""The public search through the binding (``spf``, docs/semantics.md#4-search): results, the
DAG's canonical form, the validation order on unnarrowed values, the enumerations, and the
all-or-nothing construction of the result objects."""

from __future__ import annotations

import itertools
import traceback

import numpy as np
import pytest
from helpers import (
    DET,
    L62,
    PREFER,
    RES,
    build,
    exact,
    ids,
    plain,
    sel,
    sweep_python_allocation,
    to_ids,
)
from private_hooks import needs_hooks

import netgraph_core as ngc

FWD = ngc.Direction.FORWARD
REV = ngc.Direction.REVERSE


def via(dag):
    return sorted(int(e) for e in dag.via_edges())


def opts(direction=FWD, selection=None, **kw):
    return ngc.SpfOptions(direction, selection or sel(), **kw)


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


def test_staged_single_path_with_a_shared_bottleneck_takes_the_fewest_hops():
    # s -> v (1, prio 1), v -> t (3, prio 1), v -> a -> b -> t (1 + 1 + 1, prio 2): the
    # bottleneck of every walk is 1, so the two-hop walk wins under both tie modes.
    g = build(5, [(0, 1, 1, 1), (1, 4, 1, 3), (1, 2, 2, 1), (2, 3, 2, 1), (3, 4, 2, 1)])
    for tie in (DET, PREFER):
        s = sel(tie, mp=False)
        f = ngc.spf(g, 0, opts(selection=s, targets=4))
        assert via(f.dag) == ids(g, {0, 1})
        assert f.target_cost == 4
        r = ngc.spf(g, 4, opts(REV, selection=s, targets=0))
        assert r.dag == f.dag


def test_residual_gating_and_priority():
    # Capacities 5; residuals by input [5, 4, 1, 3, 2]: 0->1->3 has bottleneck 4.
    g = build(4, [(0, 1, 5, 1), (1, 3, 5, 1), (0, 2, 5, 1), (2, 3, 5, 1), (0, 3, 5, 2)])
    residual = np.array(to_ids(g, [5, 4, 1, 3, 2]), dtype=np.int64)
    s = sel(PREFER, RES, mp=False)
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
    # The residual's own checks are the shared ones of test_edge_amounts.py.
    g = build(3, [(0, 1, 2, 1), (1, 2, 2, 1)])
    other = build(3, [(0, 1, 2, 1)])
    with pytest.raises(ValueError, match=exact("roots: the node set is empty")):
        ngc.spf(g, [], opts())
    with pytest.raises(
        ValueError, match=exact("roots: node 3 out of node range [0, 3)")
    ):
        ngc.spf(g, 3, opts())
    big = 2**70
    with pytest.raises(ValueError, match=rf"^roots: node {big} out of node range"):
        ngc.spf(g, [0, big], opts())
    with pytest.raises(ValueError, match=exact("targets: the node set is empty")):
        ngc.spf(g, 0, opts(targets=[]))
    with pytest.raises(ValueError, match=r"^targets: node -1 out of node range"):
        ngc.spf(g, 0, opts(targets=[-1]))
    with pytest.raises(
        ValueError, match=r"^node 1 is in both endpoint sets: roots and targets"
    ):
        ngc.spf(g, [0, 1], opts(targets=[1]))
    single = (
        "multipath = false requires targets (the single path joins a root to a target)"
    )
    with pytest.raises(ValueError, match=exact(single)):
        ngc.spf(g, 0, opts(selection=sel(mp=False)))
    with pytest.raises(
        ValueError, match=exact("min_residual must be at least 1 unit (got 0)")
    ):
        ngc.spf(g, 0, opts(min_residual=0))
    with pytest.raises(
        ValueError, match=exact(f"min_residual must be at least 1 unit (got {-big})")
    ):
        ngc.spf(g, 0, opts(min_residual=-big))
    with pytest.raises(
        ValueError, match=exact("min_residual must be an integer (got float)")
    ):
        ngc.spf(g, 0, opts(min_residual=1.5))
    with pytest.raises(ValueError, match="^foreign Eligibility"):
        ngc.spf(g, 0, opts(elig=other.eligibility([], [0])))
    # Several defects: the first in the documented order (roots, targets, overlap, single
    # path, residual, min_residual, Eligibility).
    with pytest.raises(ValueError, match=exact("roots: the node set is empty")):
        ngc.spf(g, [], opts(targets=[], residual=[1], min_residual=0))
    with pytest.raises(ValueError, match=rf"^targets: node {big} out of node range"):
        ngc.spf(g, 0, opts(targets=[0, big], residual=[]))
    with pytest.raises(ValueError, match="both endpoint sets"):
        ngc.spf(g, 0, opts(targets=0, residual=[]))
    with pytest.raises(ValueError, match=exact(single)):
        ngc.spf(
            g,
            0,
            opts(selection=sel(mp=False), residual=[]),
        )
    with pytest.raises(ValueError, match="^residual"):
        ngc.spf(g, 0, opts(residual=[1, 9], min_residual=0))
    with pytest.raises(ValueError, match="^min_residual"):
        ngc.spf(g, 0, opts(min_residual=0, elig=other.eligibility([], [0])))
    # A min_residual above INT64_MAX is no error and acts as INT64_MAX does: it gates every
    # residual edge, even one whose residual is the largest unit value 2^62 - 1, which
    # min_residual = 2^62 - 1 itself admits.
    top = build(2, [(0, 1, L62 - 1, 1)])
    for m in (L62 - 1, 2**63 - 1, big):
        r = ngc.spf(top, 0, opts(residual=[L62 - 1], min_residual=m))
        assert list(r.costs) == [0, 1 if m == L62 - 1 else ngc.UNREACHABLE], m
    huge = ngc.spf(g, 0, opts(residual=[2, 2], min_residual=big))
    assert list(huge.costs)[1:] == [ngc.UNREACHABLE, ngc.UNREACHABLE]
    assert plain(huge) == plain(
        ngc.spf(g, 0, opts(residual=[2, 2], min_residual=2**63 - 1))
    )


def test_options_types_and_enumerations():
    g = build(2, [(0, 1, 1, 1)])
    o = opts(targets=[1], min_residual=3)
    assert o.direction == FWD and o.targets == [1] and o.min_residual == 3
    assert o.residual is None and o.elig is None
    with pytest.raises(
        ValueError, match=exact("direction: 7 is not a value of Direction")
    ):
        ngc.SpfOptions(ngc.Direction(7), sel())
    with pytest.raises(TypeError):
        ngc.SpfOptions(FWD, sel(), fanout_edges=[0])
    with pytest.raises(TypeError, match=r"^spf\(\): incompatible function arguments"):
        ngc.spf(g, 0, opts(), fanout_edges=[0])
    with pytest.raises(
        TypeError, match=exact("SpfOptions expected (got PathSelection)")
    ):
        ngc.spf(g, 0, sel())
    with pytest.raises(TypeError):
        ngc.SpfOptions(FWD)
    # Roots as a NumPy array, with duplicates.
    r = ngc.spf(g, np.array([0, 0], dtype=np.int32), opts())
    assert list(r.costs) == [0, 1]
    # SpfOptions keeps targets, residual and min_residual as given (the same objects) and
    # validates only the direction: spf validates the rest at the call, in its order.
    assert opts().min_residual == 1
    targets, residual, low = [1, 1], [9], -(2**70)
    o = opts(targets=targets, residual=residual, min_residual=low)
    assert o.targets is targets and o.residual is residual and o.min_residual is low
    assert o.selection == sel()
    with pytest.raises(ValueError, match=r"^residual\[0\] = 9 is outside"):
        ngc.spf(g, 0, o)


@needs_hooks
def test_python_result_construction_is_all_or_nothing():
    # Every Python object allocation of the call fails once: the call raises (never crashes)
    # and a clean retry gives the same result, for a search with targets and a full search.
    g = build(5, ANYCAST)
    for o in (opts(targets=[2, 4]), opts()):
        failures = sweep_python_allocation(lambda o=o: (ngc.spf, (g, 0, o), list), 200)
        assert failures >= 3, o.targets


@needs_hooks
def test_the_python_allocation_sweep_rejects_a_wrong_call():
    # sweep_python_allocation rejects a wrong call, each on its own clause: one that absorbs the
    # injected failure and returns a value, one that raises another exception in its place, one
    # that changes the observable state before the failure reaches it; one whose completing call
    # raises after its result, returns another value or makes fewer allocations than failed
    # before it (a prebuilt result), one whose retry after a failure returns another value; and
    # one that has not completed within the attempt limit.
    g = build(5, ANYCAST)

    def absorbs(o):
        try:
            return ngc.spf(g, 0, o)
        except MemoryError:
            return "absorbed"

    def replaces(o):
        try:
            return ngc.spf(g, 0, o)
        except MemoryError:
            raise LookupError("replaced") from None

    def changing():
        log = []

        def changes(o):
            log.append(o)
            return ngc.spf(g, 0, o)

        return changes, (opts(),), lambda: len(log)

    def root(r):
        def search(o):
            return ngc.spf(g, r, o)

        return search

    root0, root1 = root(0), root(1)
    full = ngc.spf(g, 0, opts())

    def prebuilt(o):
        return full

    def raising():
        # The completing call raises after its result. It raises in the frame the sweep calls:
        # CPython keeps an exception whose traceback entry for an outer frame failed to allocate.
        made = itertools.count()

        def make():
            m = next(made)
            calls = itertools.count()

            def raises(o):
                first = not next(calls)
                result = ngc.spf(g, 0, o)
                if m and first:
                    raise LookupError("after its result")
                return result

            return raises, (opts(),), list

        return make

    def by_turn(turns):
        # The c-th call within the m-th make() runs turns(m, c); m = 0 is the sweep's reference.
        made = itertools.count()

        def make():
            m = next(made)
            calls = itertools.count()
            return (lambda o: turns(m, next(calls))(o)), (opts(),), list

        return make

    for make, clause in (
        (lambda: (absorbs, (opts(),), list), "assert value is None, (where, value)"),
        (
            lambda: (replaces, (opts(),), list),
            "assert isinstance(error, (MemoryError, RuntimeError)), (where, repr(error))",
        ),
        (changing, "assert snapshot() == before, where"),
        (raising(), "assert error is None, (where, repr(error))"),
        (
            by_turn(lambda m, c: root1 if m and not c else root0),
            "assert plain(value) == want, where",
        ),
        (
            by_turn(lambda m, c: prebuilt if m >= 3 and not c else root0),
            "assert seen == failures, (where, seen)  # every allocation failed once",
        ),
        (
            by_turn(lambda m, c: root1 if m and c == 1 else root0),
            "assert plain(fn(*args)) == want, where",
        ),
    ):
        with pytest.raises(AssertionError) as info:
            sweep_python_allocation(make, 200)
        assert traceback.extract_tb(info.tb)[-1].line == clause
    # A call that needs more attempts than the limit: the sweep fails, never passes.
    with pytest.raises(pytest.fail.Exception) as info:
        sweep_python_allocation(lambda: (ngc.spf, (g, 0, opts()), list), 2)
    assert str(info.value) == "spf never completed"
