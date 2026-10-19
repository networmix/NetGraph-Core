"""The public search through the binding (``spf``, docs/semantics.md#4-search): results, the
DAG's canonical form, the validation order on unnarrowed values, the enumerations, and the
all-or-nothing construction of the result objects."""

from __future__ import annotations

import itertools
import operator
import traceback

import numpy as np
import pytest
from helpers import (
    FEWER,
    HIGHER_AVAIL,
    HIGHER_CAP,
    L62,
    assert_read_only,
    build,
    exact,
    ids,
    plain,
    sweep_python_allocation,
    to_ids,
)
from private_hooks import needs_hooks

import netgraph_core as ngc

FWD = ngc.Direction.FORWARD
REV = ngc.Direction.REVERSE


def via(dag):
    return sorted(int(e) for e in dag.via_edges())


def opts(direction=FWD, **kw):
    return ngc.SpfOptions(direction, **kw)


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
    # bottleneck of every walk is 1, so the two-hop walk wins under FEWER_HOPS and
    # HIGHER_CAPACITY.
    g = build(5, [(0, 1, 1, 1), (1, 4, 1, 3), (1, 2, 2, 1), (2, 3, 2, 1), (3, 4, 2, 1)])
    for tie in (FEWER, HIGHER_CAP):
        f = ngc.spf(g, 0, opts(tie_break=tie, multipath=False, targets=4))
        assert via(f.dag) == ids(g, {0, 1})
        assert f.target_cost == 4
        r = ngc.spf(g, 4, opts(REV, tie_break=tie, multipath=False, targets=0))
        assert r.dag == f.dag


def test_residual_gating_and_tie_scores():
    # Capacities 5; residuals by input [5, 4, 1, 3, 2]: 0->1->3 has bottleneck 4.
    g = build(4, [(0, 1, 5, 1), (1, 3, 5, 1), (0, 2, 5, 1), (2, 3, 5, 1), (0, 3, 5, 2)])
    residual = np.array(to_ids(g, [5, 4, 1, 3, 2]), dtype=np.int64)
    r = ngc.spf(
        g,
        0,
        opts(tie_break=HIGHER_AVAIL, multipath=False, targets=3, residual=residual),
    )
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
        ngc.spf(g, 0, opts(multipath=False))
    ineffective = "spf tie_break requires multipath=False or parallel_edges=False"
    with pytest.raises(ValueError, match=exact(ineffective)):
        ngc.spf(g, 0, opts(tie_break=HIGHER_CAP))
    for multipath, parallel_edges in ((False, True), (True, False), (False, False)):
        ngc.spf(
            g,
            0,
            opts(
                tie_break=HIGHER_CAP,
                multipath=multipath,
                parallel_edges=parallel_edges,
                targets=None if multipath else 2,
            ),
        )
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
        ngc.spf(g, 0, opts(eligibility=other.eligibility([], [0])))
    # Several defects: the first in the documented order (roots, targets, overlap, single
    # path, the ineffective tie rule, residual, min_residual, Eligibility).
    with pytest.raises(ValueError, match=exact("roots: the node set is empty")):
        ngc.spf(g, [], opts(targets=[], residual=[1], min_residual=0))
    with pytest.raises(ValueError, match=rf"^targets: node {big} out of node range"):
        ngc.spf(g, 0, opts(targets=[0, big], residual=[]))
    with pytest.raises(ValueError, match="both endpoint sets"):
        ngc.spf(g, 0, opts(targets=0, residual=[]))
    with pytest.raises(ValueError, match=exact(single)):
        ngc.spf(g, 0, opts(multipath=False, tie_break=HIGHER_CAP, residual=[]))
    with pytest.raises(ValueError, match=exact(ineffective)):
        ngc.spf(g, 0, opts(tie_break=HIGHER_CAP, residual=[], min_residual=0))
    with pytest.raises(ValueError, match="^residual"):
        ngc.spf(g, 0, opts(residual=[1, 9], min_residual=0))
    with pytest.raises(ValueError, match="^min_residual"):
        ngc.spf(g, 0, opts(min_residual=0, eligibility=other.eligibility([], [0])))
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
    assert o.residual is None and o.eligibility is None
    with pytest.raises(
        ValueError, match=exact("direction: 7 is not a value of Direction")
    ):
        ngc.SpfOptions(ngc.Direction(7))
    with pytest.raises(TypeError):
        ngc.SpfOptions(FWD, fanout_edges=[0])
    with pytest.raises(TypeError, match=r"^spf\(\): incompatible function arguments"):
        ngc.spf(g, 0, opts(), fanout_edges=[0])
    with pytest.raises(TypeError, match=exact("SpfOptions expected (got KspOptions)")):
        ngc.spf(g, 0, ngc.KspOptions(1))
    with pytest.raises(TypeError):
        ngc.SpfOptions()
    # The route choice: keywords only, bools strictly (None and 1 are no bool).
    with pytest.raises(TypeError):
        ngc.SpfOptions(FWD, None)
    for bad in (None, 1, "yes"):
        for name in ("multipath", "parallel_edges"):
            with pytest.raises(TypeError):
                ngc.SpfOptions(FWD, **{name: bad})
    assert ngc.SpfOptions(FWD, multipath=np.bool_(False)).multipath is False
    # Roots as a NumPy array, with duplicates.
    r = ngc.spf(g, np.array([0, 0], dtype=np.int32), opts())
    assert list(r.costs) == [0, 1]
    # SpfOptions keeps targets, residual and min_residual as given (the same objects) and
    # validates only the direction: spf validates the rest at the call, in its order.
    assert opts().min_residual == 1
    targets, residual, low = [1, 1], [9], -(2**70)
    o = opts(targets=targets, residual=residual, min_residual=low)
    assert o.targets is targets and o.residual is residual and o.min_residual is low
    assert (o.tie_break, o.multipath, o.parallel_edges) == (FEWER, True, True)
    assert_read_only(o, "direction", "tie_break", "multipath", "parallel_edges")
    with pytest.raises(ValueError, match=r"^residual\[0\] = 9 is outside"):
        ngc.spf(g, 0, o)


@needs_hooks
def test_python_result_construction_is_all_or_nothing():
    # Every Python object allocation of the call fails once: the call raises (never crashes)
    # and a clean retry gives the same result, for a search with targets and a full search.
    g = build(5, ANYCAST)
    for o in (opts(targets=[2, 4]), opts()):
        failures = sweep_python_allocation(lambda o=o: (ngc.spf, (g, 0, o), list), 200)
        assert failures >= 2, o.targets


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
            raise RuntimeError("replaced") from None

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

    def fewer_allocations():
        made = itertools.count()
        # This C callable returns the existing result without allocating a Python object
        # or entering a Python frame, so a tracer cannot add allocations to that call.
        prebuilt = operator.itemgetter(0)

        def make():
            if next(made) >= 3:
                return prebuilt, ([full],), list
            return ngc.spf, (g, 0, opts()), list

        return make

    def raising():
        # The completing call raises after its result. It raises in the frame the sweep calls:
        # CPython keeps an exception whose traceback entry for an outer frame failed to allocate.
        made = itertools.count()

        def make():
            m = next(made)
            injected_call = False

            def snapshot():
                # The sweep reads this before injection and again before a clean retry.
                # A tracer can fail before raises() enters its body, so body-entry counts
                # cannot distinguish those two calls.
                nonlocal injected_call
                injected_call = not injected_call
                return []

            def raises(o):
                result = ngc.spf(g, 0, o)
                if m and injected_call:
                    raise LookupError("after its result")
                return result

            return raises, (opts(),), snapshot

        return make

    def by_turn(turns):
        # m = 0 is the reference, c = 0 the injected call, c = 1 its clean retry.
        # Advance at snapshot boundaries, even if a tracer failed before entering the call.
        made = itertools.count()

        def make():
            m = next(made)
            snapshots = itertools.count()
            c = 0

            def snapshot():
                nonlocal c
                c = next(snapshots)
                return []

            return (lambda o: turns(m, c)(o)), (opts(),), snapshot

        return make

    for make, clause in (
        (lambda: (absorbs, (opts(),), list), "assert value is None, (where, value)"),
        (
            lambda: (replaces, (opts(),), list),
            "assert type(error) is MemoryError, (where, repr(error))",
        ),
        (changing, "assert snapshot() == before, where"),
        (raising(), "assert error is None, (where, repr(error))"),
        (
            by_turn(lambda m, c: root1 if m and not c else root0),
            "assert plain(value) == want, where",
        ),
        (
            fewer_allocations(),
            "assert seen == failures, (where, seen)  # every allocation failed once",
        ),
        (
            by_turn(lambda m, c: root1 if m and c == 1 else root0),
            "assert plain(fn(*args)) == want, where",
        ),
    ):
        try:
            sweep_python_allocation(make, 200)
        except AssertionError as error:
            actual = traceback.extract_tb(error.__traceback__)[-1].line
            assert actual == clause, (actual, clause)
        else:
            pytest.fail(f"the planted bug passed its clause: {clause}")
    # A call that needs more attempts than the limit: the sweep fails, never passes.
    with pytest.raises(pytest.fail.Exception) as info:
        sweep_python_allocation(lambda: (ngc.spf, (g, 0, opts()), list), 2)
    assert str(info.value) == "spf never completed"


# The route choice. TRI: three parallel cost-1 links 0 -> 1 of capacities 4, 10,
# 100 and a cost-2 link of 50; availability (the residual) 4, 8, 5, 50.
TRI = [(0, 1, 4, 1), (0, 1, 10, 1), (0, 1, 100, 1), (0, 1, 50, 2)]
TRI_AVAILABLE = [4, 8, 5, 50]
TIE = ngc.RouteTieBreak
# The scores: HIGHER_CAPACITY 4, 10, 100; HIGHER_AVAILABLE_CAPACITY 4, 8, 5; HIGHER_HEADROOM
# 4/4 = 1, 8/10 = 4/5, 5/100 = 1/20. The costlier e3 never wins.
TRI_WINNER = {
    TIE.FEWER_HOPS: 0,
    TIE.HIGHER_CAPACITY: 2,
    TIE.HIGHER_AVAILABLE_CAPACITY: 1,
    TIE.HIGHER_HEADROOM: 0,
}


def test_tie_parallel():
    # TIE-parallel: the single path takes one link per rule.
    g = build(2, TRI)
    for tie, want in TRI_WINNER.items():
        options = opts(
            targets=[1], multipath=False, tie_break=tie, residual=TRI_AVAILABLE
        )
        assert via(ngc.spf(g, 0, options).dag) == [want], tie


def test_parallel_collapse():
    # Multipath with one member per (u, v) pair: the DAG keeps one link of the cost-1 class,
    # chosen by the same rule.
    g = build(2, TRI)
    for tie, want in TRI_WINNER.items():
        options = opts(parallel_edges=False, tie_break=tie, residual=TRI_AVAILABLE)
        assert via(ngc.spf(g, 0, options).dag) == [want], tie


def test_tie_routes():
    # TIE-routes: a direct link of capacity 3 (e1, cost 2) against two hops of capacity 10 (e0,
    # e2, cost 1 each). Without a residual every score reads static capacity, and HEADROOM reads
    # 1 on every positive link, so its bottlenecks tie and hops decide. TIE-bottleneck: with e2
    # at capacity 1 the two-hop bottleneck is 1 and every rule takes the direct link, although
    # the two-hop route's first link is the larger one.
    for e2, want in (
        (10, {TIE.FEWER_HOPS: [1], TIE.HIGHER_CAPACITY: [0, 2]}),
        (1, {TIE.FEWER_HOPS: [1], TIE.HIGHER_CAPACITY: [1]}),
    ):
        g = build(3, [(0, 1, 10, 1), (0, 2, 3, 2), (1, 2, e2, 1)])
        want[TIE.HIGHER_AVAILABLE_CAPACITY] = want[TIE.HIGHER_CAPACITY]
        want[TIE.HIGHER_HEADROOM] = [1]
        for tie in TIE.__members__.values():
            r = ngc.spf(g, 0, opts(targets=[2], multipath=False, tie_break=tie))
            assert via(r.dag) == want[tie], (e2, tie)


def test_headroom_zero_capacity():
    # Without a residual HIGHER_HEADROOM reads 0 on a zero-capacity link and 1 on a positive
    # one, so it takes e1; FEWER_HOPS takes the lower EdgeId e0.
    g = build(2, [(0, 1, 0, 1), (0, 1, 5, 1)])
    for tie, want in ((TIE.HIGHER_HEADROOM, [1]), (TIE.FEWER_HOPS, [0])):
        r = ngc.spf(g, 0, opts(targets=[1], multipath=False, tie_break=tie))
        assert via(r.dag) == want, tie


def test_tie_rule_validation():
    # spf's entry: the ineffective tie rule after the single-path rule and before the residual.
    g = build(3, [(0, 1, 2, 1), (1, 2, 2, 1)])
    ineffective = "spf tie_break requires multipath=False or parallel_edges=False"
    single = (
        "multipath = false requires targets (the single path joins a root to a target)"
    )
    with pytest.raises(ValueError, match=exact(ineffective)):
        ngc.spf(g, 0, opts(tie_break=TIE.HIGHER_CAPACITY))
    with pytest.raises(ValueError, match=exact(ineffective)):
        ngc.spf(g, 0, opts(tie_break=TIE.HIGHER_CAPACITY, residual=[1]))
    with pytest.raises(ValueError, match=exact(single)):
        ngc.spf(g, 0, opts(tie_break=TIE.HIGHER_CAPACITY, multipath=False))
    assert ngc.spf(
        g, 0, opts(tie_break=TIE.HIGHER_CAPACITY, parallel_edges=False)
    ).costs.tolist() == [0, 1, 2]
    assert ngc.spf(g, 0, opts()).costs.tolist() == [0, 1, 2]
