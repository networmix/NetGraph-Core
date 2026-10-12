"""The DagCache's modes through the binding (docs/semantics.md#9-the-dag-cache): the
transition table, DagCacheMiss as a dedicated exception propagated unchanged by the evaluators,
the evaluators' sharing rules (a frozen cache shared by every worker, a
building one rejected with two or more), the reuse of prepopulated Eligibility uids, frozen
lookups from several Python threads, and a frozen cache's statistics under failed calls."""

from __future__ import annotations

import itertools
import threading
import traceback

import _netgraph_core as ext
import pytest
from helpers import ALL_OUT, COSTONLY, GLOBAL, build, stats_of
from private_hooks import needs_hooks

import netgraph_core as ngc


def fabric():
    """Leaves 0..5 and spines 6..8, every leaf joined to every spine both ways (capacity 4,
    cost 1)."""
    rows = [(leaf, s) for leaf in range(6) for s in range(6, 9)]
    rows += [(s, leaf) for leaf, s in rows]
    return build(9, [(u, v, 4, 1) for u, v in rows])


def prepopulate(g, keys):
    """A building cache holding the CostOnly support of each destination set (a capacity call
    looks up the destinations' reverse support; the source is spine 6), frozen afterwards."""
    cache = ngc.DagCache(g)
    for dst in keys:
        ngc.capacity(g, 6, dst, ngc.CapacityOptions(COSTONLY, cache=cache))
    cache.freeze()
    return cache


def stats(cache):
    """Hits, misses and entries (the three statistics the transition table speaks of)."""
    return stats_of(cache)[:3]


def test_transition_table():
    g = fabric()
    cache = ngc.DagCache(g)
    opts = ngc.CapacityOptions(COSTONLY, cache=cache)
    want = ngc.capacity(g, 0, 3, ngc.CapacityOptions(COSTONLY)).placed
    assert want == 12  # three spines of 4: the comparisons below are not of zeros
    assert ngc.capacity(g, 0, 3, opts).placed == want
    assert stats(cache) == (0, 1, 1)
    cache.freeze()
    cache.freeze()  # no effect
    assert ngc.capacity(g, 1, 3, opts).placed == want  # a hit
    assert stats(cache) == (1, 1, 1)
    with pytest.raises(
        ngc.DagCacheMiss, match=r"no entry for roots \{4\}, direction Reverse"
    ):
        ngc.capacity(g, 0, 4, opts)
    assert stats(cache) == (1, 2, 1)  # one miss, nothing inserted
    with pytest.raises(ValueError, match="rejected while the cache is frozen"):
        cache.clear()
    assert stats(cache) == (1, 2, 1)
    cache.thaw()
    cache.thaw()  # no effect
    assert (
        ngc.capacity(g, 0, 4, opts).placed == want
    )  # building again: searches, inserts
    assert stats(cache) == (1, 3, 2)
    cache.clear()
    assert stats(cache) == (1, 3, 0)


def test_miss_is_a_dedicated_exception():
    assert issubclass(ngc.DagCacheMiss, Exception)
    assert not issubclass(ngc.DagCacheMiss, (RuntimeError, ValueError, TypeError))
    assert ngc.DagCacheMiss.__doc__


def test_frozen_cache_shared_by_every_worker():
    g = fabric()
    pairs = [(0, 3), (1, 4), (2, 5), ([0, 1], [3, 4])]
    want = [
        r.placed for r in ngc.batch_capacity(g, pairs, ngc.CapacityOptions(COSTONLY))
    ]
    building = ngc.DagCache(g)
    with pytest.raises(ValueError, match="workers == 1"):
        ngc.batch_capacity(
            g, pairs, ngc.CapacityOptions(COSTONLY, cache=building), workers=2
        )
    assert all(w > 0 for w in want)  # every pair places: the comparison is not of zeros
    cache = prepopulate(g, [[3], [4], [5], [3, 4]])
    for workers in (1, 2, 4):
        before = stats(cache)
        got = ngc.batch_capacity(
            g, pairs, ngc.CapacityOptions(COSTONLY, cache=cache), workers=workers
        )
        assert [r.placed for r in got] == want
        assert stats(cache) == (before[0] + 4, before[1], 4), workers  # a hit per pair
    # The lowest failing pair index's miss: pairs 0 and 2 miss, pair 0's key is reported, and
    # the call counts that one miss.
    partial = prepopulate(g, [[4]])
    for workers in (1, 2, 4):
        before = stats(partial)
        with pytest.raises(ngc.DagCacheMiss, match=r"roots \{3\}"):
            ngc.batch_capacity(
                g,
                pairs[:3],
                ngc.CapacityOptions(COSTONLY, cache=partial),
                workers=workers,
            )
        assert stats(partial) == (before[0], before[1] + 1, 1), workers


def test_prepopulated_uids_and_intersections():
    g = fabric()
    common = g.eligibility(excluded_edges=[0])
    narrow = g.eligibility(excluded_edges=[0, 1])  # contained in common
    other = g.eligibility(excluded_edges=[1])  # the intersection with common is new
    pairs = [(1, 4), (2, 5)]
    cache = ngc.DagCache(g)
    ngc.capacity(g, 1, 4, ngc.CapacityOptions(COSTONLY, elig=narrow, cache=cache))
    ngc.capacity(g, 2, 5, ngc.CapacityOptions(COSTONLY, elig=common, cache=cache))
    cache.freeze()
    opts = ngc.CapacityOptions(COSTONLY, elig=common, cache=cache)
    ngc.batch_capacity(g, pairs, opts, [narrow, None], workers=2)  # returns the operand
    ngc.batch_capacity(
        g,
        pairs,
        ngc.CapacityOptions(COSTONLY, cache=cache),
        [narrow, common],
        workers=2,
    )
    assert stats(cache) == (4, 2, 2)
    with pytest.raises(ngc.DagCacheMiss, match=r"roots \{4\}"):
        ngc.batch_capacity(g, pairs, opts, [other, None], workers=2)


def test_sensitivity_propagates_the_miss_and_global_reads_no_cache():
    g = fabric()
    cache = prepopulate(g, [[4]])
    opts = ngc.SensitivityOptions(COSTONLY, cache=cache)
    for workers in (1, 2, 4):
        before = stats(cache)
        with pytest.raises(ngc.DagCacheMiss, match=r"roots \{4\}"):
            ngc.sensitivity(g, 1, 4, opts, workers)
        assert stats(cache) == (before[0], before[1] + 1, 1), workers  # the miss, once
    want = ngc.sensitivity(g, [0, 1], [4, 5], ngc.SensitivityOptions(GLOBAL))
    assert want.baseline.placed > 0
    assert any(e.loss > 0 for e in want.entries)  # the comparison below is not of zeros
    empty = ngc.DagCache(g)
    empty.freeze()
    for workers in (1, 2, 4):
        got = ngc.sensitivity(
            g, [0, 1], [4, 5], ngc.SensitivityOptions(GLOBAL, cache=empty), workers
        )
        assert got.baseline.placed == want.baseline.placed
        assert [(e.edge, e.loss) for e in got.entries] == [
            (e.edge, e.loss) for e in want.entries
        ]
    assert stats(empty) == (0, 0, 0)


def test_frozen_lookups_from_python_threads():
    g = fabric()
    cache = prepopulate(g, [[d] for d in range(6)])
    want = {
        d: ngc.capacity(g, 6, d, ngc.CapacityOptions(COSTONLY)).placed for d in range(6)
    }
    # The spine's one edge of 4 into each leaf: the comparisons below are not of zeros.
    assert want == dict.fromkeys(range(6), 4)
    errors = []

    def reader(t):
        try:
            for i in range(60):
                d = (i + t) % 6
                r = ngc.capacity(g, 6, d, ngc.CapacityOptions(COSTONLY, cache=cache))
                assert r.placed == want[d]
        except BaseException as exc:  # reported below; pytest's own failures included
            errors.append(exc)

    before = stats(cache)
    threads = [threading.Thread(target=reader, args=(t,)) for t in range(4)]
    for t in threads:
        t.start()
    for t in threads:
        t.join()
    assert errors == []
    assert stats(cache) == (before[0] + 240, before[1], 6)


def plus(stats, hits=0, misses=0):
    return (stats[0] + hits, stats[1] + misses, *stats[2:])


def test_a_failed_call_changes_no_statistic_but_its_miss():
    # A public call publishes a frozen cache's hits only when it succeeds: a batch whose pair 0
    # hits and pair 1 misses raises the miss, keeps no hit and counts the miss once, with one
    # worker and with several; CostOnly sensitivity drops its baseline's hit the same way.
    g = fabric()
    cache = prepopulate(g, [[3]])
    opts = ngc.CapacityOptions(COSTONLY, cache=cache)
    for workers in (1, 2, 4):
        before = stats_of(cache)
        with pytest.raises(ngc.DagCacheMiss, match=r"roots \{4\}"):
            ngc.batch_capacity(
                g, [(0, 3), (1, 4), (2, 5), (0, 3)], opts, workers=workers
            )
        assert stats_of(cache) == plus(before, misses=1)
        before = stats_of(cache)
        with pytest.raises(ngc.DagCacheMiss, match=r"roots \{3\}"):
            ngc.sensitivity(
                g, 0, 3, ngc.SensitivityOptions(COSTONLY, cache=cache), workers
            )
        assert stats_of(cache) == plus(before, misses=1)
    before = stats_of(cache)
    assert ngc.batch_capacity(g, [(0, 3), (1, 3)], opts, workers=2)
    assert stats_of(cache) == plus(before, hits=2)


def test_concurrent_calls_keep_every_successful_hit():
    # Python threads (the evaluators release the GIL) on one frozen cache, each alternating a
    # successful capacity (one hit) and a failing batch (its hit dropped, one miss).
    g = fabric()
    cache = prepopulate(g, [[3], [4]])
    opts = ngc.CapacityOptions(COSTONLY, cache=cache)
    errors = []

    def caller():
        try:
            for _ in range(40):
                assert ngc.capacity(g, 0, 3, opts).placed > 0
                with pytest.raises(ngc.DagCacheMiss):
                    ngc.batch_capacity(g, [(1, 4), (2, 5)], opts, workers=2)
        except BaseException as exc:  # reported below; pytest's own failures included
            errors.append(exc)

    before = stats_of(cache)
    threads = [threading.Thread(target=caller) for _ in range(4)]
    for t in threads:
        t.start()
    for t in threads:
        t.join()
    assert errors == []
    assert stats_of(cache) == plus(before, hits=160, misses=160)


def sweep_frozen(cache, fn, args, hits=0, miss=False):
    """fn(*args) with each Python object allocation failing in turn: a failed call leaves all
    five statistics as they were, the result object's allocation included (for example
    capacity's allocation 8); the call that completes adds exactly `hits` hits, or, when its
    outcome is the DagCacheMiss, one miss. One case keeps the miss: when the call has ended in
    the miss and pybind11's translation of DagCacheMiss into its Python exception fails
    (MemoryError instead), the miss is counted (pybind11's own machinery is outside the
    all-or-nothing guarantee); every such failure comes after every allocation of the call
    itself. Returns the number of failures inside the call (the translation's excluded)."""
    inside = 0
    translated = (
        None  # the first allocation index past the call (translation of the miss)
    )
    for k in range(1000):
        before = stats_of(cache)
        hit, _, _, error = ext._fail_python_allocation(fn, args, k)
        after = stats_of(cache)
        if not hit:
            if miss:
                assert isinstance(error, ngc.DagCacheMiss), repr(error)
                assert after == plus(before, misses=1), (k, before, after)
            else:
                assert error is None, repr(error)
                assert after == plus(before, hits=hits), (k, before, after)
            return inside
        assert isinstance(error, (MemoryError, RuntimeError)), (k, repr(error))
        if miss and after == plus(before, misses=1):
            assert isinstance(error, MemoryError), (k, repr(error))
            translated = k if translated is None else translated
            continue
        assert translated is None, (
            k,
            translated,
            "a failure inside the call after its translation",
        )
        assert after == before, (k, before, after)
        inside += 1
    pytest.fail(f"{fn.__name__} never completed")


@needs_hooks
def test_frozen_statistics_survive_python_allocation_failures():
    g = fabric()
    cache = prepopulate(g, [[3], [4], [3, 4]])
    opts = ngc.CapacityOptions(
        COSTONLY,
        cache=cache,
        outputs=ALL_OUT,
    )
    pairs = [(0, 3), (1, 4), ([0, 1], [3, 4])]
    assert sweep_frozen(cache, ngc.capacity, (g, 0, 3, opts), hits=1) >= 3
    assert sweep_frozen(cache, ngc.capacity, (g, 0, 5, opts), miss=True) >= 3
    for workers in (1, 2):
        assert (
            sweep_frozen(
                cache, ngc.batch_capacity, (g, pairs, opts, None, workers), hits=3
            )
            >= 3
        )
        assert (
            sweep_frozen(
                cache,
                ngc.batch_capacity,
                (g, [(0, 3), (1, 5)], opts, None, workers),
                miss=True,
            )
            >= 3
        )
    sens = ngc.SensitivityOptions(COSTONLY, cache=cache)
    assert sweep_frozen(cache, ngc.sensitivity, (g, 0, 3, sens, 1), miss=True) >= 3
    # CostOnly sensitivity from a source that reaches no destination: one hit (the baseline's)
    # and no placed deletion.
    island = build(3, [(0, 1, 4, 1), (1, 0, 4, 1)])
    small = ngc.DagCache(island)
    ngc.capacity(island, 0, 1, ngc.CapacityOptions(COSTONLY, cache=small))
    small.freeze()
    for workers in (1, 2):
        assert (
            sweep_frozen(
                small,
                ngc.sensitivity,
                (
                    island,
                    2,
                    1,
                    ngc.SensitivityOptions(COSTONLY, cache=small),
                    workers,
                ),
                hits=1,
            )
            >= 3
        )


@needs_hooks
def test_the_frozen_sweep_rejects_a_wrong_call():
    # sweep_frozen rejects a wrong call, each on its own clause: one that absorbs the injected
    # failure; one that publishes a hit before the call the failure reaches; a hit where a miss
    # is expected; a miss counted by another cache; a completing call that raises; a hit count
    # other than expected; a failed translation of the miss reported as RuntimeError; a failure
    # inside a call after a translation; and a call that never completes.
    g = fabric()
    cache = prepopulate(g, [[3]])
    other = prepopulate(g, [[3]])
    opts = ngc.CapacityOptions(COSTONLY, cache=cache, outputs=ALL_OUT)
    elsewhere = ngc.CapacityOptions(COSTONLY, cache=other, outputs=ALL_OUT)
    hit = (g, 0, 3, opts)
    miss = (g, 0, 5, opts)

    def absorbs(*args):
        try:
            return ngc.capacity(*args)
        except MemoryError:
            return None

    def hits_first(*args):
        ngc.capacity(*args)
        return ngc.capacity(*args)

    def runtime_error(*args):
        try:
            return ngc.capacity(*args)
        except MemoryError:
            raise RuntimeError("in place of MemoryError") from None

    # The miss call until its allocations are used up (its translation fails last), then a
    # longer call of hits that fails inside.
    turn = ext._fail_python_allocation(ngc.capacity, miss, 10**6)[1]
    longer = (g, [(0, 3), (1, 3)], opts, None, 1)
    assert ext._fail_python_allocation(ngc.batch_capacity, longer, 10**6)[1] > turn + 1
    attempts = itertools.count()

    def inside_after_translation(*args):
        if next(attempts) < turn:
            return ngc.capacity(*args)
        return ngc.batch_capacity(*longer)

    for fn, args, expect, clause in (
        (
            absorbs,
            hit,
            {"hits": 2},
            "assert isinstance(error, (MemoryError, RuntimeError)), (k, repr(error))",
        ),
        (hits_first, hit, {"hits": 2}, "assert after == before, (k, before, after)"),
        (
            ngc.capacity,
            hit,
            {"miss": True},
            "assert isinstance(error, ngc.DagCacheMiss), repr(error)",
        ),
        (
            ngc.capacity,
            (g, 0, 5, elsewhere),
            {"miss": True},
            "assert after == plus(before, misses=1), (k, before, after)",
        ),
        (ngc.capacity, (g, 0, 99, opts), {}, "assert error is None, repr(error)"),
        (
            ngc.capacity,
            hit,
            {"hits": 2},
            "assert after == plus(before, hits=hits), (k, before, after)",
        ),
        (
            runtime_error,
            miss,
            {"miss": True},
            "assert isinstance(error, MemoryError), (k, repr(error))",
        ),
        (
            inside_after_translation,
            miss,
            {"miss": True},
            "assert translated is None, (",
        ),
    ):
        with pytest.raises(AssertionError) as info:
            sweep_frozen(cache, fn, args, **expect)
        assert traceback.extract_tb(info.tb)[-1].line == clause, clause

    # Made outside the sweep, as every bound object.
    uncached = ngc.CapacityOptions(COSTONLY)

    def endless(*args):
        while True:
            ngc.capacity(g, 0, 3, uncached)

    with pytest.raises(pytest.fail.Exception) as info:
        sweep_frozen(cache, endless, hit)
    assert str(info.value) == "endless never completed"
