"""The ShortestPathCache's modes through the binding (docs/semantics.md#9-the-shortest-path-cache): the
transition table, ShortestPathCacheMiss as a dedicated exception propagated unchanged by the evaluators,
the evaluators' sharing rules (a frozen cache shared by every worker, a
building one rejected with two or more), the reuse of prepopulated Eligibility uids, frozen
lookups from several Python threads, and a frozen cache's statistics under failed calls."""

from __future__ import annotations

import threading

import _netgraph_core as ext
import pytest
from concurrency_support import run_test_in_subprocess
from helpers import ALL_OUT, GLOBAL, SHORTEST_PATHS, build, exact, stats_of
from private_hooks import needs_hooks

import netgraph_core as ngc


def fabric():
    """Leaves 0..5 and spines 6..8, every leaf joined to every spine both ways (capacity 4,
    cost 1)."""
    rows = [(leaf, s) for leaf in range(6) for s in range(6, 9)]
    rows += [(s, leaf) for leaf, s in rows]
    return build(9, [(u, v, 4, 1) for u, v in rows])


def prepopulate(g, keys):
    """A building cache holding the ShortestPaths support of each destination set (a capacity call
    looks up the destinations' reverse support; the source is spine 6), frozen afterwards."""
    cache = ngc.ShortestPathCache(g)
    for dst in keys:
        ngc.evaluate_demand(
            g, 6, dst, ngc.EvaluationOptions(SHORTEST_PATHS, cache=cache)
        )
    cache.freeze()
    return cache


def stats(cache):
    """Hits, misses and entries (the three statistics the transition table speaks of)."""
    return stats_of(cache)[:3]


def test_transition_table():
    g = fabric()
    cache = ngc.ShortestPathCache(g)
    opts = ngc.EvaluationOptions(SHORTEST_PATHS, cache=cache)
    want = ngc.evaluate_demand(
        g, 0, 3, ngc.EvaluationOptions(SHORTEST_PATHS)
    ).call_delivered
    assert want == 12  # three spines of 4: the comparisons below are not of zeros
    assert ngc.evaluate_demand(g, 0, 3, opts).call_delivered == want
    assert stats(cache) == (0, 1, 1)
    cache.freeze()
    cache.freeze()  # no effect
    assert ngc.evaluate_demand(g, 1, 3, opts).call_delivered == want  # a hit
    assert stats(cache) == (1, 1, 1)
    with pytest.raises(
        ngc.ShortestPathCacheMiss, match=r"no entry for roots \{4\}, direction Reverse"
    ):
        ngc.evaluate_demand(g, 0, 4, opts)
    assert stats(cache) == (1, 2, 1)  # one miss, nothing inserted
    with pytest.raises(ValueError) as error:
        cache.clear()
    assert (
        str(error.value)
        == "ShortestPathCache.clear is rejected while the cache is frozen (thaw() it, then clear())"
    )
    assert stats(cache) == (1, 2, 1)
    cache.thaw()
    cache.thaw()  # no effect
    assert (
        ngc.evaluate_demand(g, 0, 4, opts).call_delivered == want
    )  # building again: searches, inserts
    assert stats(cache) == (1, 3, 2)
    cache.clear()
    assert stats(cache) == (1, 3, 0)


def test_miss_is_a_dedicated_exception():
    assert issubclass(ngc.ShortestPathCacheMiss, Exception)
    assert not issubclass(
        ngc.ShortestPathCacheMiss, (RuntimeError, ValueError, TypeError)
    )


def test_frozen_cache_shared_by_every_worker():
    g = fabric()
    pairs = [(0, 3), (1, 4), (2, 5), ([0, 1], [3, 4])]
    want = [
        r.call_delivered
        for r in ngc.evaluate_pairs(g, pairs, ngc.EvaluationOptions(SHORTEST_PATHS))
    ]
    building = ngc.ShortestPathCache(g)
    with pytest.raises(ValueError, match="workers == 1"):
        ngc.evaluate_pairs(
            g, pairs, ngc.EvaluationOptions(SHORTEST_PATHS, cache=building), workers=2
        )
    assert all(w > 0 for w in want)  # every pair places: the comparison is not of zeros
    cache = prepopulate(g, [[3], [4], [5], [3, 4]])
    for workers in (1, 2, 4):
        before = stats(cache)
        got = ngc.evaluate_pairs(
            g,
            pairs,
            ngc.EvaluationOptions(SHORTEST_PATHS, cache=cache),
            workers=workers,
        )
        assert [r.call_delivered for r in got] == want
        assert stats(cache) == (before[0] + 4, before[1], 4), workers  # a hit per pair
    # The lowest failing pair index's miss: pairs 0 and 2 miss, pair 0's key is reported, and
    # the call counts that one miss.
    partial = prepopulate(g, [[4]])
    for workers in (1, 2, 4):
        before = stats(partial)
        with pytest.raises(ngc.ShortestPathCacheMiss, match=r"roots \{3\}"):
            ngc.evaluate_pairs(
                g,
                pairs[:3],
                ngc.EvaluationOptions(SHORTEST_PATHS, cache=partial),
                workers=workers,
            )
        assert stats(partial) == (before[0], before[1] + 1, 1), workers


def test_prepopulated_uids_and_intersections():
    g = fabric()
    common = g.eligibility(excluded_edges=[0])
    narrow = g.eligibility(excluded_edges=[0, 1])  # contained in common
    other = g.eligibility(excluded_edges=[1])  # the intersection with common is new
    pairs = [(1, 4), (2, 5)]
    cache = ngc.ShortestPathCache(g)
    ngc.evaluate_demand(
        g, 1, 4, ngc.EvaluationOptions(SHORTEST_PATHS, eligibility=narrow, cache=cache)
    )
    ngc.evaluate_demand(
        g, 2, 5, ngc.EvaluationOptions(SHORTEST_PATHS, eligibility=common, cache=cache)
    )
    cache.freeze()
    opts = ngc.EvaluationOptions(SHORTEST_PATHS, eligibility=common, cache=cache)
    ngc.evaluate_pairs(g, pairs, opts, [narrow, None], workers=2)  # returns the operand
    ngc.evaluate_pairs(
        g,
        pairs,
        ngc.EvaluationOptions(SHORTEST_PATHS, cache=cache),
        [narrow, common],
        workers=2,
    )
    assert stats(cache) == (4, 2, 2)
    with pytest.raises(ngc.ShortestPathCacheMiss, match=r"roots \{4\}"):
        ngc.evaluate_pairs(g, pairs, opts, [other, None], workers=2)


def test_edge_failures_propagates_the_miss_and_global_rejects_a_cache():
    g = fabric()
    cache = prepopulate(g, [[4]])
    opts = ngc.EvaluationOptions(SHORTEST_PATHS, cache=cache)
    for workers in (1, 2, 4):
        before = stats(cache)
        with pytest.raises(ngc.ShortestPathCacheMiss, match=r"roots \{4\}"):
            ngc.evaluate_edge_failures(g, 1, 4, opts, workers)
        assert stats(cache) == (before[0], before[1] + 1, 1), workers  # the miss, once
    # Only ShortestPaths reads a ShortestPathCache: CentralizedTE given one is a configuration error at every worker
    # count, and the cache is not read.
    empty = ngc.ShortestPathCache(g)
    empty.freeze()
    for workers in (1, 2, 4):
        with pytest.raises(
            ValueError, match=exact("CentralizedTE does not use ShortestPathCache")
        ):
            ngc.evaluate_edge_failures(
                g, [0, 1], [4, 5], ngc.EvaluationOptions(GLOBAL, cache=empty), workers
            )
    assert stats(empty) == (0, 0, 0)


def test_frozen_lookups_from_python_threads():
    run_test_in_subprocess(__file__, "_frozen_lookups_from_python_threads")


def _frozen_lookups_from_python_threads():
    g = fabric()
    cache = prepopulate(g, [[d] for d in range(6)])
    want = {
        d: ngc.evaluate_demand(
            g, 6, d, ngc.EvaluationOptions(SHORTEST_PATHS)
        ).call_delivered
        for d in range(6)
    }
    # The spine's one edge of 4 into each leaf: the comparisons below are not of zeros.
    assert want == dict.fromkeys(range(6), 4)
    errors = []

    def reader(t):
        try:
            for i in range(60):
                d = (i + t) % 6
                r = ngc.evaluate_demand(
                    g, 6, d, ngc.EvaluationOptions(SHORTEST_PATHS, cache=cache)
                )
                assert r.call_delivered == want[d]
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
    # worker and with several; ShortestPaths' evaluate_edge_failures drops its baseline's hit the same way.
    g = fabric()
    cache = prepopulate(g, [[3]])
    opts = ngc.EvaluationOptions(SHORTEST_PATHS, cache=cache)
    for workers in (1, 2, 4):
        before = stats_of(cache)
        with pytest.raises(ngc.ShortestPathCacheMiss, match=r"roots \{4\}"):
            ngc.evaluate_pairs(
                g, [(0, 3), (1, 4), (2, 5), (0, 3)], opts, workers=workers
            )
        assert stats_of(cache) == plus(before, misses=1)
        before = stats_of(cache)
        with pytest.raises(ngc.ShortestPathCacheMiss, match=r"roots \{3\}"):
            ngc.evaluate_edge_failures(
                g, 0, 3, ngc.EvaluationOptions(SHORTEST_PATHS, cache=cache), workers
            )
        assert stats_of(cache) == plus(before, misses=1)
    before = stats_of(cache)
    assert ngc.evaluate_pairs(g, [(0, 3), (1, 3)], opts, workers=2)
    assert stats_of(cache) == plus(before, hits=2)


def test_concurrent_calls_keep_every_successful_hit():
    run_test_in_subprocess(__file__, "_concurrent_calls_keep_every_successful_hit")


def _concurrent_calls_keep_every_successful_hit():
    # Python threads (the evaluators release the GIL) on one frozen cache, each alternating a
    # successful capacity (one hit) and a failing batch (its hit dropped, one miss).
    g = fabric()
    cache = prepopulate(g, [[3], [4]])
    opts = ngc.EvaluationOptions(SHORTEST_PATHS, cache=cache)
    errors = []

    def caller():
        try:
            for _ in range(40):
                assert ngc.evaluate_demand(g, 0, 3, opts).call_delivered > 0
                with pytest.raises(ngc.ShortestPathCacheMiss):
                    ngc.evaluate_pairs(g, [(1, 4), (2, 5)], opts, workers=2)
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
    five statistics as they were, the result object's allocation included; the call that
    completes adds exactly `hits` hits, or, when its outcome is the ShortestPathCacheMiss, one
    miss. One case keeps the miss: when the call has ended in
    the miss and pybind11's translation of ShortestPathCacheMiss into its Python exception fails
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
                assert isinstance(error, ngc.ShortestPathCacheMiss), repr(error)
                assert after == plus(before, misses=1), (k, before, after)
            else:
                assert error is None, repr(error)
                assert after == plus(before, hits=hits), (k, before, after)
            return inside
        assert type(error) is MemoryError, (k, repr(error))
        if miss and after == plus(before, misses=1):
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
    opts = ngc.EvaluationOptions(
        SHORTEST_PATHS,
        cache=cache,
        outputs=ALL_OUT,
    )
    pairs = [(0, 3), (1, 4), ([0, 1], [3, 4])]
    # Each result shell exists before the lease, including calls that end in a
    # miss. Failing a shell leaves statistics unchanged; failing translation
    # after the native miss keeps exactly that one miss. List shell allocation
    # may reuse CPython's free list; each result instance is still a new object.
    assert sweep_frozen(cache, ngc.evaluate_demand, (g, 0, 3, opts), hits=1) >= 1
    assert sweep_frozen(cache, ngc.evaluate_demand, (g, 0, 5, opts), miss=True) >= 1
    for workers in (1, 2):
        assert (
            sweep_frozen(
                cache, ngc.evaluate_pairs, (g, pairs, opts, None, workers), hits=3
            )
            >= 3
        )
        assert (
            sweep_frozen(
                cache,
                ngc.evaluate_pairs,
                (g, [(0, 3), (1, 5)], opts, None, workers),
                miss=True,
            )
            >= 2
        )
    sens = ngc.EvaluationOptions(SHORTEST_PATHS, cache=cache)
    assert (
        sweep_frozen(cache, ngc.evaluate_edge_failures, (g, 0, 3, sens, 1), miss=True)
        >= 1
    )
    # ShortestPaths' evaluate_edge_failures from a source that reaches no destination: one hit (the baseline's)
    # and no placed deletion.
    island = build(3, [(0, 1, 4, 1), (1, 0, 4, 1)])
    small = ngc.ShortestPathCache(island)
    ngc.evaluate_demand(
        island, 0, 1, ngc.EvaluationOptions(SHORTEST_PATHS, cache=small)
    )
    small.freeze()
    for workers in (1, 2):
        assert (
            sweep_frozen(
                small,
                ngc.evaluate_edge_failures,
                (
                    island,
                    2,
                    1,
                    ngc.EvaluationOptions(SHORTEST_PATHS, cache=small),
                    workers,
                ),
                hits=1,
            )
            >= 1
        )
