"""The DagCache's modes through the binding (APPROACH 10, 11): the transition table, DagCacheMiss
as a dedicated exception (not a RuntimeError, so never a sensitivity status), the evaluators'
sharing rules (a frozen cache shared by every worker, a building one rejected with two or more),
the reuse of prepopulated Eligibility uids, frozen lookups from several Python threads, and a
frozen cache's statistics under failed calls (review-1d F1, NOTES R27)."""

from __future__ import annotations

import threading

import _netgraph_core as ext
import pytest
from private_hooks import needs_hooks

import netgraph_core as ngc


def fabric():
    """Leaves 0..5 and spines 6..8, every leaf joined to every spine both ways (capacity 4,
    cost 1)."""
    rows = [(leaf, s) for leaf in range(6) for s in range(6, 9)]
    rows += [(s, leaf) for leaf, s in rows]
    return ngc.Graph(
        9,
        [r[0] for r in rows],
        [r[1] for r in rows],
        [4] * len(rows),
        [1] * len(rows),
        1.0,
    )


def costonly():
    return ngc.FlowPolicyConfig(
        routing=ngc.CostOnly(),
        placement=ngc.FlowPlacement.EQUAL_BALANCED,
        admission=ngc.Admission.LOSSLESS,
        selection=ngc.PathSelection(
            ngc.EdgeTieBreak.DETERMINISTIC, ngc.EdgePriority.CAPACITY
        ),
    )


def prepopulate(g, keys, elig=None):
    """A building cache holding the CostOnly support of each destination set (a capacity call
    looks up the destinations' reverse support; the source is spine 6), frozen afterwards."""
    cache = ngc.DagCache(g)
    for dst in keys:
        ngc.capacity(g, 6, dst, ngc.CapacityOptions(costonly(), elig=elig, cache=cache))
    cache.freeze()
    return cache


def stats(cache):
    s = cache.stats()
    return (s.hits, s.misses, s.entries)


def test_transition_table():
    g = fabric()
    cache = ngc.DagCache(g)
    opts = ngc.CapacityOptions(costonly(), cache=cache)
    want = ngc.capacity(g, 0, 3, ngc.CapacityOptions(costonly())).placed
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
        r.placed for r in ngc.batch_capacity(g, pairs, ngc.CapacityOptions(costonly()))
    ]
    building = ngc.DagCache(g)
    with pytest.raises(ValueError, match="workers == 1"):
        ngc.batch_capacity(
            g, pairs, ngc.CapacityOptions(costonly(), cache=building), workers=2
        )
    cache = prepopulate(g, [[3], [4], [5], [3, 4]])
    before = stats(cache)
    for workers in (1, 2, 4):
        got = ngc.batch_capacity(
            g, pairs, ngc.CapacityOptions(costonly(), cache=cache), workers=workers
        )
        assert [r.placed for r in got] == want
    assert stats(cache) == (before[0] + 12, before[1], 4)
    # The lowest failing pair index's miss: pairs 0 and 2 miss, pair 0's key is reported.
    partial = prepopulate(g, [[4]])
    for workers in (1, 2, 4):
        with pytest.raises(ngc.DagCacheMiss, match=r"roots \{3\}"):
            ngc.batch_capacity(
                g,
                pairs[:3],
                ngc.CapacityOptions(costonly(), cache=partial),
                workers=workers,
            )


def test_prepopulated_uids_and_intersections():
    g = fabric()
    common = g.eligibility(excluded_edges=[0])
    narrow = g.eligibility(excluded_edges=[0, 1])  # contained in common
    other = g.eligibility(excluded_edges=[1])  # the intersection with common is new
    pairs = [(1, 4), (2, 5)]
    cache = ngc.DagCache(g)
    ngc.capacity(g, 1, 4, ngc.CapacityOptions(costonly(), elig=narrow, cache=cache))
    ngc.capacity(g, 2, 5, ngc.CapacityOptions(costonly(), elig=common, cache=cache))
    cache.freeze()
    opts = ngc.CapacityOptions(costonly(), elig=common, cache=cache)
    ngc.batch_capacity(g, pairs, opts, [narrow, None], workers=2)  # returns the operand
    ngc.batch_capacity(
        g,
        pairs,
        ngc.CapacityOptions(costonly(), cache=cache),
        [narrow, common],
        workers=2,
    )
    assert stats(cache) == (4, 2, 2)
    with pytest.raises(ngc.DagCacheMiss, match=r"roots \{4\}"):
        ngc.batch_capacity(g, pairs, opts, [other, None], workers=2)


def test_sensitivity_propagates_the_miss_and_global_reads_no_cache():
    g = fabric()
    cache = prepopulate(g, [[4]])
    opts = ngc.SensitivityOptions(costonly(), cache=cache)
    for workers in (1, 2, 4):
        with pytest.raises(ngc.DagCacheMiss, match=r"roots \{4\}"):
            ngc.sensitivity(g, 1, 4, opts, workers)
    glob = ngc.FlowPolicyConfig(
        routing=ngc.Global(),
        placement=ngc.FlowPlacement.MAX_FLOW,
        admission=ngc.Admission.LOSSLESS,
        selection=ngc.PathSelection(
            ngc.EdgeTieBreak.DETERMINISTIC, ngc.EdgePriority.CAPACITY
        ),
    )
    want = ngc.sensitivity(g, [0, 1], [4, 5], ngc.SensitivityOptions(glob))
    empty = ngc.DagCache(g)
    empty.freeze()
    for workers in (1, 2, 4):
        got = ngc.sensitivity(
            g, [0, 1], [4, 5], ngc.SensitivityOptions(glob, cache=empty), workers
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
        d: ngc.capacity(g, 6, d, ngc.CapacityOptions(costonly())).placed
        for d in range(6)
    }
    errors = []

    def reader(t):
        try:
            for i in range(60):
                d = (i + t) % 6
                r = ngc.capacity(g, 6, d, ngc.CapacityOptions(costonly(), cache=cache))
                assert r.placed == want[d]
        except Exception as exc:  # noqa: BLE001 - reported below
            errors.append(exc)

    before = stats(cache)
    threads = [threading.Thread(target=reader, args=(t,)) for t in range(4)]
    for t in threads:
        t.start()
    for t in threads:
        t.join()
    assert errors == []
    assert stats(cache) == (before[0] + 240, before[1], 6)


def all_stats(cache):
    s = cache.stats()
    return (s.hits, s.misses, s.entries, s.bytes, s.live_values)


def plus(stats, hits=0, misses=0):
    return (stats[0] + hits, stats[1] + misses, *stats[2:])


def test_a_failed_call_changes_no_statistic_but_its_miss():
    # A public call publishes a frozen cache's hits only when it succeeds: a batch whose pair 0
    # hits and pair 1 misses raises the miss, keeps no hit and counts the miss once, with one
    # worker and with several; CostOnly sensitivity drops its baseline's hit the same way.
    g = fabric()
    cache = prepopulate(g, [[3]])
    opts = ngc.CapacityOptions(costonly(), cache=cache)
    for workers in (1, 2, 4):
        before = all_stats(cache)
        with pytest.raises(ngc.DagCacheMiss, match=r"roots \{4\}"):
            ngc.batch_capacity(
                g, [(0, 3), (1, 4), (2, 5), (0, 3)], opts, workers=workers
            )
        assert all_stats(cache) == plus(before, misses=1)
        before = all_stats(cache)
        with pytest.raises(ngc.DagCacheMiss, match=r"roots \{3\}"):
            ngc.sensitivity(
                g, 0, 3, ngc.SensitivityOptions(costonly(), cache=cache), workers
            )
        assert all_stats(cache) == plus(before, misses=1)
    before = all_stats(cache)
    assert ngc.batch_capacity(g, [(0, 3), (1, 3)], opts, workers=2)
    assert all_stats(cache) == plus(before, hits=2)


def test_concurrent_calls_keep_every_successful_hit():
    # Python threads (the evaluators release the GIL) on one frozen cache, each alternating a
    # successful capacity (one hit) and a failing batch (its hit dropped, one miss).
    g = fabric()
    cache = prepopulate(g, [[3], [4]])
    opts = ngc.CapacityOptions(costonly(), cache=cache)
    errors = []

    def caller():
        try:
            for _ in range(40):
                assert ngc.capacity(g, 0, 3, opts).placed > 0
                with pytest.raises(ngc.DagCacheMiss):
                    ngc.batch_capacity(g, [(1, 4), (2, 5)], opts, workers=2)
        except Exception as exc:  # noqa: BLE001 - reported below
            errors.append(exc)

    before = all_stats(cache)
    threads = [threading.Thread(target=caller) for _ in range(4)]
    for t in threads:
        t.start()
    for t in threads:
        t.join()
    assert errors == []
    assert all_stats(cache) == plus(before, hits=160, misses=160)


def sweep_frozen(cache, fn, args, hits=0, miss=False):
    """fn(*args) with each Python object allocation failing in turn: a failed call leaves all
    five statistics as they were, the result object's allocation included (the reviewer's
    witness: capacity's allocation 8); the call that completes adds exactly `hits` hits, or, when
    its outcome is the DagCacheMiss, one miss. One case keeps the miss: when the call has ended
    in the miss and pybind11's translation of DagCacheMiss into its Python exception fails
    (MemoryError instead), the miss is counted (NOTES R27, APPROACH 13 limit 1); every such
    failure comes after every allocation of the call itself. Returns the number of failures."""
    failures = 0
    translated = (
        None  # the first allocation index past the call (translation of the miss)
    )
    for k in range(1000):
        before = all_stats(cache)
        hit, _, _, error = ext._fail_python_allocation(fn, args, k)
        after = all_stats(cache)
        if not hit:
            if miss:
                assert isinstance(error, ngc.DagCacheMiss), repr(error)
                assert after == plus(before, misses=1), (k, before, after)
            else:
                assert error is None, repr(error)
                assert after == plus(before, hits=hits), (k, before, after)
            return failures
        failures += 1
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
    pytest.fail(f"{fn.__name__} never completed")


@needs_hooks
def test_frozen_statistics_survive_python_allocation_failures():
    g = fabric()
    cache = prepopulate(g, [[3], [4], [3, 4]])
    opts = ngc.CapacityOptions(
        costonly(),
        cache=cache,
        outputs=ngc.OutputOptions(edge_flows=True, residual=True, bottleneck=True),
    )
    pairs = [(0, 3), (1, 4), ([0, 1], [3, 4])]
    assert sweep_frozen(cache, ngc.capacity, (g, 0, 3, opts), hits=1) >= 3
    assert sweep_frozen(cache, ngc.capacity, (g, 0, 5, opts), miss=True) >= 1
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
            >= 1
        )
    sens = ngc.SensitivityOptions(costonly(), cache=cache)
    assert sweep_frozen(cache, ngc.sensitivity, (g, 0, 3, sens, 1), miss=True) >= 1
    # CostOnly sensitivity from a source that reaches no destination: one hit (the baseline's),
    # no placed deletion, and a result with an entry per edge.
    island = ngc.Graph(3, [0, 1], [1, 0], [4, 4], [1, 1], 1.0)
    small = ngc.DagCache(island)
    ngc.capacity(island, 0, 1, ngc.CapacityOptions(costonly(), cache=small))
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
                    ngc.SensitivityOptions(costonly(), cache=small),
                    workers,
                ),
                hits=1,
            )
            >= 3
        )
