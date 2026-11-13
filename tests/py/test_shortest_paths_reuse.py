"""Baseline capture preserves fresh placement, program identities and refusals."""

from concurrent.futures import ThreadPoolExecutor
from functools import partial

import pytest
from concurrency_support import run_python, run_test_in_subprocess
from helpers import (
    build,
    full_ledger_state,
    ledger_state,
    plain,
    stats_of,
    sweep_python_allocation,
)
from private_hooks import HOOKS

import netgraph_core as c

OUT = c.OutputOptions(
    edge_flows=True,
    residual=True,
    saturated_boundary_edges=True,
    edge_reservations=True,
    min_cut=True,
    program=True,
    by_source=True,
)


def baseline(g, sources=0, destinations=1, mask=None, config=None):
    ledger = c.FlowLedger(g)
    config = config or c.FlowPolicyConfig(routing=c.ShortestPaths())
    p = c.FlowPolicy(ledger, sources, destinations, config, mask)
    p.replace(100, outputs=OUT)
    return ledger, p, c.ShortestPathsReuse([p])


def test_precommit_programs_are_fresh_and_commit_retains_only_its_program():
    g = build(2, [(0, 1, 10, 1), (0, 1, 10, 1)])
    _, base, reuse = baseline(g)
    mask = g.eligibility(excluded_edges=[0])
    hint = reuse.prepare(mask)[0]
    assert hint is not None
    ledger = c.FlowLedger(g)
    p = c.FlowPolicy(ledger, 0, 1, base.config, mask, preparation=hint)
    if HOOKS:
        assert p._reuse_consumed()
    assert p.program is None and p.place(0, OUT).program is None
    a, b = p.evaluate(7, outputs=OUT), p.evaluate(7, outputs=OUT)
    assert a.program is not b.program and p.program is None
    assert ledger.carried().tolist() == [0, 0]
    committed = p.replace(7, outputs=OUT)
    assert committed.program is p.program
    assert committed.program is not a.program and committed.program is not b.program
    assert p.evaluate(9, outputs=OUT).program is committed.program
    p.remove()
    assert p.program is None
    assert p.replace(8, outputs=OUT).program is committed.program
    if HOOKS:
        p = c.FlowPolicy(c.FlowLedger(g), 0, 1, base.config, mask, preparation=hint)
        fresh = c.FlowPolicy(c.FlowLedger(g), 0, 1, base.config, mask)
        result, trace = p._place_traced(7, OUT)
        expected, expected_trace = fresh._place_traced(7, OUT)
        assert plain(result) == plain(expected)
        assert trace.as_dict() == expected_trace.as_dict()


def test_effective_mask_alias_restoration_and_surviving_roster():
    g = build(3, [(0, 1, 2, 1), (0, 2, 2, 1)])
    masked = g.eligibility(excluded_nodes=[1], excluded_edges=[0])
    _, base, reuse = baseline(g, destinations=[1, 2], mask=masked)
    alias = g.eligibility(excluded_nodes=[1])
    assert masked.uid() != alias.uid()
    hint = reuse.prepare(alias)[0]
    assert hint is not None and reuse.prepare(None) == (None,)
    p = c.FlowPolicy(c.FlowLedger(g), 0, [1, 2], base.config, alias, preparation=hint)
    assert p.replace(100, outputs=OUT).delivered == 2
    assert p.records()[0].destinations == [2]
    empty = g.eligibility(excluded_nodes=[0, 1])
    hint = reuse.prepare(empty)[0]
    assert hint is not None
    p = c.FlowPolicy(c.FlowLedger(g), 0, [1, 2], base.config, empty, preparation=hint)
    r = p.replace(100, outputs=OUT)
    assert r.delivered == 0 and r.termination == c.Termination.NO_ROUTE
    assert not r.routes_changed and not r.weights_changed


def test_internal_distance_change_accepted_but_changed_source_rebuilt():
    edges = [(0, 1, 10, 1), (0, 2, 10, 1), (1, 3, 10, 1), (1, 2, 10, 3), (2, 3, 10, 1)]
    g = build(4, edges)
    _, base, reuse = baseline(g, destinations=3)
    mask = g.eligibility(excluded_edges=[2])
    hint = reuse.prepare(mask)[0]
    assert hint is not None
    p = c.FlowPolicy(c.FlowLedger(g), 0, 3, base.config, mask, preparation=hint)
    fresh = c.FlowPolicy(c.FlowLedger(g), 0, 3, base.config, mask)
    assert plain(p.replace(100, outputs=OUT)) == plain(fresh.replace(100, outputs=OUT))
    all_paths_lost = g.eligibility(excluded_edges=[2, 4])
    assert reuse.prepare(all_paths_lost) == (None,)


@pytest.mark.parametrize(
    "options",
    [
        {"multipath": False},
        {"parallel_edges": False},
        {"allocation": c.RouteAllocation.WEIGHTED},
    ],
)
def test_unsupported_choice_stays_fresh(options):
    g = build(2, [(0, 1, 10, 1), (0, 1, 10, 1)])
    cfg = c.FlowPolicyConfig(routing=c.ShortestPaths(**options))
    _, _, reuse = baseline(g, config=cfg)
    assert reuse.prepare(g.eligibility(excluded_edges=[0])) == (None,)


def test_validation_frozen_cache_and_simultaneous_guard_leave_state():
    g = build(2, [(0, 1, 10, 1), (0, 1, 10, 1)])
    _, base, reuse = baseline(g)
    mask = g.eligibility(excluded_edges=[0])
    hint = reuse.prepare(mask)[0]
    ledger = c.FlowLedger(g)
    before = ledger_state(ledger)
    with pytest.raises(ValueError, match="preparation Eligibility does not match"):
        c.FlowPolicy(ledger, 0, 1, base.config, preparation=hint)
    cache = c.ShortestPathCache(g)
    cache.freeze()
    with pytest.raises(c.ShortestPathCacheMiss):
        c.FlowPolicy(ledger, 0, 1, base.config, mask, cache, preparation=hint)
    assert ledger_state(ledger) == before
    p = c.FlowPolicy(ledger, 0, 1, base.config, mask, preparation=hint)
    with pytest.raises(ValueError, match="SIMULTANEOUS requires admission=LOSSY"):
        c.apply_actions(ledger, [c.Replace(p, 7)], mode=c.BatchMode.SIMULTANEOUS)
    assert ledger_state(ledger) == before
    other = build(2, [(0, 1, 10, 1)])
    with pytest.raises(ValueError, match="preparation belongs to another graph"):
        c.FlowPolicy(c.FlowLedger(other), 0, 1, base.config, preparation=hint)


def test_capture_outlives_original_owners_and_is_shared_across_workers():
    g = build(2, [(0, 1, 10, 1), (0, 1, 10, 1)])
    ledger, base, reuse = baseline(g)
    cfg = base.config
    base.remove()
    del base, ledger

    def case(i):
        mask = g.eligibility(excluded_edges=[i % 2])
        hint = reuse.prepare(mask)[0]
        p = c.FlowPolicy(c.FlowLedger(g), 0, 1, cfg, mask, preparation=hint)
        if HOOKS:
            assert p._reuse_consumed()
        return p.replace(100, outputs=OUT).delivered

    with ThreadPoolExecutor(max_workers=8) as workers:
        assert list(workers.map(case, range(64))) == [10] * 64
    if HOOKS:
        run_python("""
import weakref
network = graph()
ledger = ngc.FlowLedger(network)
owner = policy(ledger)
owner.replace(7)
owners = [owner]
weak = weakref.ref(owner)
del owner
holder, waiter = new_probe(pause='ledger_acquired'), new_probe()
a = start(lambda: ledger.carried(), holder)
assert holder.wait('ledger_acquired')
b = start(lambda: ngc.ShortestPathsReuse(owners), waiter)
assert waiter.wait('ledger_wait')
cleared = threading.Event()
def clear():
    owners.clear()
    cleared.set()
drop = start(clear)
assert cleared.wait(0.5), 'capture must own copied policy arguments while its gate is queued'
assert weak() is not None
holder.release()
finish(a)
reuse = finish(b)
finish(drop)
gc.collect()
assert weak() is None, 'completed capture must release the original policies'
assert reuse.prepare()[0] is not None
""")


def test_factory_and_binding_allocation_failures_are_atomic_and_retryable():
    if HOOKS:
        run_test_in_subprocess(__file__, "_factory_and_binding_allocation_failures")
    else:
        _factory_and_binding_allocation_failures()


def _factory_and_binding_allocation_failures():
    g = build(4, [(0, 1, 10, 1), (0, 2, 10, 1), (1, 3, 10, 1), (2, 3, 10, 1)])
    base_ledger, base, reuse = baseline(g, destinations=3)
    mask = g.eligibility(excluded_edges=[2])
    hint = reuse.prepare(mask)[0]
    ledger = c.FlowLedger(g)
    config = base.config
    factories = [
        lambda: c.ShortestPathsReuse([base]).prepare(mask)[0] is not None,
        lambda: reuse.prepare(mask)[0] is not None,
        lambda: bool(c.FlowPolicy(ledger, 0, 3, config, mask, preparation=hint)),
    ]
    for factory in factories:
        assert factory()
    if HOOKS:
        import _netgraph_core as ext

        def snapshot():
            return [
                full_ledger_state(base_ledger),
                full_ledger_state(ledger),
                base._state(),
            ]

        def unfailed(fn, args):
            hit, seen, value, error, _ = ext._fail_native_allocation(fn, args, 1000000)
            assert not hit and error is None
            return seen, value

        def trailing(index, seen, matched_base):
            assert matched_base and seen >= 2 and index == seen - 1

        base_seen, base_value = unfailed(c.ShortestPaths, ())
        base_repr = repr(base_value)
        del base_value
        assert base_seen == 3
        with pytest.raises(AssertionError):
            trailing(base_seen - 2, base_seen, True)
        print(
            f"REJECTED_SYNTHETIC_ESCAPE ShortestPaths k={base_seen - 2} seen={base_seen}"
        )
        before = snapshot()
        with pytest.raises(MemoryError):
            ext._fail_native_allocation(c.ShortestPaths, (), base_seen - 1)
        assert snapshot() == before
        assert repr(c.ShortestPaths()) == base_repr
        clean_seen, clean_value = unfailed(c.ShortestPaths, ())
        assert clean_seen == base_seen and repr(clean_value) == base_repr
        del clean_value
        print(
            f"ACCEPTED_BASE_INIT_INSTANCE ShortestPaths k={base_seen - 1} seen={base_seen}"
        )

        owners = [base]
        calls = [
            (
                "capture",
                c.ShortestPathsReuse,
                (owners,),
                lambda value: tuple(h is not None for h in value.prepare(mask)),
                True,
            ),
            (
                "prepare",
                reuse.prepare,
                (mask,),
                lambda value: tuple(h is not None for h in value),
                False,
            ),
            (
                "policy",
                partial(c.FlowPolicy, ledger, 0, 3, config, mask, preparation=hint),
                (),
                lambda value: (value.program is None, value._reuse_consumed()),
                False,
            ),
        ]
        for name, fn, args, observe, allow_trailing in calls:
            seen, value = unfailed(fn, args)
            expected = observe(value)
            del value
            assert seen > 0
            if name == "capture":
                assert seen - 1 == 22, f"capture escape 22 is not last: seen={seen}"
                with pytest.raises(AssertionError):
                    trailing(seen - 2, seen, True)
                print(f"REJECTED_SYNTHETIC_ESCAPE capture k={seen - 2} seen={seen}")
            for at in range(seen + 1):
                before = snapshot()
                escaped = False
                try:
                    hit, _, value, error, live = ext._fail_native_allocation(
                        fn, args, at
                    )
                except MemoryError:
                    trailing(at, seen, allow_trailing)
                    escaped = True
                    assert snapshot() == before
                    print(f"ACCEPTED_BASE_INIT_INSTANCE {name} k={at} seen={seen}")
                if escaped:
                    retry = fn(*args)
                    assert observe(retry) == expected
                    del retry
                    clean_seen, clean_value = unfailed(fn, args)
                    assert clean_seen == seen and observe(clean_value) == expected
                    del clean_value
                elif hit:
                    assert type(error) is MemoryError and value is None
                    assert snapshot() == before and live == 0
                    del error
                    retry = fn(*args)
                    assert observe(retry) == expected
                    del retry
                else:
                    assert at == seen and error is None and observe(value) == expected
                    del value
                    break
                print(
                    f"FAULT_ROW {name} k={at} seen={seen} escaped={escaped} state_retry=PASS"
                )

        def make():
            return (
                lambda: tuple(h is not None for h in reuse.prepare(mask)),
                (),
                snapshot,
            )

        assert sweep_python_allocation(make, limit=1000) > 0


def test_heap_costs_default_quantum_wide_requests_and_infinite_retry():
    g = build(3, [(0, 1, 1, 1), (0, 2, 1, 65537), (1, 2, 1, 65536)], q=2**-18)
    _, base, reuse = baseline(g, destinations=2)
    mask = g.eligibility(excluded_edges=[0])
    hint = reuse.prepare(mask)[0]
    assert hint is not None
    ledger, reference = c.FlowLedger(g), c.FlowLedger(g)
    p = c.FlowPolicy(ledger, 0, 2, base.config, mask, preparation=hint)
    fresh = c.FlowPolicy(reference, 0, 2, base.config, mask)
    offer = 2**61 - 1
    for count in range(1, 10):
        got, want = p.place(offer, OUT), fresh.place(offer, OUT)
        assert got.requested == count * offer
        assert plain(got) == plain(want)
        assert ledger_state(ledger) == ledger_state(reference)
    assert got.requested > 2**64 and got.delivered == 1
    assert plain(p.replace(float("inf"), outputs=OUT)) == plain(
        fresh.replace(float("inf"), outputs=OUT)
    )
    before = ledger_state(ledger)
    with pytest.raises(ValueError) as error:
        p.place(1)
    with pytest.raises(ValueError) as expected:
        fresh.place(1)
    assert str(error.value) == str(expected.value)
    assert ledger_state(ledger) == before == ledger_state(reference)


def test_failed_first_source_and_case_replay_basis_match_fresh():
    g = build(3, [(0, 2, 3, 2), (1, 2, 2, 1)])
    _, base, reuse = baseline(g, sources=[0, 1], destinations=2)
    mask = g.eligibility(excluded_nodes=[0])
    hint = reuse.prepare(mask)[0]
    assert hint is not None
    p = c.FlowPolicy(c.FlowLedger(g), [0, 1], 2, base.config, mask, preparation=hint)
    fresh = c.FlowPolicy(c.FlowLedger(g), [0, 1], 2, base.config, mask)
    got, want = p.replace(100, outputs=OUT), fresh.replace(100, outputs=OUT)
    assert plain(got) == plain(want)
    assert [(s.source, s.delivered) for s in got.by_source] == [(0, 0), (1, 2)]
    assert [r.sources for r in p.records()] == [[1]]
    replay = c.FlowPolicy(
        c.FlowLedger(g), [0, 1], 2, c.FlowPolicyConfig(routing=c.Pinned(got.program))
    )
    control = c.FlowPolicy(
        c.FlowLedger(g), [0, 1], 2, c.FlowPolicyConfig(routing=c.Pinned(want.program))
    )
    assert plain(replay.replace(100, outputs=OUT)) == plain(
        control.replace(100, outputs=OUT)
    )
    assert replay.delivered() == 2


def test_ordered_earlier_demand_recomputes_case_admission():
    state = full_ledger_state if HOOKS else ledger_state
    g = build(4, [(0, 2, 10, 1), (1, 2, 10, 1), (1, 3, 10, 2), (2, 3, 10, 1)])
    cfg = c.FlowPolicyConfig(routing=c.ShortestPaths())
    baseline_ledger = c.FlowLedger(g)
    base = [c.FlowPolicy(baseline_ledger, s, 3, cfg) for s in (1, 0)]
    original = c.apply_actions(
        baseline_ledger, [c.Replace(p, 10) for p in base], outputs=OUT
    )
    assert [r.delivered for r in original] == [10, 5]
    hints = c.ShortestPathsReuse(base).prepare(g.eligibility(excluded_edges=[2]))
    assert all(h is not None for h in hints)
    mask = g.eligibility(excluded_edges=[2])
    ledger, control = c.FlowLedger(g), c.FlowLedger(g)
    case = [
        c.FlowPolicy(ledger, s, 3, cfg, mask, preparation=h)
        for s, h in zip((1, 0), hints, strict=True)
    ]
    fresh = [c.FlowPolicy(control, s, 3, cfg, mask) for s in (1, 0)]
    if HOOKS:
        assert all(p._reuse_consumed() for p in case)
    got = c.apply_actions(ledger, [c.Replace(p, 10) for p in case], outputs=OUT)
    want = c.apply_actions(control, [c.Replace(p, 10) for p in fresh], outputs=OUT)
    assert [r.delivered for r in got] == [10, 0]
    assert plain(got) == plain(want)
    assert state(ledger) == state(control)
    case[0].remove()
    fresh[0].remove()
    assert plain(case[1].replace(10, outputs=OUT)) == plain(
        fresh[1].replace(10, outputs=OUT)
    )
    assert case[1].delivered() == 10
    assert state(ledger) == state(control)


@pytest.mark.parametrize("frozen", [False, True])
def test_supplied_cache_keeps_all_five_statistics_and_uid_misses(frozen):
    g = build(2, [(0, 1, 10, 1), (0, 1, 10, 1)])
    _, base, reuse = baseline(g)
    mask = g.eligibility(excluded_edges=[0])
    hint = reuse.prepare(mask)[0]
    caches = [c.ShortestPathCache(g), c.ShortestPathCache(g)]

    def construct(current):
        actual = c.FlowPolicy(
            c.FlowLedger(g), 0, 1, base.config, current, caches[0], preparation=hint
        )
        expected = c.FlowPolicy(c.FlowLedger(g), 0, 1, base.config, current, caches[1])
        if HOOKS:
            assert not actual._reuse_consumed()
        assert stats_of(caches[0]) == stats_of(caches[1])
        assert plain(actual.replace(7, outputs=OUT)) == plain(
            expected.replace(7, outputs=OUT)
        )
        return actual, expected

    cold = construct(mask)
    if frozen:
        for cache in caches:
            cache.freeze()
    warm = construct(mask)
    alias = g.eligibility(excluded_edges=[0])
    assert alias.uid() != mask.uid()
    if frozen:
        errors = []
        for index, cache in enumerate(caches):
            kwargs = {"preparation": hint} if index == 0 else {}
            with pytest.raises(c.ShortestPathCacheMiss) as error:
                c.FlowPolicy(c.FlowLedger(g), 0, 1, base.config, alias, cache, **kwargs)
            errors.append(str(error.value))
        assert errors[0] == errors[1]
    else:
        alias_owners = construct(alias)
        del alias_owners
    assert stats_of(caches[0]) == stats_of(caches[1])
    for cache in caches:
        cache.thaw()
        cache.clear()
    assert stats_of(caches[0]) == stats_of(caches[1])
    del cold, warm
    assert stats_of(caches[0]) == stats_of(caches[1])


def test_case_program_replay_does_not_restore_failed_parallel_basis():
    state = full_ledger_state if HOOKS else ledger_state
    g = build(2, [(0, 1, 10, 1), (0, 1, 10, 1)])
    _, base, reuse = baseline(g)
    mask = g.eligibility(excluded_edges=[0])
    hint = reuse.prepare(mask)[0]
    case = c.FlowPolicy(c.FlowLedger(g), 0, 1, base.config, mask, preparation=hint)
    fresh = c.FlowPolicy(c.FlowLedger(g), 0, 1, base.config, mask)
    program = case.replace(7, outputs=OUT).program
    control_program = fresh.replace(7, outputs=OUT).program
    ledger, control = c.FlowLedger(g), c.FlowLedger(g)
    replay = c.FlowPolicy(ledger, 0, 1, c.FlowPolicyConfig(routing=c.Pinned(program)))
    expected = c.FlowPolicy(
        control, 0, 1, c.FlowPolicyConfig(routing=c.Pinned(control_program))
    )
    # No mask is supplied on replay: the case program itself owns the restricted basis.
    assert plain(replay.replace(7, outputs=OUT)) == plain(
        expected.replace(7, outputs=OUT)
    )
    assert replay.records()[0].route.via_edges().tolist() == [1]
    assert state(ledger) == state(control)


def test_isolated_node_restoration_is_not_mask_containment():
    g = build(3, [(0, 2, 2, 1)])
    mask = g.eligibility(excluded_nodes=[1])
    _, _, reuse = baseline(g, destinations=[1, 2], mask=mask)
    assert mask.edge_allowed().tolist() == [True]
    # Restoring isolated node 1 changes containment while every edge stays eligible.
    assert reuse.prepare() == (None,)


def test_every_surviving_source_must_keep_its_distance():
    state = full_ledger_state if HOOKS else ledger_state
    edges = [(0, 1, 4, 1), (1, 2, 4, 1), (5, 2, 4, 1), (5, 3, 4, 2), (3, 2, 4, 1)]
    g = build(6, edges)
    _, base, reuse = baseline(g, sources=[0, 5], destinations=2)
    # Source 0 keeps cost 2. Source 5 loses cost 1 and takes the cost-3 detour.
    # Each independent route carries 4, so equal-source lossless admission is 8.
    failed = g.input_positions().tolist().index(2)
    mask = g.eligibility(excluded_edges=[failed])
    assert reuse.prepare(mask) == (None,)
    ledger, control = c.FlowLedger(g), c.FlowLedger(g)
    actual = c.FlowPolicy(ledger, [0, 5], 2, base.config, mask)
    fresh = c.FlowPolicy(control, [0, 5], 2, base.config, mask)
    got, want = actual.replace(100, outputs=OUT), fresh.replace(100, outputs=OUT)
    assert got.admitted == got.delivered == 8
    assert [(row.source, row.delivered) for row in got.by_source] == [(0, 4), (5, 4)]
    assert [(r.sources, r.min_route_cost) for r in actual.records()] == [
        ([0], 2),
        ([5], 3),
    ]
    assert plain(got) == plain(want)
    assert state(ledger) == state(control)


def test_case_replay_basis_keeps_every_surviving_source():
    state = full_ledger_state if HOOKS else ledger_state
    g = build(5, [(0, 3, 3, 5), (2, 3, 3, 1), (4, 3, 3, 2)])
    _, base, reuse = baseline(g, sources=[0, 2, 4], destinations=3)
    mask = g.eligibility(excluded_nodes=[2])
    hint = reuse.prepare(mask)[0]
    assert hint is not None
    ledger, control = c.FlowLedger(g), c.FlowLedger(g)
    actual = c.FlowPolicy(ledger, [0, 2, 4], 3, base.config, mask, preparation=hint)
    fresh = c.FlowPolicy(control, [0, 2, 4], 3, base.config, mask)
    got, want = actual.replace(100, outputs=OUT), fresh.replace(100, outputs=OUT)
    # Two surviving sources each carry 3; their distinct costs are 5 and 2.
    assert got.admitted == got.delivered == 6
    positions = g.input_positions().tolist()
    expected = [([0], 5, [positions.index(0)]), ([4], 2, [positions.index(2)])]
    assert [
        (r.sources, r.min_route_cost, r.route.via_edges().tolist())
        for r in actual.records()
    ] == expected
    assert plain(got) == plain(want)
    assert state(ledger) == state(control)
    replay_ledger, replay_control = c.FlowLedger(g), c.FlowLedger(g)
    replay = c.FlowPolicy(
        replay_ledger, [0, 2, 4], 3, c.FlowPolicyConfig(routing=c.Pinned(got.program))
    )
    reference = c.FlowPolicy(
        replay_control, [0, 2, 4], 3, c.FlowPolicyConfig(routing=c.Pinned(want.program))
    )
    # No mask: the program's own basis must retain both surviving roots only.
    assert plain(replay.replace(100, outputs=OUT)) == plain(
        reference.replace(100, outputs=OUT)
    )
    assert replay.delivered() == 6
    assert [(r.sources, r.min_route_cost) for r in replay.records()] == [
        ([0], 5),
        ([4], 2),
    ]
    assert state(replay_ledger) == state(replay_control)


def test_changed_middle_source_cannot_be_overwritten_by_the_last_source():
    state = full_ledger_state if HOOKS else ledger_state
    edges = [
        (0, 1, 4, 1),
        (1, 2, 4, 1),
        (5, 2, 4, 1),
        (5, 3, 4, 2),
        (3, 2, 4, 1),
        (6, 2, 4, 1),
    ]
    g = build(7, edges)
    _, base, reuse = baseline(g, sources=[0, 5, 6], destinations=2)
    positions = g.input_positions().tolist()
    mask = g.eligibility(excluded_edges=[positions.index(2)])
    # Only the middle source moves: costs [2, 1, 1] become [2, 3, 1].
    assert reuse.prepare(mask) == (None,)
    ledger, control = c.FlowLedger(g), c.FlowLedger(g)
    actual = c.FlowPolicy(ledger, [0, 5, 6], 2, base.config, mask)
    fresh = c.FlowPolicy(control, [0, 5, 6], 2, base.config, mask)
    got, want = actual.replace(100, outputs=OUT), fresh.replace(100, outputs=OUT)
    # Three independent cap-4 routes admit 12 with equal-source lossless admission.
    assert got.admitted == got.delivered == 12
    assert [(row.source, row.delivered) for row in got.by_source] == [
        (0, 4),
        (5, 4),
        (6, 4),
    ]
    assert [(r.sources, r.min_route_cost) for r in actual.records()] == [
        ([0], 2),
        ([5], 3),
        ([6], 1),
    ]
    assert plain(got) == plain(want)
    assert state(ledger) == state(control)
    # Removing the off-route detour edge changes no baseline source distance.
    unchanged = g.eligibility(excluded_edges=[positions.index(4)])
    assert reuse.prepare(unchanged)[0] is not None
