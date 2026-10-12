"""Native allocation failures in the staged result constructions and in FlowLedger.reset: the
k-th C++ allocation of spf, ksp, resolve_to_paths, of the calls whose result objects are made
before anything is published (place, capacity, batch_capacity, sensitivity) and of reset fails,
for every k in turn (the private hook _fail_native_allocation). Every unabsorbed failure
raises MemoryError, never aborts, leaves nothing it allocated behind and the observable state as
it was. Standard-library optional buffers may absorb a failure; every successful call, including
such a fallback, must return the same value and leave the same state as a call without failure.

Each sweep runs in a child interpreter, so a crash (a double free) fails its test instead of
the whole session: python -B test_native_allocation.py SCENARIO EXTENSION loads the parent's
extension file and prints one JSON row per attempt."""

from __future__ import annotations

import gc
import importlib.util
import itertools
import json
import os
import subprocess
import sys

if (
    __name__ == "__main__"
):  # the child: the parent's extension, whatever the import path finds
    _spec = importlib.util.spec_from_file_location("_netgraph_core", sys.argv[2])
    sys.modules["_netgraph_core"] = importlib.util.module_from_spec(_spec)
    _spec.loader.exec_module(sys.modules["_netgraph_core"])

import _netgraph_core as ext
from helpers import (
    ALL_OUT,
    COSTONLY,
    GLOBAL,
    build,
    full_ledger_state,
    plain,
    sel,
    stats_of,
)

import netgraph_core as ngc

EDGES = [  # (u, v, capacity units, cost): unequal costs, parallel members, several routes
    (0, 1, 5, 1),
    (1, 3, 5, 1),
    (0, 2, 3, 2),
    (2, 3, 3, 2),
    (1, 2, 2, 1),
    (0, 1, 2, 1),
    (3, 4, 9, 1),
    (2, 4, 2, 4),
]


def graph():
    return build(5, EDGES)


def caps(g):
    return [int(c) for c in g.capacity_units()]


def epoch_state(policy):
    """The policy's epoch (the requested total, the +inf flag, its flows), with each FlowId's
    local id only (the ledger uid differs between the scenario's fresh ledgers)."""
    state = policy._state()
    for flows in state["cohorts"].values():
        for flow in flows:
            flow["flow_id"] = flow["flow_id"][1]
    return state


# A scenario makes fresh state and returns (fn, args, snapshot): snapshot() is the observable
# state a failed call must leave as it was.
def stateless(call):
    def make():
        g = graph()
        fn, args = call(g)
        return fn, args, lambda: [caps(g), g.cost().tolist()]

    return make


def placement(cfg, traced):
    def make():
        g = graph()
        ledger = ngc.FlowLedger(g)
        policy = ngc.FlowPolicy(ledger, 0, 4, cfg)
        policy.place(
            2
        )  # an earlier call's flows in the same epoch: a failed call keeps them

        def snap():
            records = [(r.flow_id[1], r.delivered) for r in policy.records()]
            return [
                full_ledger_state(ledger),
                policy.placed_demand(),
                records,
                epoch_state(policy),
            ]

        return (policy._place_traced if traced else policy.place), (6, ALL_OUT), snap

    return make


def frozen_cache_capacity():
    g = graph()
    cache = ngc.DagCache(g)
    ngc.capacity(g, 0, 4, ngc.CapacityOptions(COSTONLY, cache=cache))
    cache.freeze()
    opts = ngc.CapacityOptions(COSTONLY, cache=cache, outputs=ALL_OUT)
    return ngc.capacity, (g, 0, 4, opts), lambda: stats_of(cache)


def with_orphan(call):
    """A ledger holding an orphaned record (its policy destroyed without remove()); call(ledger)
    gives (fn, args), and the snapshot is the ledger's state."""

    def make():
        g = graph()
        ledger = ngc.FlowLedger(g)
        ngc.FlowPolicy(ledger, 0, 4, COSTONLY).place(2)
        gc.collect()  # the policy is gone: its record is an orphan
        fn, args = call(ledger)
        return fn, args, lambda: full_ledger_state(ledger)

    return make


def half(g):
    return [c // 2 for c in caps(g)]


def spf_opts(direction=ngc.Direction.FORWARD, multipath=True, **kw):
    return ngc.SpfOptions(direction, sel(mp=multipath), **kw)


SCENARIOS = {
    "spf.full": stateless(lambda g: (ngc.spf, (g, 0, spf_opts()))),
    "spf.targets_residual": stateless(
        lambda g: (
            ngc.spf,
            (g, [0, 1], spf_opts(ngc.Direction.REVERSE, targets=[4], residual=caps(g))),
        )
    ),
    "spf.single_path": stateless(
        lambda g: (ngc.spf, (g, 0, spf_opts(multipath=False, targets=4)))
    ),
    "ksp": stateless(
        lambda g: (ngc.ksp, (g, 0, 4, ngc.KspOptions(4, sel(), residual=caps(g))))
    ),
    "FlowLedger.reset": with_orphan(lambda fg: (fg.reset, (half(fg.graph()),))),
    "resolve_to_paths": stateless(
        lambda g: (ngc.resolve_to_paths, (g, ngc.spf(g, 0, spf_opts()).dag, 0, 4, True))
    ),
    "place.costonly": placement(COSTONLY, traced=False),
    "place.global": placement(GLOBAL, traced=False),
    "place_traced.global": placement(GLOBAL, traced=True),
    "capacity.global": stateless(
        lambda g: (
            ngc.capacity,
            (g, 0, 4, ngc.CapacityOptions(GLOBAL, outputs=ALL_OUT)),
        )
    ),
    "capacity.frozen_cache": frozen_cache_capacity,
    "batch_capacity.global": stateless(
        lambda g: (
            ngc.batch_capacity,
            (g, [(0, 4), (1, 3)], ngc.CapacityOptions(GLOBAL), None, 1),
        )
    ),
    "sensitivity.global": stateless(
        lambda g: (ngc.sensitivity, (g, 0, 4, ngc.SensitivityOptions(GLOBAL), 1))
    ),
}


def absorbing(g):
    """A wrong call for the row checks' controls: it absorbs the injected failure."""

    def absorbs(o):
        try:
            return ngc.spf(g, 0, o)
        except MemoryError:
            return "absorbed"

    return absorbs, (spf_opts(),)


def changing():
    """A wrong call for the row checks' controls: it changes the observable state before the
    failure reaches it."""
    g = graph()
    log = []

    def changes(o):
        log.append(None)
        return ngc.spf(g, 0, o)

    return changes, (spf_opts(),), lambda: len(log)


def leaking():
    """A wrong call for the row checks' controls: it keeps a result it made before the failure
    reaches it (outside the snapshot), so that result's native allocations stay live."""
    g = graph()
    kept = []

    def leaks(o):
        kept.append(ngc.spf(g, 0, o))
        return ngc.spf(g, 0, o)

    return leaks, (spf_opts(),), lambda: [caps(g), g.cost().tolist()]


def replacing(g):
    """A wrong call for the row checks' controls: it raises another exception in place of the
    injected failure."""

    def replaces(o):
        try:
            return ngc.spf(g, 0, o)
        except MemoryError:
            raise LookupError("replaced") from None

    return replaces, (spf_opts(),)


def exiting(status, after):
    """Wrong calls for the row checks' controls: the process ends with `status` when the failure
    of an attempt after the first `after` reaches the call (a crash, or an exit with status 0)."""

    def make():
        m = next(MADE)  # 0: the reference call
        g = graph()

        def exits(o):
            try:
                return ngc.spf(g, 0, o)
            except MemoryError:
                if m > after:
                    os._exit(status)
                raise

        return exits, (spf_opts(),), lambda: [caps(g), g.cost().tolist()]

    return make


def walks():
    """A call for the row checks' controls that needs more attempts than the sweep makes: the
    walks of ten diamonds in a row (1,024 walks, far more than 5,000 native allocations)."""
    edges = []
    for i in range(10):
        a = 3 * i
        edges += [
            (a, a + 1, 1, 1),
            (a, a + 2, 1, 1),
            (a + 1, a + 3, 1, 1),
            (a + 2, a + 3, 1, 1),
        ]
    g = build(31, edges)
    dag = ngc.spf(g, 0, spf_opts()).dag
    return (
        ngc.resolve_to_paths,
        (g, dag, 0, 30, True),
        lambda: [caps(g), g.cost().tolist()],
    )


def by_turn(call):
    """A wrong call for the row checks' controls whose completing call differs from the
    reference: call(m, g) gives the call of the m-th make() (m = 0: the reference call)."""

    def make():
        g = graph()
        return call(next(MADE), g), (spf_opts(),), lambda: [caps(g), g.cost().tolist()]

    return make


def raises_after(m, g):
    """The completing call raises after its result (the reference call returns it)."""

    def raises(o):
        result = ngc.spf(g, 0, o)
        if m:
            raise LookupError("after its result")
        return result

    return raises


def appending():
    """A wrong call for the row checks' controls: the completing call adds two entries to the
    observable log where the reference call added one."""
    m = next(MADE)
    g = graph()
    log = []

    def appends(o):
        result = ngc.spf(g, 0, o)
        log.extend([None] * (2 if m else 1))
        return result

    return appends, (spf_opts(),), lambda: len(log)


def prebuilt(m, g):
    """After ten failures the call returns a result built beforehand: it completes with no native
    allocation, fewer than the failures before it."""
    built = ngc.spf(g, 0, spf_opts())

    def call(o):
        return built if m > 10 else ngc.spf(g, 0, o)

    return call


def discarding(m, g):
    """Keep the reference result, then discard it in each swept call. Comparing the returned
    value detects the discard even if thread-local scratch keeps native allocations live."""

    def discards(o):
        result = ngc.spf(g, 0, o)
        return 7 if m else result

    return discards


MADE = itertools.count()  # make() calls of a control in the child (each child runs one)

# Wrong calls that the row checks must reject (run by the controls, not by the sweep test), one
# per check. Their options are made before the call: a bound object constructed inside it would
# meet the failure in its constructor.
CONTROLS = {
    "control.absorbs": stateless(absorbing),
    "control.replaces": stateless(replacing),
    "control.leaks": leaking,
    "control.changes_state": changing,
    "control.crashes": exiting(3, 0),
    "control.never_completes": walks,
    "control.exits": exiting(0, 1),
    # reset under the floor of ten (its own scenario has a floor of its own, FEW)
    "control.few": with_orphan(lambda fg: (fg.reset, (half(fg.graph()),))),
    "control.raises_when_done": by_turn(raises_after),
    "control.other_value": by_turn(lambda m, g: lambda o: ngc.spf(g, 1 if m else 0, o)),
    "control.other_state": appending,
    "control.fewer_allocations": by_turn(prebuilt),
    "control.discards": by_turn(discarding),
}


def sweep(name: str) -> None:
    """Attempts k = 0, 1, ... of one scenario until a call completes without reaching k; one
    ROW line per attempt, written before the next attempt starts."""
    make = {**SCENARIOS, **CONTROLS}[name]
    fn, args, snap = make()
    reference = plain(fn(*args))  # also the warm-up of lazily sized scratch
    reference_after = snap()
    gc.collect()
    gc.disable()  # no unrelated object is released inside an attempt
    for k in range(5000):
        fn, args, snap = make()
        before = snap()
        print(f"START {k}", flush=True)
        hit, seen, value, error, live = ext._fail_native_allocation(fn, args, k)
        row = {"k": k, "hit": hit, "seen": seen, "live": live}
        row["error"] = None if error is None else f"{type(error).__name__}: {error}"
        if error is None:
            row["returns"] = value is not None
            row["value_ok"] = plain(value) == reference
            row["state_ok"] = snap() == reference_after
        else:
            row["state_ok"] = snap() == before
        del value, error
        print("ROW " + json.dumps(row), flush=True)
        if not hit:
            return
    print("ROW " + json.dumps({"k": None}), flush=True)


def run_child(name):
    proc = subprocess.run(
        [sys.executable, "-B", __file__, name, ext.__file__],
        capture_output=True,
        text=True,
        timeout=600,
    )
    rows, started = [], None
    for line in proc.stdout.splitlines():
        if line.startswith("START "):
            started = int(line[6:])
        elif line.startswith("ROW "):
            rows.append(json.loads(line[4:]))
    return proc, rows, started


if __name__ != "__main__":
    import traceback

    import pytest
    from private_hooks import needs_hooks

    # The scenarios with fewer than ten native allocations, and at least how many they make:
    # reset converts the residual from Python, then builds the new baseline.
    FEW = {"FlowLedger.reset": 3}

    def check(name, proc, rows, started):
        """The checks of one scenario's rows."""
        reached = [r for r in rows if r["k"] is not None and r["hit"]]
        for row in reached:
            if row["error"] is None:
                # std::stable_partition's optional nothrow buffer may fail and fall back.
                # Absorption is valid only with exactly the reference result and state.
                assert row["value_ok"], row
                assert row["state_ok"], row
            else:
                assert row["error"].startswith("MemoryError"), row
                assert row["live"] == 0, row
                assert row["state_ok"], row
        assert proc.returncode == 0, (
            f"{name}: the child died (rc {proc.returncode}) in attempt {started}\n"
            f"{proc.stderr[-2000:]}"
        )
        assert rows and rows[-1]["k"] is not None, f"{name}: the call never completed"
        # Stop only when k is not reached; an absorbed failure still advances the sweep.
        # A child that ended early with status 0 leaves a reached failure last.
        assert len(reached) == len(rows) - 1, f"{name}: {len(reached)} attempts"
        assert len(reached) >= FEW.get(name, 10), f"{name}: {len(reached)} attempts"
        done = rows[-1]
        assert done["error"] is None, done
        assert done["value_ok"], done
        assert done["state_ok"], done
        assert done["seen"] == len(reached), done
        for row in reached:
            if row["error"] is None:
                # A successful fallback holds exactly the result's native allocations,
                # not an extra temporary buffer or a leaked staged result.
                assert row["live"] == done["live"], row
        if done["returns"]:
            # The probe counts what the result object holds.
            assert done["live"] > 0, done

    def test_absorbed_rows_require_identical_success():
        # Exercise fallback acceptance on every platform, even when its standard library
        # does not request an optional buffer for this small graph.
        proc = subprocess.CompletedProcess([], 0, stderr="")
        rows = [
            {
                "k": k,
                "hit": True,
                "seen": k + 1,
                "live": 0,
                "error": "MemoryError: std::bad_alloc",
                "state_ok": True,
            }
            for k in range(10)
        ]
        rows[0].update(error=None, live=3, value_ok=True)
        rows.append(
            {
                "k": 10,
                "hit": False,
                "seen": 10,
                "live": 3,
                "error": None,
                "value_ok": True,
                "state_ok": True,
                "returns": True,
            }
        )
        check("fallback", proc, rows, 10)
        for field, wrong in (("value_ok", False), ("state_ok", False), ("live", 4)):
            corrupted = [dict(row) for row in rows]
            corrupted[0][field] = wrong
            with pytest.raises(AssertionError) as info:
                check("fallback", proc, corrupted, 10)
            clause = (
                'assert row["live"] == done["live"], row'
                if field == "live"
                else f'assert row["{field}"], row'
            )
            assert traceback.extract_tb(info.tb)[-1].line == clause
        corrupted = [dict(row) for row in rows]
        corrupted[-1]["live"] = 0
        corrupted[0]["live"] = 0
        with pytest.raises(AssertionError) as info:
            check("fallback", proc, corrupted, 10)
        assert traceback.extract_tb(info.tb)[-1].line == 'assert done["live"] > 0, done'

    @needs_hooks
    def test_allocation_hooks_preserve_exception_identity():
        class MarkedError(Exception):
            pass

        for hook in (ext._fail_python_allocation, ext._fail_native_allocation):
            cause = LookupError("cause")
            expected = MarkedError("instance")
            expected.marker = object()

            def raises(expected=expected, cause=cause):
                raise expected from cause

            hit, _, value, error, *_ = hook(raises, (), -1)
            assert not hit
            assert value is None
            assert error is expected
            assert error.marker is expected.marker
            assert error.__cause__ is cause
            assert error.__traceback__.tb_frame.f_code is raises.__code__
            # The pending exception was cleared; another call returns normally.
            assert hook(lambda: 7, (), -1)[2:4] == (7, None)

    @needs_hooks
    @pytest.mark.parametrize("name", list(SCENARIOS))
    def test_native_allocation_failure_raises_and_changes_nothing(name):
        check(name, *run_child(name))

    @needs_hooks
    def test_the_row_checks_reject_a_wrong_call():
        # Each row check rejects its wrong call: one that absorbs a failure with a wrong
        # result, raises another exception, leaks native allocations or changes the state;
        # one that crashes, outlasts its sweep, or exits with status 0 after a failure;
        # one with fewer allocations than the floor; and one whose completing call raises,
        # returns another value or state, makes fewer allocations than failed before it, or returns
        # a plain value in place of its result.
        attempts = (
            'assert len(reached) == len(rows) - 1, f"{name}: {len(reached)} attempts"'
        )
        for name, clause in (
            ("control.absorbs", 'assert row["value_ok"], row'),
            ("control.replaces", 'assert row["error"].startswith("MemoryError"), row'),
            ("control.leaks", 'assert row["live"] == 0, row'),
            ("control.changes_state", 'assert row["state_ok"], row'),
            ("control.crashes", "assert proc.returncode == 0, ("),
            (
                "control.never_completes",
                'assert rows and rows[-1]["k"] is not None, f"{name}: the call never completed"',
            ),
            ("control.exits", attempts),
            (
                "control.few",
                'assert len(reached) >= FEW.get(name, 10), f"{name}: {len(reached)} attempts"',
            ),
            ("control.raises_when_done", 'assert done["error"] is None, done'),
            ("control.other_value", 'assert done["value_ok"], done'),
            ("control.other_state", 'assert done["state_ok"], done'),
            ("control.fewer_allocations", 'assert done["seen"] == len(reached), done'),
            ("control.discards", 'assert done["value_ok"], done'),
        ):
            proc, rows, started = run_child(name)
            try:
                check(name, proc, rows, started)
            except AssertionError as error:
                actual = traceback.extract_tb(error.__traceback__)[-1].line
                assert actual == clause, name
            else:
                pytest.fail(f"{name}: the planted bug passed; last row: {rows[-1:]}")

else:
    if ngc.spf is not ext.spf:  # the package binds to the loaded extension
        sys.exit("netgraph_core does not use the extension the parent passed")
    if not hasattr(ext, "_fail_native_allocation"):
        sys.exit("the extension lacks the private hooks")
    sweep(sys.argv[1])
