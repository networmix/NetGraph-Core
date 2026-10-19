"""The 31 fixtures of the independent reference replayed through the bindings (the ledger and
the policies).

``tests/data/policy_fixtures.json`` was exported from the independent reference
(its frozen copy, the exporter, the regeneration check and the meaning of the data's historical
``source`` string: ``tests/reference/policy/README.md``): every ledger, policy and call of
every fixture in order, with the reference's values after each operation. The
replay asserts every one of them on the build: each call's result (placed, remaining,
termination, flags, drops, the cost curve and the requested outputs), its trace (the
reference's keys: CostOnly's prepared routes, offers and admitted totals; the controllers'
view, cap, flow-core payload and shares), the ledger after it (allocations, residuals, every
record field, next FlowId, owners) and the policy state (roster, requested, +inf, placed
demand, per source the flows' ids, routes, T, deliveries and allocations). An operation the
reference rejects must raise here too (the messages differ), and a rejected placement leaves
the policy's state as it was. Normalised on both sides: the flow core's operation counters and
ledger uids. The replay reads the trace through the private hooks, so it is skipped without
them.
"""

from __future__ import annotations

import copy
import json
import math
import traceback
from pathlib import Path

import pytest
from helpers import build, dag_of
from private_hooks import needs_hooks

import netgraph_core as ngc

DATA = Path(__file__).resolve().parents[1] / "data" / "policy_fixtures.json"
FIXTURES = json.loads(DATA.read_text())["fixtures"]

PLACEMENT = {
    "MaxFlow": ngc.FlowPlacement.MAX_FLOW,
    "Weighted": ngc.FlowPlacement.WEIGHTED,
    "EqualBalanced": ngc.FlowPlacement.EQUAL_BALANCED,
}
ADMISSION = {"Lossless": ngc.Admission.LOSSLESS, "Lossy": ngc.Admission.LOSSY}
ADMISSION_NAME = {v: k for k, v in ADMISSION.items()}
TIE = {
    "Deterministic": ngc.EdgeTieBreak.DETERMINISTIC,
    "PreferHigherPriority": ngc.EdgeTieBreak.PREFER_HIGHER_PRIORITY,
}
PRIORITY = {
    "Residual": ngc.EdgePriority.RESIDUAL,
    "Capacity": ngc.EdgePriority.CAPACITY,
    "Headroom": ngc.EdgePriority.HEADROOM,
}
WEIGHTS = {
    "Capacity": ngc.WeightBasis.CAPACITY,
    "Explicit": ngc.WeightBasis.EXPLICIT,
    "ResidualSnapshot": ngc.WeightBasis.RESIDUAL_SNAPSHOT,
}


def norm(x):
    """String map keys; probe counts, counters and ledger uids dropped; FlowIds as ids."""
    if isinstance(x, dict):
        out = {}
        for k, v in x.items():
            if k in ("ledger_uid", "counters"):
                continue
            if (
                k in ("flow_id", "flow")
                and isinstance(v, (list, tuple))
                and len(v) == 2
            ):
                out[str(k)] = v[1]
                continue
            out[str(k)] = norm(v)
        return out
    if isinstance(x, (list, tuple)):
        return [norm(v) for v in x]
    return x


def exactly(a, b):
    """Whether a and b are equal with the same type at every level: a dict has the same keys
    and exactly equal values, a list the same length and exactly equal items, and every other
    value the same type and value (so True differs from 1, and 1 from 1.0)."""
    if type(a) is not type(b):
        return False
    if isinstance(a, dict):
        return a.keys() == b.keys() and all(exactly(a[k], b[k]) for k in a)
    if isinstance(a, (list, tuple)):
        return len(a) == len(b) and all(
            exactly(x, y) for x, y in zip(a, b, strict=True)
        )
    return a == b


def build_graph(spec):
    g = build(spec["n"], spec["edges"], spec["q"])
    assert exactly(g.input_index().tolist(), spec["input_index"])
    return g


def config_of(g, c):
    r = c["routing"]
    kind = r["kind"]
    if kind == "CostOnly":
        routing = ngc.CostOnly()
    elif kind == "Global":
        support = r.get("support")
        routing = ngc.Global(support=None if support is None else dag_of(g, support))
    else:
        raise AssertionError(f"the data names an unknown routing ({kind})")
    sel = c["selection"]
    return ngc.FlowPolicyConfig(
        routing=routing,
        placement=None if c["placement"] is None else PLACEMENT[c["placement"]],
        admission=ADMISSION[c["admission"]],
        selection=None
        if sel is None
        else ngc.PathSelection(
            TIE[sel["tie_break"]],
            PRIORITY[sel["priority"]],
            sel["multipath"],
            sel["multi_edge"],
        ),
        weights=None if c.get("weights") is None else WEIGHTS[c["weights"]],
    )


def result_json(r):
    out = {
        "placed": r.placed,
        "remaining": "inf" if r.remaining == math.inf else r.remaining,
        "delivery_delta": r.delivery_delta,
        "call_drops": [[e, x] for e, x in r.call_drops],
        "call_source_loss": r.call_source_loss,
        "termination": r.termination.name.lower(),
        "routes_changed": r.routes_changed,
        "weights_changed": r.weights_changed,
        "requested": "inf" if r.requested == math.inf else r.requested,
        "admitted": r.admitted,
        "delivered": r.delivered,
        "source_loss": r.source_loss,
        "drops": [[e, x] for e, x in r.drops],
        "by_source": [
            [x.source, x.offered, x.admitted, x.delivered, x.source_loss, x.dropped]
            for x in r.by_source
        ],
    }
    for name in ("offered", "unadmitted", "total_cost"):
        value = getattr(r, name)
        if value is not None:
            out[name] = value
    if r.cost_distribution is not None:
        out["cost_distribution"] = [[c, u] for c, u in r.cost_distribution]
    for name in ("edge_flows", "residual", "bottleneck", "min_cut"):
        value = getattr(r, name)
        if value is not None:
            out[name] = [int(x) for x in value]
    return out


def ledger_json(fg):
    st = fg._state()
    return {
        "allocated": [int(x) for x in fg.allocated()],
        "residual": [int(x) for x in fg.residual()],
        "records": [
            {
                "flow_id": r.flow_id[1],
                "owner": r.owner,
                "sources": list(r.sources),
                "destinations": list(r.destinations),
                "admission": ADMISSION_NAME[r.admission],
                "bundle": sorted(int(e) for e in r.bundle.via_edges()),
                "route_cost": r.route_cost,
                "offered": r.offered,
                "delivered": r.delivered,
                "edges": [[e, x] for e, x in r.edges],
                "drops": [[e, x] for e, x in r.drops],
            }
            for r in fg.records()
        ],
        "next_flow_id": st["next_flow_id"],
        "bound": list(st["bound"]),
        "orphans": list(st["orphans"]),
    }


def policy_json(pol):
    st = pol._state()
    lossy = pol.config.admission == ngc.Admission.LOSSY
    cohorts = []
    for s, flows in sorted(st["cohorts"].items()):
        if not flows:
            continue
        out = []
        for f in flows:
            out.append(
                {
                    "id": f["flow_id"][1],
                    "route": list(f["route"]),
                    "T": None if lossy else f["T"],
                    "delivered": f["delivered"],
                    "alloc": [[e, x] for e, x in f["alloc"]],
                }
            )
        cohorts.append([s, out])
    return {
        "roster": list(st["roster"]),
        "requested": st["requested"],
        "infinite": st["infinite"],
        "placed_demand": st["placed_demand"],
        "cohorts": cohorts,
    }


def replay(fx, at=None):
    """Replays fx, asserting every expectation; `at`, when given, collects the index of each
    operation as it starts (so a failure can be attributed to its operation)."""
    graphs = [build_graph(spec) for spec in fx["graphs"]]
    ledgers: dict[int, object] = {}
    policies: dict[int, object] = {}
    for i, op in enumerate(fx["ops"]):
        if at is not None:
            at.append(i)
        where = f"{fx['name']} op {i} ({op['op']})"
        kind = op["op"]
        expect = op.get("expect", {})
        if kind == "ledger":
            # The data holds no ledger the reference rejects: a data file with one fails here,
            # rather than in a branch no fixture would reach.
            assert not op["error"], where
            ledgers[op["id"]] = ngc.FlowLedger(graphs[op["graph"]], op["baseline"])
            continue
        if kind == "policy":
            fg = ledgers[op["ledger"]]
            g = fg.graph()

            def construct(op=op, fg=fg, g=g):
                elig = None
                if op["elig"] is not None:
                    elig = g.eligibility(
                        excluded_nodes=op["elig"]["nodes"],
                        excluded_edges=op["elig"]["edges"],
                    )
                cfg = config_of(g, op["config"])
                return ngc.FlowPolicy(fg, op["sources"], op["destinations"], cfg, elig)

            if op["error"]:
                with pytest.raises((ValueError, TypeError)):
                    construct()
                continue
            pol = construct()
            policies[op["id"]] = (pol, fg)
            assert exactly(policy_json(pol), expect["policy"]), where
            continue
        pol, fg = policies[op["policy"]]
        # A removal cannot fail, so the reference must accept it.
        assert kind != "remove" or not op["error"], where
        if kind == "remove":
            pol.remove()
            assert exactly(ledger_json(fg), expect["ledger"]), where
            assert exactly(policy_json(pol), expect["policy"]), where
            continue
        assert kind == "place", where
        volume = math.inf if op["volume"] == "inf" else op["volume"]
        outputs = None
        if op["outputs"] is not None:
            e, r, b, m = op["outputs"]
            outputs = ngc.OutputOptions(
                edge_flows=e, residual=r, bottleneck=b, min_cut=m
            )
        if op["error"]:
            before = policy_json(pol)
            with pytest.raises((ValueError, TypeError)):
                pol.place(volume, outputs)
            assert exactly(ledger_json(fg), expect["ledger"]), where
            # The reference records the ledger only; the policy's own state (requested, +inf,
            # the flows) must be as before the call too.
            assert exactly(policy_json(pol), before), where
            continue
        r, trace = pol._place_traced(volume, outputs)
        assert exactly(result_json(r), expect["result"]), where
        assert exactly(norm(trace), expect["trace"]), where
        assert exactly(ledger_json(fg), expect["ledger"]), where
        assert exactly(policy_json(pol), expect["policy"]), where


def test_the_data_holds_every_fixture():
    assert len(FIXTURES) == 31
    assert len({fx["name"] for fx in FIXTURES}) == 31


@needs_hooks
@pytest.mark.parametrize("fx", FIXTURES, ids=[fx["name"] for fx in FIXTURES])
def test_fixture(fx):
    replay(fx)


def _first_leaf(x):
    """(container, key) of the first leaf of x in document order; None when x has none."""
    items = (
        x.items()
        if isinstance(x, dict)
        else enumerate(x)
        if isinstance(x, list)
        else ()
    )
    for k, v in items:
        if not isinstance(v, (dict, list)):
            return x, k
        found = _first_leaf(v)
        if found is not None:
            return found
    return None


def _corrupt(container, key):
    """Changes one leaf: None becomes 0, a bool flips, a number grows by one, a string gains a
    character."""
    v = container[key]
    if v is None:
        container[key] = 0
    elif isinstance(v, bool):
        container[key] = not v
    elif isinstance(v, (int, float)):
        container[key] = v + 1
    else:
        container[key] = v + "#"


def _failure(fx, at=None):
    """Where the replay of fx fails: (the index of the operation being replayed, None before
    the first; the failing source line), or None when the replay passes. `at`, when given, is
    the list the replay appends each operation's index to (so a control can read it)."""
    at = [] if at is None else at
    try:
        replay(fx, at)
    except (AssertionError, ValueError, TypeError, pytest.fail.Exception) as e:
        return (at[-1] if at else None), traceback.extract_tb(e.__traceback__)[-1].line
    return None


@needs_hooks
@pytest.mark.parametrize("key", ["result", "trace", "ledger", "policy"])
def test_a_corrupted_expectation_fails(key):
    # The replay is not vacuous (the twin of PolicyFixtures.ACorruptedExpectationFails): in
    # every fixture, the first integer of the expectation `key` of the first successful
    # placement that has one is corrupted; the replay must fail at that operation, on that
    # expectation's assertion.
    corrupted = 0
    for fx in FIXTURES:
        bad = copy.deepcopy(fx)
        n = next(
            (
                n
                for n, o in enumerate(bad["ops"])
                if o["op"] == "place"
                and not o["error"]
                and _corrupt_first_int(o["expect"][key])
            ),
            None,
        )
        if n is None:
            continue
        index, line = _failure(bad)
        assert index == n, fx["name"]
        assert f'expect["{key}"]' in line, fx["name"]
        corrupted += 1
    assert corrupted == len(FIXTURES)  # every fixture has such a placement


def test_exactly_compares_types_at_every_level():
    assert exactly(
        {"a": [1, [True, None]], "b": "x"}, {"a": [1, [True, None]], "b": "x"}
    )
    for a, b in (
        (True, 1),
        (1, 1.0),
        (0, False),
        ([1], (1,)),
        ({"a": [1]}, {"a": [1.0]}),
        ({"a": [[0]]}, {"a": [[False]]}),
        ({"a": 1}, {"a": 1, "b": 1}),
        ([1, 2], [1]),
    ):
        assert not exactly(a, b), (a, b)
        assert not exactly(b, a), (b, a)


@needs_hooks
@pytest.mark.parametrize("change", ["integer as float", "bool as integer"])
def test_an_expected_value_of_another_type_fails(change):
    # The replay compares types too: in every fixture, the first successful placement's
    # expected `placed` becomes the same number as a float, or its `routes_changed` the same
    # truth value as an integer; the replay must fail at that operation, on the result.
    changed = 0
    for fx in FIXTURES:
        bad = copy.deepcopy(fx)
        n = next(
            (
                n
                for n, o in enumerate(bad["ops"])
                if o["op"] == "place" and not o["error"]
            ),
            None,
        )
        if n is None:
            continue
        result = bad["ops"][n]["expect"]["result"]
        if change == "integer as float":
            result["placed"] = float(result["placed"])
        else:
            result["routes_changed"] = int(result["routes_changed"])
        assert result == fx["ops"][n]["expect"]["result"]  # equal under ==
        index, line = _failure(bad)
        assert index == n, fx["name"]
        assert 'expect["result"]' in line, fx["name"]
        changed += 1
    assert changed == len(FIXTURES)


@needs_hooks
def test_every_other_replay_clause_reports_its_operation(monkeypatch):
    # The replay's other assertions, each failing at the operation it concerns (the twin of the
    # other-clause controls inside PolicyFixtures.ACorruptedExpectationFails): the first leaf of every section
    # of every other operation's expectations, the error flag of the first operation of every
    # kind and outcome in each fixture, the input order of the first graph, an unknown
    # operation, a policy changed by a call the data rejects, and the curve a Global zero call
    # reports against the remembered one.
    sections = set()
    flags = set()
    for fx in FIXTURES:
        flipped = set()
        for n, op in enumerate(fx["ops"]):
            combination = f"{op['op']} {'failed' if op['error'] else 'ok'}"
            if "expect" in op and (op["op"] != "place" or op["error"]):
                for key in op["expect"]:
                    bad = copy.deepcopy(fx)
                    leaf = _first_leaf(bad["ops"][n]["expect"][key])
                    assert leaf is not None, (fx["name"], n, key)
                    _corrupt(*leaf)
                    index, line = _failure(bad)
                    assert index == n, (fx["name"], n, key)
                    assert f'expect["{key}"]' in line, (fx["name"], n, key)
                    sections.add(f"{combination} {key}")
            if combination not in flipped:
                flipped.add(combination)
                bad = copy.deepcopy(fx)
                bad["ops"][n]["error"] = not op["error"]
                failure = _failure(bad)
                assert failure is not None and failure[0] == n, (fx["name"], n)
                flags.add(combination)
        bad = copy.deepcopy(fx)
        _corrupt(bad["graphs"][0]["input_index"], 0)
        assert _failure(bad) == (
            None,
            'assert exactly(g.input_index().tolist(), spec["input_index"])',
        )
    assert sections == {
        "place failed ledger",
        "policy ok policy",
        "remove ok ledger",
        "remove ok policy",
    }
    assert flags == {
        "ledger ok",
        "place failed",
        "place ok",
        "policy failed",
        "policy ok",
        "remove ok",
    }
    # An unknown operation in place of op 3 of global_epochs (a placement).
    epochs = next(fx for fx in FIXTURES if fx["name"] == "global_epochs")
    renamed = copy.deepcopy(epochs)
    renamed["ops"][3]["op"] = "move"
    assert _failure(renamed) == (3, 'assert kind == "place", where')
    # Every call the data rejects, with the policy's +inf flag read flipped after that call only
    # (the binding gives no way to change a policy's state directly, so its state read is
    # changed): the comparison with the policy before the call must fail there. The ledger is
    # left as it is, so the ledger comparison before it still passes.
    place, read_state = ngc.FlowPolicy.place, ngc.FlowPolicy._state
    at, changed, target = [], set(), None

    def rejected_then_changed(self, *args):
        try:
            return place(self, *args)
        except (ValueError, TypeError):
            if at[-1] == target:
                changed.add(id(self))
            raise

    def changed_state(self):
        st = read_state(self)
        if id(self) in changed:
            st["infinite"] = not st["infinite"]
        return st

    monkeypatch.setattr(ngc.FlowPolicy, "place", rejected_then_changed)
    monkeypatch.setattr(ngc.FlowPolicy, "_state", changed_state)
    rejected_calls = 0
    for fx in FIXTURES:
        for n, op in enumerate(fx["ops"]):
            if op["op"] != "place" or not op["error"]:
                continue
            rejected_calls += 1
            at.clear()
            changed.clear()
            target = n
            assert _failure(fx, at) == (
                n,
                "assert exactly(policy_json(pol), before), where",
            ), (fx["name"], n)
    assert rejected_calls == 3
    monkeypatch.undo()
    # Op 2 of global_epochs is a Global zero call on a fresh policy (nothing prepared: no curve),
    # op 3 a positive call with the curve [[1, 1]]. With op 3 moved before op 2, the zero call
    # reports the remembered curve [[1, 1]] where the data, recorded before any placement, has
    # none: the result line of that operation fails.
    assert epochs["ops"][2]["volume"] == 0
    assert "cost_distribution" not in epochs["ops"][2]["expect"]["result"]
    swapped = copy.deepcopy(epochs)
    swapped["ops"][2], swapped["ops"][3] = swapped["ops"][3], swapped["ops"][2]
    assert _failure(swapped) == (
        3,
        'assert exactly(result_json(r), expect["result"]), where',
    )


def _corrupt_first_int(x):
    """Adds 1 to the first integer (not a bool) of x in document order; whether one was found."""
    items = (
        x.items()
        if isinstance(x, dict)
        else enumerate(x)
        if isinstance(x, list)
        else ()
    )
    for k, v in items:
        if isinstance(v, int) and not isinstance(v, bool):
            x[k] = v + 1
            return True
        if _corrupt_first_int(v):
            return True
    return False
