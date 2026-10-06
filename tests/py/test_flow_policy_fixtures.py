"""Prototype 2's 59 fixtures replayed through the bindings (APPROACH 8, 9).

``tests/data/policy_fixtures.json`` was exported from the reference, prototype 2 of the design
tree (its ``source`` names it; the design tree is untracked by design; its maintained exporter
is ``equivalence/export_policy_fixtures.py`` there): every ledger, policy
and call of every fixture in order, with the prototype's values after each operation. The
replay asserts every one of them on the build: each call's result (placed, remaining,
termination, local outcomes, flags, drops, the cost curve and the requested outputs), its
trace (prototype 2's keys: offers, created, per source the creation reservations, reroute
decisions, amounts and growth steps; Global's view, cap, flow-core payload and shares), the
ledger after it (allocations, residuals, every record field, next FlowId, owners) and the
policy state (roster, requested, +inf, placed demand, usable bundles, per source the flows'
ids, routes, T, deliveries, allocations, members and weights). An operation the prototype
rejects must raise here too (the messages differ). Normalised on both sides: the T_max probe
counts, the flow core's operation counters and ledger uids. The replay reads the trace and sets
the creation cap through the private hooks, so it is skipped without them.
"""

from __future__ import annotations

import gc
import json
import math
from pathlib import Path

import _netgraph_core as ext
import pytest
from private_hooks import needs_hooks

import netgraph_core as ngc

DATA = Path(__file__).resolve().parents[1] / "data" / "policy_fixtures.json"
FIXTURES = json.loads(DATA.read_text())["fixtures"]

PLACEMENT = {
    "MaxFlow": ngc.FlowPlacement.MAX_FLOW,
    "Weighted": ngc.FlowPlacement.WEIGHTED,
    "EqualBalanced": ngc.FlowPlacement.EQUAL_BALANCED,
}
PLACEMENT_NAME = {v: k for k, v in PLACEMENT.items()}
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


def norm(x):
    """String map keys; probe counts, counters and ledger uids dropped; FlowIds as ids."""
    if isinstance(x, dict):
        out = {}
        for k, v in x.items():
            if k in ("evals", "ledger_uid", "counters", "evaluations"):
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


def build_graph(spec):
    edges = spec["edges"]
    g = ngc.Graph(
        spec["n"],
        [e[0] for e in edges],
        [e[1] for e in edges],
        [e[2] for e in edges],
        [e[3] for e in edges],
        spec["q"],
    )
    assert list(g.input_index()) == spec["input_index"]
    return g


def dag_of(g, edges):
    per = [[] for _ in range(g.num_nodes())]
    for e in sorted(edges):
        per[int(g.edge_dst()[e])].append(e)
    offsets, parents, via = [0], [], []
    for v in range(g.num_nodes()):
        for e in per[v]:
            parents.append(int(g.edge_src()[e]))
            via.append(e)
        offsets.append(len(via))
    return ngc.PredDAG.from_csr(g, offsets, parents, via)


def config_of(g, c):
    r = c["routing"]
    kind = r["kind"]
    if kind == "CostOnly":
        routing = ngc.CostOnly()
    elif kind == "Global":
        routing = ngc.Global()
    elif kind == "Pinned":
        routing = ngc.Pinned([dag_of(g, b) for b in r["bundles"]])
    else:
        ceiling = None
        if r["max_path_cost"] is not None:
            name, value = r["max_path_cost"]
            ceiling = ngc.Absolute(value) if name == "Absolute" else ngc.Stretch(value)
        dim = (
            None
            if r["diminishing"] is None
            else ngc.DiminishingReturns(*r["diminishing"])
        )
        routing = ngc.CapacityAware(
            flow_count=r["flow_count"],
            max_flow_count=r["max_flow_count"],
            max_path_cost=ceiling,
            min_flow_bandwidth=r["min_flow_bandwidth"],
            reoptimize=r["reoptimize"],
            diminishing=dim,
        )
    sel = c["selection"]
    return ngc.FlowPolicyConfig(
        routing=routing,
        placement=PLACEMENT[c["placement"]],
        admission=ADMISSION[c["admission"]],
        selection=ngc.PathSelection(
            TIE[sel["tie_break"]],
            PRIORITY[sel["priority"]],
            sel["multipath"],
            sel["multi_edge"],
        ),
    )


def result_json(r):
    out = {
        "placed": r.placed,
        "remaining": "inf" if r.remaining == math.inf else r.remaining,
        "termination": r.termination.name.lower(),
        "routes_changed": r.routes_changed,
        "weights_changed": r.weights_changed,
        "drops": [[e, x] for e, x in r.drops],
        "local_outcomes": [[s, o.name.lower()] for s, o in r.local_outcomes],
    }
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
                "placement": PLACEMENT_NAME[r.placement],
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
            members = list(f["route"] if f["members"] is None else f["members"])
            weights = None
            if f["weights"] is not None:
                weights = [[e, w] for e, w in zip(members, f["weights"], strict=True)]
            out.append(
                {
                    "id": f["flow_id"][1],
                    "route": list(f["route"]),
                    "T": None if lossy else f["T"],
                    "delivered": f["delivered"],
                    "alloc": [[e, x] for e, x in f["alloc"]],
                    "members": members,
                    "weights": weights,
                }
            )
        cohorts.append([s, out])
    return {
        "roster": list(st["roster"]),
        "requested": st["requested"],
        "infinite": st["infinite"],
        "placed_demand": st["placed_demand"],
        "usable": pol.usable_bundles(),
        "cohorts": cohorts,
    }


def replay(fx):
    graphs = [build_graph(spec) for spec in fx["graphs"]]
    ledgers: dict[int, object] = {}
    policies: dict[int, object] = {}
    alive = []
    curves: dict[
        int, list
    ] = {}  # per policy: the cumulative Global curve of the epoch (F5)
    for i, op in enumerate(fx["ops"]):
        where = f"{fx['name']} op {i} ({op['op']})"
        kind = op["op"]
        expect = op.get("expect", {})
        if kind == "ledger":
            g = graphs[op["graph"]]
            if op["error"]:
                with pytest.raises((ValueError, TypeError)):
                    ngc.FlowLedger(g, op["baseline"])
                continue
            ledgers[op["id"]] = ngc.FlowLedger(g, op["baseline"])
            continue
        if kind == "copy":
            src = ledgers[op["ledger"]]
            if op["error"]:
                with pytest.raises((ValueError, TypeError)):
                    src.copy()
                assert ledger_json(src) == expect["ledger"], where
                continue
            ledgers[op["id"]] = src.copy()
            assert ledger_json(ledgers[op["id"]]) == expect["ledger"], where
            continue
        if kind == "reset":
            fg = ledgers[op["ledger"]]
            if op["error"]:
                with pytest.raises((ValueError, TypeError)):
                    fg.reset(op["baseline"])
            else:
                fg.reset(op["baseline"])
            assert ledger_json(fg) == expect["ledger"], where
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
            alive.append(pol)
            policies[op["id"]] = (pol, fg)
            curves[op["id"]] = []
            assert policy_json(pol) == expect["policy"], where
            continue
        pol, fg = policies[op["policy"]]
        if kind == "remove":
            pol.remove()
            curves[op["policy"]] = []
            assert ledger_json(fg) == expect["ledger"], where
            assert policy_json(pol) == expect["policy"], where
            continue
        if kind == "destroy":
            # The destruction law: dropping the last reference orphans the records.
            alive[:] = [x for x in alive if x is not pol]
            policies[op["policy"]] = None
            del pol
            gc.collect()
            assert ledger_json(fg) == expect["ledger"], where
            continue
        assert kind == "place", where
        volume = math.inf if op["volume"] == "inf" else op["volume"]
        outputs = None
        if op["outputs"] is not None:
            e, r, b = op["outputs"]
            outputs = ngc.OutputOptions(edge_flows=e, residual=r, bottleneck=b)
        old = ext._max_flows_per_call(op["cap"])
        try:
            if op["error"]:
                with pytest.raises((ValueError, TypeError)):
                    pol.place(volume, outputs)
                assert ledger_json(fg) == expect["ledger"], where
                continue
            r, trace = pol._place_traced(volume, outputs)
        finally:
            ext._max_flows_per_call(old)
        got = result_json(r)
        if isinstance(pol.config.routing, ngc.Global):
            # Contract over prototype (1c NOTES P3, review-1c F5): a Global zero call reports
            # the cumulative curve (9.1), which prototype 2 omits; it must equal the curve of the
            # epoch's last positive call ([] before any), and only then is it set aside.
            if volume == 0 and "cost_distribution" not in expect["result"]:
                assert got.pop("cost_distribution") == curves[op["policy"]], where
            if "cost_distribution" in expect["result"]:
                curves[op["policy"]] = expect["result"]["cost_distribution"]
        assert got == expect["result"], where
        assert norm(trace) == expect["trace"], where
        assert ledger_json(fg) == expect["ledger"], where
        assert policy_json(pol) == expect["policy"], where


def test_the_data_holds_every_fixture():
    assert len(FIXTURES) == 59
    assert len({fx["name"] for fx in FIXTURES}) == 59


@needs_hooks
@pytest.mark.parametrize("fx", FIXTURES, ids=[fx["name"] for fx in FIXTURES])
def test_fixture(fx):
    replay(fx)
