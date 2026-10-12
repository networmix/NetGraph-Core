"""Export prototype 2's 59 fixtures as replay data for the Core tests (sub-PR 1c).

    venv/bin/python -I tests/reference/policy/export_policy_fixtures.py [--check] [--out PATH]

Relocated copy (cleanup batch K2, decision D1): this file and the prototype modules beside it in
``p2_transaction/`` and ``p1_amount/`` are a frozen copy of the design tree's exporter and of
its import closure; ``README.md`` here lists the original paths, hashes and every relocation
change. The text below describes the original layout.

The maintained exporter of the policy fixture data (task 1c-grooming, brief item 8): the saved
exporter ``gates/1c/export_fixture_data.py`` wrote ``tests/data/p2_fixtures.json`` and stays as
the 1c gate's history; this copy writes the renamed file ``tests/data/policy_fixtures.json``
with the same derivation and the same bytes. The data's own ``source`` string still names the
saved exporter that first wrote it: it is part of the content, which the rename keeps byte for
byte (a later change of that string is a deliberate data change, not a rename).

Every fixture of ``prototypes/p2_transaction/fixtures.py`` runs on the PROTOTYPE (the
reference; never the build) with a recording shim installed before ``fixtures`` is
imported: ``ledger.FlowGraph`` and ``policy.FlowPolicy`` are replaced by subclasses that log
every construction and call, and ``policy.MAX_FLOWS_PER_CALL`` (the fixtures' test hook) is
read at every placement. The output ``tests/data/policy_fixtures.json`` (compact, deterministic)
holds per fixture:

* ``graphs``: ``n``, ``q`` (the resolution), the edges in input coordinates
  ``[u, v, capacity units, cost]`` and the B' ``input_index`` (the replay checks that the
  build ranks the edges the same way);
* ``ops`` in call order: ``ledger`` (graph, baseline in B' id order or null), ``copy``,
  ``reset``, ``policy`` (ledger, sources and destinations as given, the full configuration
  with Pinned bundles as ascending B' edge lists, the Eligibility as excluded node and edge
  ids, or ``error``), ``place`` (volume as an int or ``"inf"``, the outputs flags, the
  creation cap ``cap``), ``remove`` and ``destroy``;
* after every operation the prototype's expected values (``expect``): the call's result
  (``to_json``: placed, remaining, termination, local outcomes, flags, drops, the cost curve
  and the requested outputs), the ledger state (allocated and residual in id order, every
  record field with the bundle as an ascending edge list, next FlowId, bound and orphaned
  owners), the policy state (roster, requested, +inf, placed demand, usable bundles, per
  source the flows' id, route, T (null for Lossy: prototype 2 keeps none), delivered, alloc,
  members and weights) and the normalised trace (prototype 2's keys; the probe counts
  ``evals``, the flow core's ``counters`` and ledger uids dropped; FlowIds as ids).

An operation the prototype rejects is recorded with ``"error": true`` (the build must raise
too; the messages differ) and the ledger state after it (unchanged). ``--check`` recomputes
the data and compares it with the file.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
P2 = HERE / "p2_transaction"
OUT = REPO / "tests" / "data" / "policy_fixtures.json"
# The modules of the reference: each must come from this bundle (p2_transaction/ or
# p1_amount/), never from an installed package of the same name.
CLOSURE = (
    ("_p1", "fixtures", "policy", "recorder", "run_fixtures"),
    (
        "apportion",
        "certify",
        "flowcore",
        "graph",
        "kernels",
        "ledger",
        "search",
        "snap",
        "units",
    ),
)

# The original directory held the old-build recorder.py, and prototype 2 has its own recorder
# module (its fixtures import PSnap from it): the directory leaves sys.path before the prototypes
# are imported, so no module of the tooling can shadow one of the reference's (kept in the copy).
sys.path[:] = [p for p in sys.path if Path(p or ".").resolve() != HERE]
sys.path.insert(0, str(P2))

import _p1  # noqa: E402, F401  (prototype 1 appended to sys.path)
import ledger  # noqa: E402
from snap import jsonable  # noqa: E402
from units import ContractError  # noqa: E402

BaseFlowGraph = ledger.FlowGraph


class Recorder:
    def __init__(self) -> None:
        self.graphs: list[dict] = []
        self.graph_ids: dict[int, int] = {}
        self.ledgers: list[object] = []
        self.policies: list[object] = []
        self.ops: list[dict] = []

    def graph_id(self, g) -> int:
        got = self.graph_ids.get(g.uid)
        if got is not None:
            return got
        pos = g.input_index()
        edges = [None] * g.num_edges
        for e in range(g.num_edges):
            edges[pos[e]] = [g.src[e], g.dst[e], g.capacity[e], g.cost[e]]
        self.graphs.append(
            {
                "n": g.num_nodes,
                "q": float(g.q),
                "edges": edges,
                "input_index": list(pos),
            }
        )
        self.graph_ids[g.uid] = len(self.graphs) - 1
        return len(self.graphs) - 1

    def ledger_id(self, fg) -> int:
        for i, x in enumerate(self.ledgers):
            if x is fg:
                return i
        raise AssertionError("a ledger the recorder did not see")

    def policy_id(self, pol) -> int:
        for i, x in enumerate(self.policies):
            if x is pol:
                return i
        raise AssertionError("a policy the recorder did not see")


REC: Recorder | None = None


# ------------------------------------------------------------------ normalisation


def norm(x):
    """Drop process-dependent and probe-count fields; FlowIds to their id."""
    if isinstance(x, dict):
        out = {}
        for k, v in x.items():
            if k in ("evals", "ledger_uid", "counters", "evaluations"):
                continue
            if k in ("flow_id", "flow") and isinstance(v, list) and len(v) == 2:
                out[k] = v[1]
                continue
            out[k] = norm(v)
        return out
    if isinstance(x, list):
        return [norm(v) for v in x]
    return x


def ledger_state(fg) -> dict:
    return {
        "allocated": list(fg.allocated()),
        "residual": list(fg.residual()),
        "records": [
            {
                "flow_id": r.flow_id[1],
                "owner": r.owner,
                "sources": list(r.sources),
                "destinations": list(r.destinations),
                "placement": r.placement,
                "admission": r.admission,
                "bundle": sorted(r.bundle),
                "route_cost": None
                if r.route_cost == sys.modules["units"].UNREACHABLE
                else r.route_cost,
                "offered": r.offered,
                "delivered": r.delivered,
                "edges": [[e, x] for e, x in sorted(r.edges.items()) if x],
                "drops": [[e, x] for e, x in sorted(r.drops.items()) if x],
            }
            for r in fg.records()
        ],
        "next_flow_id": fg.next_flow_id(),
        "bound": sorted(fg._bound),
        "orphans": sorted(fg._orphans),
    }


def policy_state(pol) -> dict:
    lossy = pol.admission == "Lossy"
    cohorts = []
    for s in sorted(pol.cohorts):
        flows = pol.cohorts[s]
        if not flows:
            continue
        out = []
        for f in flows:
            w = f.plan.weights
            out.append(
                {
                    "id": f.flow_id[1],
                    "route": list(f.route),
                    "T": None if lossy else f.T,
                    "delivered": f.delivered,
                    "alloc": [[e, x] for e, x in sorted(f.alloc.items()) if x],
                    "members": list(f.plan.members),
                    "weights": None
                    if w is None
                    else [[e, x] for e, x in sorted(w.items())],
                }
            )
        cohorts.append([s, out])
    return {
        "roster": list(pol.roster),
        "requested": pol.requested,
        "infinite": pol.infinite,
        "placed_demand": pol.placed_demand(),
        "usable": pol.usable_bundles(),
        "cohorts": cohorts,
    }


def config_json(cfg) -> dict:
    r = cfg.routing
    kind = type(r).__name__
    routing: dict = {"kind": kind}
    if kind == "Pinned":
        routing["bundles"] = [sorted(b) for b in r.bundles]
    elif kind == "CapacityAware":
        routing.update(
            flow_count=r.flow_count,
            max_flow_count=r.max_flow_count,
            max_path_cost=None
            if r.max_path_cost is None
            else [
                r.max_path_cost[0],
                float(r.max_path_cost[1])
                if r.max_path_cost[0] == "Stretch"
                else r.max_path_cost[1],
            ],
            min_flow_bandwidth=r.min_flow_bandwidth,
            reoptimize=r.reoptimize,
            diminishing=None
            if r.diminishing is None
            else [r.diminishing.window, float(r.diminishing.fraction)],
        )
    sel = cfg.selection
    return {
        "routing": routing,
        "placement": cfg.placement,
        "admission": cfg.admission,
        "selection": {
            "multipath": sel.multipath,
            "multi_edge": sel.multi_edge,
            "tie_break": sel.tie_break,
            "priority": sel.priority,
        },
    }


def elig_json(elig):
    if elig is None:
        return None
    return {
        "nodes": [v for v, ok in enumerate(elig.node_allowed) if not ok],
        "edges": [e for e, ok in enumerate(elig.edge_allowed) if not ok],
    }


# ------------------------------------------------------------------ the shim


class RecFlowGraph(BaseFlowGraph):
    def __init__(self, g, residual=None) -> None:
        op = {
            "op": "ledger",
            "graph": REC.graph_id(g),
            "baseline": None if residual is None else list(residual),
        }
        try:
            super().__init__(g, residual)
        except ContractError:
            REC.ops.append({**op, "error": True})
            raise
        REC.ledgers.append(self)
        REC.ops.append({**op, "id": len(REC.ledgers) - 1, "error": False})

    def copy(self):
        op = {"op": "copy", "ledger": REC.ledger_id(self)}
        try:
            other = super().copy()
        except ContractError:
            REC.ops.append(
                {**op, "error": True, "expect": {"ledger": ledger_state(self)}}
            )
            raise
        REC.ledgers.append(other)
        REC.ops.append(
            {
                **op,
                "id": len(REC.ledgers) - 1,
                "error": False,
                "expect": {"ledger": ledger_state(other)},
            }
        )
        return other

    def reset(self, residual=None) -> None:
        op = {
            "op": "reset",
            "ledger": REC.ledger_id(self),
            "baseline": None if residual is None else list(residual),
        }
        try:
            super().reset(residual)
        except ContractError:
            REC.ops.append(
                {**op, "error": True, "expect": {"ledger": ledger_state(self)}}
            )
            raise
        REC.ops.append({**op, "error": False, "expect": {"ledger": ledger_state(self)}})


def make_rec_policy(P):
    Base = P.FlowPolicy

    class RecPolicy(Base):
        def __init__(
            self, fg, sources, destinations, cfg, elig=None, cache=None
        ) -> None:
            op = {
                "op": "policy",
                "ledger": REC.ledger_id(fg),
                "sources": jsonable(sources),
                "destinations": jsonable(destinations),
                "config": config_json(cfg),
                "elig": elig_json(elig),
            }
            try:
                super().__init__(fg, sources, destinations, cfg, elig, cache)
            except ContractError:
                REC.ops.append({**op, "error": True})
                raise
            REC.policies.append(self)
            REC.ops.append(
                {
                    **op,
                    "id": len(REC.policies) - 1,
                    "error": False,
                    "expect": {"policy": policy_state(self)},
                }
            )

        def place(self, volume, outputs=None):
            op = {
                "op": "place",
                "policy": REC.policy_id(self),
                "volume": "inf"
                if isinstance(volume, float) and math.isinf(volume) and volume > 0
                else jsonable(volume),
                "outputs": None
                if outputs is None
                else [outputs.edge_flows, outputs.residual, outputs.bottleneck],
                "cap": P.MAX_FLOWS_PER_CALL,
            }
            try:
                r = super().place(volume, outputs)
            except ContractError:
                REC.ops.append(
                    {**op, "error": True, "expect": {"ledger": ledger_state(self.fg)}}
                )
                raise
            REC.ops.append(
                {
                    **op,
                    "error": False,
                    "expect": {
                        "result": jsonable(r.to_json()),
                        "trace": norm(jsonable(r.trace)),
                        "ledger": ledger_state(self.fg),
                        "policy": policy_state(self),
                    },
                }
            )
            return r

        def remove(self) -> None:
            super().remove()
            REC.ops.append(
                {
                    "op": "remove",
                    "policy": REC.policy_id(self),
                    "error": False,
                    "expect": {
                        "ledger": ledger_state(self.fg),
                        "policy": policy_state(self),
                    },
                }
            )

        def destroy(self) -> None:
            super().destroy()
            REC.ops.append(
                {
                    "op": "destroy",
                    "policy": REC.policy_id(self),
                    "error": False,
                    "expect": {"ledger": ledger_state(self.fg)},
                }
            )

    return RecPolicy


def export() -> dict:
    global REC
    ledger.FlowGraph = RecFlowGraph
    import policy as P

    P.FlowGraph = RecFlowGraph
    P.FlowPolicy = make_rec_policy(P)
    from fixtures import FIXTURES
    from run_fixtures import run_one

    for directory, names in zip((P2, HERE / "p1_amount"), CLOSURE, strict=True):
        for name in names:
            origin = Path(sys.modules[name].__file__).resolve()
            if origin != directory / f"{name}.py":
                raise AssertionError(f"{name} came from {origin}, not the prototypes")
    if any(m.split(".")[0] in ("netgraph_core", "_netgraph_core") for m in sys.modules):
        raise AssertionError("the reference imported the build")
    out = []
    for fx in FIXTURES:
        REC = Recorder()
        data = run_one(fx)
        if not data["passed"]:
            raise AssertionError(f"{fx.name}: the prototype fails its own fixture")
        out.append({"name": fx.name, "graphs": REC.graphs, "ops": REC.ops})
    return {
        "format": "p2-fixtures/1",
        "source": "prototype 2 of the design tree dev/design/2026-10-streamline (untracked by "
        "design): prototypes/p2_transaction in integer units (APPROACH 9), exported by "
        "gates/1c/export_fixture_data.py",
        "fixtures": out,
    }


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--check", action="store_true")
    ap.add_argument("--out", type=Path, default=OUT)
    args = ap.parse_args(argv)
    text = json.dumps(export(), separators=(",", ":"), sort_keys=False) + "\n"
    if args.check:
        same = args.out.exists() and args.out.read_text() == text
        print("policy fixture data:", "identical" if same else "DIFFERS")
        return 0 if same else 1
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(text)
    data = json.loads(text)
    n_ops = sum(len(f["ops"]) for f in data["fixtures"])
    print(
        f"wrote {args.out}: {len(data['fixtures'])} fixtures, {n_ops} operations, {len(text)} bytes"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
