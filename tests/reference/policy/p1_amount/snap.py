"""Snapshot recording for prototype 1 fixtures.

Every fixture records its graphs (input coordinates, the B' order and
``input_index``), every operation with its arguments, its result, its trace (views,
plans, ``T`` values, ``T_max`` searches, flow-core payloads with tiers) and the
ledger state after it (allocations and residuals in input coordinates, records with
ids and input positions), and every expectation with the expected and actual values.
"""

from __future__ import annotations

import math
from fractions import Fraction

from graph import Graph, graph_from_input, input_to_id
from units import ContractError


def jsonable(x):
    if isinstance(x, bool) or x is None or isinstance(x, (int, str)):
        return x
    if isinstance(x, float):
        if math.isinf(x):
            return "inf" if x > 0 else "-inf"
        if math.isnan(x):
            return "nan"
        return x
    if isinstance(x, Fraction):
        return str(x)
    if isinstance(x, dict):
        return {str(k): jsonable(v) for k, v in x.items()}
    if isinstance(x, (list, tuple, set, frozenset)):
        seq = sorted(x) if isinstance(x, (set, frozenset)) else x
        return [jsonable(v) for v in seq]
    if hasattr(x, "to_json"):
        return jsonable(x.to_json())
    return repr(x)


class Snap:
    def __init__(self, name: str, family: str, description: str, rules) -> None:
        self.data: dict = {
            "fixture": name,
            "family": family,
            "description": description,
            "rules": list(rules),
            "graphs": {},
            "steps": [],
            "checks": [],
        }
        self.failures: list[str] = []
        self.ledgers: dict[str, object] = {}
        self.labels: dict[int, str] = {}

    # ------------------------------------------------------------ inputs

    def graph(self, label: str, n: int, edges, resolution=1) -> tuple[Graph, list[int]]:
        g = graph_from_input(n, edges, resolution)
        inv = input_to_id(g)
        self.labels[g.uid] = label
        self.data["graphs"][label] = {
            "num_nodes": n,
            "resolution": str(g.q),
            "edges_input": [list(e) for e in edges],
            "input_index": list(g.input_index()),
            "edges_by_id": [list(g.edge_tuple(e)) for e in range(g.num_edges)],
        }
        return g, inv

    def note(self, text: str, **data) -> None:
        self.data["steps"].append({"op": "note", "text": text, **jsonable(data)})

    def step(self, op: str, **data) -> None:
        self.data["steps"].append({"op": op, **jsonable(data)})

    # ------------------------------------------------------------ operations

    def place(self, label: str, policy, volume, fg=None):
        entry: dict = {"op": "place", "policy": label, "request": jsonable(volume)}
        try:
            r = policy.place(volume)
        except ContractError as exc:
            entry["error"] = f"{type(exc).__name__}: {exc}"
            self.data["steps"].append(entry)
            return exc
        entry["result"] = jsonable(r.to_json())
        entry["trace"] = jsonable(getattr(r, "trace", {}))
        fg = fg if fg is not None else getattr(policy, "fg", None)
        if fg is not None:
            entry["graph"] = self.labels.get(fg.graph.uid)
            entry["state"] = jsonable(fg.state())
        self.data["steps"].append(entry)
        return r

    def call(self, label: str, fn, *args, fg=None, **kw):
        entry: dict = {"op": label}
        try:
            out = fn(*args, **kw)
        except ContractError as exc:
            entry["error"] = f"{type(exc).__name__}: {exc}"
            self.data["steps"].append(entry)
            return exc
        entry["result"] = jsonable(out)
        if fg is not None:
            entry["state"] = jsonable(fg.state())
        self.data["steps"].append(entry)
        return out

    # ------------------------------------------------------------ expectations

    def expect(self, name: str, actual, expected) -> bool:
        ok = jsonable(actual) == jsonable(expected)
        self.data["checks"].append(
            {
                "name": name,
                "expected": jsonable(expected),
                "actual": jsonable(actual),
                "ok": ok,
            }
        )
        if not ok:
            self.failures.append(f"{name}: expected {expected!r}, got {actual!r}")
        return ok

    def check(self, name: str, ok: bool, detail=None) -> bool:
        self.data["checks"].append(
            {"name": name, "ok": bool(ok), "detail": jsonable(detail)}
        )
        if not ok:
            self.failures.append(f"{name}: {detail!r}")
        return bool(ok)

    def raises(self, name: str, fn, *args, **kw) -> bool:
        try:
            fn(*args, **kw)
        except ContractError as exc:
            self.data["checks"].append(
                {"name": name, "ok": True, "error": f"{type(exc).__name__}: {exc}"}
            )
            return True
        self.failures.append(f"{name}: no ContractError")
        self.data["checks"].append({"name": name, "ok": False, "error": None})
        return False

    def finish(self) -> dict:
        self.data["passed"] = not self.failures
        self.data["failures"] = list(self.failures)
        self.data["counts"] = {
            "steps": len(self.data["steps"]),
            "checks": len(self.data["checks"]),
        }
        return self.data


def to_input(g: Graph, per_id) -> list:
    """A per-edge array in id coordinates, re-indexed by input position."""
    out = [0] * g.num_edges
    pos = g.input_index()
    for e, x in enumerate(per_id):
        out[pos[e]] = x
    return out


def dict_to_input(g: Graph, d: dict) -> list:
    out = [0] * g.num_edges
    pos = g.input_index()
    for e, x in d.items():
        out[pos[e]] = x
    return out


def ids(inv: list[int], positions) -> list[int]:
    """Input positions to B' ids (through the caller's inverse of input_index)."""
    return sorted(inv[p] for p in positions)
