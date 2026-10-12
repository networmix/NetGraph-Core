"""Snapshot recording for prototype 2 (prototype 1's ``Snap`` plus policy state)."""

from __future__ import annotations

import _p1  # noqa: F401
from snap import Snap, jsonable
from units import ContractError


class PSnap(Snap):
    def place(self, label: str, policy, volume, fg=None, outputs=None):
        entry: dict = {"op": "place", "policy": label, "request": jsonable(volume)}
        try:
            r = (
                policy.place(volume)
                if outputs is None
                else policy.place(volume, outputs)
            )
        except ContractError as exc:
            entry["error"] = f"{type(exc).__name__}: {exc}"
            self.data["steps"].append(entry)
            return exc
        entry["result"] = jsonable(r.to_json())
        entry["trace"] = jsonable(r.trace)
        fg = fg if fg is not None else policy.fg
        entry["graph"] = self.labels.get(fg.graph.uid)
        entry["state"] = jsonable(fg.state())
        if hasattr(policy, "state"):
            entry["policy_state"] = jsonable(policy.state())
        self.data["steps"].append(entry)
        return r
