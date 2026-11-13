"""Conservative semantic call fingerprints, including masks and ledger state.

Graph identity is retained; policy/ledger allocation identities are normalized.
An already prepared program is retained by identity, because its split-DAG state
is not completely exposed. Thus repeats can be undercounted, never merged just
because public program fields happen to agree. Fingerprinting is outside timers.
"""

import hashlib
import json
import threading
import time
from collections import Counter, defaultdict

import numpy as np


def digest(value):
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


class Tracer:
    def __init__(self, core, out, stage):
        self.core, self.out, self.stage = core, out, stage
        self.policies = {}
        self.programs = {}
        self.calls = []
        self.lock = threading.Lock()

    def value(self, x):
        if x is None or isinstance(x, (str, int, float, bool)):
            return x
        if isinstance(x, np.ndarray):
            return {
                "dtype": str(x.dtype),
                "shape": list(x.shape),
                "sha256": hashlib.sha256(x.tobytes()).hexdigest(),
            }
        if isinstance(x, (tuple, list)):
            return [self.value(y) for y in x]
        if isinstance(x, dict):
            return {k: self.value(v) for k, v in x.items()}
        if isinstance(x, self.core.Graph):
            return {"graph_uid": x.uid()}
        if isinstance(x, self.core.FlowLedger):
            return {
                "graph": self.value(x.graph()),
                "residual": self.value(x.residual()),
                "reserved": self.value(x.reserved()),
                "unreserved": self.value(x.unreserved()),
                "carried": self.value(x.carried()),
                "live_records": len(x.records()),
            }
        if isinstance(x, self.core.Eligibility):
            return {
                "graph_uid": x.graph_uid(),
                "nodes": self.value(x.node_allowed()),
                "edges": self.value(x.edge_allowed()),
            }
        if isinstance(x, self.core.ForwardingProgram):
            self.programs[id(x)] = x
            return {"prepared_program_identity": id(x)}
        if isinstance(x, self.core.PredecessorDAG):
            return {
                "graph_uid": x.graph_uid(),
                "offsets": self.value(x.parent_offsets()),
                "parents": self.value(x.parents()),
                "via": self.value(x.via_edges()),
            }
        if isinstance(x, self.core.Path):
            return {
                "graph_uid": x.graph_uid(),
                "hops": self.value(x.hops()),
                "source": x.source(),
                "destination": x.destination(),
            }
        if hasattr(x, "name") and hasattr(x, "value"):
            return {"enum": type(x).__name__, "name": x.name}
        fields = {
            a: self.value(getattr(x, a))
            for a in dir(x)
            if not a.startswith("_") and not callable(getattr(x, a))
        }
        if not fields:
            raise TypeError(f"uncovered native input {type(x).__name__}")
        return {"type": type(x).__name__, **fields}

    def note(self, operation, full, roots, masks, seconds, units=1):
        record = {
            "operation": operation,
            "stage": self.stage(),
            "full": digest(full),
            "roots": digest(roots),
            "masks": digest(masks),
            "native_s": seconds,
            "units": units,
            "input": full,
        }
        with self.lock:
            self.calls.append(record)

    def install(self):
        import ngraph.analysis.msd as msd

        initialize = self.core.FlowPolicy.__init__

        def policy_init(
            policy, ledger, sources, destinations, config, eligibility=None, cache=None
        ):
            full = {
                "ledger": self.value(ledger),
                "sources": self.value(sources),
                "destinations": self.value(destinations),
                "config": self.value(config),
                "eligibility": self.value(eligibility),
                "cache": self.value(cache),
            }
            roots = [
                self.value(ledger.graph()),
                self.value(sources),
                self.value(destinations),
            ]
            masks = [self.value(ledger.graph()), self.value(eligibility)]
            start = time.perf_counter()
            initialize(
                policy, ledger, sources, destinations, config, eligibility, cache
            )
            elapsed = time.perf_counter() - start
            self.note("FlowPolicy.__init__", full, roots, masks, elapsed)

        self.core.FlowPolicy.__init__ = policy_init
        import ngraph.analysis.placement as placement

        roster = placement.policy_roster

        def policy_roster(topology, demands, mask, entries=None):
            ledger, policies = roster(topology, demands, mask, entries)
            if entries is None:
                entries = [(d.config, mask, None) for d in demands]
            for demand, policy, (config, eligible, _) in zip(
                demands, policies, entries, strict=True
            ):
                sources, targets = placement._endpoints(topology, demand)
                self.policies[policy] = {
                    "sources": sources,
                    "destinations": targets,
                    "config": self.value(config),
                    "eligibility": self.value(eligible),
                }
            return ledger, policies

        placement.policy_roster = msd.policy_roster = policy_roster
        apply = self.core.apply_actions

        def apply_actions(ledger, actions, *args, **kwargs):
            entries = []
            for action in actions:
                meta = self.policies[action.policy]
                entries.append(
                    {
                        **meta,
                        "action": type(action).__name__,
                        "offer": self.value(action.offer),
                        "design_demand": self.value(action.design_demand),
                        "program": self.value(action.policy.program),
                        "delivered_before": action.policy.delivered(),
                    }
                )
            full = {
                "graph": self.value(ledger.graph()),
                "actions": entries,
                "residual": self.value(ledger.residual()),
                "reserved": self.value(ledger.reserved()),
                "unreserved": self.value(ledger.unreserved()),
                "carried": self.value(ledger.carried()),
                "args": self.value(args),
                "kwargs": self.value(kwargs),
            }
            roots = [
                full["graph"],
                [(e["sources"], e["destinations"]) for e in entries],
            ]
            masks = [full["graph"], [e["eligibility"] for e in entries]]
            start = time.perf_counter()
            result = apply(ledger, actions, *args, **kwargs)
            elapsed = time.perf_counter() - start
            self.note("apply_actions", full, roots, masks, elapsed, len(actions))
            return result

        self.core.apply_actions = apply_actions
        for name in (
            "evaluate_demand",
            "evaluate_pairs",
            "evaluate_edge_failures",
            "spf",
            "ksp",
            "resolve_to_paths",
        ):
            original = getattr(self.core, name)

            def wrapper(*args, _name=name, _original=original, **kwargs):
                full = {"args": self.value(args), "kwargs": self.value(kwargs)}
                # Group the graph and endpoint arguments independently of the options.
                endpoint_args = (
                    args[1:2] if _name in ("evaluate_pairs", "spf") else args[1:3]
                )
                roots = [self.value(args[0]), self.value(endpoint_args)]
                opts = next(
                    (a for a in reversed(args) if hasattr(a, "eligibility")), None
                )
                masks = [
                    self.value(args[0]),
                    self.value(getattr(opts, "eligibility", None)),
                ]
                start = time.perf_counter()
                result = _original(*args, **kwargs)
                elapsed = time.perf_counter() - start
                self.note(
                    _name,
                    full,
                    roots,
                    masks,
                    elapsed,
                    len(args[1]) if _name == "evaluate_pairs" else 1,
                )
                return result

            setattr(self.core, name, wrapper)

    def finish(self):
        groups = defaultdict(list)
        for call in self.calls:
            groups[call["operation"]].append(call)
        summary = {}
        for operation, calls in groups.items():
            seen = set()
            repeated_s = 0.0
            roots_masks = defaultdict(set)
            masks_roots = defaultdict(set)
            for c in calls:
                if c["full"] in seen:
                    repeated_s += c["native_s"]
                seen.add(c["full"])
                roots_masks[c["roots"]].add(c["masks"])
                masks_roots[c["masks"]].add(c["roots"])
            total = sum(c["native_s"] for c in calls)
            summary[operation] = {
                "calls": len(calls),
                "pair_or_action_units": sum(c["units"] for c in calls),
                "distinct_inputs": len(seen),
                "repeat_calls": len(calls) - len(seen),
                "native_s": total,
                "repeated_native_s": repeated_s,
                "repeat_native_time_share": repeated_s / total if total else 0,
                "same_roots_multiple_masks_groups": sum(
                    len(s) > 1 for s in roots_masks.values()
                ),
                "same_mask_multiple_roots_groups": sum(
                    len(s) > 1 for s in masks_roots.values()
                ),
                "calls_in_same_roots_multiple_masks_groups": sum(
                    len(roots_masks[c["roots"]]) > 1 for c in calls
                ),
                "calls_in_same_mask_multiple_roots_groups": sum(
                    len(masks_roots[c["masks"]]) > 1 for c in calls
                ),
                "stages": dict(Counter(c["stage"] for c in calls)),
            }
        with (self.out / "calls.jsonl").open("x") as stream:
            for call in self.calls:
                stream.write(json.dumps(call) + "\n")
            stream.write(
                json.dumps(
                    {"summary": "native trace completed", "calls": len(self.calls)}
                )
                + "\n"
            )
        (self.out / "repeats.json").write_text(json.dumps(summary, indent=2) + "\n")
