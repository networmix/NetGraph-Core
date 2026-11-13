"""Time native API spans without cProfile or fingerprinting (one worker).

Core/NetGraph are imported before this diagnostic: its import phase is not a
cold import measurement. All user run/load/prepare/export phases are fresh.
"""

import json
import runpy
import sys
import time
from collections import defaultdict
from pathlib import Path

base = Path(__file__).resolve().parents[1]
sys.path.insert(
    0,
    str(
        "/Users/networmix/ws/work/2026-11-core-routing-research.noindex/worktrees/netgraph-orient"
    ),
)
import netgraph_core as core  # noqa: E402

import ngraph  # noqa: E402

stage = "load"
timers = defaultdict(lambda: [0, 0.0])


def timed(name, original):
    def wrapper(*args, **kwargs):
        start = time.perf_counter()
        try:
            return original(*args, **kwargs)
        finally:
            elapsed = time.perf_counter() - start
            row = timers[(stage, name)]
            row[0] += 1
            row[1] += elapsed

    return wrapper


for name in (
    "apply_actions",
    "evaluate_pairs",
    "evaluate_demand",
    "evaluate_edge_failures",
    "spf",
    "ksp",
    "resolve_to_paths",
    "to_units",
):
    setattr(core, name, timed(name, getattr(core, name)))
for name in (
    "Graph",
    "FlowLedger",
    "FlowPolicy",
    "OutputOptions",
    "EvaluationOptions",
    "SpfOptions",
    "FlowPolicyConfig",
    "ShortestPaths",
    "DistributedTE",
    "CentralizedTE",
    "Replace",
    "Forward",
    "Pinned",
    "TunnelCount",
    "AbsoluteCost",
    "CostStretch",
):
    kind = getattr(core, name)
    kind.__init__ = timed(name + ".__init__", kind.__init__)
core.Graph.eligibility = timed("Graph.eligibility", core.Graph.eligibility)
core.FlowPolicy.remove = timed("FlowPolicy.remove", core.FlowPolicy.remove)
core.Eligibility.intersect = staticmethod(
    timed("Eligibility.intersect", core.Eligibility.intersect)
)


def stage_call(name, original):
    def wrapper(*args, **kwargs):
        global stage
        stage = name
        try:
            return original(*args, **kwargs)
        finally:
            stage = "export"

    return wrapper


ngraph.prepare = stage_call("prepare", ngraph.prepare)
ngraph.run = stage_call("run", ngraph.run)
out = Path(sys.argv[3])
assert sys.argv[2] == "1", "aggregate native spans must not overlap Python workers"
runpy.run_path(str(base / "scripts/run_wall.py"), run_name="__main__")
rows = [
    {"phase": phase, "operation": name, "calls": row[0], "native_s": row[1]}
    for (phase, name), row in sorted(timers.items())
]
(out / "native-timers.json").write_text(json.dumps(rows, indent=2) + "\n")
print(json.dumps({"summary": "native API timers completed", "operations": len(rows)}))
