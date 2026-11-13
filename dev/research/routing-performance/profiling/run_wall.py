"""Unprofiled fresh user-path phase control; no graph inventory in this process."""

import argparse
import gc
import json
import sys
import time
from contextlib import contextmanager
from pathlib import Path

base = Path(__file__).resolve().parents[1]
sys.path.insert(
    0,
    str(
        "/Users/networmix/ws/work/2026-11-core-routing-research.noindex/worktrees/netgraph-orient"
    ),
)
parser = argparse.ArgumentParser()
parser.add_argument("workload")
parser.add_argument("workers", choices=("1", "2", "8", "default"))
parser.add_argument("out", type=Path)
args = parser.parse_args()
args.out.mkdir()
phases, steps = {}, []


@contextmanager
def phase(name):
    start = time.perf_counter()
    yield
    phases[name] = time.perf_counter() - start


with phase("import"):
    import yaml

    import ngraph
    from ngraph.scenario import Scenario

text = (base / f"workloads/{args.workload}.yml").read_text()
if args.workers != "default":
    document = yaml.safe_load(text)
    for step in document["workflow"]:
        if "workers" in step or step["type"] in (
            "TrafficMatrixPlacement",
            "MaxFlow",
            "TunnelCountAnalysis",
        ):
            step["workers"] = int(args.workers)
    text = yaml.safe_dump(document, sort_keys=False)
with phase("load"):
    scenario = Scenario.from_yaml(text)
with phase("prepare"):
    prepared = ngraph.prepare(scenario)


@contextmanager
def hook(step):
    start = time.perf_counter()
    yield
    steps.append({"name": step.name, "wall_s": time.perf_counter() - start})


with phase("run"):
    result = ngraph.run(prepared, step_hook=hook)
with phase("export"):
    with (args.out / "results.json").open("x") as stream:
        ngraph.write_json(result, stream)
with phase("cleanup"):
    del result, prepared, scenario
    gc.collect()
summary = {
    "exploratory": True,
    "profiled": False,
    "workload": args.workload,
    "workers_override": args.workers,
    "phases_s": phases,
    "steps": steps,
    "output_bytes": (args.out / "results.json").stat().st_size,
}
(args.out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
print(json.dumps({"summary": "unprofiled control completed", **summary}), flush=True)
