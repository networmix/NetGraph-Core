"""Fresh-process cProfile and phase wall clocks (CPython 3.14 all-thread capture)."""

import argparse
import cProfile
import gc
import hashlib
import json
import os
import platform
import resource
import sys
import time
from collections import Counter
from contextlib import contextmanager
from importlib.metadata import version
from pathlib import Path

BASE = Path(__file__).resolve().parents[1]
sys.path.insert(
    0,
    str(
        "/Users/networmix/ws/work/2026-11-core-routing-research.noindex/worktrees/netgraph-orient"
    ),
)
parser = argparse.ArgumentParser()
parser.add_argument("workload")
parser.add_argument("workers", choices=("default", "1"))
parser.add_argument("out", type=Path)
parser.add_argument("--trace", action="store_true")
args = parser.parse_args()
args.out.mkdir()
phases = {}
profiles = {}
worker_profiles = {}
stage = "import"


@contextmanager
def phase(name):
    global stage
    stage = name
    profile = cProfile.Profile()
    start = time.perf_counter()
    profile.enable()
    try:
        yield
    finally:
        profile.disable()
        phases[name] = time.perf_counter() - start
        profile.dump_stats(str(args.out / f"{name}.pstats"))
        profiles[name] = profile


started = time.perf_counter()
with phase("import"):
    import netgraph_core as core
    import yaml

    import ngraph
    from ngraph.scenario import Scenario

assert Path(ngraph.__file__).is_relative_to(
    "/Users/networmix/ws/work/2026-11-core-routing-research.noindex/worktrees/netgraph-orient"
)
manifest = json.loads((BASE / "workloads/manifest.json").read_text())
input_path = Path(manifest[args.workload]["path"])
if args.workers != "default":
    # YAML override is research harness work, outside the measured load phase.
    document = yaml.safe_load(input_path.read_text())
    for step in document["workflow"]:
        if "workers" in step or step["type"] in (
            "TrafficMatrixPlacement",
            "MaxFlow",
            "TunnelCountAnalysis",
        ):
            step["workers"] = int(args.workers)
    text = yaml.safe_dump(document, sort_keys=False)
else:
    text = input_path.read_text()

tracer = None
if args.trace:
    from trace_calls import Tracer

    tracer = Tracer(core, args.out, lambda: stage)
    tracer.install()

with phase("load"):
    scenario = Scenario.from_yaml(text)
del text
with phase("prepare"):
    prepared = ngraph.prepare(scenario)

steps = []


@contextmanager
def step_hook(step):
    global stage
    stage = step.name
    start = time.perf_counter()
    yield
    steps.append({"name": step.name, "wall_s": time.perf_counter() - start})
    stage = "run"


with phase("run"):
    result = ngraph.run(prepared, step_hook=step_hook)
with phase("export"):
    with (args.out / "results.json").open("x") as stream:
        ngraph.write_json(result, stream)

# Inventory is intentionally outside the user-path clocks/profiles.
topology = prepared.topology
graph = topology.graph
links = scenario.network.links
pairs = Counter(tuple(sorted((link.source, link.target))) for link in links.values())
costs = Counter(link.cost for link in links.values())
capacities = [link.capacity for link in links.values()]
import math  # noqa: E402 -- inventory imports are outside the user-path clocks

import networkx as nx  # noqa: E402

physical = nx.Graph()
physical.add_nodes_from(scenario.network.nodes)
physical.add_edges_from(pairs)
diameters = [
    nx.diameter(physical.subgraph(nodes)) for nodes in nx.connected_components(physical)
]
inventory = {
    "nodes": len(scenario.network.nodes),
    "links": len(links),
    "native_nodes": graph.num_nodes(),
    "native_edges": graph.num_edges(),
    "parallel_multiplicity_histogram": dict(Counter(pairs.values())),
    "cost_histogram": dict(costs),
    "cost_min": min(costs),
    "cost_max": max(costs),
    "cost_gcd": math.gcd(*costs),
    "unit_cost_link_share": costs[1] / len(links),
    "capacity_min": min(capacities),
    "capacity_max": max(capacities),
    "zero_capacity_links": capacities.count(0),
    "components": len(diameters),
    "hop_diameter": max(diameters),
    "resolution": topology.q,
    "demands": {},
    "steps": [],
}
for name, demand_set in prepared.demand_sets.items():
    inventory["demands"][name] = [
        {
            "id": d.id,
            "base_id": d.base_id,
            "sources": list(d.endpoints.sources),
            "destinations": list(d.endpoints.destinations),
            "nominal_units": d.nominal_units,
            "model": type(d.config.routing).__name__,
            "admission": str(d.config.admission),
        }
        for d in demand_set.demands
    ]
for step in result.steps:
    value = step.result
    cases = getattr(value, "cases", None)
    inventory["steps"].append(
        {
            "name": step.name,
            "kind": type(value).__name__,
            "workers": getattr(getattr(value, "metadata", None), "workers", None),
            "cases": None if cases is None else len(cases),
            "occurrences": None
            if cases is None
            else sum(c.case.occurrences for c in cases),
            "rows": None
            if cases is None
            else sum(len(s.rows) for c in cases for s in c.states.values()),
            "msd_probes": len(getattr(value, "probes", getattr(value, "history", ()))),
        }
    )
(args.out / "inventory.json").write_text(json.dumps(inventory, indent=2) + "\n")
if tracer:
    tracer.finish()
for key, profile in worker_profiles.items():
    profile.dump_stats(str(args.out / f"worker-{key}.pstats"))
export_path = args.out / "results.json"
summary = {
    "exploratory": True,
    "workload": args.workload,
    "workers_override": args.workers,
    "phases_s": phases,
    "steps": steps,
    "profile_scope": "CPython 3.14 monitoring includes worker threads; aggregate times overlap",
    "output_bytes": export_path.stat().st_size,
    "output_sha256": hashlib.sha256(export_path.read_bytes()).hexdigest(),
    "python": sys.version,
    "platform": platform.platform(),
    "cpu_count": os.cpu_count(),
    "ngraph": version("ngraph"),
    "netgraph_core": version("netgraph-core"),
    "ngraph_path": ngraph.__file__,
    "core_path": core.__file__,
    "peak_rss_bytes": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
    "process_elapsed_before_cleanup_s": time.perf_counter() - started,
}
with phase("cleanup"):
    del result, prepared, scenario
    gc.collect()
(args.out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
print(
    json.dumps(
        {
            "summary": "completed exploratory profile",
            "workload": args.workload,
            "workers": args.workers,
            "phases_s": phases,
        }
    ),
    flush=True,
)
