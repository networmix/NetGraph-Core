#!/usr/bin/env python3
"""Realistic NetGraph workloads for the Core old/new A/B/A.

Run with NetGraph's interpreter and PYTHONPATH pointing at the Core build under
test (see run_netgraph_aba.sh). Prints CSV rows:
    workload,label,median_s,samples,checksum,samples_s
Every workload returns a checksum of its results that must be identical across
builds; the driver compares them.
"""

from __future__ import annotations

import hashlib
import json
import random
import statistics
import sys
import time
from pathlib import Path

NETGRAPH = Path.home() / "ws" / "NetGraph"
sys.path.insert(0, str(NETGRAPH))

import numpy as np  # noqa: E402
from ngraph.analysis import FailureManager  # noqa: E402
from ngraph.analysis.context import AnalysisContext  # noqa: E402
from ngraph.analysis.functions import demand_placement_analysis  # noqa: E402
from ngraph.model.network import Link, Network, Node  # noqa: E402
from ngraph.scenario import Scenario  # noqa: E402
from ngraph.types.base import Mode  # noqa: E402

import netgraph_core  # noqa: E402


def norm(obj):
    """JSON-safe, order-stable view of nested results (tuple keys become strings)."""
    if isinstance(obj, dict):
        skip = {
            "metadata",
            "execution_time",
            "duration",
            "timestamp",
            "started_at",
            "finished_at",
            "run_id",
            "id",
        }
        return {
            str(k): norm(v)
            for k, v in sorted(obj.items(), key=lambda kv: str(kv[0]))
            if str(k) not in skip
        }
    if isinstance(obj, (list, tuple, set, frozenset)):
        items = [norm(v) for v in obj]
        return sorted(items, key=str) if isinstance(obj, (set, frozenset)) else items
    if isinstance(obj, float):
        return repr(obj)
    if isinstance(obj, (str, int, bool)) or obj is None:
        return obj
    if hasattr(obj, "to_dict"):
        return norm(obj.to_dict())
    return str(obj)


def digest(obj) -> str:
    return hashlib.sha256(json.dumps(norm(obj), sort_keys=True).encode()).hexdigest()[
        :16
    ]


def clos(leaves: int, spines: int, cap: float = 100.0) -> Network:
    net = Network()
    for i in range(leaves):
        net.add_node(Node(f"leaf/leaf{i:03d}"))
    for j in range(spines):
        net.add_node(Node(f"spine/spine{j:03d}"))
    for i in range(leaves):
        for j in range(spines):
            net.add_link(
                Link(f"leaf/leaf{i:03d}", f"spine/spine{j:03d}", capacity=cap, cost=1.0)
            )
    return net


def grid(rows: int, cols: int, seed: int = 3) -> Network:
    rng = random.Random(seed)
    net = Network()
    for r in range(rows):
        for c in range(cols):
            net.add_node(Node(f"g/r{r:03d}c{c:03d}"))
    for r in range(rows):
        for c in range(cols):
            if c + 1 < cols:
                net.add_link(
                    Link(
                        f"g/r{r:03d}c{c:03d}",
                        f"g/r{r:03d}c{c + 1:03d}",
                        capacity=10.0,
                        cost=float(rng.randint(1, 20)),
                    )
                )
            if r + 1 < rows:
                net.add_link(
                    Link(
                        f"g/r{r:03d}c{c:03d}",
                        f"g/r{r + 1:03d}c{c:03d}",
                        capacity=10.0,
                        cost=float(rng.randint(1, 20)),
                    )
                )
    return net


def clos3(pods: int, leaves: int, spines: int, supers: int) -> Network:
    """Three-tier Clos: per pod a leaf/spine mesh, spines meshed to super-spines."""
    net = Network()
    for k in range(supers):
        net.add_node(Node(f"ss/ss{k:03d}"))
    for p in range(pods):
        for i in range(leaves):
            net.add_node(Node(f"pod{p:02d}/leaf/leaf{i:03d}"))
        for j in range(spines):
            net.add_node(Node(f"pod{p:02d}/spine/spine{j:03d}"))
        for i in range(leaves):
            for j in range(spines):
                net.add_link(
                    Link(
                        f"pod{p:02d}/leaf/leaf{i:03d}",
                        f"pod{p:02d}/spine/spine{j:03d}",
                        capacity=100.0,
                        cost=1.0,
                    )
                )
        for j in range(spines):
            for k in range(supers):
                net.add_link(
                    Link(
                        f"pod{p:02d}/spine/spine{j:03d}",
                        f"ss/ss{k:03d}",
                        capacity=400.0,
                        cost=1.0,
                    )
                )
    return net


def scenario_text(name: str, iterations: int) -> str:
    path = NETGRAPH / "scenarios" / name
    text = path.read_text()
    text = text.replace("parallelism: auto", "parallelism: 1").replace(
        "parallelism: 8", "parallelism: 1"
    )
    text = text.replace("iterations: 1000", f"iterations: {iterations}")
    return text


def results_digest(scenario: Scenario) -> str:
    out = {}
    for step in scenario.results.get_steps_by_execution_order():
        data = scenario.results.get_step(step)
        out[step] = {k: v for k, v in data.items() if k not in ("metadata",)}
    return digest(out)


def build_workloads() -> dict[str, tuple[callable, int]]:
    """Return name -> (fn returning checksum, inner repetitions)."""
    wl: dict[str, tuple[callable, int]] = {}
    sel = netgraph_core.EdgeSelection(
        multi_edge=True,
        require_capacity=False,
        tie_break=netgraph_core.EdgeTieBreak.DETERMINISTIC,
    )

    clos200 = clos(200, 200)
    ctx_clos200 = AnalysisContext.from_network(clos200)
    n_clos = ctx_clos200.multidigraph.num_nodes()

    def spf_clos():
        acc = []
        for s in range(n_clos):
            dist, dag = ctx_clos200.algorithms.spf(ctx_clos200.handle, s, selection=sel)
            acc.append((int(np.asarray(dist).sum()), int(len(dag.parents))))
        return digest(acc)

    wl["spf_full_clos_200x200"] = (spf_clos, 1)

    grid100 = grid(100, 100)
    ctx_grid = AnalysisContext.from_network(grid100)

    def spf_grid():
        acc = []
        for s in range(0, 10000, 200):
            dist, dag = ctx_grid.algorithms.spf(ctx_grid.handle, s, selection=sel)
            d = np.asarray(dist)
            acc.append((int(d[d < 2**62].sum()), int(len(dag.parents))))
        return digest(acc)

    wl["spf_full_grid_100x100_weighted"] = (spf_grid, 5)

    ctx_mf_clos = AnalysisContext.from_network(
        clos200,
        source="^leaf/leaf00[0-9]$",
        sink="^leaf/leaf19[0-9]$",
        mode=Mode.COMBINE,
    )
    wl["maxflow_clos_200x200_combine"] = (lambda: digest(ctx_mf_clos.max_flow()), 50)

    clos100 = clos(100, 100)
    ctx_fail = AnalysisContext.from_network(
        clos100,
        source="^leaf/leaf00[0-9]$",
        sink="^leaf/leaf09[0-9]$",
        mode=Mode.COMBINE,
    )
    link_ids = sorted(clos100.links)
    rng = random.Random(11)
    exclusions = [set(rng.sample(link_ids, 12)) for _ in range(100)]

    def maxflow_failures():
        return digest([ctx_fail.max_flow(excluded_links=ex) for ex in exclusions])

    wl["maxflow_failures_clos_100x100_x100"] = (maxflow_failures, 1)

    ctx_mf_grid = AnalysisContext.from_network(
        grid100, source="^g/r00[0-9]c000$", sink="^g/r09[0-9]c099$", mode=Mode.COMBINE
    )
    wl["maxflow_grid_100x100_weighted_combine"] = (
        lambda: digest(ctx_mf_grid.max_flow()),
        20,
    )
    wl["maxflow_grid_100x100_weighted_pairwise"] = (
        lambda: digest(
            AnalysisContext.from_network(
                grid100, source="^g/r000c000$", sink="^g/r099c099$", mode=Mode.PAIRWISE
            ).max_flow()
        ),
        5,
    )

    backbone = Scenario.from_yaml(scenario_text("backbone_clos.yml", 1000))
    demands = [
        td.to_dict() for td in backbone.demand_set.get_set("baseline_traffic_matrix")
    ]
    wl["placement_backbone_clos"] = (
        lambda: digest(
            demand_placement_analysis(backbone.network, set(), set(), demands)
        ),
        50,
    )

    fm = FailureManager(
        network=backbone.network,
        failure_policy_set=backbone.failure_policy_set,
        policy_name="weighted_modes",
    )
    wl["placement_mc_backbone_clos_x300"] = (
        lambda: digest(
            fm.run_demand_placement_monte_carlo(
                demands_config=demands, iterations=300, parallelism=1, seed=42
            )
        ),
        1,
    )
    wl["maxflow_mc_backbone_clos_x300"] = (
        lambda: digest(
            fm.run_max_flow_monte_carlo(
                source="^metro1/dc1/.*",
                target="^metro2/dc1/.*",
                mode="combine",
                iterations=300,
                parallelism=1,
                seed=7,
            )
        ),
        1,
    )

    fabric = clos3(pods=8, leaves=16, spines=8, supers=16)
    ctx_fabric = AnalysisContext.from_network(fabric)
    wl["ksp_clos3_pairwise_k8"] = (
        lambda: digest(
            {
                f"{k[0]}->{k[1]}": [str(p) for p in v]
                for k, v in ctx_fabric.k_shortest_paths(
                    "^pod00/leaf/leaf00[0-3]$",
                    "^pod0[1-7]/leaf/leaf00[0-1]$",
                    mode=Mode.PAIRWISE,
                    max_k=8,
                ).items()
            }
        ),
        1,
    )
    wl["maxflow_clos3_pod_pairs_combine"] = (
        lambda: digest(
            [
                AnalysisContext.from_network(
                    fabric,
                    source=f"^pod{a:02d}/leaf/.*",
                    sink=f"^pod{b:02d}/leaf/.*",
                    mode=Mode.COMBINE,
                ).max_flow()
                for a in range(8)
                for b in range(8)
                if a != b
            ]
        ),
        1,
    )
    fabric_demands = [
        {
            "source": f"^pod{a:02d}/leaf/.*",
            "target": f"^pod{b:02d}/leaf/.*",
            "mode": "pairwise",
            "priority": 0,
            "volume": 2000.0,
            "flow_policy": policy,
        }
        for a in range(8)
        for b in range(8)
        if a != b
        for policy in ("TE_WCMP_UNLIM",)
    ] + [
        {
            "source": f"^pod{a:02d}/leaf/.*",
            "target": f"^pod{(a + 1) % 8:02d}/leaf/.*",
            "mode": "pairwise",
            "priority": 1,
            "volume": 500.0,
            "flow_policy": "SHORTEST_PATHS_ECMP",
        }
        for a in range(8)
    ]
    wl["placement_clos3_wcmp_and_ecmp"] = (
        lambda: digest(demand_placement_analysis(fabric, set(), set(), fabric_demands)),
        1,
    )
    fabric_links = sorted(fabric.links)
    rng2 = random.Random(5)
    fabric_exclusions = [set(rng2.sample(fabric_links, 6)) for _ in range(20)]
    wl["placement_clos3_failures_x20"] = (
        lambda: digest(
            [
                demand_placement_analysis(fabric, set(), ex, fabric_demands)
                for ex in fabric_exclusions
            ]
        ),
        1,
    )

    def run_scenario(name: str, iterations: int):
        def fn():
            sc = Scenario.from_yaml(scenario_text(name, iterations))
            sc.run()
            return results_digest(sc)

        return fn

    wl["scenario_square_mesh_iter100"] = (run_scenario("square_mesh.yaml", 100), 1)
    wl["scenario_nsfnet_iter200"] = (run_scenario("nsfnet.yaml", 200), 1)
    wl["scenario_backbone_clos_iter100"] = (run_scenario("backbone_clos.yml", 100), 1)
    return wl


def main() -> int:
    label = sys.argv[1] if len(sys.argv) > 1 else "run"
    samples = int(sys.argv[2]) if len(sys.argv) > 2 else 5
    only = sys.argv[3].split(",") if len(sys.argv) > 3 else None
    wl = build_workloads()
    print("workload,label,median_s,samples,checksum,samples_s", flush=True)
    for name, (fn, reps) in wl.items():
        if only and name not in only:
            continue
        checksum = fn()  # warm-up and reference checksum
        times = []
        for _ in range(samples):
            t0 = time.perf_counter()
            for _ in range(reps):
                c = fn()
            times.append(time.perf_counter() - t0)
            if c != checksum:
                raise SystemExit(
                    f"nondeterministic result in {name}: {c} != {checksum}"
                )
        print(
            f"{name},{label},{statistics.median(times):.6f},{samples},{checksum},{';'.join(f'{t:.6f}' for t in times)}",
            flush=True,
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
