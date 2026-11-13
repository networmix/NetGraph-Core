"""Generate a synthetic TopoGen-style backbone scenario (exploratory workload, not a real network).

    python gen_backbone.py OUT.yml [--metros 25] [--k 3] [--partners 5] [--draws 1000] [--preset TE_MAXFLOW_UNLIM]

Each metro has two PoPs, each a Clos of 16 leaves and 4 spines (leaf-spine mesh, cost 1), the
PoPs joined leaf to leaf (one to one, cost 163), and one DC node linked to every leaf of both
PoPs (cost 1). Metros sit on a circle; each is joined to its k nearest neighbours on each side
by corridors joining leaf i of pop p to leaf i of pop p (cost from the chord length). Demands:
each DC to its `partners` nearest DCs (both directions), pairwise, the given preset, as
backbone_clos.yml declares them. Workflow: MSD, then a placement under `draws` draws of a
weighted multi-mode failure policy with cost reports, as TopoGen's design_analysis_brief.
"""

import argparse
import math

import yaml

p = argparse.ArgumentParser()
p.add_argument("out")
p.add_argument("--metros", type=int, default=25)
p.add_argument("--k", type=int, default=2)
p.add_argument("--partners", type=int, default=5)
p.add_argument("--draws", type=int, default=1000)
p.add_argument("--preset", default="TE_MAXFLOW_UNLIM")
p.add_argument("--leaves", type=int, default=16)
p.add_argument("--spines", type=int, default=4)
a = p.parse_args()

M, L, S = a.metros, a.leaves, a.spines
nodes, links = {}, []


def link(src, dst, cap, cost, kind, adj):
    links.append(
        {
            "source": f"^{src}$",
            "target": f"^{dst}$",
            "pattern": "one_to_one",
            "capacity": float(cap),
            "cost": int(cost),
            "attrs": {"link_type": kind, "adjacency_id": adj},
        }
    )


for m in range(1, M + 1):
    for pop in (1, 2):
        for i in range(1, L + 1):
            nodes[f"metro{m}/pop{pop}/leaf{i}"] = {"attrs": {"role": "leaf", "node_type": "pop"}}
        for j in range(1, S + 1):
            nodes[f"metro{m}/pop{pop}/spine{j}"] = {"attrs": {"role": "spine", "node_type": "pop"}}
        for i in range(1, L + 1):
            for j in range(1, S + 1):
                link(f"metro{m}/pop{pop}/leaf{i}", f"metro{m}/pop{pop}/spine{j}", 3200, 1, "leaf_spine", f"pop:{m}:{pop}")
    nodes[f"metro{m}/dc1"] = {"attrs": {"role": "dc", "node_type": "dc_region"}}
    for i in range(1, L + 1):
        link(f"metro{m}/pop1/leaf{i}", f"metro{m}/pop2/leaf{i}", 1400, 163, "intra_metro", f"intra:{m}")
        for pop in (1, 2):
            link(f"metro{m}/pop{pop}/leaf{i}", f"metro{m}/dc1", 2000, 1, "dc_to_pop", f"dc:{m}")

pos = {m: (math.cos(2 * math.pi * m / M), math.sin(2 * math.pi * m / M)) for m in range(1, M + 1)}
pairs = set()
for m in range(1, M + 1):
    for d in range(1, a.k + 1):
        n = (m - 1 + d) % M + 1
        pairs.add(tuple(sorted((m, n))))
    # one long chord per metro for path diversity
    n = (m - 1 + M // 2) % M + 1
    if m % 3 == 0:
        pairs.add(tuple(sorted((m, n))))
for m, n in sorted(pairs):
    dist = math.dist(pos[m], pos[n])
    cost = max(1, int(round(400 * dist)))
    for pop in (1, 2):
        for i in range(1, L + 1):
            link(f"metro{m}/pop{pop}/leaf{i}", f"metro{n}/pop{pop}/leaf{i}", 800, cost, "inter_metro_corridor", f"inter:{m}-{n}")

demands = []
for m in range(1, M + 1):
    near = sorted((n for n in range(1, M + 1) if n != m), key=lambda n: (math.dist(pos[m], pos[n]), n))
    for n in near[: a.partners]:
        demands.append(
            {"source": f"^metro{m}/dc1$", "target": f"^metro{n}/dc1$", "mode": "pairwise", "volume": 1000.0, "preset": a.preset}
        )

scenario = {
    "seed": 42,
    "network": {"nodes": nodes, "links": links},
    "risk_groups": [],
    "failures": {
        "weighted_modes": {
            "modes": [
                {"weight": 0.35, "rules": [{"scope": "link", "mode": "choice", "count": 3,
                    "match": {"conditions": [{"attr": "link_type", "op": "==", "value": "dc_to_pop"}], "logic": "and"}}]},
                {"weight": 0.3, "rules": [{"scope": "link", "mode": "choice", "count": 1,
                    "match": {"conditions": [{"attr": "link_type", "op": "==", "value": "inter_metro_corridor"}], "logic": "and"}}]},
                {"weight": 0.25, "rules": [{"scope": "node", "mode": "choice", "count": 1,
                    "match": {"conditions": [{"attr": "node_type", "op": "!=", "value": "dc_region"}], "logic": "and"}}]},
                {"weight": 0.1, "rules": [{"scope": "link", "mode": "choice", "count": 4,
                    "match": {"conditions": [{"attr": "link_type", "op": "==", "value": "leaf_spine"}], "logic": "and"}}]},
            ]
        }
    },
    "demands": {"baseline_traffic_matrix": demands},
    "workflow": [
        {"type": "MaximumSupportedDemand", "name": "msd_baseline", "demand_set": "baseline_traffic_matrix",
         "alpha_min": 0.001, "alpha_max": 1000000.0, "alpha_tolerance": 0.05, "max_bracket_iters": 16},
        {"type": "TrafficMatrixPlacement", "name": "tm_placement", "seed": 42, "demand_set": "baseline_traffic_matrix",
         "failure_policy": "weighted_modes", "iterations": a.draws, "outputs": {"cost_distribution": True},
         "alpha_from_step": "msd_baseline"},
    ],
}
del scenario["risk_groups"]
with open(a.out, "w") as f:
    yaml.safe_dump(scenario, f, sort_keys=False)
print(f"summary: {len(nodes)} nodes, {len(links)} links, {len(demands)} demands, {len(pairs)} corridors")
