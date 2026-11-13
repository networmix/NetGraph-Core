"""Make fresh research inputs without writing in a consumer checkout."""

import gzip
import hashlib
import json
import shutil
from pathlib import Path

import networkx as nx
from netlab_topogen.config import TopologyConfig
from netlab_topogen.library import load_library
from netlab_topogen.scenario import build_scenario
from netlab_topogen.validation import validate_scenario

BASE = Path(__file__).resolve().parents[1]
REPO = BASE.parents[1] / "worktrees/netgraph-orient"
out = BASE / "workloads"
out.mkdir()
manifest = {}
for name, filename in (
    ("readme", "readme_example.yml"),
    ("square", "square_mesh.yaml"),
    ("nsfnet", "nsfnet.yaml"),
):
    source = REPO / "scenarios" / filename
    target = out / f"{name}.yml"
    shutil.copyfile(source, target)
    manifest[name] = {"source": str(source), "path": str(target), "changes": "none"}
source = REPO / "docs/design/perf-0.24-vs-0.26/inputs/dcbb96-0.26.0.yml.gz"
target = out / "dcbb96.yml"
target.write_bytes(gzip.decompress(source.read_bytes()))
manifest["dcbb96"] = {
    "source": str(source),
    "path": str(target),
    "changes": "decompressed only; full 200 draws per step, historical replay risk groups",
}
# Use the checked-in offline corridor shape, not a claimed US reconstruction.
graph = nx.MultiGraph()
for i, name in enumerate(("alpha", "beta", "gamma")):
    graph.add_node(
        i,
        node_type="metro",
        name=name,
        name_orig=name,
        metro_id=str(i),
        x=i * 300000.0,
        y=0.0,
        radius_km=20.0,
    )
for a, b in ((0, 1), (1, 2), (0, 2)):
    graph.add_edge(a, b, edge_type="corridor", length_km=300.0)
config = TopologyConfig()
config.build.build_defaults.pop_per_metro = 1
config.build.build_defaults.dc_regions_per_metro = 1
config.traffic.gbps_per_mw = 1
config.traffic.mw_per_dc_region = 10
config.traffic.priority_ratios = {0: 1}
config.output.scenario_seed = 42
result = build_scenario(graph, config, load_library(None))
validation = validate_scenario(result.yaml, config.components)
assert not validation.issues, validation.issues
target = out / "topogen.yml"
target.write_text(result.yaml)
manifest["topogen"] = {
    "source": "/Users/networmix/ws/NetLab-TopoGen/tests/conftest.py:10-47; public build_scenario API",
    "path": str(target),
    "changes": "fresh offline build from three synthetic metros; seed explicitly 42",
}
for item in manifest.values():
    item["sha256"] = hashlib.sha256(Path(item["path"]).read_bytes()).hexdigest()
(out / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
print(
    json.dumps(
        {"summary": "generated and validated inputs", "workloads": list(manifest)}
    )
)
