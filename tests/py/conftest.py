"""Shared pytest fixtures for the test suite.

Provides small, reusable builders to reduce repetition across tests.

Fixtures (logical order):
- build_graph: generic builder from edge tuples.
- line1_graph: 3-node line with tiered B->C parallels.
- square1_graph: shortest vs longer alternative path.
- square2_graph: two equal-cost shortest paths.
- two_disjoint_shortest_graph: two disjoint equal-cost routes with bottlenecks.
- graph3: mixed 6-node topology with parallels and a longer route.
- square4_graph: 4-node square with mixed tiers.
- fully_connected_graph: factory for K_n directed (no self-loops).
"""

from __future__ import annotations

import numpy as np
import pytest

import netgraph_core as ngc

# Global Algorithms fixture (CPU backend) and helper to build graph handles


@pytest.fixture
def algs():
    be = ngc.Backend.cpu()
    return ngc.Algorithms(be)


@pytest.fixture
def to_handle(algs):
    def _h(g: ngc.StrictMultiDiGraph):
        return algs.build_graph(g)

    return _h


# Centralized helper fixtures and shared assertions for tests


@pytest.fixture
def t_cost():
    """Return a function to extract distance to a target node as float."""

    def _t_cost(dist: np.ndarray, node: int) -> float:
        return float(dist[int(node)])

    return _t_cost


@pytest.fixture
def flows_by_eid():
    """Return a function mapping EdgeId -> edge flow value for a graph."""

    def _flows_by_eid(
        g: ngc.StrictMultiDiGraph, edge_flows: np.ndarray
    ) -> dict[int, float]:
        arr = np.asarray(edge_flows, dtype=float)
        return {int(eid): float(arr[eid]) for eid in range(g.num_edges())}

    return _flows_by_eid


@pytest.fixture
def cost_dist_dict():
    """Return a function converting FlowSummary (costs, flows) to a dict."""

    def _to_dict(summary) -> dict[float, float]:
        costs = np.asarray(getattr(summary, "costs", []), dtype=float)
        flows = np.asarray(getattr(summary, "flows", []), dtype=float)
        return {float(c): float(f) for c, f in zip(costs, flows, strict=False)}

    return _to_dict


@pytest.fixture
def assert_pred_dag_integrity():
    """Return an assertion helper to validate PredDAG shapes and id ranges."""

    def _assert(g: ngc.StrictMultiDiGraph, dag: ngc.PredDAG) -> None:
        off = np.asarray(dag.parent_offsets)
        par = np.asarray(dag.parents)
        via = np.asarray(dag.via_edges)
        # Offsets length and monotonicity
        assert off.shape == (g.num_nodes() + 1,)
        assert np.all(off[:-1] <= off[1:])
        total = int(off[-1])
        # Parents/via sizes match total entries
        assert par.shape == (total,)
        assert via.shape == (total,)
        if total > 0:
            # Id ranges
            assert int(par.min()) >= 0 and int(par.max()) < g.num_nodes()
            assert int(via.min()) >= 0 and int(via.max()) < g.num_edges()

    return _assert


@pytest.fixture
def assert_edge_flows_shape():
    """Return an assertion helper to validate edge_flows presence/shape.

    expected_present=True means summary.edge_flows must exist and be length N.
    If False, allow either empty list or omitted; if present, allow 0 or N.
    """

    def _assert(
        g: ngc.StrictMultiDiGraph, summary, expected_present: bool = True
    ) -> None:
        edge_flows = getattr(summary, "edge_flows", None)
        if expected_present:
            assert edge_flows is not None
            arr = np.asarray(edge_flows)
            assert arr.shape == (g.num_edges(),)
        else:
            if edge_flows is not None:
                arr = np.asarray(edge_flows)
                assert arr.size in (0, g.num_edges())

    return _assert


@pytest.fixture
def assert_valid_min_cut():
    """Return an assertion helper to validate MinCut edge ids are unique and valid."""

    def _assert(g: ngc.StrictMultiDiGraph, min_cut) -> None:
        edges = [int(e) for e in getattr(min_cut, "edges", [])]
        assert len(edges) == len(set(edges))
        for e in edges:
            assert 0 <= e < g.num_edges()

    return _assert


@pytest.fixture
def build_graph():
    """Build a StrictMultiDiGraph from edge tuples.

    Args:
        num_nodes: Number of nodes in the graph
        edges: List of tuples (src, dst, cost, cap[, ext_id])

    Returns:
        StrictMultiDiGraph instance
    """

    def _builder(num_nodes: int, edges: list[tuple]) -> ngc.StrictMultiDiGraph:
        if not edges:
            return ngc.StrictMultiDiGraph.from_arrays(
                num_nodes,
                np.empty((0,), dtype=np.int32),
                np.empty((0,), dtype=np.int32),
                np.empty((0,), dtype=np.float64),
                np.empty((0,), dtype=np.int64),
                np.empty((0,), dtype=np.int64),
            )
        src = np.array([e[0] for e in edges], dtype=np.int32)
        dst = np.array([e[1] for e in edges], dtype=np.int32)
        cost = np.array([int(e[2]) for e in edges], dtype=np.int64)
        cap = np.array([float(e[3]) for e in edges], dtype=np.float64)
        ext_edge_ids = np.arange(len(edges), dtype=np.int64)
        return ngc.StrictMultiDiGraph.from_arrays(
            num_nodes, src, dst, cap, cost, ext_edge_ids
        )

    return _builder


@pytest.fixture
def line1_graph(build_graph):
    """Line graph A->B->C with unit costs; caps as in typical tests.

    Nodes: 0->1 (cap 5, cost 1), 1->2 has parallel edges with costs {1,1,2} and
    capacities {1,3,7} to test tiering.
    """

    edges = [
        (0, 1, 1, 5, 0),  # A->B
        (1, 2, 1, 1, 2),  # B->C (min-cost)
        (1, 2, 1, 3, 4),  # B->C (min-cost)
        (1, 2, 2, 7, 6),  # B->C (higher cost)
    ]
    return build_graph(3, edges)


@pytest.fixture
def square1_graph(build_graph):
    """Square1: one shortest route and one longer alternative.

    0->1->2 with total cost 2; 0->3->2 with total cost 4.
    """

    edges = [
        (0, 1, 1, 1, 0),
        (1, 2, 1, 1, 1),
        (0, 3, 2, 2, 2),
        (3, 2, 2, 2, 3),
    ]
    return build_graph(4, edges)


@pytest.fixture
def square2_graph(build_graph):
    """Square2: two equal-cost shortest routes from 0->2 via 1 and 3.

    All edge costs are 1; both paths have total cost 2.
    """

    edges = [
        (0, 1, 1, 1, 0),
        (1, 2, 1, 1, 1),
        (0, 3, 1, 2, 2),
        (3, 2, 1, 2, 3),
    ]
    return build_graph(4, edges)


@pytest.fixture
def graph3(build_graph):
    """Six-node mixed topology with parallel edges and multiple paths."""

    edges = [
        (0, 1, 1, 2, 0),  # A->B parallels
        (0, 1, 1, 4, 1),
        (0, 1, 1, 6, 2),
        (1, 2, 1, 1, 3),  # B->C parallels
        (1, 2, 1, 2, 4),
        (1, 2, 1, 3, 5),
        (2, 3, 2, 3, 6),  # C->D
        (0, 4, 1, 5, 7),  # A->E
        (4, 2, 1, 4, 8),  # E->C
        (0, 3, 4, 2, 9),  # A->D
        (2, 5, 1, 1, 10),  # C->F
        (5, 3, 1, 2, 11),  # F->D
    ]
    return build_graph(6, edges)


@pytest.fixture
def square4_graph(build_graph):
    """A 4-node square with parallel/higher-cost tiers used in max-flow tests.

    Nodes: A=0, B=1, C=2, D=3
    """

    edges = [
        (0, 1, 1, 100, 0),  # A->B
        (1, 2, 1, 125, 1),  # B->C
        (0, 3, 1, 75, 2),  # A->D
        (3, 2, 1, 50, 3),  # D->C
        (1, 3, 1, 50, 4),  # B->D
        (3, 1, 1, 50, 5),  # D->B
        (0, 1, 2, 200, 6),  # A->B (higher cost)
        (1, 3, 2, 200, 7),  # B->D (higher cost)
        (3, 2, 2, 200, 8),  # D->C (higher cost)
    ]
    return build_graph(4, edges)


@pytest.fixture
def triangle1_graph(build_graph):
    """Triangle graph with higher capacities on A<->B and B<->C; lower on A<->C.

    Nodes: 0(A), 1(B), 2(C)
    Edges directed both ways to match the fixture used across tests.
    """

    edges = [
        (0, 1, 1, 15, 0),
        (1, 0, 1, 15, 1),
        (1, 2, 1, 15, 2),
        (2, 1, 1, 15, 3),
        (0, 2, 1, 5, 4),
        (2, 0, 1, 5, 5),
    ]
    return build_graph(3, edges)


@pytest.fixture
def square3_graph(build_graph):
    """Square with cross-links B<->D and mixed capacities.

    Nodes: 0(A),1(B),2(C),3(D)
    """

    edges = [
        (0, 1, 1, 100, 0),
        (1, 2, 1, 125, 1),
        (0, 3, 1, 75, 2),
        (3, 2, 1, 50, 3),
        (1, 3, 1, 50, 4),
        (3, 1, 1, 50, 5),
    ]
    return build_graph(4, edges)


@pytest.fixture
def graph1_graph(build_graph):
    """Five-node graph with bidirectional edges between B and C."""

    edges = [
        (0, 1, 1, 1, 0),
        (0, 2, 1, 1, 1),
        (1, 3, 1, 1, 2),
        (2, 3, 1, 1, 3),
        (1, 2, 1, 1, 4),
        (2, 1, 1, 1, 5),
        (3, 4, 1, 1, 6),
    ]
    return build_graph(5, edges)


@pytest.fixture
def graph2_graph(build_graph):
    """Five-node branched graph with bidirectional C-D edges."""

    edges = [
        (0, 1, 1, 1, 0),
        (1, 2, 1, 1, 1),
        (1, 3, 1, 1, 2),
        (2, 3, 1, 1, 3),
        (3, 2, 1, 1, 4),
        (2, 4, 1, 1, 5),
        (3, 4, 1, 1, 6),
    ]
    return build_graph(5, edges)


@pytest.fixture
def graph4_graph(build_graph):
    """Graph with three parallel branches from A to C of increasing cost/capacity.

    Nodes: 0(A),1(B),2(B1),3(B2),4(C)
    """

    edges = [
        (0, 1, 1, 1, 0),
        (1, 4, 1, 1, 1),
        (0, 2, 2, 2, 2),
        (2, 4, 2, 2, 3),
        (0, 3, 3, 3, 4),
        (3, 4, 3, 3, 5),
    ]
    return build_graph(5, edges)


@pytest.fixture
def make_fully_connected_graph(build_graph):
    """Return a builder for fully-connected directed graphs (no self-loops).

    Usage: make_fully_connected_graph(n, cost=1.0, cap=1.0)
    """

    def _build(
        n: int, *, cost: float = 1.0, cap: float = 1.0
    ) -> ngc.StrictMultiDiGraph:
        edges = []
        lid = 0
        for s in range(n):
            for d in range(n):
                if s == d:
                    continue
                edges.append((s, d, float(cost), float(cap), lid))
                lid += 1
        return build_graph(n, edges)

    return _build


@pytest.fixture
def dag_to_pred_map():
    """Return a converter from PredDAG to {node: {parent: [EdgeId...]}} mapping."""

    def _convert(g: ngc.StrictMultiDiGraph, dag: ngc.PredDAG):
        pred: dict[int, dict[int, list[int]]] = {}
        offsets = np.asarray(dag.parent_offsets)
        parents = np.asarray(dag.parents)
        via = np.asarray(dag.via_edges)
        n = g.num_nodes()
        for v in range(n):
            start = int(offsets[v])
            end = int(offsets[v + 1])
            if start == end:
                continue
            group: dict[int, list[int]] = {}
            for i in range(start, end):
                p = int(parents[i])
                e = int(via[i])
                group.setdefault(p, []).append(e)
            pred[v] = group
        if 0 not in pred and n > 0:
            pred[0] = {}
        return pred

    return _convert


@pytest.fixture
def make_pred_map(dag_to_pred_map):
    """Alias fixture for converting PredDAG to {node: {parent: [EdgeId]}}.

    Provided to make intent clearer at call sites.
    """

    return dag_to_pred_map


@pytest.fixture
def assert_paths_concrete():
    """Validate path tuples returned by resolve_to_paths.

    Ensures structure is ((node, (edge_ids...)), ..., (dst, ())).
    If expect_split=True, each hop (except src/dst) must have exactly one edge id.
    If False, hops must have >=1 edge id.
    """

    def _assert(paths, src: int, dst: int, expect_split: bool) -> None:
        for path in paths:
            # path should be a tuple/list of elements
            assert isinstance(path, (tuple, list)) and len(path) >= 2
            # First element is (src, ())
            first = path[0]
            assert isinstance(first, (tuple, list)) and len(first) == 2
            assert int(first[0]) == int(src)
            assert isinstance(first[1], (tuple, list)) and len(first[1]) >= 1
            if expect_split:
                assert len(first[1]) == 1
            # Last element is (dst, ())
            last = path[-1]
            assert isinstance(last, (tuple, list)) and len(last) == 2
            assert int(last[0]) == int(dst)
            assert isinstance(last[1], (tuple, list)) and len(last[1]) == 0
            # Intermediate hops
            for elem in path[1:-1]:
                assert isinstance(elem, (tuple, list)) and len(elem) == 2
                node, edges = elem
                # node id is int
                assert isinstance(int(node), int)
                # edges is tuple/list of ints
                assert isinstance(edges, (tuple, list)) and len(edges) >= 1
                if expect_split:
                    assert len(edges) == 1
                for e in edges:
                    assert isinstance(int(e), int)

    return _assert


@pytest.fixture
def two_disjoint_shortest_graph(build_graph):
    """Two disjoint shortest routes 0->1->3 and 0->2->3 with equal cost.

    Capacities configured so that total is limited by per-path bottlenecks.
    """

    edges = [
        (0, 1, 1, 3, 0),  # S->A cap 3
        (1, 3, 1, 2, 1),  # A->T cap 2
        (0, 2, 1, 4, 2),  # S->B cap 4
        (2, 3, 1, 1, 3),  # B->T cap 1
    ]
    return build_graph(4, edges)


@pytest.fixture
def graph5(build_graph):
    """Fully connected 5-node directed graph with unit costs/caps."""

    edges = []
    lid = 0
    for s in range(5):
        for d in range(5):
            if s == d:
                continue
            edges.append((s, d, 1.0, 1.0, lid))
            lid += 1
    return build_graph(5, edges)


@pytest.fixture
def square5_graph(build_graph):
    """Five-node 'square5' topology used in KSP tests.

    Nodes: 0(A),1(B),2(C),3(D),4(E)
    Edges: A->B, A->C, B->D, C->D, B->C, C->B; E is isolated for negative tests.
    """

    edges = [
        (0, 1, 1, 1, 0),  # A->B
        (0, 2, 1, 1, 1),  # A->C
        (1, 3, 1, 1, 2),  # B->D
        (2, 3, 1, 1, 3),  # C->D
        (1, 2, 1, 1, 4),  # B->C
        (2, 1, 1, 1, 5),  # C->B
    ]
    return build_graph(5, edges)


# FlowPolicy canonical configurations
# These represent common use cases for flow routing and traffic engineering


@pytest.fixture
def make_flow_policy_config():
    """Factory for canonical FlowPolicy configurations.

    Canonical configurations:
    - SHORTEST_PATHS_ECMP: IP/IGP-style ECMP (hash-based multipath, shortest paths only)
    - SHORTEST_PATHS_WCMP: IP/IGP-style WCMP (proportional multipath, shortest paths only)
    - TE_WCMP_UNLIM: Traffic engineering with unlimited proportional flows
    - TE_ECMP_UP_TO_256_LSP: TE with up to 256 ECMP LSPs (tunnel-based)
    - TE_ECMP_16_LSP: TE with exactly 16 ECMP LSPs (fixed allocation)

    Usage:
        config = make_flow_policy_config("SHORTEST_PATHS_ECMP")
        policy = ngc.FlowPolicy(algs, graph, config)
    """

    def _make(config_name: str) -> ngc.FlowPolicyConfig:
        if config_name == "SHORTEST_PATHS_ECMP":
            # IP/IGP ECMP: Hash-based equal splitting across shortest paths
            # - multipath=True: Each flow splits across all equal-cost next-hops
            # - require_capacity=False: Routes based on costs only (IGP behavior)
            # - shortest_path=True: Only use lowest-cost paths
            # - max_flow_count=1: Single flow (models single demand stream)
            config = ngc.FlowPolicyConfig()
            config.path_alg = ngc.PathAlg.SPF
            config.flow_placement = ngc.FlowPlacement.EQUAL_BALANCED
            config.selection = ngc.EdgeSelection(
                multi_edge=True,
                require_capacity=False,
                tie_break=ngc.EdgeTieBreak.DETERMINISTIC,
            )
            config.require_capacity = False
            config.multipath = True
            config.shortest_path = True
            config.min_flow_count = 1
            config.max_flow_count = 1
            return config

        elif config_name == "SHORTEST_PATHS_WCMP":
            # IP/IGP WCMP: Proportional splitting across shortest paths
            # - Similar to ECMP but uses PROPORTIONAL placement
            # - Flow splits proportionally to available capacity
            config = ngc.FlowPolicyConfig()
            config.path_alg = ngc.PathAlg.SPF
            config.flow_placement = ngc.FlowPlacement.PROPORTIONAL
            config.selection = ngc.EdgeSelection(
                multi_edge=True,
                require_capacity=False,
                tie_break=ngc.EdgeTieBreak.DETERMINISTIC,
            )
            config.require_capacity = False
            config.multipath = True
            config.shortest_path = True
            config.min_flow_count = 1
            config.max_flow_count = 1
            return config

        elif config_name == "TE_WCMP_UNLIM":
            # Traffic Engineering: Unlimited proportional flows
            # - require_capacity=True: Respects capacity, progressive fill
            # - No flow count limits: Creates flows as needed
            # - Uses all available paths (not just shortest)
            config = ngc.FlowPolicyConfig()
            config.path_alg = ngc.PathAlg.SPF
            config.flow_placement = ngc.FlowPlacement.PROPORTIONAL
            config.selection = ngc.EdgeSelection(
                multi_edge=True,
                require_capacity=True,
                tie_break=ngc.EdgeTieBreak.PREFER_HIGHER_RESIDUAL,
            )
            config.require_capacity = True
            config.multipath = True
            config.min_flow_count = 1
            # max_flow_count not set: unlimited
            return config

        elif config_name == "TE_ECMP_UP_TO_256_LSP":
            # Traffic Engineering: Up to 256 ECMP LSPs (tunnels)
            # - multipath=False: Each flow is a single-path LSP/tunnel
            # - max_flow_count=256: Up to 256 parallel LSPs
            # - ECMP across LSPs (not within each LSP)
            config = ngc.FlowPolicyConfig()
            config.path_alg = ngc.PathAlg.SPF
            config.flow_placement = ngc.FlowPlacement.EQUAL_BALANCED
            config.selection = ngc.EdgeSelection(
                multi_edge=False,
                require_capacity=True,
                tie_break=ngc.EdgeTieBreak.PREFER_HIGHER_RESIDUAL,
            )
            config.require_capacity = True
            config.multipath = False
            config.min_flow_count = 1
            config.max_flow_count = 256
            config.reoptimize_flows_on_each_placement = True
            return config

        elif config_name == "TE_ECMP_16_LSP":
            # Traffic Engineering: Exactly 16 ECMP LSPs (fixed allocation)
            # - Similar to TE_ECMP_UP_TO_256_LSP but fixed at 16 flows
            # - Models hardware with fixed LSP resources
            config = ngc.FlowPolicyConfig()
            config.path_alg = ngc.PathAlg.SPF
            config.flow_placement = ngc.FlowPlacement.EQUAL_BALANCED
            config.selection = ngc.EdgeSelection(
                multi_edge=False,
                require_capacity=True,
                tie_break=ngc.EdgeTieBreak.PREFER_HIGHER_RESIDUAL,
            )
            config.require_capacity = True
            config.multipath = False
            config.min_flow_count = 16
            config.max_flow_count = 16
            config.reoptimize_flows_on_each_placement = True
            return config

        else:
            raise ValueError(f"Unknown configuration: {config_name}")

    return _make


@pytest.fixture
def make_flow_policy(make_flow_policy_config):
    """Factory for creating FlowPolicy with canonical configurations and optional masks.

    Usage:
        policy = make_flow_policy(
            "SHORTEST_PATHS_ECMP",
            algs,
            graph_handle,
            node_mask=mask_array,
            edge_mask=mask_array
        )
    """

    def _make(
        config_name: str,
        algs: ngc.Algorithms,
        graph_handle,
        node_mask=None,
        edge_mask=None,
    ) -> ngc.FlowPolicy:
        # Reuse the single-source preset factory to avoid drift.
        base_cfg = make_flow_policy_config(config_name)
        # Pass masks directly to FlowPolicy constructor (not stored on config).
        return ngc.FlowPolicy(
            algs, graph_handle, base_cfg, node_mask=node_mask, edge_mask=edge_mask
        )

    return _make


# ============================================================================
# Validation helpers for path distribution and flow placement
# ============================================================================


def analyze_path_usage(
    fg: ngc.FlowGraph, graph: ngc.StrictMultiDiGraph, source: int, target: int
) -> dict[int, float]:
    """Analyze which parallel paths are used and their volumes.

    For topologies with parallel paths between source and target,
    identifies which middle nodes are used and how much flow goes through each.

    Args:
        fg: FlowGraph with placed flows
        graph: Topology graph
        source: Source node
        target: Target node

    Returns:
        Dict mapping middle_node_id -> total_flow_volume through that path
        For 2-hop parallel paths: {middle_node: flow_volume}
    """
    edge_flows = fg.edge_flow_view()
    edge_src = graph.edge_src_view()
    edge_dst = graph.edge_dst_view()

    path_usage = {}

    # Find all edges from source and sum their flows
    for edge_id in range(graph.num_edges()):
        if edge_src[edge_id] == source:
            middle_node = edge_dst[edge_id]
            flow = edge_flows[edge_id]
            if flow > 1e-9:  # Only count non-zero flows
                path_usage[middle_node] = flow

    return path_usage


def count_lsps_per_path(
    policy: ngc.FlowPolicy,
    fg: ngc.FlowGraph,
    graph: ngc.StrictMultiDiGraph,
    source: int,
    target: int,
) -> dict[int, list[int]]:
    """Count how many LSPs use each parallel path.

    For topologies with parallel paths, determines which LSPs (flows)
    use which paths.

    Args:
        policy: FlowPolicy with placed flows
        fg: FlowGraph with placed flows
        graph: Topology graph
        source: Source node
        target: Target node

    Returns:
        Dict mapping middle_node_id -> list of flow_ids using that path
    """
    edge_src = graph.edge_src_view()
    edge_dst = graph.edge_dst_view()

    path_lsps: dict[int, list[int]] = {}

    # For each flow, determine which path it uses
    for flow_id in range(policy.flow_count()):
        try:
            flow_index = ngc.FlowIndex(source, target, 0, flow_id)
            flow_edges = fg.get_flow_edges(flow_index)

            # Find the middle node this flow uses (for 2-hop paths)
            # flow_edges is list of (edge_id, volume) tuples
            for edge_id, _ in flow_edges:
                src = edge_src[edge_id]
                dst = edge_dst[edge_id]
                if src == source:
                    middle_node = dst
                    if middle_node not in path_lsps:
                        path_lsps[middle_node] = []
                    path_lsps[middle_node].append(flow_id)
                    break
        except (KeyError, RuntimeError):
            pass  # Flow might not exist or have no edges

    return path_lsps


def check_distribution_balanced(
    lsps_per_path: dict[int, list[int]], tolerance: int = 1
) -> tuple[bool, dict]:
    """Check if LSPs are evenly distributed across paths.

    Args:
        lsps_per_path: Dict mapping path_id -> list of LSP ids
        tolerance: Maximum allowed difference in LSP count between paths

    Returns:
        (is_balanced, stats)
        is_balanced: True if max_count - min_count <= tolerance
        stats: Dict with min, max, counts for debugging
    """
    if not lsps_per_path:
        return (True, {"counts": [], "min": 0, "max": 0, "variance": 0})

    counts = [len(lsps) for lsps in lsps_per_path.values()]
    min_count = min(counts)
    max_count = max(counts)
    variance = max_count - min_count

    stats = {
        "counts": counts,
        "min": min_count,
        "max": max_count,
        "variance": variance,
        "per_path": {path_id: len(lsps) for path_id, lsps in lsps_per_path.items()},
    }

    is_balanced = variance <= tolerance

    return (is_balanced, stats)


def calculate_theoretical_capacity_parallel_paths(
    capacity_per_path: float, num_paths: int, num_lsps: int
) -> float:
    """Calculate theoretical maximum capacity with optimal LSP distribution.

    For parallel equal-capacity paths with single-path LSPs.

    Args:
        capacity_per_path: Capacity of each path
        num_paths: Number of available paths
        num_lsps: Number of LSPs to distribute

    Returns:
        Theoretical maximum total capacity with optimal distribution
    """
    if num_lsps <= num_paths:
        # Each LSP gets its own path
        return num_lsps * capacity_per_path
    else:
        # LSPs must share paths
        lsps_per_path = num_lsps / num_paths
        capacity_per_lsp = capacity_per_path / lsps_per_path
        return num_lsps * capacity_per_lsp
