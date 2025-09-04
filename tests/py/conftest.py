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


@pytest.fixture
def build_graph():
    """Return a function that builds a StrictMultiDiGraph from edges.

    The returned function accepts:
      - num_nodes: int
      - edges: list of tuples (src, dst, cost, cap[, ext_id])
      - add_reverse: bool (default False)
    """

    def _builder(
        num_nodes: int, edges: list[tuple], *, add_reverse: bool = False
    ) -> ngc.StrictMultiDiGraph:
        if not edges:
            # Build an empty graph quickly
            return ngc.StrictMultiDiGraph.from_arrays(
                num_nodes,
                np.empty((0,), dtype=np.int32),
                np.empty((0,), dtype=np.int32),
                np.empty((0,), dtype=np.float64),
                np.empty((0,), dtype=np.float64),
                add_reverse=add_reverse,
            )
        src = np.array([e[0] for e in edges], dtype=np.int32)
        dst = np.array([e[1] for e in edges], dtype=np.int32)
        cost = np.array([float(e[2]) for e in edges], dtype=np.float64)
        cap = np.array([float(e[3]) for e in edges], dtype=np.float64)
        # External ids at this layer are ignored by core and used only in tests
        return ngc.StrictMultiDiGraph.from_arrays(
            num_nodes, src, dst, cap, cost, add_reverse=add_reverse
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
    """Graph3: six-node mixed topology used across SPF/max-flow tests.

    Nodes: 0..5; includes parallel min-cost edges and a longer route via node 5.
    """

    edges = [
        (0, 1, 1, 2, 0),  # A->B parallels
        (0, 1, 1, 4, 1),
        (0, 1, 1, 6, 2),
        (1, 2, 1, 1, 3),  # B->C parallels
        (1, 2, 1, 2, 4),
        (1, 2, 1, 3, 5),
        (2, 3, 2, 3, 6),  # C->D cost 2
        (0, 4, 1, 5, 7),  # A->E cost 1
        (4, 2, 1, 4, 8),  # E->C cost 1
        (0, 3, 4, 2, 9),  # A->D cost 4
        (2, 5, 1, 1, 10),  # C->F cost 1
        (5, 3, 1, 2, 11),  # F->D cost 1
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
    """Small five-node graph with cross edges between B and C.

    Nodes: 0(A),1(B),2(C),3(D),4(E)
    """

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
    """Branched five-node graph with C/D to E.

    Nodes: 0(A),1(B),2(C),3(D),4(E)
    """

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
