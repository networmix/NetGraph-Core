"""Oracle-free correctness certificates for max-flow.

The suite historically checked max-flow results against hand-computed values on a
handful of fixed topologies. That misses defects whose trigger is a graph nobody
thought to write down -- the 0.7.2 under-reporting bug (fixed in 0.8.0) was not
exposed by any fixture in the suite, at any src/dst pair.

These tests take the opposite approach: generate graphs in the region where such
defects live (dispersed edge costs, full-duplex links) and check a *certificate of
optimality* rather than a recorded number. See the ``certify_max_flow`` fixture for
the four properties; together they prove maximality, so no external oracle and no
stored expected value is required.
"""

from __future__ import annotations

import numpy as np
import pytest

import netgraph_core as ngc

# A fixed seed keeps these deterministic, so a given defect is either caught on every
# run or never -- which means the corpus has to be shown to have teeth, not assumed to.
# This seed was picked because its corpus trips the 0.7.2 under-reporting bug within the
# first few graphs; both sweeps below were confirmed to fail against 0.7.2 and pass
# against 0.8.0. An arbitrary seed needed ~5750 graphs to hit the same defect, so this
# is purely a cheaper corpus, not a narrower one.
SEED = 198
QUICK_GRAPHS = 2000
THOROUGH_GRAPHS = 20000


def _random_graph(rng: np.random.Generator) -> tuple[ngc.StrictMultiDiGraph, int, int]:
    """Build a random full-duplex graph with dispersed costs.

    Cost dispersion is the essential ingredient: with a single cost tier there is
    nothing for shortest-path tier ordering to get wrong, and randomized checks over
    equal-cost graphs find nothing.
    """
    n = int(rng.integers(4, 9))
    m = int(rng.integers(n, 2 * n))
    links = []
    for _ in range(m):
        u, v = rng.choice(n, 2, replace=False)
        cap = float(rng.integers(1, 6))
        cost = int(rng.integers(1, 21))
        links.append((int(u), int(v), cap, cost))

    edges = [e for u, v, c, k in links for e in ((u, v, c, k), (v, u, c, k))]
    src = np.array([e[0] for e in edges], dtype=np.int32)
    dst = np.array([e[1] for e in edges], dtype=np.int32)
    cap = np.array([e[2] for e in edges], dtype=np.float64)
    cost = np.array([e[3] for e in edges], dtype=np.int64)
    g = ngc.StrictMultiDiGraph.from_arrays(
        n, src, dst, cap, cost, np.arange(len(edges), dtype=np.int64)
    )
    return g, 0, n - 1


def _sweep(count: int, algs, certify_max_flow) -> None:
    rng = np.random.default_rng(SEED)
    for i in range(count):
        g, s, t = _random_graph(rng)
        total, summary = algs.max_flow(algs.build_graph(g), s, t, with_edge_flows=True)
        try:
            certify_max_flow(g, summary, total, s, t)
        except AssertionError as exc:  # pragma: no cover - only on a real defect
            pytest.fail(f"certificate failed on random graph #{i} (seed {SEED}): {exc}")


def test_random_graphs_certify(algs, certify_max_flow):
    """Every generated max-flow result must carry a valid optimality certificate."""
    _sweep(QUICK_GRAPHS, algs, certify_max_flow)


@pytest.mark.slow
def test_random_graphs_certify_thorough(algs, certify_max_flow):
    """Wider sweep of the same generator."""
    _sweep(THOROUGH_GRAPHS, algs, certify_max_flow)


FIXTURES = [
    "line1_graph",
    "square1_graph",
    "square2_graph",
    "graph3",
    "square4_graph",
    "triangle1_graph",
]


@pytest.mark.parametrize("fixture_name", FIXTURES)
def test_named_fixtures_certify_all_pairs(
    fixture_name, request, algs, certify_max_flow
):
    """Certify the shared fixtures for every src/dst pair, not just the tested one.

    A fixture is usually queried for a single pair; the remaining pairs are free
    coverage of the same topology.
    """
    g = request.getfixturevalue(fixture_name)
    handle = algs.build_graph(g)
    n = g.num_nodes()
    for s in range(n):
        for t in range(n):
            if s == t:
                continue
            total, summary = algs.max_flow(handle, s, t, with_edge_flows=True)
            certify_max_flow(g, summary, total, s, t)


def test_certificate_rejects_a_non_maximal_flow(square1_graph, algs, certify_max_flow):
    """The certificate must have teeth.

    ``shortest_path=True`` performs a single augmentation and is deliberately not
    maximal, so it is a genuine non-maximal result reachable through the public API.
    Certifying it as maximal must fail. Without this guard the certificate could
    silently degrade into an assertion that cannot fail -- which is exactly how the
    previous ``assert_valid_min_cut`` helper missed the 0.7.2 bug.
    """
    g = square1_graph
    total, summary = algs.max_flow(
        algs.build_graph(g), 0, 2, shortest_path=True, with_edge_flows=True
    )
    assert total == pytest.approx(1.0)

    # Feasibility and conservation still hold for a partial placement.
    certify_max_flow(g, summary, total, 0, 2, maximal=False)

    with pytest.raises(AssertionError):
        certify_max_flow(g, summary, total, 0, 2, maximal=True)


def test_equal_balanced_is_feasible_but_not_certified_maximal(
    line1_graph, algs, certify_max_flow
):
    """EQUAL_BALANCED is an admission model, so only feasibility/conservation apply."""
    g = line1_graph
    total, summary = algs.max_flow(
        algs.build_graph(g),
        0,
        2,
        flow_placement=ngc.FlowPlacement.EQUAL_BALANCED,
        with_edge_flows=True,
    )
    certify_max_flow(g, summary, total, 0, 2, maximal=False)


def test_masked_runs_certify(algs, certify_max_flow):
    """Masking must not break the certificate, and must not leak flow onto masked edges.

    Masks are a good place for a future defect of this kind to hide: no external
    solver models them, so an oracle-based test cannot reach them at all.

    Node and edge masks are checked separately and together, because they are
    separate branches in ``calc_max_flow`` -- including inside the residual
    completion phase added in 0.8.0, whose forward and reverse arc loops each
    consult ``edge_mask`` independently. A node-mask-only sweep leaves those
    branches, and the certificate's own ``edge_mask`` handling, unexercised.
    """
    rng = np.random.default_rng(SEED + 1)
    saw_masked_node = False
    saw_masked_edge = False

    for i in range(300):
        g, s, t = _random_graph(rng)
        n = g.num_nodes()

        node_mask = None
        edge_mask = None
        mode = i % 3  # 0: nodes only, 1: edges only, 2: both
        if mode in (0, 2):
            node_mask = np.ones(n, dtype=bool)
            victim = int(rng.integers(0, n))
            if victim not in (s, t):
                node_mask[victim] = False
                saw_masked_node = True
        if mode in (1, 2):
            edge_mask = rng.random(g.num_edges()) > 0.25
            if not edge_mask.all():
                saw_masked_edge = True

        total, summary = algs.max_flow(
            algs.build_graph(g),
            s,
            t,
            node_mask=node_mask,
            edge_mask=edge_mask,
            with_edge_flows=True,
        )
        certify_max_flow(
            g, summary, total, s, t, node_mask=node_mask, edge_mask=edge_mask
        )

    # Guard against the sweep silently degrading into an unmasked one.
    assert saw_masked_node, "node-mask branch was never exercised"
    assert saw_masked_edge, "edge-mask branch was never exercised"
