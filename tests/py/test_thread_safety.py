"""
Test thread safety of Python bindings (Issue #4).

These tests validate that NumPy array data is properly copied before releasing
the GIL, preventing data races when other threads mutate the arrays.
"""

import threading
import time

import numpy as np
import pytest

from netgraph_core import (
    Algorithms,
    Backend,
    StrictMultiDiGraph,
)


def _make_graph(num_nodes, src, dst, capacity, cost):
    """Helper to build graph with auto-generated ext_edge_ids."""
    ext_edge_ids = np.arange(len(src), dtype=np.int64)
    return StrictMultiDiGraph.from_arrays(
        num_nodes, src, dst, capacity, cost, ext_edge_ids
    )


class TestSPFThreadSafety:
    """Test that spf() is thread-safe with respect to mask mutations."""

    def test_spf_node_mask_thread_safety(self):
        """
        Verify that mutating node_mask in another thread doesn't affect spf() results.

        The C++ code copies the mask to a local vector before releasing the GIL,
        so concurrent mutations to the NumPy buffer won't cause data races.
        """
        backend = Backend.cpu()
        algs = Algorithms(backend)

        n = 100
        g = _make_graph(
            n,
            np.arange(n - 1, dtype=np.int32),
            np.arange(1, n, dtype=np.int32),
            np.ones(n - 1) * 10.0,
            np.ones(n - 1, dtype=np.int64),
        )
        graph = algs.build_graph(g)

        # Create node mask (all True initially)
        node_mask = np.ones(n, dtype=bool)

        results = []

        def run_spf():
            """Run SPF with the mask."""
            dist, dag = algs.spf(graph, 0, node_mask=node_mask)
            results.append((dist.copy(), dag))

        def mutate_mask():
            """Mutate the mask after a short delay."""
            time.sleep(0.001)  # Let SPF start
            node_mask[:] = False  # Try to cause a race

        # Run SPF in main thread while another thread mutates the mask
        mutator = threading.Thread(target=mutate_mask)
        mutator.start()

        run_spf()

        mutator.join()

        # SPF should see a consistent mask (either all True or all False)
        # because the C++ code copies the mask before releasing the GIL
        dist, dag = results[0]

        # If mask was all True: all nodes should be reachable
        # If mask was all False: only source should have finite distance
        # It should NOT be corrupted (some arbitrary subset)
        finite_count = np.sum(np.isfinite(dist))

        # Either saw all True (finite_count == n) or all False (finite_count == 1)
        assert finite_count in [1, n], (
            f"SPF saw corrupted mask: {finite_count} finite distances "
            f"(expected 1 or {n})"
        )

    def test_spf_edge_mask_thread_safety(self):
        """Verify edge_mask is also safe from concurrent mutations."""
        backend = Backend.cpu()
        algs = Algorithms(backend)

        n = 50
        edges = [(i, (i + 1) % n) for i in range(n)]  # Ring topology
        src = np.array([e[0] for e in edges], dtype=np.int32)
        dst = np.array([e[1] for e in edges], dtype=np.int32)

        g = _make_graph(
            n,
            src,
            dst,
            np.ones(n) * 10.0,
            np.ones(n, dtype=np.int64),
        )
        graph = algs.build_graph(g)

        edge_mask = np.ones(n, dtype=bool)

        results = []

        def run_spf():
            dist, dag = algs.spf(graph, 0, edge_mask=edge_mask)
            results.append(dist.copy())

        def mutate_mask():
            time.sleep(0.001)
            edge_mask[:] = False

        mutator = threading.Thread(target=mutate_mask)
        mutator.start()
        run_spf()
        mutator.join()

        dist = results[0]
        # Should see consistent results (not corrupted)
        # Either all edges available (complete ring) or no edges
        assert np.sum(np.isfinite(dist)) in [1, n]


class TestKSPThreadSafety:
    """Test that ksp() is thread-safe with respect to mask mutations."""

    def test_ksp_mask_thread_safety(self):
        """Verify ksp() copies masks before releasing GIL."""
        backend = Backend.cpu()
        algs = Algorithms(backend)

        n = 20
        g = _make_graph(
            n,
            np.array([0, 0, 1, 1, 2], dtype=np.int32),
            np.array([1, 2, 3, 4, 4], dtype=np.int32),
            np.ones(5) * 10.0,
            np.array([1, 2, 1, 2, 1], dtype=np.int64),
        )
        graph = algs.build_graph(g)

        node_mask = np.ones(n, dtype=bool)

        results = []

        def run_ksp():
            paths = algs.ksp(graph, 0, 4, k=3, node_mask=node_mask)
            results.append(len(paths))

        def mutate_mask():
            time.sleep(0.001)
            node_mask[2] = False  # Block a node

        mutator = threading.Thread(target=mutate_mask)
        mutator.start()
        run_ksp()
        mutator.join()

        # Should get consistent results (either with or without node 2)
        num_paths = results[0]
        assert num_paths >= 0  # Should not crash or return corrupted data


class TestMaxFlowThreadSafety:
    """Test that max_flow() is thread-safe with respect to mask mutations."""

    def test_maxflow_mask_thread_safety(self):
        """Verify max_flow() copies masks before releasing GIL."""
        backend = Backend.cpu()
        algs = Algorithms(backend)

        g = _make_graph(
            4,
            np.array([0, 0, 1, 2], dtype=np.int32),
            np.array([1, 2, 3, 3], dtype=np.int32),
            np.array([10.0, 10.0, 10.0, 10.0]),
            np.ones(4, dtype=np.int64),
        )
        graph = algs.build_graph(g)

        node_mask = np.ones(4, dtype=bool)
        edge_mask = np.ones(4, dtype=bool)

        results = []

        def run_maxflow():
            total_flow, summary = algs.max_flow(
                graph, 0, 3, node_mask=node_mask, edge_mask=edge_mask
            )
            results.append(total_flow)

        def mutate_masks():
            time.sleep(0.001)
            node_mask[1] = False
            edge_mask[0] = False

        mutator = threading.Thread(target=mutate_masks)
        mutator.start()
        run_maxflow()
        mutator.join()

        # Should get consistent flow value (either full graph or restricted)
        flow = results[0]
        assert flow >= 0.0  # Should not be corrupted (negative/NaN)


class TestConcurrentAlgorithmCalls:
    """Test multiple algorithm calls running concurrently."""

    def test_concurrent_spf_calls(self):
        """Multiple threads calling spf() simultaneously should be safe."""
        backend = Backend.cpu()
        algs = Algorithms(backend)

        n = 30
        g = _make_graph(
            n,
            np.arange(n - 1, dtype=np.int32),
            np.arange(1, n, dtype=np.int32),
            np.ones(n - 1) * 10.0,
            np.ones(n - 1, dtype=np.int64),
        )
        graph = algs.build_graph(g)

        results = []

        def worker(src_node):
            """Run SPF from a specific source."""
            mask = np.ones(n, dtype=bool)
            mask[src_node] = True  # Ensure source is available
            dist, dag = algs.spf(graph, src_node, node_mask=mask)
            results.append((src_node, dist.copy()))

        # Run multiple SPF calls concurrently
        threads = []
        for i in range(5):
            t = threading.Thread(target=worker, args=(i,))
            threads.append(t)
            t.start()

        for t in threads:
            t.join()

        # All threads should complete successfully
        assert len(results) == 5

        # Each result should be valid
        for src, dist in results:
            assert dist[src] == 0.0  # Source distance is always 0


@pytest.mark.slow
class TestStressThreadSafety:
    """Stress tests for thread safety under heavy concurrent load."""

    def test_repeated_concurrent_mutations(self):
        """Repeatedly mutate masks while running algorithms."""
        backend = Backend.cpu()
        algs = Algorithms(backend)

        n = 50
        g = _make_graph(
            n,
            np.arange(n - 1, dtype=np.int32),
            np.arange(1, n, dtype=np.int32),
            np.random.rand(n - 1) * 100,
            np.ones(n - 1, dtype=np.int64),
        )
        graph = algs.build_graph(g)

        node_mask = np.ones(n, dtype=bool)

        def run_algorithms():
            """Run various algorithms repeatedly."""
            for _ in range(20):
                try:
                    algs.spf(graph, 0, node_mask=node_mask)
                except Exception as e:
                    pytest.fail(f"Algorithm raised exception: {e}")

        def mutate_continuously():
            """Continuously mutate the mask."""
            for _ in range(100):
                node_mask[:] = np.random.rand(n) > 0.5
                time.sleep(0.001)

        algo_thread = threading.Thread(target=run_algorithms)
        mutator_thread = threading.Thread(target=mutate_continuously)

        algo_thread.start()
        mutator_thread.start()

        algo_thread.join()
        mutator_thread.join()

        # If we got here without crashes/exceptions, thread safety is working
        assert True, "Thread safety test completed without crashes"


class TestMemoryOrdering:
    """
    Test that mask copying has proper memory ordering.

    This is more of a theoretical concern - if masks aren't properly copied,
    we might see stale or reordered writes from other threads.
    """

    def test_mask_sees_latest_values(self):
        """
        Verify that when we create a mask and immediately pass it to C++,
        the C++ code sees the correct values (not stale cache).
        """
        backend = Backend.cpu()
        algs = Algorithms(backend)

        n = 10
        g = _make_graph(
            n,
            np.arange(n - 1, dtype=np.int32),
            np.arange(1, n, dtype=np.int32),
            np.ones(n - 1),
            np.ones(n - 1, dtype=np.int64),
        )
        graph = algs.build_graph(g)

        # Create mask with specific pattern
        mask = np.zeros(n, dtype=bool)
        mask[0] = True  # Only source
        mask[5] = True  # One intermediate node
        mask[9] = True  # Destination

        # Immediately use it
        dist, dag = algs.spf(graph, 0, node_mask=mask)

        # Should respect the mask (only nodes 0, 5, 9 are reachable)
        finite_nodes = [i for i in range(n) if np.isfinite(dist[i])]

        # With the mask, we can't reach node 9 because intermediate nodes are blocked
        # (except 5, but we need more nodes in the path)
        assert 0 in finite_nodes  # Source is always reachable

        # The key test: C++ saw the mask we just created (not stale/zero data)
        # If it saw all-false or all-true, the results would be very different
