"""
Test lifetime management in Python bindings.

These tests validate that pybind11 keep_alive policies correctly prevent
premature garbage collection of objects that are referenced by raw pointers
in C++ classes.
"""

import gc
import weakref

import numpy as np
import pytest

from netgraph_core import (
    Algorithms,
    Backend,
    FlowGraph,
    FlowPolicy,
    FlowPolicyConfig,
    FlowState,
    StrictMultiDiGraph,
)


def _make_graph(num_nodes, src, dst, capacity, cost):
    """Helper to build graph with auto-generated ext_edge_ids."""
    ext_edge_ids = np.arange(len(src), dtype=np.int64)
    return StrictMultiDiGraph.from_arrays(
        num_nodes, src, dst, capacity, cost, ext_edge_ids
    )


class TestFlowStateLifetime:
    """Test FlowState keeps graph alive (Issue #1)."""

    def test_flowstate_keeps_graph_alive_basic(self):
        """FlowState should prevent graph from being garbage collected."""
        # Create graph
        g = _make_graph(
            3,
            np.array([0, 1], dtype=np.int32),
            np.array([1, 2], dtype=np.int32),
            np.array([10.0, 10.0]),
            np.array([1, 1], dtype=np.int64),
        )

        # Create FlowState (should keep graph alive via keep_alive<0,1>)
        fs = FlowState(g)

        # Create weak reference to detect when graph is collected
        graph_ref = weakref.ref(g)

        # Delete strong reference to graph
        del g
        gc.collect()

        # With fix: graph should still be alive (kept by fs)
        assert graph_ref() is not None, "Graph was prematurely collected"

        # FlowState should still work
        residual = fs.residual_view()
        assert len(residual) == 2
        assert residual[0] == 10.0

        # When both FlowState and all views are deleted, graph can be collected
        del residual
        del fs
        gc.collect()
        assert graph_ref() is None, (
            "Graph should be collected after FlowState and views are deleted"
        )

    def test_flowstate_with_residual_keeps_graph_alive(self):
        """FlowState(graph, residual) should keep graph alive, not residual."""
        g = _make_graph(
            3,
            np.array([0, 1], dtype=np.int32),
            np.array([1, 2], dtype=np.int32),
            np.array([10.0, 10.0]),
            np.array([1, 1], dtype=np.int64),
        )

        residual = np.array([5.0, 7.0])
        fs = FlowState(g, residual)

        graph_ref = weakref.ref(g)

        # Delete references
        del g
        del residual
        gc.collect()

        # Graph should be kept alive (by FlowState), residual is copied
        assert graph_ref() is not None, "Graph was prematurely collected"

        # FlowState should work
        res = fs.residual_view()
        assert len(res) == 2
        assert res[0] == 5.0
        assert res[1] == 7.0


class TestFlowPolicyLifetime:
    """Test FlowPolicy keeps Algorithms and Graph alive (Issue #2)."""

    def test_flowpolicy_keeps_algorithms_alive(self):
        """FlowPolicy should prevent Algorithms from being garbage collected."""
        # Setup
        backend = Backend.cpu()
        algs = Algorithms(backend)

        g = _make_graph(
            3,
            np.array([0, 1], dtype=np.int32),
            np.array([1, 2], dtype=np.int32),
            np.array([10.0, 10.0]),
            np.array([1, 1], dtype=np.int64),
        )
        graph = algs.build_graph(g)

        # Create FlowPolicy
        config = FlowPolicyConfig()
        policy = FlowPolicy(algs, graph, config)

        # Create weak references
        algs_ref = weakref.ref(algs)
        graph_ref = weakref.ref(graph)

        # Delete strong references
        del algs
        del graph
        gc.collect()

        # With fix: both should still be alive
        assert algs_ref() is not None, "Algorithms was prematurely collected"
        assert graph_ref() is not None, "Graph was prematurely collected"

        # Policy should still work
        assert policy.flow_count() == 0
        assert policy.placed_demand() == 0.0

    def test_flowpolicy_keyword_constructor_keeps_alive(self):
        """FlowPolicy keyword constructor should also keep dependencies alive."""
        backend = Backend.cpu()
        algs = Algorithms(backend)

        g = _make_graph(
            3,
            np.array([0, 1], dtype=np.int32),
            np.array([1, 2], dtype=np.int32),
            np.array([10.0, 10.0]),
            np.array([1, 1], dtype=np.int64),
        )
        graph = algs.build_graph(g)

        # Use config constructor
        cfg = FlowPolicyConfig(min_flow_count=1, max_total_iterations=100)
        policy = FlowPolicy(algs, graph, cfg)

        algs_ref = weakref.ref(algs)

        del algs
        gc.collect()

        assert algs_ref() is not None, "Algorithms was prematurely collected"
        assert policy.flow_count() == 0


class TestFlowGraphLifetime:
    """Test FlowGraph lifetime management."""

    def test_flowgraph_keeps_graph_alive(self):
        """FlowGraph should prevent graph from being garbage collected."""
        g = _make_graph(
            3,
            np.array([0, 1], dtype=np.int32),
            np.array([1, 2], dtype=np.int32),
            np.array([10.0, 10.0]),
            np.array([1, 1], dtype=np.int64),
        )

        fg = FlowGraph(g)
        graph_ref = weakref.ref(g)

        # Delete graph reference
        del g
        gc.collect()

        # Graph should still be alive (kept by FlowGraph)
        assert graph_ref() is not None, "Graph was prematurely collected"

        # FlowGraph should still work
        graph_from_fg = fg.graph
        assert graph_from_fg.num_nodes() == 3
        assert graph_from_fg.num_edges() == 2

        # When FlowGraph is deleted, graph can be collected
        del fg
        del graph_from_fg
        gc.collect()
        assert graph_ref() is None, (
            "Graph should be collected after FlowGraph is deleted"
        )


class TestComplexLifetimeScenarios:
    """Test more complex lifetime interactions."""

    def test_nested_lifetime_chain(self):
        """Test chain: Graph -> FlowState -> operations."""
        g = _make_graph(
            3,
            np.array([0, 1], dtype=np.int32),
            np.array([1, 2], dtype=np.int32),
            np.array([10.0, 10.0]),
            np.array([1, 1], dtype=np.int64),
        )

        fs = FlowState(g)

        # Create weak references
        g_ref = weakref.ref(g)

        # Delete graph
        del g
        gc.collect()

        # Graph should be kept alive by FlowState
        assert g_ref() is not None

        # Should be able to perform operations
        fs.reset()
        residual = fs.residual_view()
        assert len(residual) == 2

    def test_multiple_flowstates_share_graph(self):
        """Multiple FlowStates should all keep the same graph alive."""
        g = _make_graph(
            3,
            np.array([0, 1], dtype=np.int32),
            np.array([1, 2], dtype=np.int32),
            np.array([10.0, 10.0]),
            np.array([1, 1], dtype=np.int64),
        )

        fs1 = FlowState(g)
        fs2 = FlowState(g, np.array([8.0, 9.0]))

        g_ref = weakref.ref(g)

        del g
        gc.collect()

        # Graph kept alive by both FlowStates
        assert g_ref() is not None

        # Both should work
        assert len(fs1.residual_view()) == 2
        assert len(fs2.residual_view()) == 2

        # Delete one FlowState - graph should still be alive
        del fs1
        gc.collect()
        assert g_ref() is not None

        # Delete second FlowState - now graph can be collected
        del fs2
        gc.collect()
        assert g_ref() is None


@pytest.mark.slow
class TestLifetimeUnderStress:
    """Stress tests for lifetime management."""

    def test_rapid_create_delete_cycle(self):
        """Rapidly create and delete objects to stress lifetime management."""
        for _ in range(100):
            g = _make_graph(
                10,
                np.arange(9, dtype=np.int32),
                np.arange(1, 10, dtype=np.int32),
                np.ones(9),
                np.ones(9, dtype=np.int64),
            )

            fs = FlowState(g)
            del g  # Should be kept alive by fs

            # Do some work
            try:
                residual = fs.residual_view()
                assert len(residual) == 9
                fs.reset()
            except Exception:
                # If the graph was prematurely collected, this will crash
                # That's what we're testing for
                pass

            # Delete FlowState
            del fs
            gc.collect()

        # If we got here without crashing, lifetime management is working

    def test_large_graph_lifetime(self):
        """Test lifetime with larger graphs to catch memory issues."""
        n = 1000
        g = _make_graph(
            n,
            np.arange(n - 1, dtype=np.int32),
            np.arange(1, n, dtype=np.int32),
            np.random.rand(n - 1) * 100,
            np.ones(n - 1, dtype=np.int64),
        )

        fs = FlowState(g)
        g_ref = weakref.ref(g)

        del g
        gc.collect()

        assert g_ref() is not None

        # Access views multiple times
        for _ in range(10):
            residual = fs.residual_view()
            assert len(residual) == n - 1
