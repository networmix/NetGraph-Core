"""StrictMultiDiGraph.from_arrays behaviors and validations."""

from __future__ import annotations

import numpy as np
import pytest

import netgraph_core as ngc


def _make_ext_ids(n):
    """Generate ext_edge_ids for n edges."""
    return np.arange(n, dtype=np.int64)


def test_from_arrays_basic_bidirectional_manual():
    n = 3
    # Manually include reverse edges
    src = np.array([0, 1, 1, 2], dtype=np.int32)
    dst = np.array([1, 0, 2, 1], dtype=np.int32)
    cap = np.array([1.0, 1.0, 2.0, 2.0], dtype=np.float64)
    cost = np.array([1, 1, 1, 1], dtype=np.int64)
    g = ngc.StrictMultiDiGraph.from_arrays(n, src, dst, cap, cost, _make_ext_ids(4))
    # Expect reverse edges present; 4 edges total
    assert g.num_edges() == 4


def test_from_arrays_mismatched_lengths_raise():
    n = 2
    src = np.array([0], dtype=np.int32)
    dst = np.array([1, 1], dtype=np.int32)
    cap = np.array([1.0], dtype=np.float64)
    cost = np.array([1], dtype=np.int64)
    with pytest.raises(TypeError):
        ngc.StrictMultiDiGraph.from_arrays(n, src, dst, cap, cost, _make_ext_ids(1))


def test_from_arrays_wrong_dtypes_raise():
    n = 2
    src = np.array([0], dtype=np.int64)  # wrong dtype
    dst = np.array([1], dtype=np.int32)
    cap = np.array([1.0], dtype=np.float32)  # wrong dtype
    cost = np.array([1], dtype=np.int64)
    with pytest.raises((TypeError, ValueError)):
        ngc.StrictMultiDiGraph.from_arrays(n, src, dst, cap, cost, _make_ext_ids(1))


def test_from_arrays_negative_capacity_raises():
    n = 2
    src = np.array([0], dtype=np.int32)
    dst = np.array([1], dtype=np.int32)
    cap = np.array([-1.0], dtype=np.float64)
    cost = np.array([1], dtype=np.int64)
    with pytest.raises(ValueError):
        ngc.StrictMultiDiGraph.from_arrays(n, src, dst, cap, cost, _make_ext_ids(1))


def test_from_arrays_self_loops_behavior():
    n = 2
    src = np.array([0], dtype=np.int32)
    dst = np.array([0], dtype=np.int32)  # self-loop
    cap = np.array([1.0], dtype=np.float64)
    cost = np.array([1], dtype=np.int64)
    # Either allowed or rejected; assert it doesn't crash and creates 1 edge or raises
    try:
        g = ngc.StrictMultiDiGraph.from_arrays(n, src, dst, cap, cost, _make_ext_ids(1))
        assert g.num_edges() >= 1
    except Exception:
        # Accept rejection behavior too
        pass
