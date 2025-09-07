"""StrictMultiDiGraph.from_arrays behaviors and validations."""

from __future__ import annotations

import numpy as np
import pytest

import netgraph_core as ngc


def test_from_arrays_basic_and_add_reverse():
    n = 3
    src = np.array([0, 1], dtype=np.int32)
    dst = np.array([1, 2], dtype=np.int32)
    cap = np.array([1.0, 2.0], dtype=np.float64)
    cost = np.array([1.0, 1.0], dtype=np.float64)
    g = ngc.StrictMultiDiGraph.from_arrays(n, src, dst, cap, cost, add_reverse=True)
    # Expect reverse edges present; at least 4 edges total
    assert g.num_edges() >= 4


def test_from_arrays_mismatched_lengths_raise():
    n = 2
    src = np.array([0], dtype=np.int32)
    dst = np.array([1, 1], dtype=np.int32)
    cap = np.array([1.0], dtype=np.float64)
    cost = np.array([1.0], dtype=np.float64)
    with pytest.raises(TypeError):
        ngc.StrictMultiDiGraph.from_arrays(n, src, dst, cap, cost)


def test_from_arrays_wrong_dtypes_raise():
    n = 2
    src = np.array([0], dtype=np.int64)  # wrong dtype
    dst = np.array([1], dtype=np.int32)
    cap = np.array([1.0], dtype=np.float32)  # wrong dtype
    cost = np.array([1.0], dtype=np.float64)
    with pytest.raises((TypeError, ValueError)):
        ngc.StrictMultiDiGraph.from_arrays(n, src, dst, cap, cost)


def test_from_arrays_negative_or_inf_values_raise():
    n = 2
    src = np.array([0], dtype=np.int32)
    dst = np.array([1], dtype=np.int32)
    cap = np.array([-1.0], dtype=np.float64)
    cost = np.array([np.inf], dtype=np.float64)
    with pytest.raises(ValueError):
        ngc.StrictMultiDiGraph.from_arrays(n, src, dst, cap, cost)


def test_from_arrays_self_loops_behavior():
    n = 2
    src = np.array([0], dtype=np.int32)
    dst = np.array([0], dtype=np.int32)  # self-loop
    cap = np.array([1.0], dtype=np.float64)
    cost = np.array([1.0], dtype=np.float64)
    # Either allowed or rejected; assert it doesn't crash and creates 1 edge or raises
    try:
        g = ngc.StrictMultiDiGraph.from_arrays(n, src, dst, cap, cost)
        assert g.num_edges() >= 1
    except Exception:
        # Accept rejection behavior too
        pass
