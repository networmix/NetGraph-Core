"""Typing and docstrings for netgraph_core public API.

This module is intentionally light; runtime implementations live in the
compiled extension. The goal is to provide type hints and help text.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from typing import TYPE_CHECKING, Optional

if TYPE_CHECKING:  # only for typing; runtime comes from extension
    from . import StrictMultiDiGraph  # type: ignore
import numpy as np


class EdgeSelect(Enum):
    ALL_MIN_COST = 1
    SINGLE_MIN_COST = 2
    ALL_MIN_COST_WITH_CAP_REMAINING = 3


class FlowPlacement(Enum):
    PROPORTIONAL = 1
    EQUAL_BALANCED = 2


class PredDAG:
    parent_offsets: np.ndarray
    parents: np.ndarray
    via_edges: np.ndarray


@dataclass(frozen=True)
class Path:
    nodes: np.ndarray
    edges: np.ndarray
    cost: float


def spf(
    g: "StrictMultiDiGraph",
    src: int,
    dst: Optional[int] = None,
    *,
    edge_select: EdgeSelect = EdgeSelect.ALL_MIN_COST,
    multipath: bool = True,
    node_mask: Optional[np.ndarray] = None,
    edge_mask: Optional[np.ndarray] = None,
    eps: float = 1e-10,
):
    """Shortest paths wrapper.

    Returns: (dist: float64[N], pred: PredDAG)
    """
    ...


def ksp(
    g: "StrictMultiDiGraph",
    src: int,
    dst: int,
    *,
    k: int,
    max_cost_factor: Optional[float] = None,
    unique: bool = True,
    node_mask: Optional[np.ndarray] = None,
    edge_mask: Optional[np.ndarray] = None,
    eps: float = 1e-10,
): ...


def calc_max_flow(
    g: "StrictMultiDiGraph",
    src: int,
    dst: int,
    *,
    flow_placement: FlowPlacement = FlowPlacement.PROPORTIONAL,
    shortest_path: bool = False,
    eps: float = 1e-10,
    with_edge_flows: bool = False,
    node_mask: Optional[np.ndarray] = None,
    edge_mask: Optional[np.ndarray] = None,
): ...
