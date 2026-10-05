"""Typing and docstrings for netgraph_core public API.

This module is intentionally light; runtime implementations live in the
compiled extension. The goal is to provide type hints and help text.

TODO: replace this hand-written stub with generated `.pyi` (e.g. via
pybind11-stubgen) to prevent drift from the runtime bindings (Phase 5).
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Sequence, overload

if TYPE_CHECKING:  # only for typing; runtime comes from extension
    import numpy as np  # type: ignore[reportMissingImports]
    # Use forward references (strings) below; avoid importing from the package
    # to prevent circular imports during type checking.


@overload
def to_units(x: float, resolution: float) -> int: ...
@overload
def to_units(x: "Sequence[float] | np.ndarray", resolution: float) -> "np.ndarray": ...
def to_units(x, resolution):  # type: ignore[no-untyped-def]
    """Convert physical value(s) to integer units of the resolution, rounding half up.

    The only way physical numbers enter Core: ``to_units(x, q) = floor(x / q + 1/2)``.

    Args:
        x: A physical float, or a 1-D array of them.
        resolution: The unit ``q``, a power of two in ``[2**-32, 2**32]``.

    Returns:
        An ``int`` for a scalar, an ``int64`` array for an array.

    Raises:
        ValueError: an invalid resolution (checked first, also for an empty array);
            NaN, infinity or a negative value (no unit value); a result at or above
            ``2**62`` units (the domain bound).
    """
    ...


@overload
def as_float(x: int, resolution: float) -> float: ...
@overload
def as_float(x: "Sequence[int] | np.ndarray", resolution: float) -> "np.ndarray": ...
def as_float(x, resolution):  # type: ignore[no-untyped-def]
    """Convert integer units to physical value(s): ``x * q``, correctly rounded.

    Exact below ``2**53`` units. Returns a ``float`` for an int of any size and a
    ``float64`` array for an integer array or sequence. The resolution is validated first
    (also for an empty array); a result beyond the finite doubles raises ValueError.
    """
    ...


class Eligibility:
    """Immutable, graph-bound static eligibility: dense node and edge flags.

    Built by ``StrictMultiDiGraph.eligibility(excluded_nodes, excluded_edges)``. An
    edge is allowed when it is not excluded and both of its endpoints are allowed.
    Passing it to a call of another graph raises ValueError.
    """

    def graph_uid(self) -> int:
        """uid of the graph this value belongs to."""
        ...

    def uid(self) -> int:
        """0 for the all-allowed value, a fresh uid otherwise."""
        ...

    def node_allowed(self) -> "np.ndarray":
        """bool[N] read-only view."""
        ...

    def edge_allowed(self) -> "np.ndarray":
        """bool[E] read-only view."""
        ...


class StrictMultiDiGraph:
    """Immutable directed multigraph on integer units (layout B').

    Edges are given in input order. The public EdgeId of an edge is its rank in the
    stable order ``(src, cost, dst, input)``: the out-edges of node ``u`` are the ids
    ``[row_offsets()[u], row_offsets()[u + 1])``, and ``input_index()[e]`` is the input
    position of edge ``e`` (the inverse permutation is the caller's). Every other
    per-edge array is in id order. Accessors return read-only NumPy views that keep the
    graph alive.
    """

    def __init__(
        self,
        num_nodes: int,
        src: "Sequence[int] | np.ndarray",
        dst: "Sequence[int] | np.ndarray",
        capacity_units: "Sequence[int] | np.ndarray",
        cost: "Sequence[int] | np.ndarray",
        resolution: float,
    ) -> None:
        """Build a graph from parallel edge arrays in input order.

        Args:
            num_nodes: Node count; every src/dst must lie in ``[0, num_nodes)``.
            src: Source node ids (integers).
            dst: Destination node ids (integers).
            capacity_units: Capacities as integer units of ``resolution``, each in
                ``[0, 2**62)``; physical capacities go through ``to_units`` first.
            cost: Positive integer costs (``cost >= 1``).
            resolution: The unit ``q``, a power of two in ``[2**-32, 2**32]``.

        Raises:
            ValueError: On the first defect, in this order: the resolution, the unit
                values, the array shapes, zero or negative costs, the remaining bounds
                (endpoints in range, then the capacity and cost sums below ``2**62``).
                Integers are read without narrowing, so the reported defect does not
                depend on how the input is represented; floats or bools in an integer
                array are rejected, never truncated.
        """
        ...

    def uid(self) -> int:
        """Process-unique id, never reused."""
        ...

    def resolution(self) -> float: ...
    def num_nodes(self) -> int: ...
    def num_edges(self) -> int: ...
    def max_cost(self) -> int:
        """Largest edge cost (0 for a graph without edges)."""
        ...

    def cost_gcd(self) -> int:
        """Greatest common divisor of all edge costs (0 for a graph without edges)."""
        ...

    def capacity_units(self) -> "np.ndarray":
        """int64[E] capacities in units, id order."""
        ...

    def cost(self) -> "np.ndarray":
        """int64[E] costs, id order."""
        ...

    def edge_src(self) -> "np.ndarray":
        """int32[E] source node per edge."""
        ...

    def edge_dst(self) -> "np.ndarray":
        """int32[E] destination node per edge."""
        ...

    def input_index(self) -> "np.ndarray":
        """int32[E]: id -> input position."""
        ...

    def row_offsets(self) -> "np.ndarray":
        """int32[N + 1] forward CSR offsets (id ranges of the out-edges)."""
        ...

    def in_row_offsets(self) -> "np.ndarray":
        """int32[N + 1] reverse CSR offsets."""
        ...

    def in_adj_edge_index(self) -> "np.ndarray":
        """int32[E] in-edges per node in (cost, src, input) order, as EdgeIds."""
        ...

    def eligibility(
        self,
        excluded_nodes: "Sequence[int] | np.ndarray" = (),
        excluded_edges: "Sequence[int] | np.ndarray" = (),
    ) -> Eligibility:
        """Static eligibility from excluded node ids and EdgeIds.

        Duplicates and unsorted lists are accepted; out-of-range ids raise ValueError.
        Empty lists mean all allowed (uid 0).
        """
        ...
