"""Typing and docstrings for netgraph_core public API.

This module is intentionally light; runtime implementations live in the
compiled extension. The goal is to provide type hints and help text.

TODO: replace this hand-written stub with generated `.pyi` (e.g. via
pybind11-stubgen) to prevent drift from the runtime bindings (Phase 5).
"""

from __future__ import annotations

from typing import TYPE_CHECKING, ClassVar, Sequence, overload

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

    @staticmethod
    def intersect(a: "Eligibility", b: "Eligibility") -> "Eligibility":
        """The intersection of two values of one graph (ValueError for different graphs).

        Returns an operand itself when the result equals it (``a`` first, then ``b``),
        otherwise a new value with a fresh uid; equal intersections built separately get
        distinct uids, so a consumer builds each state once and reuses it.
        """
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


class EdgeTieBreak:
    """Ties between equal-cost candidates, at both tie levels (parallel members and walks).

    ``DETERMINISTIC`` takes the lowest EdgeId and ignores capacities;
    ``PREFER_HIGHER_PRIORITY`` prefers the higher ``priority`` value first.
    """

    DETERMINISTIC: ClassVar["EdgeTieBreak"]
    PREFER_HIGHER_PRIORITY: ClassVar["EdgeTieBreak"]
    __members__: ClassVar[dict[str, "EdgeTieBreak"]]

    def __init__(self, value: int) -> None: ...
    def __int__(self) -> int: ...
    def __index__(self) -> int: ...
    @property
    def name(self) -> str: ...
    @property
    def value(self) -> int: ...


class EdgePriority:
    """The per-edge value ``PREFER_HIGHER_PRIORITY`` reads, and the Weighted kernel's
    membership and weights: ``RESIDUAL`` (residual units), ``CAPACITY`` (static capacity
    units), ``HEADROOM`` (the exact rational residual / capacity; rejected with WEIGHTED).
    A search given no residual reads static capacities for every priority."""

    RESIDUAL: ClassVar["EdgePriority"]
    CAPACITY: ClassVar["EdgePriority"]
    HEADROOM: ClassVar["EdgePriority"]
    __members__: ClassVar[dict[str, "EdgePriority"]]

    def __init__(self, value: int) -> None: ...
    def __int__(self) -> int: ...
    def __index__(self) -> int: ...
    @property
    def name(self) -> str: ...
    @property
    def value(self) -> int: ...


class FlowPlacement:
    """The kernel of one flow (APPROACH 7).

    ``MAX_FLOW``: the minimum-cost flow core on the flow's split graph, a replacement on
    ``V = W + own`` (its own allocation cancellable, other owners never). ``WEIGHTED``:
    hop-by-hop WCMP, highest-averages apportionment over weights frozen at creation.
    ``EQUAL_BALANCED``: ECMP, round robin at every node.
    """

    MAX_FLOW: ClassVar["FlowPlacement"]
    WEIGHTED: ClassVar["FlowPlacement"]
    EQUAL_BALANCED: ClassVar["FlowPlacement"]
    __members__: ClassVar[dict[str, "FlowPlacement"]]

    def __init__(self, value: int) -> None: ...
    def __int__(self) -> int: ...
    def __index__(self) -> int: ...
    @property
    def name(self) -> str: ...
    @property
    def value(self) -> int: ...


class Admission:
    """``LOSSLESS`` admits only what fits; ``LOSSY`` forwards hop by hop without admission
    control and drops what does not fit on each edge (EQUAL_BALANCED and WEIGHTED only)."""

    LOSSLESS: ClassVar["Admission"]
    LOSSY: ClassVar["Admission"]
    __members__: ClassVar[dict[str, "Admission"]]

    def __init__(self, value: int) -> None: ...
    def __int__(self) -> int: ...
    def __index__(self) -> int: ...
    @property
    def name(self) -> str: ...
    @property
    def value(self) -> int: ...


class PathSelection:
    """The selection of a search (APPROACH 6).

    ``tie_break`` and ``priority`` are explicit; ``multipath = True`` keeps every equal-cost
    walk (the DAG), ``False`` the one canonical single path; ``multi_edge = True`` keeps every
    eligible parallel member, ``False`` one per ``(u, v)`` pair (minimum cost, then the
    priority under PREFER_HIGHER_PRIORITY, then the lowest EdgeId). Immutable value.
    """

    def __init__(
        self,
        tie_break: EdgeTieBreak,
        priority: EdgePriority,
        multipath: bool = True,
        multi_edge: bool = True,
    ) -> None: ...
    @property
    def tie_break(self) -> EdgeTieBreak: ...
    @property
    def priority(self) -> EdgePriority: ...
    @property
    def multipath(self) -> bool: ...
    @property
    def multi_edge(self) -> bool: ...


class PredDAG:
    """Immutable, graph-bound predecessor DAG (APPROACH 5).

    The entries of node ``v`` are ``parents()[k]``, ``via_edges()[k]`` for ``k`` in
    ``[parent_offsets()[v], parent_offsets()[v + 1])``: the edge ``via`` joins
    ``parent -> v``. Within a node the entries are in strictly ascending via EdgeId (one
    canonical form). Equality is structural. Accessors return read-only NumPy views that
    keep the DAG alive. A DAG of another graph is rejected wherever it is used.
    """

    @staticmethod
    def from_csr(
        graph: StrictMultiDiGraph,
        offsets: "Sequence[int] | np.ndarray",
        parents: "Sequence[int] | np.ndarray",
        via: "Sequence[int] | np.ndarray",
    ) -> "PredDAG":
        """Build a DAG from its CSR, validated before any traversal.

        Raises:
            ValueError: On the first defect, in this order: offsets of length N + 1,
                ``offsets[0] == 0``, monotone, ``offsets[N] == len(parents) == len(via)``;
                every parent and via in range, every via joining ``(parent, child)``, no
                duplicate ``(child, via)``, acyclic; within a node the entries in ascending
                via EdgeId.
        """
        ...

    def graph_uid(self) -> int:
        """uid of the graph this DAG belongs to."""
        ...

    def parent_offsets(self) -> "np.ndarray":
        """int32[N + 1] CSR offsets."""
        ...

    def parents(self) -> "np.ndarray":
        """int32 parent NodeId per entry."""
        ...

    def via_edges(self) -> "np.ndarray":
        """int32 via EdgeId per entry."""
        ...


class Path:
    """Immutable, graph-bound explicit walk (APPROACH 5).

    ``hops`` lists, per hop, one or more parallel members of one ``(u, v)`` pair (the
    caller's member order is kept); ``cost`` is the sum over hops of the minimum member
    cost. ``to_dag()`` gives the equivalent ``PredDAG``.
    """

    def __init__(
        self, graph: StrictMultiDiGraph, hops: "Sequence[Sequence[int]]"
    ) -> None:
        """Validate and build the walk.

        Raises:
            ValueError: On the first defect, in this order: every member id in range,
                no duplicate member, a non-empty hop list, every hop non-empty, the
                members of a hop sharing ``(u, v)``, contiguous hops, no repeated node.
        """
        ...

    def graph_uid(self) -> int: ...
    def hops(self) -> list[list[int]]: ...
    def src(self) -> int: ...
    def dst(self) -> int: ...
    def nodes(self) -> list[int]:
        """The visited nodes, ``src`` first."""
        ...

    def cost(self) -> int: ...
    def to_dag(self) -> PredDAG: ...


class FlowRecord:
    """One flow of a ledger, as a snapshot (APPROACH 8).

    ``flow_id`` is ``(ledger uid, id)``; ``owner`` the id of the policy that placed it;
    ``bundle`` the flow's route; ``route_cost`` its minimum walk cost (``None`` for an empty
    route); ``offered`` exists on Lossy records only (``None`` otherwise); ``edges`` and
    ``drops`` are ``(EdgeId, units)`` pairs in ascending EdgeId (drops cumulative over
    calls). Counters are exact ints and may exceed 64 bits.
    """

    @property
    def flow_id(self) -> tuple[int, int]: ...
    @property
    def owner(self) -> int: ...
    @property
    def sources(self) -> list[int]: ...
    @property
    def destinations(self) -> list[int]: ...
    @property
    def placement(self) -> FlowPlacement: ...
    @property
    def admission(self) -> Admission: ...
    @property
    def bundle(self) -> PredDAG: ...
    @property
    def route_cost(self) -> int | None: ...
    @property
    def offered(self) -> int | None: ...
    @property
    def delivered(self) -> int: ...
    @property
    def edges(self) -> list[tuple[int, int]]: ...
    @property
    def drops(self) -> list[tuple[int, int]]: ...


class FlowGraph:
    """The integer ledger of one graph (APPROACH 8).

    An immutable baseline residual in units and the records of the flows placed on it by
    policies; the ledger has no public placement, removal or reservation. ``allocated()``
    is the sum of every record's edges and ``residual() = baseline - allocated()``, both
    ``int64[E]`` snapshots in EdgeId order. ``reset()`` and ``copy()`` are rejected while a
    policy is bound. A policy destroyed without ``remove()`` leaves its records as
    orphans (readable, never changed again) until ``reset()``.
    """

    def __init__(
        self,
        graph: StrictMultiDiGraph,
        residual: "Sequence[int] | np.ndarray | None" = None,
    ) -> None:
        """Bind a ledger to ``graph``.

        Args:
            graph: The graph (kept alive by the ledger).
            residual: The baseline as integer units, one entry per edge in
                ``[0, capacity_units[e]]`` (a physical baseline goes through ``to_units``);
                ``None``: the capacities.

        Raises:
            ValueError: a baseline of the wrong length, a non-integer entry, or an entry
                outside ``[0, capacity_units[e]]``.
        """
        ...

    def uid(self) -> int:
        """Process-unique ledger uid (the first half of every FlowId)."""
        ...

    def graph(self) -> StrictMultiDiGraph: ...
    def capacity_units(self) -> "np.ndarray":
        """int64[E] the graph's capacities in units (read-only view)."""
        ...

    def residual(self) -> "np.ndarray":
        """int64[E] snapshot: baseline - allocated."""
        ...

    def allocated(self) -> "np.ndarray":
        """int64[E] snapshot: the sum of every record's edges."""
        ...

    def records(self) -> list[FlowRecord]:
        """Live records (bound and orphaned) as snapshots in FlowId order."""
        ...

    def record(self, flow_id: tuple[int, int]) -> FlowRecord:
        """One record; ValueError for a FlowId of another ledger, a retired or an unknown id."""
        ...

    def bottleneck(
        self, sources: "int | Sequence[int]", elig: Eligibility | None = None
    ) -> list[int]:
        """Saturation diagnostic, ascending EdgeId (not a cut).

        The edges (statically eligible, capacity >= 1 unit) leaving the set reachable
        from the allowed sources over eligible forward arcs with residual >= 1 unit.
        """
        ...

    def reset(self, residual: "Sequence[int] | np.ndarray | None" = None) -> None:
        """A new baseline (``None``: the capacities); orphans cleared, ids retired.

        Raises:
            ValueError: while a policy is bound; an invalid baseline.
        """
        ...

    def copy(self) -> "FlowGraph":
        """A new ledger (new uid) with the same baseline and records (ids kept under the
        new uid); ValueError while a policy is bound."""
        ...
