"""Typing and docstrings for netgraph_core public API.

This module is intentionally light; runtime implementations live in the
compiled extension. The goal is to provide type hints and help text.

Every public call is all-or-nothing up to its caller (APPROACH 13), with two limits: (1) the
pybind11 object machinery is outside the guarantee: constructing a bound object from Python,
or reading an attribute or method value converted into a new bound object or container outside
the staged calls, allocates through pybind11's ``make_new_instance``, which does not check the
result of ``tp_alloc``, so a Python allocation failure at that point crashes the process, with
Core state unchanged; (2) a supplied building DagCache is a memo: it keeps the complete entries
and statistics of the lookups that succeeded before a later failure of the call.

Generated ``.pyi`` stubs (pybind11-stubgen, APPROACH 13) arrive with the public search (PLAN,
sub-PR 1d); this hand-written shim is deleted in Phase 5.
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

    Built by ``Graph.eligibility(excluded_nodes, excluded_edges)``. An
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


class Graph:
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
    priority under PREFER_HIGHER_PRIORITY, then the lowest EdgeId). Immutable value; a
    ``tie_break`` or ``priority`` outside its enumeration raises ValueError.
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
    canonical form). Equality is structural. ``parents()`` and ``via_edges()`` are read-only
    NumPy views that keep the DAG alive; ``parent_offsets()`` is a fresh read-only array that
    owns its data, built on demand (the DAG is stored over its member nodes only). A DAG of
    another graph is rejected wherever it is used.
    """

    @staticmethod
    def from_csr(
        graph: Graph,
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
        """int32[N + 1] CSR offsets, a fresh read-only array built on demand."""
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

    def __init__(self, graph: Graph, hops: "Sequence[Sequence[int]]") -> None:
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


class FlowLedger:
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
        graph: Graph,
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

    def graph(self) -> Graph: ...
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

    def copy(self) -> "FlowLedger":
        """A new ledger (new uid) with the same baseline and records (ids kept under the
        new uid); ValueError while a policy is bound."""
        ...


# ---------------------------------------------------------------------------- policies

MAX_FLOWS_PER_CALL: int
"""The call-wide CapacityAware creation cap, 10,000: the flows one ``place`` call may create
(creation and accepted growth, across sources) before it stops with ``LIMIT_FLOWS``. An engine
constant (C++ ``kMaxFlowsPerCall``); not a limit on any mode's existing flow count."""


class Absolute:
    """CostCeiling: an absolute, inclusive path-cost budget (``cost >= 0``)."""

    def __init__(self, cost: int) -> None: ...
    @property
    def cost(self) -> int: ...


class Stretch:
    """CostCeiling: a factor in ``[1, inf)`` of the static shortest cost under the policy's
    Eligibility, read as the exact dyadic rational it is (compared without rounding)."""

    def __init__(self, factor: float) -> None: ...
    @property
    def factor(self) -> float: ...


class DiminishingReturns:
    """The diminishing-returns guard of CapacityAware growth.

    Growth stops once at least ``window`` steps were accepted in the call and the last
    ``window`` gains sum below ``fraction`` times the source's (finite) offer
    (``limit_diminishing``). ``window >= 1``, ``fraction`` in ``(0, 1)``.
    """

    def __init__(self, window: int, fraction: float) -> None: ...
    @property
    def window(self) -> int: ...
    @property
    def fraction(self) -> float: ...


class CostOnly:
    """Routing: IGP forwarding. One flow per roster source on the static shortest-path DAG
    (the DagCache value of the destination set; selection projected at creation on the
    working residual), never refreshed, grown or reoptimised."""

    def __init__(self) -> None: ...


class CapacityAware:
    """Routing: distributed CSPF tunnels (APPROACH 9.1, 9.2).

    ``flow_count`` flows are created by the canonical search of section 6 on the selection
    view (residual >= ``max(1, min_flow_bandwidth)`` units, the ``max_path_cost`` ceiling),
    then grown up to ``max_flow_count`` (``None``: no ceiling); ``reoptimize`` reroutes
    retained flows (adopted only when the candidate delivers the flow's amount);
    ``diminishing`` guards growth (``None``: off). Fields are keyword-only.
    """

    def __init__(
        self,
        *,
        flow_count: int = 1,
        max_flow_count: int | None = None,
        max_path_cost: "Absolute | Stretch | None" = None,
        min_flow_bandwidth: int = 0,
        reoptimize: bool = False,
        diminishing: DiminishingReturns | None = None,
    ) -> None: ...
    @property
    def flow_count(self) -> int: ...
    @property
    def max_flow_count(self) -> int | None: ...
    @property
    def max_path_cost(self) -> "Absolute | Stretch | None": ...
    @property
    def min_flow_bandwidth(self) -> int: ...
    @property
    def reoptimize(self) -> bool: ...
    @property
    def diminishing(self) -> DiminishingReturns | None: ...


class Global:
    """Routing: the omniscient controller. Every call re-solves the minimum-cost flow of the
    cumulative value over the whole eligible graph on ``V = W + own`` (placement MaxFlow
    only); records are the canonical per-source decomposition."""

    def __init__(self) -> None: ...


class Pinned:
    """Routing: explicit route structures, one flow per usable bundle in supply order
    (single endpoints). A ``Path`` converts with ``to_dag()``."""

    def __init__(self, bundles: "Sequence[PredDAG]") -> None: ...
    @property
    def bundles(self) -> list[PredDAG]: ...


class FlowPolicyConfig:
    """The four explicit axes of a policy (keyword-only, no defaults).

    ``Global`` and ``Pinned`` read ``selection.priority`` only for the Weighted kernel and
    ignore the rest; ``CostOnly`` and ``CapacityAware`` read every field. A ``placement`` or
    ``admission`` outside its enumeration raises ValueError; a routing of another type,
    TypeError.
    """

    def __init__(
        self,
        *,
        routing: "CostOnly | CapacityAware | Global | Pinned",
        placement: FlowPlacement,
        admission: Admission,
        selection: PathSelection,
    ) -> None: ...
    @property
    def routing(self) -> "CostOnly | CapacityAware | Global | Pinned": ...
    @property
    def placement(self) -> FlowPlacement: ...
    @property
    def admission(self) -> Admission: ...
    @property
    def selection(self) -> PathSelection: ...


class OutputOptions:
    """Optional outputs of a placement (absent outputs are ``None``): ``edge_flows`` (the
    policy's own allocation per edge after the call), ``residual`` and ``bottleneck`` (the
    ledger's after the call); on Global, ``bottleneck`` adds ``min_cut``."""

    def __init__(
        self,
        *,
        edge_flows: bool = False,
        residual: bool = False,
        bottleneck: bool = False,
    ) -> None: ...
    @property
    def edge_flows(self) -> bool: ...
    @property
    def residual(self) -> bool: ...
    @property
    def bottleneck(self) -> bool: ...


class Termination:
    """A source's local outcome and a call's termination (APPROACH 9.4).

    The reduction over the sources whose offer is unserved: none ``COMPLETE``, else the
    highest of ``LIMIT_FLOWS > LIMIT_DIMINISHING > EXHAUSTED > NO_ROUTE``.
    """

    COMPLETE: ClassVar["Termination"]
    LIMIT_FLOWS: ClassVar["Termination"]
    LIMIT_DIMINISHING: ClassVar["Termination"]
    NO_ROUTE: ClassVar["Termination"]
    EXHAUSTED: ClassVar["Termination"]
    __members__: ClassVar[dict[str, "Termination"]]

    def __init__(self, value: int) -> None: ...
    def __int__(self) -> int: ...
    def __index__(self) -> int: ...
    @property
    def name(self) -> str: ...
    @property
    def value(self) -> int: ...


class PlacementResult:
    """One placement call (the deltas of this call).

    ``placed`` and ``remaining = request - placed`` (``float('inf')`` for ``+inf``) in
    units; ``termination`` and ``local_outcomes`` (one ``(source, Termination)`` per roster
    source, ascending); ``routes_changed`` and ``weights_changed`` (committed events of the
    call, the flag table of 9.4); ``drops`` (this call's transit drops); Global adds
    ``cost_distribution``, the marginal curve of the policy's cumulative flow (its amounts sum
    to ``placed_demand()``; a zero call reports it unchanged, ``[]`` before the first
    placement and after ``remove()``). ``edge_flows``, ``residual``, ``bottleneck`` and
    (Global) ``min_cut`` exist when requested through ``OutputOptions``; ``min_cut`` is ``None``
    when the solve was not maximal on its view, and on a zero call (nothing is solved).
    """

    @property
    def placed(self) -> int: ...
    @property
    def remaining(self) -> "int | float": ...
    @property
    def termination(self) -> Termination: ...
    @property
    def local_outcomes(self) -> list[tuple[int, Termination]]: ...
    @property
    def routes_changed(self) -> bool: ...
    @property
    def weights_changed(self) -> bool: ...
    @property
    def drops(self) -> list[tuple[int, int]]: ...
    @property
    def cost_distribution(self) -> list[tuple[int, int]] | None: ...
    @property
    def edge_flows(self) -> "np.ndarray | None": ...
    @property
    def residual(self) -> "np.ndarray | None": ...
    @property
    def bottleneck(self) -> list[int] | None: ...
    @property
    def min_cut(self) -> list[int] | None: ...


class DagCacheStats:
    """DagCache statistics, atomic per field (not a consistent snapshot): ``hits``,
    ``misses``, ``entries``, ``bytes`` (the entries' retained payload) and
    ``live_values`` (values alive anywhere, entries included)."""

    @property
    def hits(self) -> int: ...
    @property
    def misses(self) -> int: ...
    @property
    def entries(self) -> int: ...
    @property
    def bytes(self) -> int: ...
    @property
    def live_values(self) -> int: ...


class DagCache:
    """A cache of static shortest-path supports bound to one graph (APPROACH 10).

    Keyed ``(roots as a set, direction, Eligibility uid)``; the value is the full static
    support (``multipath``, ``multi_edge``, no residual) as an immutable shared value, valid
    after ``clear()``. Selection is projected at creation, never cached. Building mode: one
    writer, no concurrent readers (a supplied cache needs ``workers == 1``). A supplied
    building cache is a memo (APPROACH 13, limit 2): a public call that fails keeps the
    complete entries and statistics of the lookups that succeeded before; it never holds a
    partial entry, held values stay valid, and a retained entry never changes a result.
    """

    def __init__(self, graph: Graph) -> None: ...
    def graph_uid(self) -> int: ...
    def clear(self) -> None:
        """Empties the cache; the hit and miss counters persist; held values stay valid."""
        ...

    def stats(self) -> DagCacheStats: ...


class FlowPolicy:
    """The one placement engine (APPROACH 9), bound to its ledger at construction.

    ``place(volume, outputs)`` runs the transaction of 9.2 (per source for CostOnly and
    CapacityAware: creation, reoptimisation, amounts, growth on the shared working
    residual; Global: the cumulative minimum-cost flow on ``V = W + own``) and commits one
    atomic exchange: a failed call raises and leaves the ledger and the policy as they were,
    a failed allocation of its Python result included (the result object exists before
    anything is published). Two limits (APPROACH 13): (1) the pybind11 object machinery is
    outside the guarantee: constructing a bound object from Python, or reading an attribute
    or method value converted into a new bound object or container outside the staged calls,
    allocates through pybind11's ``make_new_instance``, which does not check ``tp_alloc``, so a
    Python allocation failure there crashes the process, with Core state unchanged; (2) a
    supplied building DagCache is a memo: it keeps the complete entries and statistics of the
    lookups that succeeded before a later failure (of this call or of the construction);
    never a partial entry, and a retained entry never changes a result. ``remove()`` reverts every
    flow, retires the ids and opens a new epoch. Dropping the policy without ``remove()``
    leaves its records as orphans. Keeps the ledger, the Eligibility and the cache alive.

    Raises:
        ValueError: invalid endpoints (empty, out of range, overlapping sets), configurations
            (Lossy with MaxFlow, Global or CapacityAware; Global without MaxFlow; Weighted
            with Headroom; counts, ceilings, guards; Pinned bundles: empty, duplicate,
            foreign, without a src -> dst walk; Pinned with endpoint sets), foreign objects,
            and requests (outside ``[0, 2**62)`` units, a finite request after ``+inf``).
        TypeError: arguments of the wrong type.
    """

    def __init__(
        self,
        ledger: FlowLedger,
        sources: "int | Sequence[int]",
        destinations: "int | Sequence[int]",
        config: FlowPolicyConfig,
        elig: Eligibility | None = None,
        cache: DagCache | None = None,
    ) -> None: ...
    def place(
        self, volume: "int | float", outputs: OutputOptions | None = None
    ) -> PlacementResult:
        """Place ``volume`` units (an int in ``[0, 2**62)``) or ``float('inf')``."""
        ...

    def remove(self) -> None: ...
    def placed_demand(self) -> int:
        """The sum of delivered over the policy's live records."""
        ...

    def records(self) -> list[FlowRecord]:
        """The policy's live records as snapshots in FlowId order."""
        ...

    def record(self, flow_id: tuple[int, int]) -> FlowRecord:
        """One of the policy's live records; ValueError for a FlowId of another ledger or
        policy, a retired or an unknown id."""
        ...

    def usable_bundles(self) -> int:
        """Pinned: the bundles with a src -> dst walk after static pruning (else 0)."""
        ...

    @property
    def config(self) -> FlowPolicyConfig:
        """The configuration (a copy taken at construction)."""
        ...


class CapacityOptions:
    """The configuration, the static Eligibility (``None``: all allowed), the DagCache
    (``None``: a private cache per worker) and the outputs of ``capacity()`` and
    ``batch_capacity()``."""

    def __init__(
        self,
        config: FlowPolicyConfig,
        *,
        elig: Eligibility | None = None,
        cache: DagCache | None = None,
        outputs: OutputOptions | None = None,
    ) -> None: ...
    @property
    def config(self) -> FlowPolicyConfig: ...
    @property
    def elig(self) -> Eligibility | None: ...
    @property
    def cache(self) -> DagCache | None: ...
    @property
    def outputs(self) -> OutputOptions: ...


class SensitivityOptions:
    """The configuration, Eligibility, DagCache and baseline outputs of ``sensitivity()``."""

    def __init__(
        self,
        config: FlowPolicyConfig,
        *,
        elig: Eligibility | None = None,
        cache: DagCache | None = None,
        outputs: OutputOptions | None = None,
    ) -> None: ...
    @property
    def config(self) -> FlowPolicyConfig: ...
    @property
    def elig(self) -> Eligibility | None: ...
    @property
    def cache(self) -> DagCache | None: ...
    @property
    def outputs(self) -> OutputOptions: ...


class SensitivityStatus:
    """``EVALUATED`` (loss = baseline minus the capacity without the edge) or
    ``UNEVALUABLE`` (the deletion raised a model exception: ``loss`` is ``None`` and
    ``reason`` the message)."""

    EVALUATED: ClassVar["SensitivityStatus"]
    UNEVALUABLE: ClassVar["SensitivityStatus"]
    __members__: ClassVar[dict[str, "SensitivityStatus"]]

    def __init__(self, value: int) -> None: ...
    def __int__(self) -> int: ...
    def __index__(self) -> int: ...
    @property
    def name(self) -> str: ...
    @property
    def value(self) -> int: ...


class SensitivityEntry:
    """One candidate edge: ``edge``, ``loss`` (signed units, never thresholded; ``None``
    when unevaluable), ``status``, ``reason`` (``None`` unless unevaluable)."""

    @property
    def edge(self) -> int: ...
    @property
    def loss(self) -> int | None: ...
    @property
    def status(self) -> SensitivityStatus: ...
    @property
    def reason(self) -> str | None: ...


class SensitivityResult:
    """``baseline`` (the PlacementResult of ``capacity()``) and ``entries`` (ascending
    EdgeId; an absent edge has loss 0)."""

    @property
    def baseline(self) -> PlacementResult: ...
    @property
    def entries(self) -> list[SensitivityEntry]: ...


def capacity(
    graph: Graph,
    sources: "int | Sequence[int]",
    destinations: "int | Sequence[int]",
    opts: CapacityOptions,
) -> PlacementResult:
    """``place(+inf)`` of a FlowPolicy with ``opts.config`` on a fresh ledger of the graph.

    The exact maximum flow is the Global + MaxFlow configuration (the minimum-cost maximum
    flow; ``cost_distribution`` its marginal curve). A failed call raises and returns no
    result; a supplied building cache keeps the entries and statistics of the lookups that
    succeeded before the failure (a memo: APPROACH 13, limit 2).
    """
    ...


def batch_capacity(
    graph: Graph,
    pairs: "Sequence[tuple[int | Sequence[int], int | Sequence[int]]]",
    opts: CapacityOptions,
    per_pair_elig: "Sequence[Eligibility | None] | None" = None,
    workers: int | None = None,
) -> list[PlacementResult]:
    """``capacity()`` per ``(sources, destinations)`` pair, in pair order.

    ``opts.elig`` is intersected with ``per_pair_elig[i]`` (``None`` entries allowed; a given
    sequence has one entry per pair, also when it is empty). Everything is validated before a
    worker starts: the arguments first (workers, options, the configuration: Pinned is
    rejected; the ``per_pair_elig`` length), then each pair completely in pair order (its
    endpoints, their overlap, its Eligibility), so the error of the lowest failing pair index
    is raised; an empty batch returns ``[]``. ``workers``: an explicit int >= 1 is honoured; ``None`` takes
    ``NGRAPH_CORE_BATCH_THREADS`` (or 1), and 1 when a cache is supplied (a supplied
    building cache is rejected unless ``workers == 1``). A failed call raises and returns
    no result; a supplied building cache keeps the entries and statistics of the lookups
    that succeeded before the failure (a memo: APPROACH 13, limit 2).
    """
    ...


def sensitivity(
    graph: Graph,
    sources: "int | Sequence[int]",
    destinations: "int | Sequence[int]",
    opts: SensitivityOptions,
    workers: int | None = None,
) -> SensitivityResult:
    """The baseline ``capacity()`` and one deletion per candidate edge.

    Candidates: under Global the baseline's positive-flow edges; under every other model
    every edge ``opts.elig`` allows with at least 1 unit of static capacity. A loss is the
    baseline's placed minus ``capacity()`` with the edge deleted (``opts.elig`` minus the
    edge, the same configuration and cache), in signed units, never thresholded. Global
    evaluates it by the flow core's deletion re-solve from the baseline solve's final state
    (that value exactly; a failure there raises); CostOnly gives an edge outside every roster
    source's reachable part of the static support loss 0 without a placement, a
    result-equivalent predicate; every other deletion is placed on a fresh ledger, where a
    model exception makes the entry ``UNEVALUABLE`` with the message as reason. ``workers`` as ``batch_capacity`` (``NGRAPH_CORE_SENSITIVITY_THREADS``).
    A failed call raises and returns no result; a supplied building cache keeps the
    entries and statistics of the lookups that succeeded before the failure (a memo:
    APPROACH 13, limit 2).
    """
    ...
