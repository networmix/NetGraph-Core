# Generated from the compiled module; do not edit.
"""NetGraph-Core: graph, flow and placement models on exact integer units.

Every public name is re-exported from the compiled ``_netgraph_core`` module, so users
import everything from here. Type checkers read the generated stub ``__init__.pyi`` next
to this file (``py.typed`` marks the package as typed).

- The graph on integer units: ``Graph``, ``Eligibility``, ``to_units``, ``as_float``.
- Route structures: ``PredDAG``, ``Path``.
- Route selection and placement settings: ``PathSelection``, ``EdgeTieBreak``,
  ``EdgePriority``, ``FlowPlacement``, ``Admission``, ``Origin``, ``WeightBasis``.
- The ledger of placed traffic: ``FlowLedger`` and its ``FlowRecord`` snapshots.
- Forwarding programs: ``ForwardingProgram`` (with its ``Tunnel`` entries and
  ``ForwardingProgram.from_routes`` over ``Route`` values), ``TunnelStatus``.
- The policy engine: ``FlowPolicy``, ``FlowPolicyConfig``, the routings ``CostOnly``,
  ``CapacityAware`` (T1, with ``Tunnels``), ``Global`` (with ``tunnels``: T2), ``Pinned``, and
  ``PlacementResult`` (with its ``SourceFlow`` and ``TunnelFlow`` rows and its ``TunnelPlan``),
  ``OutputOptions``, ``Termination``.
- The shortest-path cache: ``DagCache`` (building and frozen modes), ``DagCacheMiss``.
- The evaluators: ``capacity``, ``batch_capacity``, ``sensitivity``.
- The searches: ``spf`` (``Direction``, ``SpfOptions``, ``SpfResult``, ``UNREACHABLE``),
  ``ksp`` (``KspOptions``), ``resolve_to_paths``.

``_netgraph_core`` also holds a private entry of the flow core, never re-exported here."""

from __future__ import annotations

import collections.abc
import fractions
import typing

import numpy
import numpy.typing

__all__: list[str] = [
    "__version__",
    "Graph",
    "Eligibility",
    "to_units",
    "as_float",
    "PredDAG",
    "Path",
    "EdgeTieBreak",
    "EdgePriority",
    "PathSelection",
    "FlowPlacement",
    "Admission",
    "Origin",
    "WeightBasis",
    "FlowLedger",
    "FlowRecord",
    "TunnelStatus",
    "Tunnel",
    "Route",
    "ForwardingProgram",
    "Absolute",
    "Stretch",
    "Tunnels",
    "CostOnly",
    "CapacityAware",
    "Global",
    "Pinned",
    "FlowPolicyConfig",
    "OutputOptions",
    "Termination",
    "PlacementResult",
    "SourceFlow",
    "TunnelFlow",
    "TunnelPlan",
    "FlowPolicy",
    "DagCache",
    "DagCacheMiss",
    "DagCacheStats",
    "QueryOptions",
    "SensitivityEntry",
    "SensitivityResult",
    "capacity",
    "batch_capacity",
    "sensitivity",
    "Direction",
    "SpfOptions",
    "SpfResult",
    "UNREACHABLE",
    "spf",
    "KspOptions",
    "ksp",
    "resolve_to_paths",
]
__version__: str

class FlowPlacement:
    """
    How traffic spreads over the routes (CostOnly's kernel; the interior split of supplied PredDAG
    routes, EQUAL_BALANCED or WEIGHTED).
    MAX_FLOW: the minimum-cost flow core over the routes (a controller restricted to them), re-solved
    over the free capacity plus the policy's own allocation (it may move its own earlier placement,
    never another's). WEIGHTED (WCMP, weighted-cost multipath): every node splits in proportion to
    weights fixed when the routes are prepared (highest averages: whole-unit shares as close to
    proportional as possible; the weights come from WeightBasis). EQUAL_BALANCED (ECMP, equal-cost
    multipath): every node splits equally, by round robin.

    Members:

      MAX_FLOW

      WEIGHTED

      EQUAL_BALANCED
    """

    EQUAL_BALANCED: typing.ClassVar[FlowPlacement]
    MAX_FLOW: typing.ClassVar[FlowPlacement]
    WEIGHTED: typing.ClassVar[FlowPlacement]
    __members__: typing.ClassVar[dict[str, FlowPlacement]]

    @typing.overload
    def __eq__(self, other: FlowPlacement, /) -> bool: ...
    @typing.overload
    def __eq__(self, other: typing.Any, /) -> bool: ...
    def __getstate__(self) -> int: ...
    def __hash__(self) -> int: ...
    def __index__(self) -> int: ...
    def __init__(self, value: typing.SupportsInt | typing.SupportsIndex) -> None: ...
    def __int__(self) -> int: ...
    @typing.overload
    def __ne__(self, other: FlowPlacement, /) -> bool: ...
    @typing.overload
    def __ne__(self, other: typing.Any, /) -> bool: ...
    def __repr__(self) -> str: ...
    def __setstate__(
        self, state: typing.SupportsInt | typing.SupportsIndex, /
    ) -> None: ...
    def __str__(self) -> str: ...
    @property
    def name(self) -> str: ...
    @property
    def value(self) -> int: ...

class Admission:
    """
    LOSSLESS admits the largest total the routes carry with no loss at all; LOSSY sends the whole
    finite offer hop by hop and records where it is lost (EQUAL_BALANCED and WEIGHTED only).

    Members:

      LOSSLESS

      LOSSY
    """

    LOSSLESS: typing.ClassVar[Admission]
    LOSSY: typing.ClassVar[Admission]
    __members__: typing.ClassVar[dict[str, Admission]]

    @typing.overload
    def __eq__(self, other: Admission, /) -> bool: ...
    @typing.overload
    def __eq__(self, other: typing.Any, /) -> bool: ...
    def __getstate__(self) -> int: ...
    def __hash__(self) -> int: ...
    def __index__(self) -> int: ...
    def __init__(self, value: typing.SupportsInt | typing.SupportsIndex) -> None: ...
    def __int__(self) -> int: ...
    @typing.overload
    def __ne__(self, other: Admission, /) -> bool: ...
    @typing.overload
    def __ne__(self, other: typing.Any, /) -> bool: ...
    def __repr__(self) -> str: ...
    def __setstate__(
        self, state: typing.SupportsInt | typing.SupportsIndex, /
    ) -> None: ...
    def __str__(self) -> str: ...
    @property
    def name(self) -> str: ...
    @property
    def value(self) -> int: ...

class Origin:
    """
    How a demand's traffic enters at its sources, separately from how it is routed. EQUAL: round
    robin over the roster (ascending NodeId, the first sources one unit more); a source whose share
    cannot go anywhere keeps it and loses it (LOSSLESS then admits only totals that give it
    nothing). OPTIMIZED: the controller chooses the injections (one virtual source).

    Members:

      EQUAL

      OPTIMIZED
    """

    EQUAL: typing.ClassVar[Origin]
    OPTIMIZED: typing.ClassVar[Origin]
    __members__: typing.ClassVar[dict[str, Origin]]

    @typing.overload
    def __eq__(self, other: Origin, /) -> bool: ...
    @typing.overload
    def __eq__(self, other: typing.Any, /) -> bool: ...
    def __getstate__(self) -> int: ...
    def __hash__(self) -> int: ...
    def __index__(self) -> int: ...
    def __init__(self, value: typing.SupportsInt | typing.SupportsIndex) -> None: ...
    def __int__(self) -> int: ...
    @typing.overload
    def __ne__(self, other: Origin, /) -> bool: ...
    @typing.overload
    def __ne__(self, other: typing.Any, /) -> bool: ...
    def __repr__(self) -> str: ...
    def __setstate__(
        self, state: typing.SupportsInt | typing.SupportsIndex, /
    ) -> None: ...
    def __str__(self) -> str: ...
    @property
    def name(self) -> str: ...
    @property
    def value(self) -> int: ...

class WeightBasis:
    """
    Where WEIGHTED member weights come from (CostOnly WEIGHTED and supplied PredDAG routes).
    CAPACITY (the default): each member's static capacity. EXPLICIT: edge_weights, one integer per
    edge (a positive weight on an edge of capacity 0 becomes 0). RESIDUAL_SNAPSHOT: each member's free
    capacity read once when the routes are prepared, against one common view, and frozen for the
    epoch. A member of weight 0 is pruned with the branches it leaves without a destination, inside
    the route only. The selection's priority never chooses weights.

    Members:

      CAPACITY

      EXPLICIT

      RESIDUAL_SNAPSHOT
    """

    CAPACITY: typing.ClassVar[WeightBasis]
    EXPLICIT: typing.ClassVar[WeightBasis]
    RESIDUAL_SNAPSHOT: typing.ClassVar[WeightBasis]
    __members__: typing.ClassVar[dict[str, WeightBasis]]

    @typing.overload
    def __eq__(self, other: WeightBasis, /) -> bool: ...
    @typing.overload
    def __eq__(self, other: typing.Any, /) -> bool: ...
    def __getstate__(self) -> int: ...
    def __hash__(self) -> int: ...
    def __index__(self) -> int: ...
    def __init__(self, value: typing.SupportsInt | typing.SupportsIndex) -> None: ...
    def __int__(self) -> int: ...
    @typing.overload
    def __ne__(self, other: WeightBasis, /) -> bool: ...
    @typing.overload
    def __ne__(self, other: typing.Any, /) -> bool: ...
    def __repr__(self) -> str: ...
    def __setstate__(
        self, state: typing.SupportsInt | typing.SupportsIndex, /
    ) -> None: ...
    def __str__(self) -> str: ...
    @property
    def name(self) -> str: ...
    @property
    def value(self) -> int: ...

class EdgeTieBreak:
    """
    How a search breaks cost ties. Between parallel members of one (u, v) pair
    (multi_edge = False): DETERMINISTIC keeps the lowest EdgeId, PREFER_HIGHER_PRIORITY the
    higher priority value, then the lowest EdgeId. Between equal-cost walks (multipath =
    False): PREFER_HIGHER_PRIORITY first takes the highest bottleneck (the smallest priority
    value along the walk), DETERMINISTIC skips that step; then the fewest hops, then the
    lexicographically smallest EdgeId sequence.

    Members:

      DETERMINISTIC

      PREFER_HIGHER_PRIORITY
    """

    DETERMINISTIC: typing.ClassVar[EdgeTieBreak]
    PREFER_HIGHER_PRIORITY: typing.ClassVar[EdgeTieBreak]
    __members__: typing.ClassVar[dict[str, EdgeTieBreak]]

    @typing.overload
    def __eq__(self, other: EdgeTieBreak, /) -> bool: ...
    @typing.overload
    def __eq__(self, other: typing.Any, /) -> bool: ...
    def __getstate__(self) -> int: ...
    def __hash__(self) -> int: ...
    def __index__(self) -> int: ...
    def __init__(self, value: typing.SupportsInt | typing.SupportsIndex) -> None: ...
    def __int__(self) -> int: ...
    @typing.overload
    def __ne__(self, other: EdgeTieBreak, /) -> bool: ...
    @typing.overload
    def __ne__(self, other: typing.Any, /) -> bool: ...
    def __repr__(self) -> str: ...
    def __setstate__(
        self, state: typing.SupportsInt | typing.SupportsIndex, /
    ) -> None: ...
    def __str__(self) -> str: ...
    @property
    def name(self) -> str: ...
    @property
    def value(self) -> int: ...

class EdgePriority:
    """
    The per-edge value PREFER_HIGHER_PRIORITY reads (single-path and parallel-member choices):
    RESIDUAL (residual units), CAPACITY (static capacity units), HEADROOM (residual / capacity,
    compared exactly as a fraction, 0 on an edge of capacity 0). A search without a residual reads
    the static capacity as the residual: RESIDUAL and CAPACITY read the capacity, HEADROOM 1 (0 on
    an edge of capacity 0). T1's CSPF reads its unreserved view U instead (RESIDUAL U, HEADROOM
    U / capacity). A search and T1 comparison only: it never enables a controller objective
    (T2's is Global(prefer_headroom)) and never chooses weights (WeightBasis does).

    Members:

      RESIDUAL

      CAPACITY

      HEADROOM
    """

    CAPACITY: typing.ClassVar[EdgePriority]
    HEADROOM: typing.ClassVar[EdgePriority]
    RESIDUAL: typing.ClassVar[EdgePriority]
    __members__: typing.ClassVar[dict[str, EdgePriority]]

    @typing.overload
    def __eq__(self, other: EdgePriority, /) -> bool: ...
    @typing.overload
    def __eq__(self, other: typing.Any, /) -> bool: ...
    def __getstate__(self) -> int: ...
    def __hash__(self) -> int: ...
    def __index__(self) -> int: ...
    def __init__(self, value: typing.SupportsInt | typing.SupportsIndex) -> None: ...
    def __int__(self) -> int: ...
    @typing.overload
    def __ne__(self, other: EdgePriority, /) -> bool: ...
    @typing.overload
    def __ne__(self, other: typing.Any, /) -> bool: ...
    def __repr__(self) -> str: ...
    def __setstate__(
        self, state: typing.SupportsInt | typing.SupportsIndex, /
    ) -> None: ...
    def __str__(self) -> str: ...
    @property
    def name(self) -> str: ...
    @property
    def value(self) -> int: ...

class Direction:
    """
    The orientation of a search: FORWARD from the roots along the edges, or REVERSE toward the roots
    against the edges (the DAG it returns keeps the edges' own direction, so its walks lead to the
    roots).

    Members:

      FORWARD

      REVERSE
    """

    FORWARD: typing.ClassVar[Direction]
    REVERSE: typing.ClassVar[Direction]
    __members__: typing.ClassVar[dict[str, Direction]]

    @typing.overload
    def __eq__(self, other: Direction, /) -> bool: ...
    @typing.overload
    def __eq__(self, other: typing.Any, /) -> bool: ...
    def __getstate__(self) -> int: ...
    def __hash__(self) -> int: ...
    def __index__(self) -> int: ...
    def __init__(self, value: typing.SupportsInt | typing.SupportsIndex) -> None: ...
    def __int__(self) -> int: ...
    @typing.overload
    def __ne__(self, other: Direction, /) -> bool: ...
    @typing.overload
    def __ne__(self, other: typing.Any, /) -> bool: ...
    def __repr__(self) -> str: ...
    def __setstate__(
        self, state: typing.SupportsInt | typing.SupportsIndex, /
    ) -> None: ...
    def __str__(self) -> str: ...
    @property
    def name(self) -> str: ...
    @property
    def value(self) -> int: ...

class PathSelection:
    """
    The selection of a search or a policy: tie_break (DETERMINISTIC by default), priority
    (CAPACITY by default); multipath None (the default) is the model's natural shape (every
    equal-cost walk for searches and CostOnly, one path for a tunnel model), True every equal-cost
    walk (the DAG; tunnel models reject it), False the one canonical single path; multi_edge = True
    keeps every eligible parallel member, False one per (u, v) pair. A tie_break or priority
    outside its enumeration raises ValueError.
    """

    def __eq__(self, arg0: typing.Any) -> bool: ...
    def __hash__(self) -> int: ...
    def __init__(
        self,
        tie_break: EdgeTieBreak = ...,
        priority: EdgePriority = ...,
        multipath: bool | None = None,
        multi_edge: bool = True,
    ) -> None: ...
    def __repr__(self) -> str: ...
    @property
    def multi_edge(self) -> bool:
        """
        True: every eligible parallel member; False: one member per (u, v) pair
        """

    @property
    def multipath(self) -> bool | None:
        """
        None: the model's natural shape; True: every equal-cost walk (a DAG); False: the canonical
        single path
        """

    @property
    def priority(self) -> EdgePriority:
        """
        the EdgePriority
        """

    @property
    def tie_break(self) -> EdgeTieBreak:
        """
        the EdgeTieBreak
        """

class Eligibility:
    """
    Which nodes and edges of one graph may be used: one flag per node and per edge, fixed when built
    (it does not depend on traffic). Built by Graph.eligibility(excluded_nodes, excluded_edges) and
    Eligibility.intersect.
    """

    def __init__(self, _no_public_constructor: typing.Never, /) -> None: ...
    @staticmethod
    def intersect(a: Eligibility, b: Eligibility) -> Eligibility:
        """
        The intersection of two values of one graph (ValueError for values of different graphs):
        an operand itself when the result equals it (a first, then b), otherwise a new value with a
        fresh uid; equal intersections built separately get distinct uids.
        """

    def edge_allowed(self) -> numpy.typing.NDArray[numpy.bool_]:
        """
        bool[E] read-only view (not excluded and both endpoints allowed)
        """

    def graph_uid(self) -> int:
        """
        uid of the graph this value belongs to
        """

    def node_allowed(self) -> numpy.typing.NDArray[numpy.bool_]:
        """
        bool[N] read-only view
        """

    def uid(self) -> int:
        """
        0 for the all-allowed value, a fresh uid otherwise
        """

class Graph:
    """
    Immutable directed multigraph on integer units (parallel edges and self-loops allowed). Edges are
    given in input order and numbered by sorting on (src, cost, dst), ties kept in input order: that
    number is the public EdgeId. input_index()[e] is the input position of edge e. Per-edge values
    use EdgeId order; in_adj_edge_index() instead lists EdgeIds in incoming-slot order. Costs are
    positive integers; capacities are int64 units of the explicit resolution q (a power of two in
    [2^-32, 2^32]); physical capacities enter through to_units.
    """

    def __init__(
        self,
        num_nodes: typing.SupportsIndex,
        src: collections.abc.Sequence[typing.SupportsIndex]
        | numpy.typing.NDArray[numpy.integer],
        dst: collections.abc.Sequence[typing.SupportsIndex]
        | numpy.typing.NDArray[numpy.integer],
        capacity_units: collections.abc.Sequence[typing.SupportsIndex]
        | numpy.typing.NDArray[numpy.integer],
        cost: collections.abc.Sequence[typing.SupportsIndex]
        | numpy.typing.NDArray[numpy.integer],
        resolution: typing.SupportsFloat | typing.SupportsIndex,
    ) -> None:
        """
        Validates, in order, and reports the first defect as ValueError: the resolution (a power of
        two in [2^-32, 2^32]), the unit values (every capacity in [0, 2^62)), the array shapes,
        zero or negative costs, the remaining bounds (endpoints in range, then the
        capacity and cost sums below 2^62). Integers are read without narrowing, so the error
        does not depend on how the input is represented.
        """

    def capacity_units(self) -> numpy.typing.NDArray[numpy.int64]:
        """
        int64[E] capacities in units, id order (read-only view)
        """

    def cost(self) -> numpy.typing.NDArray[numpy.int64]:
        """
        int64[E] costs, id order (read-only view)
        """

    def cost_gcd(self) -> int:
        """
        gcd of the edge costs (0 without edges)
        """

    def edge_dst(self) -> numpy.typing.NDArray[numpy.int32]:
        """
        int32[E] destination node per edge (read-only view)
        """

    def edge_src(self) -> numpy.typing.NDArray[numpy.int32]:
        """
        int32[E] source node per edge (read-only view)
        """

    def eligibility(
        self,
        excluded_nodes: collections.abc.Sequence[typing.SupportsIndex]
        | numpy.typing.NDArray[numpy.integer] = (),
        excluded_edges: collections.abc.Sequence[typing.SupportsIndex]
        | numpy.typing.NDArray[numpy.integer] = (),
    ) -> Eligibility:
        """
        Static eligibility from excluded node ids and EdgeIds (duplicates and unsorted lists
        accepted, out-of-range ids rejected); an edge is allowed when it is not excluded and
        both endpoints are allowed. Empty lists mean all allowed (uid 0).
        """

    def in_adj_edge_index(self) -> numpy.typing.NDArray[numpy.int32]:
        """
        int32[E]: the in-edges of v in (cost, src, input) order, as EdgeIds
        """

    def in_row_offsets(self) -> numpy.typing.NDArray[numpy.int32]:
        """
        int32[N + 1]: the in-edges of v are in_adj_edge_index()[in_row_offsets[v] : in_row_offsets[v + 1]]
        """

    def input_index(self) -> numpy.typing.NDArray[numpy.int32]:
        """
        int32[E]: EdgeId -> input position (read-only view); build the inverse yourself if needed
        """

    def max_cost(self) -> int:
        """
        largest edge cost (0 without edges)
        """

    def num_edges(self) -> int:
        """
        E: EdgeIds are 0 .. E - 1
        """

    def num_nodes(self) -> int:
        """
        N: node ids are 0 .. N - 1
        """

    def resolution(self) -> float:
        """
        the resolution q (units of capacity)
        """

    def row_offsets(self) -> numpy.typing.NDArray[numpy.int32]:
        """
        int32[N + 1]: the out-edges of u are the ids [row_offsets[u], row_offsets[u + 1])
        """

    def uid(self) -> int:
        """
        process-unique id (never reused)
        """

class PredDAG:
    """
    Immutable route structure of one graph, a predecessor DAG (directed, no cycle): for each node v,
    the edges that reach v on the route and the node each comes from (its parent). The entries of v
    are parents()[k], via_edges()[k] for k in [parent_offsets()[v], parent_offsets()[v + 1]), in
    strictly ascending EdgeId, so two DAGs of the same graph with the same entries have equal arrays
    and compare equal; DAGs of different graphs never compare equal.
    parents() and via_edges() are read-only views that keep the DAG alive; parent_offsets() is a new
    read-only array on each call (the DAG stores only the nodes that have entries). Built by
    PredDAG.from_csr, Path.to_dag(), the searches and the policies (FlowRecord.bundle).
    """

    def __init__(self, _no_public_constructor: typing.Never, /) -> None: ...
    @staticmethod
    def from_csr(
        graph: Graph,
        offsets: collections.abc.Sequence[typing.SupportsIndex]
        | numpy.typing.NDArray[numpy.integer],
        parents: collections.abc.Sequence[typing.SupportsIndex]
        | numpy.typing.NDArray[numpy.integer],
        via: collections.abc.Sequence[typing.SupportsIndex]
        | numpy.typing.NDArray[numpy.integer],
    ) -> PredDAG:
        """
        Validates, in order, before any traversal (ValueError on the first defect): offsets of
        length N + 1, offsets[0] == 0, monotone, offsets[N] == len(parents) == len(via); then every
        parent and via in range, every via joining (parent, child), no duplicate (child, via),
        acyclic; within a node the entries in ascending via EdgeId.
        """

    def __eq__(self, arg0: typing.Any) -> bool: ...
    def __hash__(self) -> int: ...
    def __repr__(self) -> str: ...
    def graph_uid(self) -> int:
        """
        uid of the graph this DAG belongs to
        """

    def parent_offsets(self) -> numpy.typing.NDArray[numpy.int32]:
        """
        int32[N + 1]: the entries of node v are [parent_offsets[v], parent_offsets[v + 1]); a new read-only array per call
        """

    def parents(self) -> numpy.typing.NDArray[numpy.int32]:
        """
        int32 parent NodeId per entry (read-only view)
        """

    def via_edges(self) -> numpy.typing.NDArray[numpy.int32]:
        """
        int32 via EdgeId per entry (read-only view)
        """

class Path:
    """
    Immutable route of one graph, hop by hop: each hop lists one or more parallel edges (members)
    from the same u to the same v, in the caller's order. cost = the sum over hops of the cost of
    each hop's cheapest member.
    """

    def __eq__(self, arg0: typing.Any) -> bool: ...
    def __hash__(self) -> int: ...
    def __init__(
        self,
        graph: Graph,
        hops: collections.abc.Sequence[
            collections.abc.Sequence[typing.SupportsIndex]
            | numpy.typing.NDArray[numpy.integer]
        ]
        | numpy.typing.NDArray[numpy.integer],
    ) -> None:
        """
        Validates, in order (ValueError on the first defect): every member id in range, no duplicate
        member, a non-empty hop list, every hop non-empty, the members of a hop sharing (u, v),
        contiguous hops, no repeated node.
        """

    def __repr__(self) -> str: ...
    def cost(self) -> int:
        """
        sum over hops of the minimum member cost
        """

    def dst(self) -> int:
        """
        the last node of the walk
        """

    def graph_uid(self) -> int:
        """
        uid of the graph this path belongs to
        """

    def hops(self) -> list[list[int]]:
        """
        the hops as lists of EdgeIds (the caller's member order)
        """

    def src(self) -> int:
        """
        the first node of the walk
        """

    def to_dag(self) -> PredDAG:
        """
        the equivalent PredDAG (every member of every hop)
        """

class FlowRecord:
    """
    One flow of a ledger, as a snapshot: flow_id (ledger uid, id), owner (the policy's id),
    sources, destinations, admission, entry (a T1, T2 or Pinned policy's program entry, else None), bundle
    (PredDAG: the flow's route),
    route_cost (its minimum walk cost; None when empty), offered (Lossy only, else None),
    delivered, edges [(EdgeId, units)] (carried load) and drops [(EdgeId, units)] (cumulative,
    ascending EdgeId), reserved_edges [(EdgeId, units)] (T1's strict holds, ascending EdgeId). offered
    and the drops units are exact Python ints (they may exceed 64 bits).
    """

    def __init__(self, _no_public_constructor: typing.Never, /) -> None: ...
    def __repr__(self) -> str: ...
    @property
    def admission(self) -> Admission:
        """
        LOSSLESS or LOSSY
        """

    @property
    def bundle(self) -> PredDAG:
        """
        the flow's route (an immutable PredDAG shared with the ledger)
        """

    @property
    def delivered(self) -> int:
        """
        the units the flow delivers to its destinations (its total)
        """

    @property
    def destinations(self) -> list[int]:
        """
        the flow's destination node ids, ascending
        """

    @property
    def drops(self) -> list[tuple[int, int]]:
        """
        Lossy: the units dropped per edge [(EdgeId, units)], summed over calls, positive entries in
        ascending EdgeId (exact ints); empty for a Lossless flow
        """

    @property
    def edges(self) -> list[tuple[int, int]]:
        """
        the flow's allocation [(EdgeId, units)], positive entries in ascending EdgeId
        """

    @property
    def entry(self) -> int | None:
        """
        the program entry the record carries (T1, T2 and Pinned: the entry identifier); None for the
        records kept per roster source (CostOnly, Global)
        """

    @property
    def flow_id(self) -> tuple[int, int]:
        """
        (ledger uid, id): ids are assigned in creation order and never reused in the ledger
        """

    @property
    def offered(self) -> int | None:
        """
        Lossy: the units offered to the flow, summed over calls (an exact int, may exceed 64 bits);
        None for a Lossless flow
        """

    @property
    def owner(self) -> int:
        """
        the ledger-local id of the policy that created the flow
        """

    @property
    def reserved_edges(self) -> list[tuple[int, int]]:
        """
        T1: the strict holds of an up tunnel [(EdgeId, units)], its request on every link of its path,
        ascending EdgeId; empty for every other record (holds never police traffic)
        """

    @property
    def route_cost(self) -> int | None:
        """
        the cost of the cheapest walk of the bundle from the sources to the destinations; None when
        the bundle has no such walk
        """

    @property
    def sources(self) -> list[int]:
        """
        the flow's source node ids, ascending
        """

class FlowLedger:
    """
    The integer ledger of one graph, in two projections that are never added: carried traffic (the
    service baseline; allocated() is the sum of every record's edges and residual() = baseline -
    allocated()) and the strict holds of T1 tunnels (the reservation baseline; reserved() is the sum of
    every record's reserved_edges and unreserved() = reservation baseline - reserved()), all int64[E]
    snapshots. Traffic may use bandwidth that is held: holds never police. The baselines are replaced
    only by reset(); the records are the flows placed on it by policies (the ledger has no public
    placement, removal or reservation). reset() and copy() are rejected while a policy is bound (a
    FlowPolicy built on this ledger still exists).
    """

    def __copy__(self) -> FlowLedger: ...
    def __deepcopy__(self, memo: typing.Any) -> FlowLedger: ...
    def __init__(
        self,
        graph: Graph,
        residual: collections.abc.Sequence[typing.SupportsIndex]
        | numpy.typing.NDArray[numpy.integer]
        | None = None,
        *,
        unreserved: collections.abc.Sequence[typing.SupportsIndex]
        | numpy.typing.NDArray[numpy.integer]
        | None = None,
    ) -> None:
        """
        residual and unreserved: the service and reservation baselines as integer units, one entry per
        edge in [0, capacity_units[e]] (a physical baseline goes through to_units first; residual checked
        first); None: the capacities. ValueError otherwise.
        """

    def allocated(self) -> numpy.typing.NDArray[numpy.int64]:
        """
        int64[E] snapshot: the sum of every record's edges
        """

    def bottleneck(
        self,
        sources: typing.SupportsIndex
        | collections.abc.Sequence[typing.SupportsIndex]
        | numpy.typing.NDArray[numpy.integer],
        elig: Eligibility | None = None,
    ) -> list[int]:
        """
        Saturation report: the eligible edges with capacity >= 1 unit that leave the set of nodes the
        allowed sources can reach (following eligible edges forward with residual >= 1 unit), in
        ascending EdgeId. It ignores the destinations, so it is not a source-destination cut.
        """

    def capacity_units(self) -> numpy.typing.NDArray[numpy.int64]:
        """
        int64[E] the graph's capacities in units (read-only view)
        """

    def copy(self) -> FlowLedger:
        """
        a new ledger (new uid) with the same baselines, records and holds, ids kept under the new uid;
        ValueError while a policy is bound
        """

    def graph(self) -> Graph:
        """
        the ledger's graph
        """

    def record(
        self, flow_id: tuple[typing.SupportsIndex, typing.SupportsIndex]
    ) -> FlowRecord:
        """
        one record by FlowId (ValueError for a FlowId of another ledger, a retired or an unknown id)
        """

    def records(self) -> list[FlowRecord]:
        """
        live records as snapshots in FlowId order: those of bound policies and orphans (records left by a
        policy destroyed without remove())
        """

    def reserved(self) -> numpy.typing.NDArray[numpy.int64]:
        """
        int64[E] snapshot: strict holds per edge (T1), the sum of every record's reserved_edges
        """

    def reset(
        self,
        residual: collections.abc.Sequence[typing.SupportsIndex]
        | numpy.typing.NDArray[numpy.integer]
        | None = None,
        *,
        unreserved: collections.abc.Sequence[typing.SupportsIndex]
        | numpy.typing.NDArray[numpy.integer]
        | None = None,
    ) -> None:
        """
        New baselines (None: the capacities); orphaned records with their allocations and holds are
        cleared and their ids retired (ids are never reused). ValueError while a policy is bound.
        """

    def residual(self) -> numpy.typing.NDArray[numpy.int64]:
        """
        int64[E] snapshot: baseline - allocated
        """

    def uid(self) -> int:
        """
        process-unique ledger uid (the first half of every FlowId)
        """

    def unreserved(self) -> numpy.typing.NDArray[numpy.int64]:
        """
        int64[E] snapshot: reservation baseline minus strict holds, per edge
        """

class TunnelStatus:
    """
    A tunnel's state. UP; NO_BANDWIDTH and COST_CEILING: T1's setup left it down (no path passes
    the bandwidth filter, or the cheapest exceeds the cost ceiling): no path, weight 0, and it
    stays down in every forwarding of its program; and in one forwarding only (a by_tunnel row)
    FAILED: its path touches a node or an edge the forwarding policy's Eligibility excludes, so it
    gets nothing and its head splits over its other entries (replay), and BLACKHOLED: up, but its
    path reaches a black hole of stale forwarding, where its traffic is dropped.

    Members:

      UP

      NO_BANDWIDTH

      COST_CEILING

      FAILED

      BLACKHOLED
    """

    BLACKHOLED: typing.ClassVar[TunnelStatus]
    COST_CEILING: typing.ClassVar[TunnelStatus]
    FAILED: typing.ClassVar[TunnelStatus]
    NO_BANDWIDTH: typing.ClassVar[TunnelStatus]
    UP: typing.ClassVar[TunnelStatus]
    __members__: typing.ClassVar[dict[str, TunnelStatus]]

    @typing.overload
    def __eq__(self, other: TunnelStatus, /) -> bool: ...
    @typing.overload
    def __eq__(self, other: typing.Any, /) -> bool: ...
    def __getstate__(self) -> int: ...
    def __hash__(self) -> int: ...
    def __index__(self) -> int: ...
    def __init__(self, value: typing.SupportsInt | typing.SupportsIndex) -> None: ...
    def __int__(self) -> int: ...
    @typing.overload
    def __ne__(self, other: TunnelStatus, /) -> bool: ...
    @typing.overload
    def __ne__(self, other: typing.Any, /) -> bool: ...
    def __repr__(self) -> str: ...
    def __setstate__(
        self, state: typing.SupportsInt | typing.SupportsIndex, /
    ) -> None: ...
    def __str__(self) -> str: ...
    @property
    def name(self) -> str: ...
    @property
    def value(self) -> int: ...

class Tunnel:
    """
    One tunnel (a path entry) of a ForwardingProgram, read-only. index: the entry identifier (T1 and
    T2: j*H + i for round j and roster position i, T2's dead heads leaving holes; Global and CostOnly
    MAX_FLOW: the path's position in the decomposition order; supplied routes: the route's position in
    the list), never a position in program.tunnels; head; local_index: its position among the head's entries of both kinds;
    destination: where its path ends; edges: its EdgeIds in path order; route_cost: the sum of their
    costs (all three None for a tunnel that is down); status: UP, NO_BANDWIDTH or COST_CEILING;
    requested and reserved: T1's strict bandwidth (its request, and what it holds on every link of
    its path: the request when UP, else 0; 0 for the other models); planned: T2's planned bandwidth
    u (metadata, never held; 0 for the other models); weight: its integer share among the head's
    entries (T1: 1 up, 0 down; T2: 1; a controller path: its solved amount).
    """

    def __init__(self, _no_public_constructor: typing.Never, /) -> None: ...
    def __repr__(self) -> str: ...
    @property
    def destination(self) -> int | None: ...
    @property
    def edges(self) -> tuple[int, ...] | None:
        """
        the path's EdgeIds in path order; None for a tunnel that is down
        """

    @property
    def head(self) -> int: ...
    @property
    def index(self) -> int: ...
    @property
    def local_index(self) -> int: ...
    @property
    def planned(self) -> int: ...
    @property
    def requested(self) -> int: ...
    @property
    def reserved(self) -> int: ...
    @property
    def route_cost(self) -> int | None: ...
    @property
    def status(self) -> TunnelStatus: ...
    @property
    def weight(self) -> int: ...

class Route:
    """
    A supplied route for ForwardingProgram.from_routes: a Path with one edge per hop (a tunnel) or a
    PredDAG (a split route) starting at `head`; weight is its integer share among the head's routes
    (0 allowed: it carries nothing). from_routes validates it.
    """

    def __init__(
        self,
        head: typing.SupportsIndex,
        bundle: Path | PredDAG,
        *,
        weight: typing.SupportsIndex = 1,
    ) -> None:
        """
        Route(head, bundle, *, weight=1): head a NodeId, bundle a Path or a PredDAG (TypeError otherwise),
        weight an integer (checked in [0, 2^62) by from_routes).
        """

    def __repr__(self) -> str: ...
    @property
    def bundle(self) -> Path | PredDAG:
        """
        the Path (a tunnel) or the PredDAG (a split route)
        """

    @property
    def head(self) -> int:
        """
        the route's head (a source of the program)
        """

    @property
    def weight(self) -> int:
        """
        the route's integer share among the head's routes
        """

class ForwardingProgram:
    """
    An immutable value: one demand's routes and integer weights. sources and destinations (the normalized
    sets); tunnels: the path entries in ascending index (split routes are not listed); source_weights: the
    captured source supplies ((source, weight) for every roster source, ascending) of the optimized
    controllers, None for equal sources (round robin over the roster). Each head splits its share by
    highest averages over its entries' weights; identical routes are never merged. CostOnly and Global
    produce one (FlowPolicy.program, PlacementResult.program), and so do T1 and T2 (their tunnels);
    from_routes builds one from given routes; Pinned forwards one. No public constructor.
    """

    def __init__(self, _no_public_constructor: typing.Never, /) -> None: ...
    @staticmethod
    def from_routes(
        graph: Graph,
        sources: typing.SupportsIndex
        | collections.abc.Sequence[typing.SupportsIndex]
        | numpy.typing.NDArray[numpy.integer],
        destinations: typing.SupportsIndex
        | collections.abc.Sequence[typing.SupportsIndex]
        | numpy.typing.NDArray[numpy.integer],
        routes: collections.abc.Sequence[Route],
        *,
        placement: FlowPlacement | None = None,
        weights: WeightBasis | None = None,
        edge_weights: collections.abc.Sequence[typing.SupportsIndex]
        | numpy.typing.NDArray[numpy.integer]
        | None = None,
        residual: collections.abc.Sequence[typing.SupportsIndex]
        | numpy.typing.NDArray[numpy.integer]
        | None = None,
    ) -> ForwardingProgram:
        """
        Validate and freeze supplied routes -> ForwardingProgram (ValueError for the first failure, in this
        order): the endpoint sets; `routes must not be empty`; when no route is a PredDAG, any option given
        (`<placement, weights, edge_weights or residual> applies only to PredDAG routes`); placement MAX_FLOW
        (`from_routes placement must be EQUAL_BALANCED or WEIGHTED`); weights only with placement=WEIGHTED,
        EXPLICIT with edge_weights, RESIDUAL_SNAPSHOT with residual; edge_weights' length and entries in [0,
        2^62) and residual's (one entry per edge in [0, capacity]); then each route in list order: its
        bundle's graph, its head a source, a Path with one edge per hop that starts at the head, ends at a
        destination and meets none before, a PredDAG with a walk from the head to a destination, its weight
        in [0, 2^62). Order and duplicates are kept; a route of weight 0 is validated in full. placement
        (None: EQUAL_BALANCED) is the interior split of the PredDAG routes; each split route is its own
        split context (a router splits only that route's share). The roster is the sources with a route.
        """

    @property
    def destinations(self) -> tuple[int, ...]:
        """
        the destination set, ascending
        """

    @property
    def source_weights(self) -> tuple[tuple[int, int], ...] | None:
        """
        the captured source weights ((source, weight), ascending), or None for equal sources
        """

    @property
    def sources(self) -> tuple[int, ...]:
        """
        the source set, ascending
        """

    @property
    def tunnels(self) -> tuple[Tunnel, ...]:
        """
        the program's tunnels (Tunnel), ascending index
        """

class Absolute:
    """
    Cost ceiling: an absolute, inclusive path-cost budget (cost >= 0).
    """

    def __eq__(self, arg0: typing.Any) -> bool: ...
    def __hash__(self) -> int: ...
    def __init__(self, cost: typing.SupportsIndex) -> None:
        """
        Absolute(cost): a 64-bit integer budget in cost units (checked >= 0 where it is used)
        """

    def __repr__(self) -> str: ...
    @property
    def cost(self) -> int:
        """
        the budget: a path qualifies when its cost is at most this
        """

class Stretch:
    """
    Cost ceiling: a route qualifies when its cost is at most factor times the static shortest cost
    under the Eligibility; factor is in [1, inf), and the float is compared exactly, without rounding.
    """

    def __eq__(self, arg0: typing.Any) -> bool: ...
    def __hash__(self) -> int: ...
    def __init__(self, factor: float | numpy.floating | typing.SupportsIndex) -> None:
        """
        Stretch(factor): a number, checked finite and >= 1 where it is used; a path qualifies when its
        cost is at most factor times the static shortest cost. Here, an integer beyond the largest finite
        float raises ValueError.
        """

    def __repr__(self) -> str: ...
    @property
    def factor(self) -> float:
        """
        the factor as a float
        """

class Tunnels:
    """
    The tunnel count of T1 (CapacityAware) and T2, per head: count in [1, 10000] (checked when a
    policy or an evaluator uses it). up_to=True (T2 only; T1 rejects it): one common count k in
    1..count for every head, the smallest that delivers the design.
    """

    def __eq__(self, arg0: typing.Any) -> bool: ...
    def __hash__(self) -> int: ...
    def __init__(self, count: typing.SupportsIndex, *, up_to: bool = False) -> None:
        """
        Tunnels(count, *, up_to=False): a 64-bit integer count (a bool, a float or an integer outside 64 bits
        raises ValueError here; [1, 10000] is checked where it is used) and a bool
        """

    def __repr__(self) -> str: ...
    @property
    def count(self) -> int:
        """
        the tunnels per head (the largest count for up_to)
        """

    @property
    def up_to(self) -> bool:
        """
        True: the smallest sufficient count up to `count` (T2 only)
        """

class CostOnly:
    """
    Routing like IGP forwarding (OSPF, IS-IS routers) on each roster source's static shortest-path
    DAG (read from the DagCache; the selection is applied at the epoch's first positive call),
    fixed for the epoch; one record per roster source. EQUAL_BALANCED and WEIGHTED pool the sources
    (origin EQUAL): where their routes meet, a router splits their traffic once. MAX_FLOW is a
    controller on the union of those routes (origin OPTIMIZED or EQUAL). An edge of capacity 0
    stays in the DAG: EQUAL_BALANCED still sends it its share, which it cannot carry, and WEIGHTED
    gives it weight 0.
    """

    def __init__(self) -> None:
        """
        CostOnly routing (no parameters)
        """

    def __repr__(self) -> str: ...

class CapacityAware:
    """
    T1: every head sets up its own strict tunnels by CSPF, no controller (like RSVP-TE head ends).
    The design volume D is split equally over the roster heads and each head's share equally over its
    tunnels.count tunnels; in round order (round j, heads ascending: tunnel index j*H + i) each tunnel
    keeps the allowed links whose unreserved capacity (the ledger's, minus the holds of every other
    owner and of this setup's earlier tunnels) is at least max(request, min_unreserved_bandwidth) (a
    request of 0 with a minimum of 0 has no filter) and takes the shortest path in the selection's
    single-path order; it comes up holding exactly its request on every link (UP), or stays down
    holding nothing: NO_BANDWIDTH (no path passes the filter) or COST_CEILING (the cheapest exceeds
    max_path_cost: Absolute is inclusive; Stretch(f) is f times the smallest static shortest cost over
    all roster heads). Traffic enters equally at the heads and splits equally over each head's up
    tunnels; a head with no tunnel up loses its share. Holds never police traffic.
    """

    def __init__(
        self,
        tunnels: Tunnels,
        *,
        max_path_cost: Absolute | Stretch | None = None,
        min_unreserved_bandwidth: typing.SupportsIndex = 0,
    ) -> None:
        """
        CapacityAware(tunnels, *, max_path_cost=None, min_unreserved_bandwidth=0): the fields are checked
        when a policy is built with them (the count in [1, 10000] and exact, the ceiling,
        min_unreserved_bandwidth in [0, 2^62) units). Here, a min_unreserved_bandwidth of the wrong type or
        outside 64 bits raises ValueError; tunnels that is not a Tunnels, a ceiling of the wrong class and
        an unknown keyword raise TypeError.
        """

    def __repr__(self) -> str: ...
    @property
    def max_path_cost(self) -> Absolute | Stretch | None:
        """
        the inclusive route cost ceiling (Absolute or Stretch); None: no ceiling
        """

    @property
    def min_unreserved_bandwidth(self) -> int:
        """
        units in [0, 2^62): a search minimum only (each attempt keeps the links with at least
        max(request, this) unreserved units), never a reservation
        """

    @property
    def tunnels(self) -> Tunnels:
        """
        the tunnel count per head (a Tunnels)
        """

class Global:
    """
    Routing: an ideal central controller that sees the whole network (with `tunnels`: T2, below).
    Every positive call re-solves the minimum-cost flow of what the policy already delivers plus this
    call's request (placed_demand() + volume; what an earlier call could not place is not requested
    again) over the eligible graph, or exactly over `support` (a PredDAG of the policy's graph, trimmed to the
    walks from the sources to the allowed destinations under the Eligibility; every other edge
    carries nothing), on the free capacity plus its own allocation (it may move its own traffic,
    never another's). It takes no placement and no selection. Origin OPTIMIZED (the default): the
    controller chooses how much enters at each source; EQUAL: every roster source sends its
    round-robin share of the largest total whose shares all fit. Its program is the solve cut into
    paths per source, each a tunnel weighted by its amount. A positive call on a network beyond the
    flow core's int32 numbering raises ValueError (see FlowPolicy.place).
    T2, Global(tunnels=Tunnels(N)) (support must be None; no selection either): a controller (like
    a PCE) lays out, for every live head (one reaching a destination over links with free capacity
    >= 1; the other heads are dead and lose their share), N single-path tunnels with one common
    whole-unit planned bandwidth: first as large as possible (U), then for a finite design no larger
    than needed (u), then the cheapest routes, then with prefer_headroom=True the routes that leave
    the most spare capacity, then one fixed unit-by-unit tie. Tunnels(L, up_to=True): the smallest
    common count k <= L whose program carries the design (else the best, then the smallest).
    Traffic splits equally over the heads and over each head's tunnels; planned bandwidth is never
    held.
    """

    def __init__(
        self,
        *,
        support: PredDAG | None = None,
        tunnels: Tunnels | None = None,
        prefer_headroom: bool = False,
    ) -> None:
        """
        Global(*, support=None, tunnels=None, prefer_headroom=False): support a PredDAG of the policy's
        graph (checked when a policy is built with it: `support belongs to another graph`) or None for
        every eligible edge; tunnels a Tunnels (T2; checked when a policy is built with it: the count in
        [1, 10000], and `Global support and tunnels are mutually exclusive`) or None; prefer_headroom
        (T2 only; True without tunnels raises `prefer_headroom requires tunnels` when a policy is built
        with it): the headroom penalty as T2's objective after route cost
        """

    def __repr__(self) -> str: ...
    @property
    def prefer_headroom(self) -> bool:
        """
        True: T2 minimizes the headroom penalty after route cost (T2 only)
        """

    @property
    def support(self) -> PredDAG | None:
        """
        the supplied support (a PredDAG), or None
        """

    @property
    def tunnels(self) -> Tunnels | None:
        """
        T2's tunnel count per head (a Tunnels), or None
        """

class Pinned:
    """
    Routing: forwarding over a given program, computing nothing (the program's endpoints must be
    the policy's). The policy's Eligibility repairs failures: a tunnel touching an excluded element is
    FAILED, a split route loses its excluded members, and the roster is rebuilt from the program's
    basis (a head that no longer reaches a destination leaves the equal pool). black_holes: the nodes
    and edges this Eligibility excludes are silent failures; the program is forwarded unchanged and
    traffic reaching them is dropped there (stale forwarding; None: no black holes). One record per
    program entry.
    """

    def __init__(
        self, program: ForwardingProgram, *, black_holes: Eligibility | None = None
    ) -> None:
        """
        Pinned(program, *, black_holes=None): a ForwardingProgram and an Eligibility of the same graph (both
        kept by reference; checked when a policy is built with them)
        """

    def __repr__(self) -> str: ...
    @property
    def black_holes(self) -> Eligibility | None:
        """
        the black holes (an Eligibility), or None
        """

    @property
    def program(self) -> ForwardingProgram:
        """
        the forwarded program (shared, immutable)
        """

class FlowPolicyConfig:
    """
    One demand's settings: routing (CostOnly, CapacityAware (T1), Global, Global(tunnels=...) (T2) or
    Pinned), placement (FlowPlacement: required for CostOnly; the other models take none), admission
    (LOSSLESS by default; LOSSY for CostOnly EQUAL_BALANCED and WEIGHTED, T1, T2 and Pinned), selection
    (PathSelection; None: PathSelection(); Global, T2 included, and Pinned take none; T1 reads tie_break,
    priority and multi_edge and rejects multipath=True), origin (Origin; None: the model's default, EQUAL for CostOnly EQUAL_BALANCED and WEIGHTED, T1 and T2,
    OPTIMIZED for CostOnly MAX_FLOW and Global; T2 takes OPTIMIZED with one input source only, the same
    calculation; Pinned takes none, its program has it), weights (WeightBasis, CostOnly WEIGHTED only;
    None: CAPACITY) and edge_weights (one integer per edge, given exactly when weights=EXPLICIT; one
    array shared by every policy built from this configuration). A field that does not apply to the model must stay None. The combination is checked
    when a policy or an evaluator uses the configuration (ValueError): the enumeration values, then the
    model's row (Global: a selection first, then support and tunnels together, then prefer_headroom
    without tunnels; T1: an up_to count; then placement, origin,
    admission, selection, weights with edge_weights), then a DagCache (only CostOnly reads one), then
    the tunnel models' parameters (the tunnel count in [1, 10000]; T1's max_path_cost and
    min_unreserved_bandwidth), then edge_weights' length and entries, then the roster and count rules
    (CostOnly with several input sources rejects multipath=False with PREFER_HIGHER_PRIORITY; T2 with
    several input sources rejects OPTIMIZED; roster heads times the tunnel count at most 2^30 - 1).
    """

    def __init__(
        self,
        *,
        routing: CostOnly | CapacityAware | Global | Pinned,
        placement: FlowPlacement | None = None,
        admission: Admission = ...,
        selection: PathSelection | None = None,
        origin: Origin | None = None,
        weights: WeightBasis | None = None,
        edge_weights: collections.abc.Sequence[typing.SupportsIndex]
        | numpy.typing.NDArray[numpy.integer]
        | None = None,
    ) -> None: ...
    def __repr__(self) -> str: ...
    @property
    def admission(self) -> Admission:
        """
        LOSSLESS or LOSSY
        """

    @property
    def edge_weights(self) -> tuple[int, ...] | None:
        """
        the explicit weights (one per edge, a tuple), or None
        """

    @property
    def origin(self) -> Origin | None:
        """
        the Origin, or None (a policy's config: the model's default resolved; Pinned: None)
        """

    @property
    def placement(self) -> FlowPlacement | None:
        """
        the FlowPlacement, or None (Global, Pinned)
        """

    @property
    def routing(self) -> CostOnly | CapacityAware | Global | Pinned:
        """
        the routing mode: CostOnly, CapacityAware, Global or Pinned
        """

    @property
    def selection(self) -> PathSelection | None:
        """
        the PathSelection, or None (a policy's config: PathSelection() resolved; Global, Pinned: None)
        """

    @property
    def weights(self) -> WeightBasis | None:
        """
        the WeightBasis, or None (a policy's config: CAPACITY resolved for CostOnly WEIGHTED)
        """

class OutputOptions:
    """
    Optional outputs of a placement (all false by default; absent outputs are None):
    edge_flows (the policy's own allocation per edge after the call), residual and
    bottleneck (the ledger's after the call), reservation_edge_flows (the policy's own strict
    holds per edge after the call; zeros for a model without holds), and min_cut (Global with
    the optimized origin, when its solve reached the maximum flow; None for every other
    model, origin or call).
    """

    def __init__(
        self,
        *,
        edge_flows: bool = False,
        residual: bool = False,
        bottleneck: bool = False,
        reservation_edge_flows: bool = False,
        min_cut: bool = False,
    ) -> None: ...
    @property
    def bottleneck(self) -> bool:
        """
        request PlacementResult.bottleneck
        """

    @property
    def edge_flows(self) -> bool:
        """
        request PlacementResult.edge_flows
        """

    @property
    def min_cut(self) -> bool:
        """
        request PlacementResult.min_cut
        """

    @property
    def reservation_edge_flows(self) -> bool:
        """
        request PlacementResult.reservation_edge_flows
        """

    @property
    def residual(self) -> bool:
        """
        request PlacementResult.residual
        """

class Termination:
    """
    How a call ended, the first that applies: COMPLETE for a zero request; NO_ROUTE for a positive
    request with an empty roster (no eligible source statically reaches an allowed destination);
    NO_ACTIVE_WEIGHT when the program has an entry with a path (an up tunnel, a controller path or
    a split route) but no head can use one (no positive weight: WEIGHTED over edges of capacity 0
    only, routes of weight 0, every up tunnel failed); a program whose tunnels are all down is not
    this case; a LOSSY finite request: COMPLETE (delivery and losses account for the whole offer,
    which does not mean nothing was lost); otherwise COMPLETE when the whole request was delivered,
    else EXHAUSTED (an infinite request is EXHAUSTED). A controller whose maximum is 0 over existing
    routes is EXHAUSTED, not NO_ROUTE.

    Members:

      COMPLETE

      NO_ROUTE

      EXHAUSTED

      NO_ACTIVE_WEIGHT
    """

    COMPLETE: typing.ClassVar[Termination]
    EXHAUSTED: typing.ClassVar[Termination]
    NO_ACTIVE_WEIGHT: typing.ClassVar[Termination]
    NO_ROUTE: typing.ClassVar[Termination]
    __members__: typing.ClassVar[dict[str, Termination]]

    @typing.overload
    def __eq__(self, other: Termination, /) -> bool: ...
    @typing.overload
    def __eq__(self, other: typing.Any, /) -> bool: ...
    def __getstate__(self) -> int: ...
    def __hash__(self) -> int: ...
    def __index__(self) -> int: ...
    def __init__(self, value: typing.SupportsInt | typing.SupportsIndex) -> None: ...
    def __int__(self) -> int: ...
    @typing.overload
    def __ne__(self, other: Termination, /) -> bool: ...
    @typing.overload
    def __ne__(self, other: typing.Any, /) -> bool: ...
    def __repr__(self) -> str: ...
    def __setstate__(
        self, state: typing.SupportsInt | typing.SupportsIndex, /
    ) -> None: ...
    def __str__(self) -> str: ...
    @property
    def name(self) -> str: ...
    @property
    def value(self) -> int: ...

class SourceFlow:
    """
    One normalized input source's forwarding in the epoch (read-only; zero rows included: a
    masked or unreachable source): source, offered (LOSSY: its injections; None when
    LOSSLESS), admitted (LOSSLESS: what it delivers; LOSSY: offered), delivered, source_loss
    (offered units lost at the source for lack of a usable route) and dropped (its losses on
    links). LOSSY: offered = delivered + source_loss + dropped. Exact ints.
    """

    def __init__(self, _no_public_constructor: typing.Never, /) -> None: ...
    def __repr__(self) -> str: ...
    @property
    def admitted(self) -> int:
        """
        LOSSLESS: the units it delivers; LOSSY: offered
        """

    @property
    def delivered(self) -> int:
        """
        the units it delivers to the destinations
        """

    @property
    def dropped(self) -> int:
        """
        LOSSY: its units dropped on links; 0 when LOSSLESS
        """

    @property
    def offered(self) -> int | None:
        """
        LOSSY: the units this source injected; None when LOSSLESS
        """

    @property
    def source(self) -> int:
        """
        the source node id
        """

    @property
    def source_loss(self) -> int:
        """
        LOSSY: its offered units lost at the source (no usable route); 0 when LOSSLESS
        """

class TunnelFlow:
    """
    One program tunnel's forwarding in the epoch (read-only): index (its entry identifier) and
    head, status (TunnelStatus: UP; NO_BANDWIDTH or COST_CEILING for a tunnel T1's setup left
    down; FAILED, BLACKHOLED in this forwarding), assigned (its share of the injections: LOSSY the
    shares sent to it, LOSSLESS what it carries), delivered and dropped (its losses on links).
    Exact ints (assigned and dropped may exceed 64 bits).
    """

    def __init__(self, _no_public_constructor: typing.Never, /) -> None: ...
    def __repr__(self) -> str: ...
    @property
    def assigned(self) -> int:
        """
        its share of the injections in the epoch (LOSSY: sent to it; LOSSLESS: what it carries)
        """

    @property
    def delivered(self) -> int:
        """
        the units it delivers
        """

    @property
    def dropped(self) -> int:
        """
        its units dropped on links
        """

    @property
    def head(self) -> int:
        """
        the tunnel's head
        """

    @property
    def index(self) -> int:
        """
        the tunnel's entry identifier
        """

    @property
    def status(self) -> TunnelStatus:
        """
        UP, a down tunnel's setup status (NO_BANDWIDTH, COST_CEILING), FAILED or BLACKHOLED in this
        forwarding
        """

class TunnelPlan:
    """
    The tunnel models' setup summary (read-only): count (the tunnels per head; T2 up to: the
    chosen count), design_demand (the design volume D in units), head_design ((head, RR(D, H, i))
    for every roster head in ascending NodeId) and reserved_total (T1's strict holds, counted once
    per up tunnel: the sum of their requests; 0 for T2). T2 only (None for T1): live_heads,
    insufficient_resolution (live heads exist but no whole unit per tunnel fits, U == 0: the
    program has no tunnels), fluid_share (the exact fluid optimum s*, a Fraction), common_share
    (U), planned_per_tunnel (u), planned_total (K*u), route_cost (the route-cost sum, can exceed
    2^64), headroom_penalty (with prefer_headroom only) and forwarding_bound (the program's
    lossless maximum; a dead head's share blocks it). Exact ints.
    """

    def __init__(self, _no_public_constructor: typing.Never, /) -> None: ...
    def __repr__(self) -> str: ...
    @property
    def common_share(self) -> int | None:
        """
        T2: U, the largest common whole-unit size
        """

    @property
    def count(self) -> int:
        """
        the tunnels per head (T2 up to: the chosen count)
        """

    @property
    def design_demand(self) -> int | float:
        """
        the design volume D in units; float('inf') when unlimited
        """

    @property
    def fluid_share(self) -> fractions.Fraction | None:
        """
        T2: the exact fluid optimum s* (a fractions.Fraction; 0 without a live head); None for T1
        """

    @property
    def forwarding_bound(self) -> int | None:
        """
        T2: F, the program's lossless maximum on the call's view (a dead head's share blocks it)
        """

    @property
    def head_design(self) -> tuple[tuple[int, int], ...]:
        """
        ((head, share), ...): every roster head's share RR(D, H, i), ascending NodeId
        """

    @property
    def headroom_penalty(self) -> int | None:
        """
        T2 with prefer_headroom: the sum over the links' tunnel slots j of floor(2^32 (B - V + j*u) / B)
        (B: static capacity; V: free capacity); None without prefer_headroom and for T1
        """

    @property
    def insufficient_resolution(self) -> bool | None:
        """
        T2: True when live heads exist but no whole unit per tunnel fits (common_share 0, no tunnels)
        """

    @property
    def live_heads(self) -> tuple[int, ...] | None:
        """
        T2: the heads that reach a destination over links with free capacity >= 1 (ascending); None for T1
        """

    @property
    def planned_per_tunnel(self) -> int | None:
        """
        T2: u, every tunnel's planned bandwidth
        """

    @property
    def planned_total(self) -> int | None:
        """
        T2: K*u over the K tunnels
        """

    @property
    def reserved_total(self) -> int:
        """
        the strict holds of the up tunnels, each tunnel's request counted once
        """

    @property
    def route_cost(self) -> int | None:
        """
        T2: the sum of cost * tunnel count over the links (exact; can exceed 2^64)
        """

class PlacementResult:
    """
    The result of one call. This call: placed (units delivered by the call; replace, forward and
    evaluate: the full new delivery), remaining (request - placed; float('inf') for an infinite
    request), delivery_delta (new minus old delivery, may be negative), call_drops [(EdgeId, units)],
    call_source_loss, termination, routes_changed, weights_changed and reservations_changed. The
    demand's state after the call (totals of the epoch): program (the installed ForwardingProgram;
    evaluate: the computed one; None before the first preparation), plan (a TunnelPlan for T1 and T2;
    None for every other model), requested (float('inf') after an infinite request), offered (LOSSY: the traffic sent; None when LOSSLESS), admitted, delivered, unadmitted
    (LOSSLESS: requested - delivered, None after an infinite request; LOSSY: 0), source_loss, drops
    [(EdgeId, units)], by_source (one SourceFlow per normalized input source, ascending NodeId),
    by_tunnel (one TunnelFlow per program tunnel, ascending index; an alternative view, never added to
    by_source), total_cost (controllers: Global and CostOnly MAX_FLOW; else None) and cost_distribution
    (unrestricted Global with the optimized origin only). LOSSY: offered = source_loss + sum(drops) +
    delivered. The rows sum to the totals, except for LOSSY with an empty roster, where the whole offer
    is the demand's source loss and every row is zero. Optional outputs, None unless requested:
    edge_flows, residual and reservation_edge_flows (int64[E]), bottleneck, min_cut (they describe the same
    state; evaluate: the evaluated one). Amounts that accumulate are exact ints (they may exceed 64 bits).
    """

    def __init__(self, _no_public_constructor: typing.Never, /) -> None: ...
    def __repr__(self) -> str: ...
    @property
    def admitted(self) -> int:
        """
        LOSSLESS: the total carried; LOSSY: offered
        """

    @property
    def bottleneck(self) -> list[int] | None:
        """
        if requested: the ledger's saturation report after the call (FlowLedger.bottleneck of the
        policy's sources and Eligibility), else None
        """

    @property
    def by_source(self) -> tuple[SourceFlow, ...]:
        """
        one SourceFlow per normalized input source, ascending NodeId (zero rows included)
        """

    @property
    def by_tunnel(self) -> tuple[TunnelFlow, ...]:
        """
        one TunnelFlow per program tunnel, ascending index (empty without tunnels)
        """

    @property
    def call_drops(self) -> list[tuple[int, int]]:
        """
        this call's drops on links [(EdgeId, units)], ascending EdgeId (empty unless LOSSY)
        """

    @property
    def call_source_loss(self) -> int:
        """
        this call's offered units lost at sources without a usable route (LOSSY)
        """

    @property
    def cost_distribution(self) -> list[tuple[int, int]] | None:
        """
        unrestricted Global with the optimized origin: the marginal cost curve [(cost, units)] of the policy's
        cumulative flow (a zero call reports it unchanged; None before the epoch's first positive
        call); in ascending cost, delivering `units` more raised the minimum total flow cost by `cost`
        per unit (more delivery can reroute earlier flow, so these are not the costs of the final
        paths). None for the other models and the equal origin
        """

    @property
    def delivered(self) -> int:
        """
        the units the demand delivers
        """

    @property
    def delivery_delta(self) -> int:
        """
        the delivery after the call minus before it (placed for a placement; may be negative for a
        full-offer call)
        """

    @property
    def drops(self) -> list[tuple[int, int]]:
        """
        the epoch's drops on links [(EdgeId, units)], ascending EdgeId (empty unless LOSSY)
        """

    @property
    def edge_flows(self) -> numpy.typing.NDArray[numpy.int64] | None:
        """
        if requested: the policy's own allocation per edge after the call (int64[E], EdgeId order),
        else None
        """

    @property
    def min_cut(self) -> list[int] | None:
        """
        if requested (OutputOptions(min_cut=True)), unrestricted Global with the optimized origin, when the solve
        reached the maximum flow: the minimum cut closest to the sources, every allowed edge leaving
        the nodes the roster sources reach in the residual network (free capacity plus the policy's
        own allocation), an edge with nothing free (capacity 0, or held by others) included,
        ascending EdgeId; [] for an empty roster. None otherwise (another model, the equal origin, a
        support, a call below the maximum, a zero call, forward)
        """

    @property
    def offered(self) -> int | None:
        """
        LOSSY: the traffic sent in the epoch (the first-hop injection of an infinite request); None
        when LOSSLESS
        """

    @property
    def placed(self) -> int:
        """
        the units placed (delivered) by this call
        """

    @property
    def plan(self) -> TunnelPlan | None:
        """
        T1 and T2: the TunnelPlan of the program after the call; None for every other model (Pinned over
        a tunnel program too)
        """

    @property
    def program(self) -> ForwardingProgram | None:
        """
        the installed ForwardingProgram after the call (evaluate: the computed one); None before the first
        preparation
        """

    @property
    def remaining(self) -> int | float:
        """
        the request minus placed, in units (an amount, not a status); float('inf') for an infinite request
        """

    @property
    def requested(self) -> int | float:
        """
        the sum of the epoch's requests; float('inf') after an infinite one
        """

    @property
    def reservation_edge_flows(self) -> numpy.typing.NDArray[numpy.int64] | None:
        """
        if requested: the policy's own strict holds per edge after the call (int64[E], EdgeId order; zeros
        for a model without holds), else None
        """

    @property
    def reservations_changed(self) -> bool:
        """
        True when a tunnel's requested or reserved amount differs from the installed program's, by
        index (an absent tunnel holds 0): only T1 sets it; forward and Pinned never do
        """

    @property
    def residual(self) -> numpy.typing.NDArray[numpy.int64] | None:
        """
        if requested: the ledger's residual after the call (int64[E], EdgeId order), else None
        """

    @property
    def routes_changed(self) -> bool:
        """
        True when the call's program differs from the installed one in its ordered heads, paths, split
        tables or duplicates (no program counts as the empty program; two programs without an entry
        are equal). An offer-only change (forward, Pinned, CostOnly increments) and an identical
        replace set nothing
        """

    @property
    def source_loss(self) -> int:
        """
        LOSSY: the offered units lost at sources (no usable route, or an empty roster); 0 when LOSSLESS
        """

    @property
    def termination(self) -> Termination:
        """
        the call's outcome (Termination)
        """

    @property
    def total_cost(self) -> int | None:
        """
        controllers (Global, CostOnly MAX_FLOW): the exact sum of cost * flow over the policy's
        allocation; None for the other models (T1, T2, Pinned)
        """

    @property
    def unadmitted(self) -> int | None:
        """
        LOSSLESS: requested - delivered, None after an infinite request; LOSSY: 0
        """

    @property
    def weights_changed(self) -> bool:
        """
        True when the call's program differs from the installed one in its injection (source weights)
        or any entry or split-table weight (a controller's re-solve that changes a path's amount; an
        absent entry weighs 0)
        """

class DagCacheStats:
    """
    DagCache statistics (atomic per field, not a consistent snapshot): hits, misses,
    entries, bytes (the entries' retained payload), live_values (values alive anywhere).
    """

    def __init__(self, _no_public_constructor: typing.Never, /) -> None: ...
    @property
    def bytes(self) -> int:
        """
        the retained payload of the entries (their DAGs' stored arrays)
        """

    @property
    def entries(self) -> int:
        """
        the entries held
        """

    @property
    def hits(self) -> int:
        """
        lookups answered from an entry (they persist across clear())
        """

    @property
    def live_values(self) -> int:
        """
        the cached values alive anywhere, entries included
        """

    @property
    def misses(self) -> int:
        """
        lookups without an entry: a building cache searched, a frozen one raised DagCacheMiss
        """

class DagCacheMiss(Exception):
    """
    A lookup in a frozen DagCache whose key is absent: a dedicated exception (not a
    RuntimeError or ValueError); the evaluators propagate it unchanged.
    """

class DagCache:
    """
    A cache of static shortest-path DAGs (from costs and the Eligibility only) for one graph, keyed by
    (the root set, the direction, the Eligibility uid). It stores the full DAG and applies each
    caller's selection afterwards, so the selection is never cached. CostOnly policies and the
    evaluators read it.
    Building mode (a new cache): one writer, no concurrent readers; a miss searches and inserts;
    clear() empties it (the hit and miss counters persist), and values already handed out stay valid.
    Frozen mode (freeze()): read-only and shareable by every worker; a miss counts one miss and raises
    DagCacheMiss; clear() raises ValueError; thaw() returns to building once the caller has joined
    every reader.
    When a public call fails: a building cache passed to it keeps the complete entries and statistics
    of the lookups that succeeded before the failure (never a partial entry; a kept entry never
    changes a result); a frozen cache's statistics stay unchanged, except one miss when the call
    raised DagCacheMiss. A frozen cache counts its hits only when the call succeeds, Python result
    included.
    """

    def __init__(self, graph: Graph) -> None:
        """
        An empty cache in building mode, bound to graph (by its uid).
        """

    def clear(self) -> None:
        """
        empty a building cache; ValueError when frozen
        """

    def freeze(self) -> None:
        """
        enter frozen mode (no effect when frozen)
        """

    def graph_uid(self) -> int:
        """
        uid of the graph this cache is bound to
        """

    def stats(self) -> DagCacheStats:
        """
        the current DagCacheStats (each field read atomically)
        """

    def thaw(self) -> None:
        """
        return to building mode (no effect when building); join every reader first
        """

class FlowPolicy:
    """
    The placement engine for one demand, bound to its ledger at construction. Every call (place,
    evaluate, replace, forward) is all-or-nothing: a failed call raises and leaves the ledger and the
    policy as they were (with the two limits in the module docstring). remove() takes every flow of the
    policy back out, clears its program and starts a new epoch (the request total starts again from 0);
    replace and forward start one too. A policy dropped without remove()
    leaves its records in the ledger as orphans: readable, never changed again, cleared by
    FlowLedger.reset(). sources and destinations are one NodeId or a sequence; the two sets must not
    overlap. ValueError for invalid endpoints, configurations and requests (non-integer node ids and
    requests included) and for objects of another graph, ledger or policy; TypeError for a ledger,
    config, elig, cache or outputs of the wrong class (a ledger of None: ValueError) and for an
    unknown keyword.
    """

    def __init__(
        self,
        ledger: FlowLedger,
        sources: typing.SupportsIndex
        | collections.abc.Sequence[typing.SupportsIndex]
        | numpy.typing.NDArray[numpy.integer],
        destinations: typing.SupportsIndex
        | collections.abc.Sequence[typing.SupportsIndex]
        | numpy.typing.NDArray[numpy.integer],
        config: FlowPolicyConfig,
        elig: Eligibility | None = None,
        cache: DagCache | None = None,
    ) -> None:
        """
        Bind a policy to ledger. elig (None: everything allowed) and cache (None: none; only CostOnly
        routing reads it, another model rejects it) belong to the ledger's graph; the policy keeps the
        ledger, elig and cache alive (a Pinned policy shares its program and black holes). The roster
        (the sources that statically reach an allowed destination, over allowed edges whatever their
        capacity, the tunnel models' heads; Pinned: the program's heads that reach an allowed destination over its basis
        under elig) is fixed here. Everything is validated before the policy binds, in this order: the endpoints, the
        graph of elig, the program, the support, the cache and the black holes (`foreign Eligibility:
        built for another graph`, `program belongs to another graph`, `support belongs to another graph`,
        `foreign DagCache: bound to another graph`, `black_holes belongs to another graph`), the
        configuration (FlowPolicyConfig's order), then the pooled single-path rule, T2's several-sources
        rule (`T2 with several sources requires origin=EQUAL`), the tunnel limit (`tunnel limit exceeded:
        <H> heads x <N> tunnels > 1073741823`) and a Pinned program's endpoints (`program endpoints do not
        match policy endpoints`).
        """

    def evaluate(
        self,
        offer: typing.SupportsIndex | float | numpy.floating,
        *,
        design_demand: typing.SupportsIndex | float | numpy.floating | None = None,
        outputs: OutputOptions | None = None,
    ) -> PlacementResult:
        """
        A full, independent computation and forwarding of offer (int units or float('inf')) from the
        policy's own released view (the residual plus its own carried load; T1 sets up on the unreserved
        capacity with its own holds released; T2 lays out on the released residual), as replace would do it -> PlacementResult (the computed
        program and plan, the flags replace would set, the outputs of the evaluated state); publishes
        nothing, creates no record, books no hold and leaves the epoch as it was. design_demand as for
        replace. LOSSY requires a finite offer.
        """

    def forward(
        self,
        offer: typing.SupportsIndex | float | numpy.floating,
        *,
        black_holes: Eligibility | None = None,
        outputs: OutputOptions | None = None,
    ) -> PlacementResult:
        """
        Forward a full new offer over the installed program, with no route computation and no setup, on the
        released view -> PlacementResult; a new epoch like replace, the program, its plan and T1's holds
        kept (no flag; T1's records are recreated with the same reserved_edges). LOSSLESS admits the
        program's lossless bound up to the offer; LOSSY sends the offer. black_holes (an Eligibility of the
        graph): the nodes and edges it excludes are black holes for this call (stale forwarding: traffic
        reaching one is dropped there; a black-holed source loses its share). In this order: Pinned raises
        `Pinned forwards with replace`; no installed program `forward requires an installed program`;
        the offer; LOSSY +inf `LOSSY forwarding requires a finite offer`; black_holes
        of another graph `black_holes belongs to another graph`.
        """

    def place(
        self,
        volume: typing.SupportsIndex | float | numpy.floating,
        outputs: OutputOptions | None = None,
    ) -> PlacementResult:
        """
        Place volume (int units in [0, 2^62), or float('inf')) all or nothing -> PlacementResult. T1 and T2
        raise `<T1 or T2> does not support incremental placement` (also for 0).
        Requests add up within an epoch (from construction, the last remove(), replace() or forward());
        after float('inf') any finite request is rejected, 0 included. A request of 0 places nothing and
        reports the current state. CostOnly prepares its program at the epoch's first positive request;
        LOSSLESS admits the largest total in [admitted, admitted + volume] that the program carries with no
        loss (a refused unit is not offered again); LOSSY forwards the new round-robin (or captured-weight)
        injections, each head splitting its increment afresh (earlier drops are not offered again); a
        controller re-solves admitted + volume. A LOSSY float('inf') is CostOnly's first-hop experiment;
        Pinned raises `LOSSY forwarding requires a finite offer`. outputs (OutputOptions; None: none)
        requests edge_flows, residual, bottleneck, reservation_edge_flows and min_cut. Everything is
        validated before any change; a
        failed call raises with the ledger and the policy unchanged. A solve of the flow core (Global,
        CostOnly MAX_FLOW, or a MAX_FLOW flow) on a network beyond its int32 numbering raises ValueError:
        more than INT32_MAX - 2 nodes, or more than 2^30 - 1 edges, sources and destinations together. A
        LOSSY forwarding whose losses wait on each other in a cycle raises ValueError when its rounds
        alternate between two states or reach the round limit, RuntimeError when the result's check fails.
        """

    def placed_demand(self) -> int:
        """
        the sum of delivered over the policy's live records
        """

    def record(
        self, flow_id: tuple[typing.SupportsIndex, typing.SupportsIndex]
    ) -> FlowRecord:
        """
        one of the policy's live records by FlowId (ValueError for a FlowId of another ledger or policy, a
        retired or an unknown id)
        """

    def records(self) -> list[FlowRecord]:
        """
        the policy's live records as snapshots in FlowId order
        """

    def remove(self) -> None:
        """
        revert every flow and hold of the policy, retire the ids, clear the program, open a new epoch
        """

    def replace(
        self,
        offer: typing.SupportsIndex | float | numpy.floating,
        *,
        design_demand: typing.SupportsIndex | float | numpy.floating | None = None,
        outputs: OutputOptions | None = None,
    ) -> PlacementResult:
        """
        Recompute from the policy's own released view and forward the full offer, published atomically ->
        PlacementResult: the old records and holds are retired (their ids never live again), new ones
        created, and a new epoch starts at this offer (delivery may go down: delivery_delta). CostOnly
        prepares its routes again; a controller re-solves; T1 sets its tunnels up again for the design and
        books their holds; T2 lays its tunnels out again for the design (planned bandwidth, never held);
        Pinned forwards its program. design_demand (T1 and T2: int units, or float('inf'); None: equal to
        the offer): T1 needs a finite design (`T1 requires a finite design demand`), T2 takes an unlimited
        one (its largest common size); every other model raises `<model> does not accept design_demand`.
        In this order: the offer, LOSSY's finite offer, then the design demand (the model, its domain,
        finiteness).
        """

    @property
    def config(self) -> FlowPolicyConfig:
        """
        the configuration as validated at construction, its absent origin and selection
        resolved to the model's defaults
        """

    @property
    def program(self) -> ForwardingProgram | None:
        """
        the installed ForwardingProgram; None before the first preparation and after remove()
        """

class QueryOptions:
    """
    The options of every evaluator (capacity, batch_capacity, sensitivity): the configuration,
    the static Eligibility (None: all allowed), the DagCache (None: private per worker), the
    optional outputs of each result (sensitivity: of its baseline), design_demand (int units or
    float('inf'); None: absent): T1's design, always required and finite (`T1 capacity requires
    design_demand`, `T1 requires a finite design demand`); T2's design (None: unlimited when
    LOSSLESS, the offered amount when LOSSY); rejected by every other model (`<model> does not
    accept design_demand`); and offered (int units or float('inf'); None: absent), the offer of a
    T1, T2 or Pinned LOSSY evaluation: required and finite there (`<T1 or T2 or Pinned> LOSSY
    capacity requires finite offered`; T2 takes a finite design instead when offered is None); a
    LOSSLESS model rejects it (`<model> LOSSLESS capacity does not accept offered`), CostOnly LOSSY
    too (`CostOnly LOSSY capacity uses first-hop injection`).
    """

    def __init__(
        self,
        config: FlowPolicyConfig,
        *,
        elig: Eligibility | None = None,
        cache: DagCache | None = None,
        outputs: OutputOptions | None = None,
        design_demand: typing.SupportsIndex | float | numpy.floating | None = None,
        offered: typing.SupportsIndex | float | numpy.floating | None = None,
    ) -> None: ...
    @property
    def cache(self) -> DagCache | None:
        """
        the supplied DagCache; None: a private cache per worker
        """

    @property
    def config(self) -> FlowPolicyConfig:
        """
        the FlowPolicyConfig of every evaluated policy (sensitivity: the baseline and every deletion)
        """

    @property
    def design_demand(self) -> int | float | None:
        """
        the tunnel models' design volume (int units, or float('inf')); None: absent
        """

    @property
    def elig(self) -> Eligibility | None:
        """
        the static Eligibility; None: everything allowed
        """

    @property
    def offered(self) -> int | float | None:
        """
        the offer of a T1, T2 or Pinned LOSSY evaluation (int units, or float('inf')); None: absent
        """

    @property
    def outputs(self) -> OutputOptions:
        """
        the OutputOptions of every result (sensitivity: of the baseline)
        """

class SensitivityEntry:
    """
    One candidate edge: edge and loss (int units: the baseline minus the capacity without the
    edge).
    """

    def __init__(self, _no_public_constructor: typing.Never, /) -> None: ...
    @property
    def edge(self) -> int: ...
    @property
    def loss(self) -> int: ...

class SensitivityResult:
    """
    baseline (the PlacementResult of capacity()) and entries (ascending EdgeId; an
    absent edge has loss 0).
    """

    def __init__(self, _no_public_constructor: typing.Never, /) -> None: ...
    @property
    def baseline(self) -> PlacementResult: ...
    @property
    def entries(self) -> list[SensitivityEntry]: ...

class KspOptions:
    """
    The options of ksp: k (>= 1) and selection are explicit (multi_edge = True lists each parallel
    edge as a separate path, False keeps one per (u, v) pair as EdgeTieBreak describes; multipath is
    ignored); max_path_cost (None, Absolute(c) with c >= 0, or Stretch(f) with f finite in [1, inf)
    times the static shortest cost under the Eligibility, residual ignored; inclusive); residual
    (int64 units per EdgeId, each in [0, capacity_units]: edges with less than 1 unit are skipped, and
    the values feed the priority; it is the only capacity filter, and without it an edge of capacity 0
    is searched like any other); elig (None: everything allowed). The fields are validated by ksp,
    in the order its docstring gives.
    """

    def __init__(
        self,
        k: typing.SupportsIndex,
        selection: PathSelection,
        *,
        max_path_cost: Absolute | Stretch | None = None,
        residual: collections.abc.Sequence[typing.SupportsIndex]
        | numpy.typing.NDArray[numpy.integer]
        | None = None,
        elig: Eligibility | None = None,
    ) -> None: ...
    @property
    def elig(self) -> Eligibility | None:
        """
        the static Eligibility; None: everything allowed
        """

    @property
    def k(self) -> typing.SupportsIndex:
        """
        the number of paths as given (validated by ksp: an integer >= 1)
        """

    @property
    def max_path_cost(self) -> Absolute | Stretch | None:
        """
        the inclusive cost ceiling as given: None, Absolute or Stretch (validated by ksp)
        """

    @property
    def residual(
        self,
    ) -> (
        collections.abc.Sequence[typing.SupportsIndex]
        | numpy.typing.NDArray[numpy.integer]
        | None
    ):
        """
        the residual units per EdgeId as given (validated by ksp); None: no capacity filter;
        priorities read static capacity
        """

    @property
    def selection(self) -> PathSelection:
        """
        the PathSelection (multipath is ignored)
        """

class SpfOptions:
    """
    The options of spf: direction and selection are explicit; targets (one NodeId or a sequence;
    None: a full search) stops the search at the nearest targets (anycast: any nearest target will
    do); residual (int64 units per EdgeId, each in [0, capacity_units]) makes the search skip edges
    with fewer than min_residual residual units (an integer >= 1; without a residual it has no
    effect) and supplies the priority values: it is the only capacity filter, and without it an edge
    of capacity 0 is searched like any other; elig (None: everything allowed). A direction outside its
    enumeration raises ValueError here; the other fields are validated by spf, in the order its
    docstring gives.
    """

    def __init__(
        self,
        direction: Direction,
        selection: PathSelection,
        *,
        targets: typing.SupportsIndex
        | collections.abc.Sequence[typing.SupportsIndex]
        | numpy.typing.NDArray[numpy.integer]
        | None = None,
        residual: collections.abc.Sequence[typing.SupportsIndex]
        | numpy.typing.NDArray[numpy.integer]
        | None = None,
        min_residual: typing.SupportsIndex = 1,
        elig: Eligibility | None = None,
    ) -> None: ...
    @property
    def direction(self) -> Direction:
        """
        FORWARD or REVERSE
        """

    @property
    def elig(self) -> Eligibility | None:
        """
        the static Eligibility; None: everything allowed
        """

    @property
    def min_residual(self) -> typing.SupportsIndex:
        """
        the gating threshold in units as given (validated by spf: an integer >= 1)
        """

    @property
    def residual(
        self,
    ) -> (
        collections.abc.Sequence[typing.SupportsIndex]
        | numpy.typing.NDArray[numpy.integer]
        | None
    ):
        """
        the residual units per EdgeId as given (validated by spf); None: no capacity filter;
        priorities read static capacity
        """

    @property
    def selection(self) -> PathSelection:
        """
        the PathSelection the search reads
        """

    @property
    def targets(
        self,
    ) -> (
        typing.SupportsIndex
        | collections.abc.Sequence[typing.SupportsIndex]
        | numpy.typing.NDArray[numpy.integer]
        | None
    ):
        """
        the targets as given (validated by spf); None: a full search
        """

class SpfResult:
    """
    The result of spf: costs (full searches: a read-only int64 array, UNREACHABLE for the nodes not
    reached; None with targets), dag (PredDAG of the tight edges, the edges on some shortest walk,
    one per (u, v) pair with multi_edge = False, each in its own direction; with targets only the
    walks between the roots and a nearest target, for REVERSE from the target to a root; for
    multipath = False the canonical single path) and
    target_cost (with targets: the nearest targets' cost, None when no target is reached; None for a
    full search).
    """

    def __init__(self, _no_public_constructor: typing.Never, /) -> None: ...
    def __repr__(self) -> str: ...
    @property
    def costs(self) -> numpy.typing.NDArray[numpy.int64] | None:
        """
        full search: each node's distance from the roots for FORWARD, or to the roots for REVERSE
        (read-only int64[N], UNREACHABLE for nodes not reached); None with targets
        """

    @property
    def dag(self) -> PredDAG:
        """
        the PredDAG of the result, in the edges' original orientation
        """

    @property
    def target_cost(self) -> int | None:
        """
        with targets: the cost of the nearest targets, None when no target is reached; None for a
        full search
        """

@typing.overload
def as_float(
    x: collections.abc.Sequence[typing.SupportsIndex]
    | numpy.typing.NDArray[numpy.integer],
    resolution: typing.SupportsFloat | typing.SupportsIndex,
) -> numpy.typing.NDArray[numpy.float64]:
    """
    Units to physical value(s): x * q, correctly rounded (exact below 2^53 units). An int, of
    any size, gives a float, an integer array a float64 array; a result beyond the finite
    doubles raises ValueError; q must be a power of two in [2^-32, 2^32].
    """

@typing.overload
def as_float(
    x: typing.SupportsIndex, resolution: typing.SupportsFloat | typing.SupportsIndex
) -> float: ...
def batch_capacity(
    graph: Graph,
    pairs: collections.abc.Sequence[
        collections.abc.Sequence[
            typing.SupportsIndex
            | collections.abc.Sequence[typing.SupportsIndex]
            | numpy.typing.NDArray[numpy.integer]
        ]
        | numpy.typing.NDArray[numpy.integer]
    ]
    | numpy.typing.NDArray[numpy.integer],
    opts: QueryOptions,
    per_pair_elig: collections.abc.Sequence[Eligibility | None] | None = None,
    workers: typing.SupportsIndex | None = None,
) -> list[PlacementResult]:
    """
    capacity() for each (sources, destinations) pair, results in pair order. Pair i uses opts.elig
    intersected with per_pair_elig[i] (None entries allowed; a given sequence has one entry per pair,
    also when empty). workers >= 1; None: NGRAPH_CORE_BATCH_THREADS or 1, and 1 when a cache is
    supplied. Everything is validated before a worker starts, the arguments first, then each pair
    completely in pair order, so the error of the lowest failing pair is raised. Pinned is rejected
    (`batch_capacity does not accept Pinned`: a program has fixed endpoints);
    an empty batch returns [].
    Caches: a supplied building cache needs workers == 1; a frozen cache is shared by every worker
    (pass workers explicitly), and its DagCacheMiss reaches the caller unchanged (the lowest failing
    pair's). An intersection that differs from both operands gets a new Eligibility uid, which a
    frozen cache cannot hold: pass the combined per-pair values with elig None, or a common value
    that contains them.
    A failed call raises; a supplied DagCache then behaves as the DagCache docstring states.
    """

def capacity(
    graph: Graph,
    sources: typing.SupportsIndex
    | collections.abc.Sequence[typing.SupportsIndex]
    | numpy.typing.NDArray[numpy.integer],
    destinations: typing.SupportsIndex
    | collections.abc.Sequence[typing.SupportsIndex]
    | numpy.typing.NDArray[numpy.integer],
    opts: QueryOptions,
) -> PlacementResult:
    """
    place(+inf) of a FlowPolicy with opts.config on a fresh ledger of the graph (Pinned LOSSY:
    place(opts.offered); T1: evaluate(+inf or opts.offered, design_demand=opts.design_demand), the
    lossless bound of the tunnels set up for the design, or what they deliver of offered; T2: the same
    with an unlimited design when LOSSLESS and none is given, and the design as the offer when LOSSY
    without offered) -> PlacementResult.
    The exact maximum flow is Global() (the optimized origin); with origin=EQUAL, the largest total whose
    round-robin source shares are all carried. A supplied frozen cache raises DagCacheMiss for an absent
    key, unchanged.
    A failed call raises; a supplied DagCache then behaves as the DagCache docstring states.
    """

def ksp(
    graph: Graph, src: typing.SupportsIndex, dst: typing.SupportsIndex, opts: KspOptions
) -> list[tuple[int, Path]]:
    """
    The k shortest simple src -> dst paths over the eligible edges (only those with residual >= 1 when
    a residual is given; an edge of capacity 0 is searched like any other) -> [(cost, Path)], each hop
    a single edge. Order: lowest cost; then, under PREFER_HIGHER_PRIORITY only, highest bottleneck
    (the smallest priority value along the path); then fewest hops; then the smallest EdgeId
    sequence. The result is the first k paths of that order
    within the inclusive cost ceiling. Validated in this order before any work (ValueError): src and
    dst in range, src == dst, k >= 1, the ceiling, the residual (length, then every entry), a foreign
    Eligibility. An excluded endpoint, no static route under Stretch, or no path gives [].
    """

def resolve_to_paths(
    graph: Graph,
    dag: PredDAG,
    src: typing.SupportsIndex,
    dst: typing.SupportsIndex,
    split_parallel_edges: bool = False,
    max_paths: typing.SupportsIndex | None = None,
) -> list[Path]:
    """
    The src -> dst walks of a DAG -> [Path], in increasing (cost, hops, edge-id sequence):
    a hop costs its cheapest member and is keyed by its lowest member id. split_parallel_edges = False
    keeps every parallel member of a (u, v) pair in one hop (members ascending); True gives one path
    per member combination, sorted by that path's own (cost, hops, edge-id sequence). max_paths
    (>= 1) truncates; None lists every walk.
    ValueError: a DAG of another graph, src or dst out of range, src == dst, max_paths < 1. An
    endpoint outside the DAG, or no walk, gives [].
    """

def sensitivity(
    graph: Graph,
    sources: typing.SupportsIndex
    | collections.abc.Sequence[typing.SupportsIndex]
    | numpy.typing.NDArray[numpy.integer],
    destinations: typing.SupportsIndex
    | collections.abc.Sequence[typing.SupportsIndex]
    | numpy.typing.NDArray[numpy.integer],
    opts: QueryOptions,
    workers: typing.SupportsIndex | None = None,
) -> SensitivityResult:
    """
    The baseline capacity() and one SensitivityEntry per candidate edge. Candidates: with unrestricted
    Global and the optimized origin, the edges that carry flow in the baseline (removing any other edge cannot
    lower a maximum flow); otherwise every allowed edge, whatever its capacity (an edge of capacity 0
    can carry a route, so deleting it can change the result, and with the equal origin a deletion that
    removes a source from the roster can raise the total: a negative loss). loss = baseline.placed
    minus the placed of capacity() without the edge, in signed units, never thresholded. Two
    shortcuts give exactly that value: unrestricted optimized Global repairs a copy of the baseline's
    solve instead of solving again, and CostOnly
    gives loss 0 without placing to an edge outside the static shortest-path support reachable from
    the roster sources (no route can use it). Global with a support solves every deletion afresh; Pinned
    replays the same program under each deletion. workers as batch_capacity
    (NGRAPH_CORE_SENSITIVITY_THREADS). Only CostOnly reads the cache, and each placed deletion uses a
    new Eligibility uid, so a frozen cache raises DagCacheMiss there.
    A failed deletion fails the call (the lowest failing candidate's exception, no partial result); a
    supplied DagCache then behaves as the DagCache docstring states.
    """

def spf(
    graph: Graph,
    roots: typing.SupportsIndex
    | collections.abc.Sequence[typing.SupportsIndex]
    | numpy.typing.NDArray[numpy.integer],
    opts: SpfOptions,
) -> SpfResult:
    """
    Shortest paths from the roots for FORWARD, or to the roots for REVERSE (one NodeId or a sequence;
    every root at distance 0), over the eligible edges (only those with residual >= min_residual when
    a residual is given; an edge of capacity 0 is searched like any other) -> SpfResult. Without
    targets every reachable node gets its distance: costs and the full DAG of tight edges (the edges
    on some shortest walk; one per (u, v) pair with multi_edge = False). With targets the search
    stops at the nearest targets' cost:
    target_cost and the DAG of the walks between the roots and a nearest target (for REVERSE, from
    the target to a root, in the edges' own direction);
    with multipath = False, instead the canonical single path, the walk minimising (cost, -bottleneck
    under PREFER_HIGHER_PRIORITY, hops, EdgeId sequence). Validated in this order before any work
    (ValueError): the roots (non-empty, ids in range), the targets (the same), a node in both sets,
    multipath = False without targets, the residual (length, then every entry in [0, capacity_units]),
    min_residual >= 1, a foreign Eligibility. Excluded roots and targets are skipped; all roots
    excluded gives an empty result.
    """

@typing.overload
def to_units(
    x: float | numpy.floating | numpy.integer,
    resolution: typing.SupportsFloat | typing.SupportsIndex,
) -> int:
    """
    Physical value(s) to units of the resolution q, rounding half up: an int for a scalar, an
    int64 array for an array. NaN, infinity, negative values and results at or above 2^62
    units raise ValueError; q must be a power of two in [2^-32, 2^32].
    """

@typing.overload
def to_units(
    x: collections.abc.Sequence[float | numpy.floating | numpy.integer]
    | numpy.typing.NDArray[numpy.floating | numpy.integer],
    resolution: typing.SupportsFloat | typing.SupportsIndex,
) -> numpy.typing.NDArray[numpy.int64]: ...

UNREACHABLE: int = 9223372036854775807
