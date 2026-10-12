"""The one graph type in layout B' (APPROACH revision 24, sections 3 and 4).

Construction arrays are in input coordinates. The public ``EdgeId`` is the rank of
the edge in the stable order ``(src, cost, dst, input)``; ``input_index()[e]`` is
the input position of ``e``. Every other per-edge array and id is in id
coordinates. The inverse permutation is the caller's (``input_to_id``).
"""

from __future__ import annotations

import itertools
from dataclasses import dataclass
from fractions import Fraction

from units import LIMIT, ContractError, check_resolution

_UID = itertools.count(1)


class Graph:
    """StrictMultiDiGraph(num_nodes, src, dst, capacity_units, cost, resolution)."""

    def __init__(self, num_nodes, src, dst, capacity_units, cost, resolution):
        q = check_resolution(resolution)
        if not isinstance(num_nodes, int) or isinstance(num_nodes, bool):
            raise ContractError("num_nodes must be an int")
        if num_nodes < 0:
            raise ContractError("num_nodes must be non-negative")
        e = len(src)
        if not (len(dst) == len(capacity_units) == len(cost) == e):
            raise ContractError("edge arrays differ in length")
        for name, arr in (("src", src), ("dst", dst)):
            for x in arr:
                if not isinstance(x, int) or isinstance(x, bool):
                    raise ContractError(f"{name} entries must be ints")
                if not 0 <= x < num_nodes:
                    raise ContractError(f"{name} node {x} out of range")
        for c in capacity_units:
            if not isinstance(c, int) or isinstance(c, bool):
                raise ContractError("capacity_units must be integer units")
            if not 0 <= c < LIMIT:
                raise ContractError(f"capacity {c} outside [0, 2^62)")
        if sum(capacity_units) >= LIMIT:
            raise ContractError("sum of capacities must be below 2^62 units")
        for w in cost:
            if not isinstance(w, int) or isinstance(w, bool):
                raise ContractError("costs must be integers")
            if w < 1:
                raise ContractError(f"cost {w} rejected: costs are positive (rule 28)")
        if sum(cost) >= LIMIT:
            raise ContractError("sum of costs must be below 2^62")
        order = sorted(range(e), key=lambda i: (src[i], cost[i], dst[i], i))
        self.uid = next(_UID)
        self.q: Fraction = q
        self.num_nodes = num_nodes
        self.num_edges = e
        self._input_index = tuple(order)
        self.src = tuple(src[i] for i in order)
        self.dst = tuple(dst[i] for i in order)
        self.capacity = tuple(capacity_units[i] for i in order)
        self.cost = tuple(cost[i] for i in order)
        offsets = [0] * (num_nodes + 1)
        for u in self.src:
            offsets[u + 1] += 1
        for u in range(num_nodes):
            offsets[u + 1] += offsets[u]
        self.row_offsets = tuple(offsets)
        # Reverse CSR in today's (cost, src, input) order; slots map to ids.
        in_lists: list[list[int]] = [[] for _ in range(num_nodes)]
        for eid in sorted(
            range(e), key=lambda i: (self.cost[i], self.src[i], self._input_index[i])
        ):
            in_lists[self.dst[eid]].append(eid)
        in_off = [0]
        slots: list[int] = []
        for v in range(num_nodes):
            slots.extend(in_lists[v])
            in_off.append(len(slots))
        self.in_row_offsets = tuple(in_off)
        self.in_adj_edge_index = tuple(slots)
        self.max_cost = max(self.cost, default=0)

    # ------------------------------------------------------------ accessors

    def input_index(self) -> tuple[int, ...]:
        return self._input_index

    def resolution(self) -> Fraction:
        return self.q

    def out_ids(self, u: int) -> range:
        return range(self.row_offsets[u], self.row_offsets[u + 1])

    def in_ids(self, v: int) -> tuple[int, ...]:
        """In-edge ids of ``v`` in reverse-CSR order ``(cost, src, input)``."""
        return self.in_adj_edge_index[
            self.in_row_offsets[v] : self.in_row_offsets[v + 1]
        ]

    def eligibility(self, excluded_nodes=(), excluded_edges=()) -> Eligibility:
        return Eligibility.build(self, excluded_nodes, excluded_edges)

    def edge_tuple(self, e: int) -> tuple[int, int, int, int]:
        return (self.src[e], self.dst[e], self.capacity[e], self.cost[e])


def input_to_id(g: Graph) -> list[int]:
    """The caller's inverse permutation of ``input_index()`` (one O(E) scatter)."""
    inv = [0] * g.num_edges
    for eid, pos in enumerate(g.input_index()):
        inv[pos] = eid
    return inv


def graph_from_input(num_nodes, edges, resolution=1) -> Graph:
    """Convenience: ``edges`` as ``(u, v, capacity_units, cost)`` in input order."""
    return Graph(
        num_nodes,
        [e[0] for e in edges],
        [e[1] for e in edges],
        [e[2] for e in edges],
        [e[3] for e in edges],
        resolution,
    )


@dataclass(frozen=True)
class Eligibility:
    """Immutable static eligibility: dense node and edge bits, graph-bound."""

    graph_uid: int
    uid: int
    node_allowed: tuple[bool, ...]
    edge_allowed: tuple[bool, ...]

    @staticmethod
    def build(g: Graph, excluded_nodes=(), excluded_edges=()) -> Eligibility:
        xn = set()
        for v in excluded_nodes:
            if not isinstance(v, int) or not 0 <= v < g.num_nodes:
                raise ContractError(f"excluded node {v} out of range")
            xn.add(v)
        xe = set()
        for e in excluded_edges:
            if not isinstance(e, int) or not 0 <= e < g.num_edges:
                raise ContractError(f"excluded edge {e} out of range")
            xe.add(e)
        node_allowed = tuple(v not in xn for v in range(g.num_nodes))
        edge_allowed = tuple(
            e not in xe and node_allowed[g.src[e]] and node_allowed[g.dst[e]]
            for e in range(g.num_edges)
        )
        uid = 0 if not xn and not xe else next(_UID)
        return Eligibility(g.uid, uid, node_allowed, edge_allowed)

    def check(self, g: Graph) -> None:
        if self.graph_uid != g.uid:
            raise ContractError("foreign Eligibility")


def all_allowed(g: Graph, elig: Eligibility | None) -> Eligibility:
    if elig is None:
        return Eligibility.build(g)
    elig.check(g)
    return elig


def check_node_set(g: Graph, nodes, name: str) -> tuple[int, ...]:
    """Normalise a root or endpoint set: duplicates removed, ids validated."""
    if isinstance(nodes, int):
        nodes = (nodes,)
    out = set()
    for v in nodes:
        if not isinstance(v, int) or isinstance(v, bool) or not 0 <= v < g.num_nodes:
            raise ContractError(f"{name}: node {v!r} out of range")
        out.add(v)
    if not out:
        raise ContractError(f"{name}: empty set")
    return tuple(sorted(out))


def check_disjoint(sources, destinations) -> None:
    if set(sources) & set(destinations):
        raise ContractError("a vertex in both endpoint sets is rejected (rule 13)")


def check_dag_edges(g: Graph, edges) -> tuple[int, ...]:
    """A DAG given by its edge ids: range checked, no duplicates, acyclic."""
    seen = set()
    for e in edges:
        if not isinstance(e, int) or not 0 <= e < g.num_edges:
            raise ContractError(f"DAG edge {e!r} out of range")
        if e in seen:
            raise ContractError(f"duplicate DAG edge {e}")
        seen.add(e)
    out_adj: dict[int, list[int]] = {}
    indeg: dict[int, int] = {}
    for e in seen:
        u, v = g.src[e], g.dst[e]
        out_adj.setdefault(u, []).append(v)
        indeg[v] = indeg.get(v, 0) + 1
        indeg.setdefault(u, 0)
    queue = [v for v, d in indeg.items() if d == 0]
    count = 0
    while queue:
        u = queue.pop()
        count += 1
        for v in out_adj.get(u, ()):
            indeg[v] -= 1
            if indeg[v] == 0:
                queue.append(v)
    if count != len(indeg):
        raise ContractError("route structure is cyclic")
    return tuple(sorted(seen))


def path_edges(g: Graph, hops) -> tuple[int, ...]:
    """``Path(g, hops)`` validation; returns the member edge ids (a DAG)."""
    if not hops:
        raise ContractError("empty path")
    members: list[int] = []
    nodes: list[int] = []
    for hop in hops:
        if not hop:
            raise ContractError("empty hop")
        for e in hop:
            if not isinstance(e, int) or not 0 <= e < g.num_edges:
                raise ContractError(f"path member {e!r} out of range")
        uv = {(g.src[e], g.dst[e]) for e in hop}
        if len(uv) != 1:
            raise ContractError("members of a hop must share (u, v)")
        (u, v) = next(iter(uv))
        if nodes and nodes[-1] != u:
            raise ContractError("hops are not contiguous")
        if not nodes:
            nodes.append(u)
        nodes.append(v)
        members.extend(hop)
    if len(set(members)) != len(members):
        raise ContractError("duplicate path member")
    if len(set(nodes)) != len(nodes):
        raise ContractError("path repeats a node")
    return tuple(sorted(members))
