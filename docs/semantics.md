# NetGraph-Core semantics

This document says exactly what NetGraph-Core computes: how amounts are counted, how edges are
numbered, which route a search returns, how each routing model spreads traffic, what a call reports
and when it fails, how the shortest-path cache behaves, and what the evaluators mean. It is the full
reference behind the API docstrings. For installation and a first example, see the
[README](../README.md#installation).

Each section starts with a short explanation in plain words and then gives the exact rules. When
the two seem to differ, the exact rules are the authority. Names are the Python spellings
(`import netgraph_core as ngc`); [the last section](#11-c-names) lists the C++ spellings. Names that
start with an underscore are private test and verification entries, not part of the API.

Every Python example below is run by the test suite, in order and in one namespace, and must print
exactly the output block that follows it.

Contents:

- [How the engine works](#how-the-engine-works)
- [Terms](#terms)

1. [Units: the one numeric domain](#1-units-the-one-numeric-domain)
2. [The graph](#2-the-graph)
3. [Route values: PredecessorDAG and Path](#3-route-values-predecessordag-and-path)
4. [Search](#4-search)
5. [Route allocations and admission](#5-route-allocations-and-admission)
6. [The ledger](#6-the-ledger)
7. [Policies](#7-policies)
8. [Failures, threads and lifetimes](#8-failures-threads-and-lifetimes)
9. [The shortest-path cache](#9-the-shortest-path-cache)
10. [Evaluators](#10-evaluators)
11. [C++ names](#11-c-names)

## How the engine works

A demand is traffic from a set of source nodes to a set of destination nodes. The engine follows it
through six things:

1. A **graph** (`Graph`): nodes and directed edges, each edge with an integer capacity and a
   positive integer cost. Capacities are counted in whole **units**; the graph fixes how big a unit
   is.
2. A **ledger** (`FlowLedger`): what the edges of one graph carry now, and the bandwidth that
   tunnels hold on them.
3. A **policy** (`FlowPolicy`): one demand on one ledger, with a configuration: a **routing model**
   (how routes are chosen) and an **admission** (what to do with traffic that does not fit).
4. A **program** (`ForwardingProgram`): the routes the model chose and the integer weights that say
   how traffic is split over them.
5. **Forwarding**: the call's offer is split over the program's routes in whole units. `LOSSLESS`
   admits only what goes through without loss; `LOSSY` sends everything and records where it is
   lost.
6. A **result** (`PlacementResult`): what the call delivered, what was lost and where, the program,
   and the demand's totals.

One demand, routed like IGP routers with equal-cost multipath (ECMP), on a network of two branches:

```python
import netgraph_core as ngc

# Node 0 reaches node 3 over 0 -> 1 -> 3 (capacity 10) and 0 -> 2 -> 3 (capacity 2), every cost 1.
g = ngc.Graph(4, [0, 1, 0, 2], [1, 3, 2, 3], [10, 10, 2, 2], [1, 1, 1, 1], 1.0)
ledger = ngc.FlowLedger(g)
ecmp = ngc.FlowPolicyConfig(routing=ngc.ShortestPaths(), admission=ngc.Admission.LOSSY)
policy = ngc.FlowPolicy(ledger, 0, 3, ecmp)
r = policy.place(10)
print(r.call_delivered, r.call_edge_drops, r.termination.name)
print(ledger.carried().tolist(), ledger.residual().tolist())
```

```text
7 [(1, 3)] COMPLETE
[5, 2, 5, 2] [5, 0, 5, 0]
```

The edges are numbered by source, cost and destination, so 0 -> 1 is edge 0, 0 -> 2 edge 1,
1 -> 3 edge 2 and 2 -> 3 edge 3. Node 0 splits the 10 units equally over its two next hops, 5 and 5.
Edge 1 carries 2 of its 5 and drops 3, so 7 arrive. `LOSSY` sends the whole offer, so the call is
`COMPLETE` although 3 units were lost; with `LOSSLESS` the same policy would admit only 5, the
largest offer whose equal split overloads no edge.

## Terms

Each term in one sentence, in the order the later terms need them; the section named after it
gives the exact rule.

- **Unit, resolution `q`.** Traffic is counted in whole units; one unit is `q` of a physical
  quantity (for example `q = 2**-12` Gbit/s). The graph fixes `q` when it is built
  ([section 1](#1-units-the-one-numeric-domain)).
- **Amount.** A whole number of units: a capacity, a residual, an offer, a delivery or a drop.
- **EdgeId.** The public number of an edge: its rank when the edges are sorted by source, cost and
  destination. It is not the edge's position in the input arrays ([section 2](#2-the-graph)).
- **Eligibility (a mask).** Which nodes and edges of a graph a call may use. A failed link or router
  is modelled by an Eligibility that excludes it ([Eligibility](#eligibility)).
- **Static.** Depending only on the graph, its costs and the Eligibility, not on traffic.
- **Residual.** The free capacity of an edge: in a ledger, its service baseline minus the traffic
  it carries; in a search, the array the caller supplies.
- **Tight edge.** An edge `u -> v` that lies on a shortest walk of a search. Searching from the
  roots, the distance of `u` plus the edge's cost equals the distance of `v`; searching toward
  them (reverse), where a distance is the distance to the roots, the distance of `v` plus the
  cost equals the distance of `u`.
- **Shortest-path DAG.** The set of all tight edges: every shortest route at once. A DAG (directed
  acyclic graph) is a directed graph with no directed cycle.
- **Parallel members.** Several edges between the same two nodes in the same direction.
- **Route choice.** The tie rule (`RouteTieBreak`), `multipath` and `parallel_edges` of a model or
  a search ([section 4](#the-route-choice)): cost always comes first, and the tie rule only
  orders equal-cost choices.
- **Member graph.** The edges a search may use: eligible (and, when a residual is supplied, at
  least the search's residual threshold), whatever their capacity, and one member per pair with
  `parallel_edges=False`.
- **Canonical path.** The one walk a search returns when a single path is asked for, chosen by
  the fixed order of [section 4](#the-canonical-single-path).
- **Bottleneck of a walk.** Under a `HIGHER_*` tie rule, the smallest edge score along the walk:
  static capacity, availability, or their exact ratio (headroom), as the rule names
  ([the route choice](#the-route-choice)).
- **Demand.** The traffic from a set of sources to a set of destinations; a route ends at the first
  destination it reaches.
- **Policy.** A `FlowPolicy`: one demand on one ledger, with its configuration. It computes the
  demand's routes, forwards its traffic over them and removes it again.
- **Routing model.** `ShortestPaths`, `DistributedTE` (T1), `CentralizedTE` (with a tunnel count:
  T2) or `Pinned`; each holds the options it reads ([configuration](#configuration)).
- **Roster.** The sources of a policy that can reach an allowed destination over allowed nodes and
  edges, whatever the capacities, in ascending node id; fixed when the policy is built
  ([the source roster](#the-source-roster)).
- **Head.** A roster source, as the start of the routes it sends on.
- **Source allocation.** How a demand's total enters at its roster sources: in equal round-robin
  shares (`SourceAllocation.EQUAL`), highest-averages shares of declared weights (`WEIGHTED`),
  or as a controller chooses (`OPTIMIZED`).
- **Route allocation.** How a demand's traffic is split at every router of its routes
  (`RouteAllocation`): equally (`EQUAL`, ECMP), by integer weights (`WEIGHTED`, WCMP), or as a
  controller restricted to those routes chooses (`OPTIMIZED`).
- **Round robin `RR(x, k, i)`.** `x` units over `k` members: `x // k` each, and one more for the
  first `x % k` members.
- **Highest averages `H(x, w)`.** `x` units over members with integer weights `w`, given one at a
  time to the member with the largest `w[i] / (assigned[i] + 1)`, ties to the earlier member; a
  zero weight gets nothing.
- **Admission.** `LOSSLESS` admits the largest total that the program, splitting by its own rules,
  carries with no edge overloaded and no share left without a route; the rest is refused before it
  enters. `LOSSY` sends the whole offer and records where it is lost ([section
  5](#5-route-allocations-and-admission)).
- **Offer.** The traffic a call asks a policy to carry (`place(offer)`, `replace(offer)`, ...).
- **Controller.** A routing that solves a flow problem for the whole demand: `CentralizedTE`
  without a tunnel count, and `ShortestPaths` with `allocation=OPTIMIZED`.
- **Support.** The DAG a `CentralizedTE` policy may be restricted to (`support=`); every other edge
  then carries nothing.
- **Service view.** The ledger's residual plus this policy's own carried traffic: the capacity a
  solve or a full-offer call of the policy may use.
- **Setup view (T1).** The ledger's unreserved capacity plus this policy's own holds: what T1's
  tunnel setup may reserve. The service view and the setup view are never added together.
- **Program.** A `ForwardingProgram`: the routes and integer weights of one demand and the rule
  that injects its traffic at the sources. Every computing model produces one; `Pinned` forwards a
  given one without computing anything ([forwarding programs](#forwarding-programs)).
- **Entry.** One route of a program: a **tunnel** (one path, one edge per hop) or a **split
  route** (a DAG whose routers split by fixed rules). Each entry has an integer identifier, its
  `entry_id`, and a weight: its share among its head's entries.
- **Split context.** The traffic a router splits as one amount. A `ShortestPaths` `EQUAL` or
  `WEIGHTED` program and a `from_dag` program add up all their sources' traffic at a shared router
  and split it once; each split route given to `from_routes` splits only its own share. The two can
  round differently ([example](#supplied-dags-from_dag)).
- **T1.** `DistributedTE(TunnelCount(N))`: every head sets up its own N strict-bandwidth tunnels by
  CSPF (constrained shortest path first: a shortest path among the links with enough unreserved
  bandwidth), with no controller ([T1](#t1-tunnels-set-up-by-their-heads)).
- **Hold.** Bandwidth an up T1 tunnel reserves on every link of its path. Holds are the ledger's
  second projection (`reserved`, `unreserved`), kept apart from the carried traffic (`carried`,
  `residual`); they never police or queue traffic.
- **T2.** `CentralizedTE(source_allocation=..., tunnel_count=TunnelCount(N))`: a central
  controller lays out N single-path tunnels for every live head with one common whole-unit planned
  bandwidth ([T2](#t2-the-controllers-equal-tunnels)).
- **Live head, dead head (T2).** A head is live when it reaches an allowed destination over links
  with at least 1 unit of free capacity; a dead head keeps its share of the traffic and loses it.
- **Design demand.** The volume T1 and T2 size their tunnels for (`design_demand`), separately
  from the offer they forward. T2's planned bandwidth is metadata, never a hold.
- **Epoch.** The period from a policy's construction, a `remove()`, a `replace()` or a `forward()`
  to the next of these; the offers of `place` add up within it.
- **Black hole.** A node or edge that fails silently: forwarding goes on unchanged and traffic
  that reaches it disappears there (stale forwarding).
- **First-hop experiment.** A `LOSSY` `place(float('inf'))` of `ShortestPaths`: each source puts
  exactly each first link's free capacity on it, which measures how much the sources can push into
  their first links at once ([section 5](#5-route-allocations-and-admission)).
- **Record.** A `FlowRecord`: one routed amount of a policy in the ledger, kept per roster source
  or per program entry ([section 6](#6-the-ledger)).
- **Gate.** The single native lock of a ledger: calls on one ledger take effect one at a time
  ([section 8](#8-failures-threads-and-lifetimes)).
- **Lease.** A call's admission to a `ShortestPathCache`: exclusive while the cache is building,
  shared while it is frozen ([section 9](#9-the-shortest-path-cache)).
- **Minimum cut.** For the optimized controller: the edges that leave the nodes the sources can
  still reach once its solve is maximal, the minimum cut closest to the sources
  ([the result](#the-result)).
- **Marginal cost curve.** `(cost, units)` pairs in ascending cost: delivering `units` more units
  raises the minimum total cost of the flow by `cost` per unit. Delivering more can reroute flow
  already placed, so the pairs are not the costs of the paths in the final allocation.

The examples from here on share these two helpers:

```python
import netgraph_core as ngc

def graph(num_nodes, edges, q=1.0):
    """edges: (source, destination, capacity units, cost) tuples in input order."""
    sources, destinations, capacities, costs = zip(*edges)
    return ngc.Graph(num_nodes, list(sources), list(destinations), list(capacities), list(costs), q)

def config(routing, admission=ngc.Admission.LOSSLESS):
    """A configuration: the routing model and the admission."""
    return ngc.FlowPolicyConfig(routing=routing, admission=admission)
```

## 1. Units: the one numeric domain

**In plain words.** Every amount of traffic is a whole number of units, and the graph says how
big a unit is: its resolution `q`. Physical numbers enter through `to_units` and leave through
`as_float`; nothing in between is a float, so totals add up exactly and nothing is lost to rounding.
Costs are a separate kind of number: positive integers that are never converted.

`to_units(x, q)` works in two steps. It first reads `x` as the binary64 value `b` it holds: a float
is its own value, and an integer becomes the nearest double. It then rounds `b / q` to a whole
number of units, half up, exactly. Converted back, the result differs from `b` by at most `q / 2`,
and a tie (exactly half a unit) goes up. Reading an integer above `2^53` as `b` can add an error of
its own.

At `q = 1/4` (a quarter of a Gbit/s, say):

```python
print(ngc.to_units([2.6, 0.375, 0.125, 0.1], 0.25).tolist())
print(ngc.as_float([10, 2], 0.25).tolist())
```

```text
[10, 2, 1, 0]
[2.5, 0.5]
```

2.6 is about 10.4 units and rounds down to 10, which is 2.5 again. 0.375 is exactly 1.5 units, a
tie, and goes up to 2, which is 0.5: more than the input. 0.125 is the tie 0.5 and gives 1, and 0.1
(about 0.4 units) gives 0. A link of 2.6 therefore has a capacity of 2.5 in this graph.

### Splitting whole units

When `I` units are split over `k` members, every member gets whole units. The members have a
fixed order, which decides every remainder and every tie: a router's next-hop edges in ascending
EdgeId, a demand's sources in ascending node id, and a head's entries in the program's entry
order. There are two rules:

- **Round robin** (`EQUAL`): every member gets `I // k` units, and the first `I % k` members in
  that order get one more. `RR(10, 3) = (4, 3, 3)`.
- **Highest averages** (`WEIGHTED`, the D'Hondt rule): the units are handed out one at a time,
  each to the member with the largest `w / (a + 1)` (its weight `w` over its current count `a`
  plus one), ties to the earlier member in that order. Equal positive weights give round robin.

Highest averages with weights 2 and 1, unit by unit:

| Unit | `w / (a + 1)`: member 0, member 1 | Goes to | Counts after |
| --- | --- | --- | --- |
| 1 | 2/1, 1/1 | 0 | 1, 0 |
| 2 | 2/2, 1/1 (a tie) | 0 | 2, 0 |
| 3 | 2/3, 1/1 | 1 | 2, 1 |
| 4 | 2/3, 1/2 | 0 | 3, 1 |
| 5 | 2/4, 1/2 (a tie) | 0 | 4, 1 |
| 6 | 2/5, 1/2 | 1 | 4, 2 |

So `H(5, [2, 1]) = (4, 1)`, not the proportional 10/3 and 5/3, and `H(6, [2, 1]) = (4, 2)`: highest
averages gives a remainder to the larger weight first.

```python
# Three parallel links of 10 share 10 units by round robin; two links weighted 2:1 share 5, then 6.
g = graph(2, [(0, 1, 10, 1)] * 3)
ledger = ngc.FlowLedger(g)
ngc.FlowPolicy(ledger, 0, 1, config(ngc.ShortestPaths())).place(10)
print("RR(10, 3):", ledger.carried().tolist())
g = graph(2, [(0, 1, 100, 1), (0, 1, 100, 1)])
wcmp = ngc.ShortestPaths(allocation=ngc.RouteAllocation.WEIGHTED,
                         weight_basis=ngc.WeightBasis.EXPLICIT, edge_weights=[2, 1])
for x in (5, 6):
    ledger = ngc.FlowLedger(g)
    ngc.FlowPolicy(ledger, 0, 1, config(wcmp)).place(x)
    print(f"H({x}, [2, 1]):", ledger.carried().tolist())
```

```text
RR(10, 3): [4, 3, 3]
H(5, [2, 1]): [4, 1]
H(6, [2, 1]): [4, 2]
```

**Exact rules.**

- Both rules are exact: every unit goes to some member, unless every weight is 0, and then no
  unit is assigned (a proportional share is undefined). A share that rounds down to zero leaves
  its unit to a sibling.
- Both are monotone: when `I` grows, no member's share shrinks.
- **The local rounding bound.** For one split into `k` members, each member's absolute rounding
  error, relative to that split's actual integer input and effective weights, is below `k` units
  (`k * q` physically), and below one unit for equal positive weights. This is a local bound, not
  an end-to-end delivery bound.
- Highest averages is not the closest proportional rounding: weights `[100, 1, 1, 1]` and 50 units
  give `[50, 0, 0, 0]` (total absolute error `300/103` against the exact shares
  `[5000/103, 50/103, 50/103, 50/103]`, where `[49, 1, 0, 0]` would have `200/103`).

### The conversion exactly

- **Resolution.** A graph's resolution `q` is a power of two in `[2^-32, 2^32]`, given
  explicitly when the graph is built; there is no default. The resolution is the precision of
  every result: there is no finer grid. A model that needs finer splits uses a smaller `q`.
- **Amounts.** Capacities, residuals, allocations, requests, deliveries, offers and drops are
  integers in units of `q`. In Python they are `int` (exact, unbounded); per-edge arrays are
  NumPy `int64[E]`. Cumulative counters that are sums over many calls (a Lossy record's
  `offered` and its per-edge `edge_drops`) may exceed 64 bits; Python returns them as exact `int`.
  Cumulative counters must remain below `2^127`. A call that would reach that limit raises
  `ValueError('<counter> limit exceeded')` and publishes nothing; the cumulative request is
  checked before computation.
- **Domain.** Every edge capacity, and the sum of all capacities, is below `2^62` units. Edge
  costs are integers of at least 1 and their sum is below `2^62`, so no path cost reaches the
  sentinel `UNREACHABLE` (the largest `int64`). An offer is an integer in `[0, 2^62)` units or
  positive infinity (`float('inf')`, "as much as possible"). Everything outside the domain is a
  `ValueError` raised before anything changes. The limits on the sizes of a call are in
  [section 8](#limits-that-refuse-a-call).
- **`to_units(x, q)`** reads each physical input as its binary64 value `b`: Python floats and
  NumPy float64 keep their exact value, float32 widens exactly, and integers become the nearest
  double, ties to even (`2^53 + 1` becomes `2^53`, `2^53 + 3` becomes `2^53 + 4`). It returns
  `floor(b/q + 1/2)`, evaluated exactly on rational numbers, not by floating addition of `1/2`.
  A scalar or 0-d NumPy array gives an `int`; a Python sequence or 1-d NumPy array gives an
  `int64` array. Negative zero gives 0. NaN, infinity, negative values and results at or above
  `2^62` units raise `ValueError`. This is the only entry for physical numbers.
- **Ties.** The ties are exactly the accepted `b = (2*n + 1)*q/2` for integers `0 <= n < 2^52`
  (at and above `2^52 * q`, `b/q` is an integer). Every tie gives `n + 1`, never the even
  neighbour. At `q = 1`, `0.5` gives 1, `2.5` gives 3, `3.5` gives 4, and the last tie
  `2^52 - 1/2` gives `2^52`.
- **Binary values and precision.** A float is never read through its printed decimal; conversion
  cannot recover a source literal. `2^53` is the limit of consecutive integers among doubles,
  not a size below which decimal fractions are exact. A decimal is exactly a double only if
  its reduced denominator is a power of two and it fits binary64; even `0.1` is rounded on
  parsing. A decimal and its parsed double give the same units exactly when both lie in the
  same cell `[(n - 1/2)*q, (n + 1/2)*q)`. No magnitude or unit count guarantees this: a parse
  error can cross or land on a boundary at an ordinary size. At `q = 2^-21`, `5501195.1796`
  denotes the exact binary tie `11536842473288.5` units, giving `11536842473289`; the typed
  decimal would be `11536842473288.4992` units. At `q = 2^-18`, the computed `400 + 2^-19`
  is the tie `104857600.5` and gives `104857601`, although its printed decimal
  `400.00000190734863` is just below the tie. `np.float32(0.1)` is exactly
  `13421773/134217728`, giving `429496736` units at `q = 2^-32`; decimal `0.1` would give
  `429496730`. `CostStretch` deliberately uses the printed decimal instead: a stretch is a
  ratio a person writes (`1.2` means `6/5`), compared exactly with ratios of integer costs.
  A physical amount is the measured or computed double the caller holds.
- **`as_float(u, q)`** returns the correctly rounded binary64 value of the exact product `u*q`
  (nearest, ties to even). Here `u` may be any signed, unbounded integer, unlike graph and
  flow amounts. An integer scalar gives a `float`, an integer vector a `float64` array.
  An int or NumPy integer, a subclass included, is never read as a sequence.
  The result is finite exactly when `abs(u*q) < 2^1024 - 2^970`; otherwise it raises
  `ValueError('as_float: <u> units times the resolution exceed the largest finite double')`.
  Every `abs(u) <= 2^53`, including the endpoint, converts exactly. Beyond `2^53`, an integer
  of binary length `L` converts exactly iff it is divisible by `2^(L - 53)`; `2^53 + 1`
  gives `2^53*q`, and `2^53 + 3` gives `(2^53 + 4)*q`. Costs are never converted.
- **Round trip.** `to_units(as_float(u, q), q) == u` for every `u` in `[0, 2^53]` at every
  resolution, and for larger `u` in `[0, 2^62)` exactly when `u` is a double. Otherwise the
  result is the units of the rounded value. The first failure, `2^53 + 1`, comes back as
  `2^53`; `2^62 - 1` becomes `2^62*q`, which `to_units` refuses.
- **The error against `b`.** Every successful `to_units` result `u` is itself a double, so
  `as_float` adds no error: `e = as_float(to_units(b, q), q) - b` satisfies `-q/2 < e <= q/2`,
  with equality at the upper end exactly at ties. Against an argument `a` changed by promotion (a
  large integer, or a wider NumPy float where the platform has one), the error is `e + (b - a)`:
  `2^53 + 1` at `q = 1` has error 0 against `b`, but -1 against the integer passed. There is no
  exact-integer path for physical inputs.
- **Sums.** Rounding `N` values separately gives a total error in `(-N*q/2, N*q/2]`. Ties
  contribute `+q/2` each; other errors go either way, so a sum has no promised sign. At `q = 1`,
  `[0.5, 0.5]` gives `[1, 1]` (2 units for a total of 1), while `[0.25, 0.25]` gives `[0, 0]`
  (0 units for a total of 0.5). When non-tie fractional parts are symmetric about `1/2`, their
  mean error cancels and the mean error per value is `(q/2)*P(tie)`; a continuous uniform
  fractional part has none. The Core assumes no distribution. Round a demand once, then split
  whole units to preserve its total: rounding `1.0` once gives 1 unit.
- **Choosing `q`.** Choose one `q` before graph construction and use the constructed graph's
  `resolution()` afterwards. Choosing `q` and reporting rounding belong to the caller (NetGraph's
  units report); the Core converts and validates, without a selector or rounding report.
- **Routability.** An edge of capacity 0 carries nothing, but it stays in routing: searches,
  routes and rosters treat it like any other edge. Zero capacity alone never excludes an edge: an
  Eligibility excludes it, the route choice can keep another of its parallel members
  ([section 4](#the-route-choice)), and a supplied residual adds a capacity threshold. An equal
  split still sends such an edge its share, which it cannot carry
  ([section 5](#5-route-allocations-and-admission)). A physical offer below `q / 2` converts to
  0 units. A `place(0)` (in a model that has `place`) is accepted after validation and does
  nothing; its result is `COMPLETE`. A full zero offer (`evaluate(0)`, `replace(0)`, `forward(0)`)
  is a whole computation or forwarding of nothing, and `replace(0)` publishes it
  ([operations](#operations-evaluate-replace-and-forward)). Neither is a removal: `remove()` is.
- **Conservation.** At every node, the traffic that arrives plus what the node injects equals what
  it sends on plus what it delivers plus what it loses there, exactly; an edge carries what it is
  sent minus what it drops. The ledger's totals are exactly the sums of its records.

```python
print(ngc.to_units([1.5, 2.5, 0.4, 0.5], 1.0).tolist())   # half up
print(ngc.to_units(10.0, 2.0**-12), ngc.as_float(40960, 2.0**-12))
```

```text
[2, 3, 0, 1]
40960 10.0
```

## 2. The graph

**In plain words.** A `Graph` is built once from arrays and never changes. Every edge gets a
public number, its EdgeId, by sorting the edges by source, cost and destination; only edges
that agree on all three keep their input order among themselves. A map gives back each edge's
input position. All per-edge arrays the library reads or returns are in EdgeId order, and every
tie rule of the library ("the lowest EdgeId") refers to this order. A caller who keeps per-edge
values in input order (explicit weights, a residual, exclusions) maps them through
`input_positions()` first.

**Exact rules.**

- **Construction.** `Graph(num_nodes, edge_sources, edge_destinations, edge_capacity_units,
  edge_costs, resolution)` builds an immutable directed multigraph: parallel edges and self-loops
  are allowed. The arrays are in input order (the caller's). Capacities are already in units; use
  `to_units` for physical values.
- **Validation order.** When an input has several defects, the first in this order is reported (a
  wrong type as a `TypeError`, a value outside its domain as a `ValueError`): the resolution (a
  bool is a `TypeError`); the capacity values (each an integer in `[0, 2^62)`; floats and booleans
  are rejected with `TypeError`, not rounded); the array shapes; a zero or negative cost; the remaining
  bounds (node ids in range, then the capacity and cost sums below `2^62`, both added up edge by
  edge in input order: the first edge at which a running sum reaches `2^62` reports that sum, the
  capacity sum before the cost sum at the same edge). Integers are read
  without narrowing, so the reported error does not depend on how the same input is represented.
- **Size.** Node ids and EdgeIds are `int32`: a graph has at most `INT32_MAX` (`2^31 - 1`) nodes
  and at most as many edges. The constructor accepts every graph within these bounds, and the
  searches and route values work on all of them. A flow solve has slightly smaller limits
  ([section 8](#limits-that-refuse-a-call)).
- **EdgeId.** The EdgeId of an edge is its rank in the stable order `(source, cost, destination,
  input position)`. `input_positions()[e]` is the input position of edge `e`; the inverse map is
  the caller's to build (one pass over `E`).
- **EdgeId order is the coordinate of every per-edge value.** Residual arrays, explicit weights,
  Eligibility exclusions, DAG and Path edges, flows, records, `saturated_boundary_edges`, `min_cut`
  and the entries of `evaluate_edge_failures` all use EdgeIds. A caller that keeps per-edge data in
  input order translates it first: an in-range id in the wrong coordinate is silently a different
  edge, and no check can catch that.
- **Layout guarantees.** The out-edges of node `u` are exactly the ids
  `outgoing_offsets()[u]` to `outgoing_offsets()[u + 1] - 1`, in `(cost, destination, input)`
  order, and `edge_destinations()[e]` is the head of `e`. The in-edges of node `v` are
  `incoming_edges()[incoming_offsets()[v]:incoming_offsets()[v + 1]]`, as EdgeIds in
  `(cost, source, input)` order. Both offset arrays have `N + 1` entries, starting at 0 and ending
  at `E`, also for an empty graph.
- **Accessors.** `edge_capacity_units()`, `edge_costs()` (`int64[E]`), `edge_sources()`,
  `edge_destinations()`, `input_positions()`, `outgoing_offsets()`, `incoming_offsets()` and
  `incoming_edges()` (`int32`) are read-only NumPy views that keep the graph alive;
  `max_edge_cost()` and `edge_cost_gcd()` are 0 for a graph without edges; `num_nodes()`,
  `num_edges()` and `resolution()` read the sizes; `uid()` is a process-unique number, never
  reused.
- **Ids belong to one graph.** Building a graph from different arrays, or with another release of
  the library, may renumber the edges. Every value bound to a graph (Eligibility, PredecessorDAG,
  Path, program, ledger, cache) is rejected with another graph.

The three edges below, sorted by (source, cost, destination), are 0 -> 1 (cost 1), 0 -> 2 (cost 5)
and 1 -> 2: EdgeIds 0, 1 and 2, from input positions 2, 1 and 0.

```python
# Input order: edge 0 = 1->2, edge 1 = 0->2 (cost 5), edge 2 = 0->1.
g = ngc.Graph(3, [1, 0, 0], [2, 2, 1], [7, 8, 9], [1, 5, 1], 1.0)
print("sources", g.edge_sources().tolist(), "destinations", g.edge_destinations().tolist(), "costs", g.edge_costs().tolist())
print("input_positions", g.input_positions().tolist(), "outgoing_offsets", g.outgoing_offsets().tolist())
```

```text
sources [0, 0, 1] destinations [1, 2, 2] costs [1, 5, 1]
input_positions [2, 1, 0] outgoing_offsets [0, 2, 3, 3]
```

### Eligibility

**In plain words.** An `Eligibility` says which nodes and edges of one graph may be used, for
example to model a failure or a maintenance window. It is built from lists of exclusions and
never changes.

**Exact rules.**

- `graph.eligibility(excluded_nodes=(), excluded_edges=())` builds the value. Duplicates and
  unsorted lists are accepted; out-of-range ids raise `ValueError`. An edge is allowed when it is
  not excluded and both its endpoints are allowed. `node_allowed()` and `edge_allowed()` are
  read-only boolean views.
- Empty lists give the all-allowed value, whose `uid()` is 0. Every other value has a fresh uid;
  equal lists built twice give two different uids. Wherever an Eligibility is optional, `None`
  means everything is allowed.
- `Eligibility.intersect(a, b)` allows what both allow. When the result equals one operand, that
  operand is returned (same uid); otherwise the result is a new value with a new uid. Both
  operands must belong to the same graph.
- Uids matter for the [shortest-path cache](#9-the-shortest-path-cache): its key contains the
  Eligibility uid, so build each state once and reuse the value.
- An excluded endpoint takes no part in a search or a policy: an excluded root or target in a
  set is skipped, and an excluded single endpoint gives an empty result.

Excluding node 2 of the graph above also excludes the two edges that touch it:

```python
mask = g.eligibility(excluded_nodes=[2])
print(mask.node_allowed().tolist(), mask.edge_allowed().tolist(), mask.uid() != 0)
```

```text
[True, True, False] [True, False, False] True
```

## 3. Route values: PredecessorDAG and Path

**In plain words.** A route is either a DAG of predecessor edges (the shape a search returns:
possibly many equal-cost walks) or a `Path` (one explicit walk, each hop one edge or several
parallel edges). Both are immutable, bound to their graph and checked once when built.

**Exact rules.**

- **`PredecessorDAG`.** For node `v`, its entries are `parents()[k]` (the parent node) and
  `via_edges()[k]` (the edge from that parent) for `k` from `parent_offsets()[v]` to
  `parent_offsets()[v + 1] - 1`, in strictly ascending EdgeId. This is the one canonical form:
  two DAGs of the same graph are equal exactly when their entries are equal; DAGs of different
  graphs are never equal. `parents()` and `via_edges()` are read-only views; `parent_offsets()`
  is built on demand.
- **`PredecessorDAG.from_csr(graph, parent_offsets, parents, via_edges)`** validates, before
  reading any structure: `parent_offsets` has `N + 1` entries, starts at 0, never decreases and
  ends at `len(parents) == len(via_edges)`; then every id is in range, every via edge joins its
  parent to its
  node, no `(node, via)` pair repeats, the DAG has no cycle and every node's entries are in
  ascending via EdgeId. A wrong type raises `TypeError`; structural and range defects raise
  `ValueError`, in the order above.
- **`Path(graph, hops)`.** `hops` is a list of hops, each a non-empty list of EdgeIds that share
  one `(u, v)` pair (parallel members). Validation, in order: every member id in range, no
  member twice, at least one hop, no empty hop, the members of a hop share `(u, v)`, consecutive
  hops connect, no node is visited twice. The members keep the caller's order. `cost()` is the
  sum over the hops of each hop's cheapest member; `source()`, `destination()` and `hops()` read
  the walk;
  `to_dag()` gives the equivalent `PredecessorDAG`. A single-member walk is `[[e1], [e2], ...]`.
- **`resolve_to_paths(graph, dag, source, destination, split_parallel_edges=False,
  max_paths=None)`** lists the `source -> destination` walks of a DAG as `Path` values in
  increasing `(cost, hop count, edge-id
  sequence)`, where a hop costs its cheapest member and is represented in the key by its lowest
  member id. With `split_parallel_edges=False` the parallel members of a pair stay together in
  one hop (ascending); with `True` every member combination is a separate path, and the paths are
  sorted by each path's own `(cost, hop count, edge-id sequence)`. `max_paths` (at least 1) keeps
  the first paths of that order. `source == destination`, a DAG of another graph or
  `max_paths < 1` raise `ValueError`; an endpoint outside the DAG, or no walk, gives `[]`.


Two parallel edges 0 -> 1 of costs 1 and 3 (EdgeIds 0 and 1), then 1 -> 2 (EdgeId 2). A hop that
holds both parallel edges costs its cheapest member, so the path costs 1 + 1. The same route as a
DAG lists node 1's two entries and node 2's one; listing its walks with the parallel members
split gives the cost-2 and the cost-4 path:

```python
g = graph(3, [(0, 1, 5, 1), (0, 1, 5, 3), (1, 2, 5, 1)])
path = ngc.Path(g, [[0, 1], [2]])
dag = ngc.PredecessorDAG.from_csr(g, [0, 0, 2, 3], [0, 0, 1], [0, 1, 2])
print(path.cost(), path.to_dag() == dag, dag.via_edges().tolist())
print([(p.cost(), p.hops()) for p in ngc.resolve_to_paths(g, dag, 0, 2)],
      [(p.cost(), p.hops()) for p in ngc.resolve_to_paths(g, dag, 0, 2, split_parallel_edges=True)])
```

```text
2 True [0, 1, 2]
[(2, [[0, 1], [2]])] [(2, [[0], [2]]), (4, [[1], [2]])]
```

## 4. Search

**In plain words.** `spf` finds shortest paths from a set of roots, forward or in reverse,
either to every node or only to the nearest of a set of targets. It returns the shortest-path
DAG, which contains every shortest route at once, or, on request, one canonical path chosen by a
fixed rule. `ksp` lists the K shortest simple paths in a fixed order. Capacity alone never
removes an edge from a search: an eligible edge of capacity 0 is a route like any other. The route
choice can keep one of several parallel members (below), and only a supplied residual adds a
capacity threshold.

### The route choice

The searches and the models that choose routes (`ShortestPaths`, T1's CSPF) share one rule set.
Cost always comes first; the rest only orders equal-cost choices.

- **`tie_break: RouteTieBreak`**, one of four values: `FEWER_HOPS` (the default: cost, then fewer
  hops, then the smaller EdgeId sequence), `HIGHER_CAPACITY`, `HIGHER_AVAILABLE_CAPACITY` and
  `HIGHER_HEADROOM` (cost, then the larger bottleneck of the named edge score, then fewer hops,
  then the smaller EdgeId sequence).
- **Edge scores.** `HIGHER_CAPACITY` reads the static capacity; `HIGHER_AVAILABLE_CAPACITY` the
  availability of the current search; `HIGHER_HEADROOM` the availability divided by the static
  capacity, compared exactly as a fraction (an edge of capacity 0 scores 0). The availability
  depends on the caller:

  | Caller | Availability |
  | --- | --- |
  | `spf`, `ksp` with a `residual` | the supplied `residual` |
  | `spf`, `ksp` without one | the static capacity (`HIGHER_HEADROOM` then reads 1 on an edge of positive capacity and 0 on an edge of capacity 0) |
  | `ShortestPaths` | its service view when the program is prepared ([routing models](#routing-models)) |
  | T1's CSPF | its setup view, minus the requests of the tunnels already up in the same setup ([T1](#t1-tunnels-set-up-by-their-heads)) |
- **`multipath: bool`** (`ShortestPaths` and `spf` only; default `True`): `True` keeps every
  equal-cost route (a DAG), `False` the one canonical path. `None` is not a value. T1 always
  builds single-path tunnels; `ksp` always returns individual paths.
- **`parallel_edges: bool`** (default `True`): `True` keeps every eligible parallel member of a
  `(u, v)` pair; `False` keeps one per pair by the member rule.
- **The member rule** (`parallel_edges=False`). Of the eligible parallel members of a pair, keep
  the one with the lowest cost; among those, under a `HIGHER_*` rule, the one with the highest
  score; then the lowest EdgeId. This choice happens before any route is searched, so it can
  change which walks exist.
- **No silent no-op.** With `multipath=True` and `parallel_edges=True` a tie rule has nothing to
  choose, so a `HIGHER_*` value there is rejected: `spf tie_break requires multipath=False or
  parallel_edges=False` (and `ShortestPaths tie_break requires multipath=False or
  parallel_edges=False` for a model). `FEWER_HOPS` is always valid and never prunes the all-routes
  DAG by hop count. T1 and `ksp` always order paths and need no such rule.
- **The tie rule never supplies weights** (`ShortestPaths` `WEIGHTED` reads a `WeightBasis`) and
  never enables a controller objective (T2's headroom objective is its own option).
- **Static eligibility** (always): the Eligibility. Capacity is never a condition: an edge of
  capacity 0 is searched like any other.
- **Residual gating** (only when a `residual` array is supplied): `residual` is an `int64` array
  in EdgeId order with every entry in `[0, edge_capacity_units[e]]`. It excludes edges whose
  residual is below the threshold (`min_residual` units for `spf`, 1 unit for `ksp`) and supplies
  the availability of the tie scores. It is read during the call only and must not change during
  it.

### Shortest paths

`spf(graph, roots, SpfOptions(direction, *, targets=None, tie_break=FEWER_HOPS, multipath=True,
parallel_edges=True, residual=None, min_residual=1, eligibility=None))` returns an `SpfResult` with
`costs`, `dag` and `target_cost`.

- **Roots and targets** are one node id or a sequence. Every root starts at distance 0.
  Duplicates are removed. Validated in this order before any work (`ValueError`): the roots
  (non-empty, ids in range), the targets (the same), a node in both sets, `multipath=False`
  without targets (`multipath = false requires targets (the single path joins a root to a
  target)`), the ineffective tie rule above, the residual (length, then every entry),
  `min_residual` (an integer of at least 1), a foreign Eligibility. Excluded roots and targets are
  skipped; when every root is excluded the result is empty, not an error.
- **Direction.** `FORWARD` searches from the roots along the edges. `REVERSE` searches toward the
  roots against the edges, but the DAG it returns keeps the edges' own orientation: its walks
  run from a target toward a root. A forward and a reverse search with the roles of roots and
  targets exchanged describe the same shortest walks.
- **Without targets** every reachable node gets its final distance: `costs` is a read-only
  `int64[N]` array with `UNREACHABLE` for the nodes not reached, `dag` holds every tight edge,
  `target_cost` is `None`.
- **With targets** the search stops at the cost of the nearest targets (anycast: the nearest of
  a set). `target_cost` is that cost (`None` when no target is reached), `dag` keeps only the
  walks that end in a target at that cost, and `costs` is `None`.
- **The DAG is unique.** With `multipath=True`, costs are positive, so the member graph, the
  roots, the direction and the targets (or their absence) determine one DAG. The member graph
  changes through `parallel_edges`, the member rule, static eligibility and residual gating; with
  the member graph fixed, the tie rule, the capacity and residual values and the C++ queue choice
  do not change the DAG. A forward search with targets and a reverse search with the
  roots and the targets exchanged return the same DAG. Without targets the two directions can
  return different DAGs, also when the reverse search is rooted at a node the forward search
  reaches: each keeps the edges tight in its own direction.

### The canonical single path

With `multipath=False` (targets required), `spf` returns the one walk from a root to a target
that is first in this order:

1. the lowest **cost**;
2. under a `HIGHER_*` tie rule only: the highest **bottleneck**, the smallest edge score along
   the walk (`FEWER_HOPS` skips this step);
3. the fewest **hops**;
4. the lexicographically smallest **sequence of EdgeIds** along the walk.

The order is total, so the walk is unique. The same walk is returned by a forward search, by a
reverse search with the roles exchanged, and by `ksp` with `k=1`. Note the order of steps 3 and
4: a walk with fewer hops wins even when its first EdgeId is larger. These steps compare whole
walks, never their first edge alone. The member rule above is a separate, earlier choice: with
`parallel_edges=False` each pair is first reduced to one member, and the canonical path is then
chosen among the walks that remain.

```python
F = ngc.Direction.FORWARD
# Two equal-cost walks 0->1->3 (capacity 1) and 0->2->3 (capacity 5).
g = graph(4, [(0, 1, 1, 1), (0, 2, 5, 1), (1, 3, 1, 1), (2, 3, 5, 1)])
for tie in (ngc.RouteTieBreak.FEWER_HOPS, ngc.RouteTieBreak.HIGHER_CAPACITY):
    r = ngc.spf(g, 0, ngc.SpfOptions(F, targets=3, tie_break=tie, multipath=False))
    print(tie.name, r.target_cost, r.dag.via_edges().tolist())

# Equal cost 2: the direct edge (EdgeId 1) against the two-hop walk [0, 2].
g = graph(3, [(0, 1, 1, 1), (1, 2, 1, 1), (0, 2, 1, 2)])
r = ngc.spf(g, 0, ngc.SpfOptions(F, targets=2, multipath=False))
print("fewest hops:", [p.hops() for p in ngc.resolve_to_paths(g, r.dag, 0, 2)])
full = ngc.spf(g, 0, ngc.SpfOptions(F))
print("costs", full.costs.tolist(), "dag", full.dag.via_edges().tolist())
```

```text
FEWER_HOPS 2 [0, 2]
HIGHER_CAPACITY 2 [1, 3]
fewest hops: [[[1]]]
costs [0, 1, 2] dag [0, 1, 2]
```

The bottleneck is the whole walk's, not its first edge's. Below, a direct link 0 -> 2 of capacity 3
and cost 2 (EdgeId 1) competes with 0 -> 1 -> 2, cost 1 + 1 (EdgeIds 0 and 2), whose first link has
capacity 10. With the second hop at 10, `HIGHER_CAPACITY` takes the two hops (bottleneck 10 against
3); with it at 1, the two-hop bottleneck is 1, and every rule takes the direct link:

```python
for second in (10, 1):
    g = graph(3, [(0, 1, 10, 1), (1, 2, second, 1), (0, 2, 3, 2)])
    picks = []
    for tie in (ngc.RouteTieBreak.FEWER_HOPS, ngc.RouteTieBreak.HIGHER_CAPACITY):
        r = ngc.spf(g, 0, ngc.SpfOptions(F, targets=2, tie_break=tie, multipath=False))
        picks.append(r.dag.via_edges().tolist())
    print("second hop", second, picks)
```

```text
second hop 10 [[1], [0, 2]]
second hop 1 [[1], [1]]
```

### K shortest paths

`ksp(graph, source, destination, KspOptions(k, *, tie_break=FEWER_HOPS, parallel_edges=True,
cost_ceiling=None, residual=None, eligibility=None))`
returns a list of `(cost, Path)` pairs; every hop of every Path has one member.

- **Order.** The simple `source -> destination` paths of the member graph in increasing `(cost,
  -bottleneck under a HIGHER_* rule, hop count, edge-id sequence)`: the canonical order above
  applied to whole paths. The result is the first `k` paths of that order that are within the
  ceiling.
- **Parallel members.** `parallel_edges=True` makes every eligible member a separate path;
  `parallel_edges=False` keeps one member per pair by the member rule. `ksp` has no `multipath`.
- **Ceiling.** `cost_ceiling` is `None`, `AbsoluteCost(c)` (cost at most `c`; `AbsoluteCost(0)` is
  a valid budget) or `CostStretch(f)` (cost at most `f` times the static shortest cost of the
  member graph under the Eligibility, residual ignored, edges of capacity 0 included). Both are
  inclusive and compared exactly ([CostStretch](#cost-ceilings-and-the-exact-stretch)). With a
  stretch, no static route gives `[]`.
- **Residual gating** at 1 unit when a residual is supplied.
- Validated in this order before any work (`ValueError` unless stated): `source` and
  `destination` in range, `source == destination` (`node <v> is in both endpoint sets: source and
  destination must not overlap`), `k >= 1` (`k must be >= 1 (got <k>)`), the ceiling (a value of
  another type: `TypeError("cost_ceiling must be AbsoluteCost, CostStretch or None")`; then
  `cost_ceiling AbsoluteCost cost must be >= 0 (got <c>)` or `cost_ceiling CostStretch factor must
  be >= 1 (got <n>/<d>)`), the residual (length, then every entry), a foreign Eligibility.
  `KspOptions` keeps the ceiling as given; `ksp` checks it in this place. An excluded endpoint or
  no path gives `[]`.

```python
g = graph(3, [(0, 1, 3, 1), (0, 1, 5, 1), (1, 2, 9, 1)])   # parallel members 0 and 1
higher = ngc.RouteTieBreak.HIGHER_CAPACITY
for options in (ngc.KspOptions(5), ngc.KspOptions(5, tie_break=higher),
                ngc.KspOptions(5, tie_break=higher, parallel_edges=False)):
    print([(c, p.hops()) for c, p in ngc.ksp(g, 0, 2, options)])
r = ngc.spf(g, 0, ngc.SpfOptions(F, targets=2))
print([p.hops() for p in ngc.resolve_to_paths(g, r.dag, 0, 2)],
      [p.hops() for p in ngc.resolve_to_paths(g, r.dag, 0, 2, split_parallel_edges=True)])
```

```text
[(2, [[0], [2]]), (2, [[1], [2]])]
[(2, [[1], [2]]), (2, [[0], [2]])]
[(2, [[1], [2]])]
[[[0, 1], [2]]] [[[0], [2]], [[1], [2]]]
```

### Cost ceilings and the exact stretch

`AbsoluteCost(cost)` is an inclusive cost bound, an integer (`AbsoluteCost cost must be an integer
(got <type>)`, `TypeError`); a negative cost constructs and is refused where it is used.
`CostStretch(factor)` is one exact reduced fraction `n/d` with `d >= 1`, converted in this order:

1. a boolean (Python `bool` or `numpy.bool_`) raises `TypeError`;
2. a `fractions.Fraction` is taken as it is;
3. a float (a Python `float` or any NumPy floating value) is `Fraction(str(f))`, the shortest
   decimal that value's own type prints: `CostStretch(1.2)`, `CostStretch(numpy.float32(1.2))` and
   `CostStretch(Fraction(6, 5))` are all exactly 6/5, never the binary value of the float. A NaN
   or an infinity raises `ValueError` (`CostStretch factor must be finite (got <value>)`);
4. an integer (an `int`, a NumPy integer, any object with `__index__`) is `n/1`;
5. anything else raises `TypeError` (`CostStretch factor must be a Fraction, a float or an integer
   (got <type>)`).

Then the reduced components are checked: `CostStretch fraction numerator must be in [-2^63, 2^63)`,
then `CostStretch fraction denominator must be in [1, 2^63)` (`Fraction(2**64, 2**62)` reduces to
4/1 and is accepted). `factor` returns the `Fraction`; equality, `hash` and `repr` follow it. A
factor below 1 is a valid value, refused where the ceiling is used. A path of cost `c` passes
`CostStretch(n/d)` against a baseline `c0` exactly when `c * d <= n * c0`, both products computed
exactly in 128 bits. For `ksp` the baseline is the static shortest cost; for T1 it is one baseline
for the whole pool of heads ([T1](#t1-tunnels-set-up-by-their-heads)).

```python
from fractions import Fraction
print(ngc.CostStretch(1.2).factor, ngc.CostStretch(1.2) == ngc.CostStretch(Fraction(6, 5)),
      ngc.CostStretch(Fraction(2**64, 2**62)))
# Two links 0 -> 1 of cost 5 and 6: a stretch of 1.2 admits the cost-6 link at equality, 11/10 does not.
g = graph(2, [(0, 1, 1, 5), (0, 1, 1, 6)])
for factor in (1.2, Fraction(11, 10)):
    paths = ngc.ksp(g, 0, 1, ngc.KspOptions(2, cost_ceiling=ngc.CostStretch(factor)))
    print(factor, [cost for cost, _ in paths])
```

```text
6/5 True CostStretch(Fraction(4, 1))
1.2 [5, 6]
11/10 [5]
```

## 5. Route allocations and admission

**In plain words.** A route allocation decides how a demand's traffic is split at every router of
its routes. `EQUAL` splits equally (ECMP); `WEIGHTED` splits by fixed integer weights (WCMP);
`OPTIMIZED` is a controller that sends as much as the routes can carry, the cheapest way.
Admission is a separate choice. `LOSSLESS` admits the largest total that the split carries with no
edge overloaded and no share left without a route, and refuses the rest before it enters. `LOSSY`
forwards the whole offer hop by hop, as routers do, and drops on each edge what does not fit.

**Exact rules.**

- **Members.** A demand uses the part of its route (a DAG or a path) that reaches its destination
  set through the edges its allocation may use: its members.
  - `EQUAL`: every edge of the route, an edge of capacity 0 included. A router splits over every
    next hop it has, so such an edge receives its round-robin share and carries none of it
    (`LOSSLESS` then admits only totals that give it nothing; `LOSSY` drops the share there).
  - `WEIGHTED`: the edges of positive weight (below). A branch left without such an edge is
    dropped from the route, inside the existing route only: it never takes a longer bypass.
  - `OPTIMIZED`: every edge of the route, an edge of capacity 0 included, so a zero-capacity
    shortest route is never bypassed; the solve can put flow only where its service view has free
    capacity.
  - `LOSSY` forwards over the same members as `LOSSLESS`.
- **`OPTIMIZED`** solves the minimum-cost flow of the demand's cumulative amount on its route,
  over its service view (the free capacity plus the policy's own current allocation). It may move
  the policy's own earlier allocation (that is what lets it reach the maximum); it never moves
  another policy's. It fills the cheapest walks of a mixed-cost route first. Per-edge changes
  between calls can be negative. With `ShortestPaths` it solves the whole demand at once, over the
  union of its sources' routes ([section 7](#routing-models)). `OPTIMIZED` is Lossless only: a
  controller never drops. A solve beyond the flow core's size limits raises `ValueError`
  ([section 8](#limits-that-refuse-a-call)).
- **`WEIGHTED`** splits at every node by highest averages over integer weights frozen when the
  route is prepared. The **weight basis** (`WeightBasis`) says where they come from:

  | `weight_basis` | A member's weight |
  | --- | --- |
  | `CAPACITY` (the default) | the static capacity of the outgoing link, not of the branch behind it |
  | `EXPLICIT` | `edge_weights[e]`: one integer per edge in EdgeId order, each in `[0, 2^62)`, validated in full; a positive weight on an edge of capacity 0 becomes 0 |
  | `RESIDUAL_SNAPSHOT` | the member's free capacity, read once when the program is prepared: for a policy, its service view at that call's entry (the ledger's residual, plus the policy's own traffic for `evaluate` and `replace`); for `from_routes` and `from_dag`, the `residual` they are given |

  The tie rule never chooses weights. Weights are never renormalised: a member whose residual
  later shrinks limits the demand. Equal positive weights give exactly round robin:
  `H(x, [w] * k) = RR(x, k)`.
- **`EQUAL`** splits at every node by round robin over its members, which are frozen when the
  program is prepared.
- **The lossless bound** (`EQUAL`, `WEIGHTED` and every fixed program). A total fits when its
  forwarding from zero, by the program's own splits and with no drops, puts no edge above its free
  capacity and leaves no share without a route (a router or source whose members all have weight
  0, a head with no usable entry). Because the split rules are monotone, raising a total never
  takes units away from an edge, so the totals that fit are exactly `0` to some largest `T_max`. A
  `LOSSLESS` call raises the demand's delivered amount `T` to the smaller of `T` plus the call's
  offer and `T_max`. When a capacity is what stops it, that edge is exactly full at `T_max`: there
  is no safety margin. A share without a route can stop it earlier, with every edge below its
  capacity (the example of two sources below).
- **Lossy forwarding** (`ShortestPaths` with `WEIGHTED` or `EQUAL`, T1, T2 and `Pinned`). Each
  call's offer is split at every node by the allocation's rule; each edge carries the smaller of
  what reaches it and its residual, and the rest is dropped on that edge (where several routes of
  one program, or of several demands forwarded together, meet on an edge,
  [the loss rule](#losses-on-shared-edges) divides what it carries). The delivery is what reaches
  the destination set. A share that has no edge to go on (a `WEIGHTED` node whose members all have
  weight 0) is lost where it is: at a source it is a source drop (`source_drops`). Records keep
  `offered`, `delivered` and per-edge `edge_drops`. Lossy calls are independent: nothing is
  offered again and two calls are not equivalent to one. An infinite Lossy offer is refused
  (`LOSSY forwarding requires a finite offer`) except for `ShortestPaths` `place(float('inf'))`.
- **The first-hop experiment** is that one exception. It asks how much the sources can push into
  their first links at once: every source injects on each of its first edges exactly that edge's
  free capacity, all at once, and forwarding then runs as usual from the next node on. The
  injection is per edge, never one total split over the members:

```python
# Two parallel links of 4 and 10. The experiment puts 4 on one and 10 on the other; an ECMP
# offer of the same 14 splits it 7 and 7 and drops 3 on the 4-unit link.
g = graph(2, [(0, 1, 4, 1), (0, 1, 10, 1)])
lossy = config(ngc.ShortestPaths(), admission=ngc.Admission.LOSSY)
for offer in (float("inf"), 14):
    r = ngc.FlowPolicy(ngc.FlowLedger(g), 0, 1, lossy).place(offer)
    print(offer, r.offered, r.call_delivered, r.call_edge_drops, r.termination.name)
```

```text
inf 14 14 [] EXHAUSTED
14 14 11 [(0, 3)] COMPLETE
```

On a diamond with paths of capacity 10 and 2, equal-share forwarding stops at 5 units. Node 0's
members are EdgeId 0 (toward node 1, capacity 10) and EdgeId 1 (toward node 2, capacity 2).
`RR(5, 2) = (3, 2)` fits; `RR(6, 2) = (3, 3)` would put 3 units on the 2-unit path. The weighted
allocation splits by the capacities 10 and 2, and the optimized one solves the flow; both carry
12:

```python
g = graph(4, [(0, 1, 10, 1), (1, 3, 10, 1), (0, 2, 2, 1), (2, 3, 2, 1)])
allocations = ngc.RouteAllocation.EQUAL, ngc.RouteAllocation.WEIGHTED, ngc.RouteAllocation.OPTIMIZED
for allocation in allocations:
    ledger = ngc.FlowLedger(g)
    r = ngc.FlowPolicy(ledger, 0, 3, config(ngc.ShortestPaths(allocation=allocation))).place(100)
    print(allocation.name, r.call_delivered, r.call_undelivered, r.termination.name, ledger.carried().tolist())

ledger = ngc.FlowLedger(g)
pol = ngc.FlowPolicy(ledger, 0, 3, config(ngc.ShortestPaths(), admission=ngc.Admission.LOSSY))
r = pol.place(10)
print("LOSSY", r.call_delivered, r.call_undelivered, r.termination.name, "drops", r.call_edge_drops)
rec = pol.records()[0]
print("record", rec.offered, rec.delivered, rec.edge_drops)
```

```text
EQUAL 5 95 EXHAUSTED [3, 2, 3, 2]
WEIGHTED 12 88 EXHAUSTED [10, 2, 10, 2]
OPTIMIZED 12 88 EXHAUSTED [10, 2, 10, 2]
LOSSY 7 3 COMPLETE drops [(1, 3)]
record 10 7 [(1, 3)]
```

The Lossy call is `COMPLETE` although 3 units were dropped: the whole offer of 10 was forwarded,
and what did not fit is accounted as drops, not as a shortfall ([section 7](#termination)).

An edge of capacity 0 on the only shortest route is still that route. Equal splitting sends it
the whole offer: Lossless places nothing, Lossy drops all of it there. Weighted splitting gives
it weight 0 and has no other edge in its route, so nothing can be sent (`NO_ACTIVE_WEIGHT`) and
a Lossy offer is lost at the source; neither takes the costlier detour. A controller
(`CentralizedTE`) uses the detour:

```python
# Edge 0: 0 -> 1 of capacity 0 (cost 1); the detour 0 -> 2 -> 1 (5 units, cost 2).
g = graph(3, [(0, 1, 0, 1), (0, 2, 5, 1), (2, 1, 5, 1)])
r = ngc.spf(g, 0, ngc.SpfOptions(F, targets=1))
print("spf", r.target_cost, r.dag.via_edges().tolist())
for allocation in (ngc.RouteAllocation.EQUAL, ngc.RouteAllocation.WEIGHTED):
    for admission in (ngc.Admission.LOSSLESS, ngc.Admission.LOSSY):
        cfg = config(ngc.ShortestPaths(allocation=allocation), admission=admission)
        r = ngc.FlowPolicy(ngc.FlowLedger(g), 0, 1, cfg).place(5)
        print(allocation.name, admission.name, r.call_delivered, r.call_edge_drops, r.call_source_drops,
              r.termination.name)
cfg = ngc.FlowPolicyConfig(routing=ngc.CentralizedTE(source_allocation=ngc.SourceAllocation.OPTIMIZED))
print("CentralizedTE", ngc.FlowPolicy(ngc.FlowLedger(g), 0, 1, cfg).place(5).call_delivered)
```

```text
spf 1 [0]
EQUAL LOSSLESS 0 [] 0 EXHAUSTED
EQUAL LOSSY 0 [(0, 5)] 0 COMPLETE
WEIGHTED LOSSLESS 0 [] 0 NO_ACTIVE_WEIGHT
WEIGHTED LOSSY 0 [] 5 NO_ACTIVE_WEIGHT
CentralizedTE 5
```

Weights: two parallel links of 10 and an offer of 12. Capacity weights split it 6 and 6, explicit
weights 1 and 2 give 4 and 8. With only 2 units free on the first link (a service baseline of 2
and 10), equal splitting admits 4, while weights snapshot from the free capacity (2 and 10) admit
all 12:

```python
g = graph(2, [(0, 1, 10, 1), (0, 1, 10, 1)])
W = ngc.RouteAllocation.WEIGHTED
for model in (ngc.ShortestPaths(allocation=W),
              ngc.ShortestPaths(allocation=W, weight_basis=ngc.WeightBasis.EXPLICIT, edge_weights=[1, 2])):
    ledger = ngc.FlowLedger(g)
    pol = ngc.FlowPolicy(ledger, 0, 1, config(model))
    print(pol.config.routing.weight_basis.name, pol.place(12).call_delivered, ledger.carried().tolist())
snapshot = ngc.ShortestPaths(allocation=W, weight_basis=ngc.WeightBasis.RESIDUAL_SNAPSHOT)
for model in (ngc.ShortestPaths(), snapshot):
    print(model.allocation.name, ngc.FlowPolicy(ngc.FlowLedger(g, [2, 10]), 0, 1, config(model)).place(12).call_delivered)
```

```text
CAPACITY 12 [6, 6]
EXPLICIT 12 [4, 8]
EQUAL 4
WEIGHTED 12
```

Two sources share a total equally, the first unit to source 0. Source 0's only link has capacity
0, so `WEIGHTED` gives it weight 0 and source 0's share has no route. `LOSSLESS` must give source 0
nothing, which only the total 0 does, although source 1's link is entirely free:

```python
g = graph(3, [(0, 2, 0, 1), (1, 2, 10, 1)])
ledger = ngc.FlowLedger(g)
r = ngc.FlowPolicy(ledger, [0, 1], 2, config(ngc.ShortestPaths(allocation=W))).place(float("inf"))
print(r.call_delivered, r.termination.name, ledger.residual().tolist())
```

```text
0 EXHAUSTED [0, 10]
```

## 6. The ledger

**In plain words.** A `FlowLedger` holds what is placed on one graph, in two projections that are
never added together. The **carried traffic** is a service baseline (the starting residual) and
one record per flow; its residual is always the baseline minus everything the records carry. The
**strict holds** are a reservation baseline and the bandwidth that up T1 tunnels reserve on their
paths; its unreserved capacity is that baseline minus every hold. A hold never polices traffic:
traffic may use bandwidth that is held, and a tunnel may carry more than it holds. Policies
place, hold and remove; the ledger has no public placement, reservation or removal of its own.
Its own `reset()`, rejected while a policy is bound, replaces both baselines and clears the
records.

**Exact rules.**

- **`FlowLedger(graph, service_baseline=None, *, reservation_baseline=None)`.** The service
  baseline B and the reservation baseline R0 are each an integer array in EdgeId order with every
  entry in `[0, edge_capacity_units[e]]` (a physical baseline goes through `to_units` first;
  anything above the capacity is rejected); `None` means the capacities. `service_baseline` is
  checked first, then `reservation_baseline`, with the same messages naming each
  (`service_baseline has length <n>, expected <E> (one entry per edge)`, then
  `service_baseline[<e>] = <v> is outside [0, edge_capacity_units[<e>] = <c>]`). The baselines
  never change except through `reset`. A ledger that should model load from outside the library
  starts from reduced baselines.
- **Snapshots.** `carried()` L is the exact sum of every record's `edge_flows` and `residual()` is
  `B - L`; `reserved()` S is the exact sum of every record's `edge_reservations` and
  `unreserved()` is `R0 - S`; all four are `int64[E]` copies. Each projection keeps its own bound
  on every edge (`0 <= L <= B`, `0 <= S <= R0`); L and S are never added, and `max(L, S)` is no
  occupancy. T1's setup reads `unreserved` (other owners' holds included); every forwarding, every
  `CentralizedTE` solve and T2's layout read `residual` (a T2 plan books no hold). Traffic without
  a hold can use bandwidth that T1 holds: there is no priority queueing, preemption or policing.
  `edge_capacity_units()` is the graph's view and `graph()` returns the graph.
- **Records.** `records()` lists the live records (of bound policies and orphans) as snapshots
  in FlowId order; `record(flow_id)` reads one. A FlowId is the pair `(ledger uid, id)`; ids are
  never reused. A FlowId of another ledger, or a retired or unknown id, raises `ValueError`.
- **`FlowRecord`** fields: `flow_id`; `policy_id` (the ledger-local creation-order identifier of
  the policy, kept after it becomes an orphan); `sources` and `destinations` (node sets);
  `admission`; `entry_id` (a T1, T2 or `Pinned` record: the identifier of the program entry it
  carries; `None` for the records kept per roster source); `route` (the record's route as a
  `PredecessorDAG`); `min_route_cost` (the cheapest walk of the route, cached; `None` for an empty
  route); `offered` (Lossy records only, else `None`); `delivered`; `edge_flows`
  (`[(EdgeId, units)]`, ascending, the carried traffic), `edge_drops` (`[(EdgeId, units)]`,
  cumulative, ascending) and `edge_reservations` (`[(EdgeId, units)]`, ascending: an up T1
  tunnel's request on every edge of its path; `[]` for every other record). A record has no
  allocation field: the policy's configuration says how its traffic was computed and the program
  how it is forwarded.
- **`saturated_boundary_edges(sources, eligibility=None)`** is a saturation report of the carried
  traffic (holds play no part): the edges (eligible, at least 1 unit of capacity) that leave the
  set of nodes reachable from the allowed sources over edges with at least 1 unit of residual, in
  ascending EdgeId. It follows edges forward only. It ignores the destinations, so it is not a
  source-destination cut, and it can be empty while a route exists.
- **Lifetimes.** The graph outlives the ledger and the ledger outlives its policies (Python keeps
  them alive for you).

| Ledger state | Reached by | `reset(service_baseline=None, *, reservation_baseline=None)` | `copy()` |
| --- | --- | --- | --- |
| a policy is bound | building a policy on it | `ValueError` | `ValueError` |
| no policy bound | no policy built yet, or every policy destroyed | allowed: orphaned records with their allocations and holds are cleared, their ids retired, the new baselines set (`None`: the capacities) | a new ledger with a new uid and the same baselines, records, holds and ids |

A policy destroyed without `remove()` leaves its records, allocations and holds in the ledger as
**orphans**: readable, never changed again, cleared only by `reset`.

The two projections of one link of 10, where a T1 tunnel holds 6 and
carries 4; a second demand on the same ledger places the residual 6, whatever is held:

```python
g = graph(2, [(0, 1, 10, 1)])
ledger = ngc.FlowLedger(g)
t1 = ngc.FlowPolicy(ledger, 0, 1, ngc.FlowPolicyConfig(routing=ngc.DistributedTE(ngc.TunnelCount(1))))
t1.replace(4, design_demand=6)
print(ledger.carried().tolist(), ledger.residual().tolist(),
      ledger.reserved().tolist(), ledger.unreserved().tolist())
print(ngc.FlowPolicy(ledger, 0, 1, config(ngc.ShortestPaths())).place(10).call_delivered)
```

```text
[4] [6] [6] [4]
6
```

The tunnel holds 6 and carries 4, so the residual is 10 - 4 = 6 and the unreserved capacity
10 - 6 = 4. The second demand reads the residual alone and places 6.

## 7. Policies

**In plain words.** A `FlowPolicy` is one demand between a set of sources and a set of
destinations, routed by one model. `place(offer)` asks it to carry more: it computes or forwards,
commits everything at once, and reports what was delivered, what stopped it, where traffic was
lost and what each source carried. T1 and T2 size their tunnels for a design demand and take whole
offers only. `evaluate`, `replace` and `forward` work on a whole new offer: `evaluate` computes it
and publishes nothing, `replace` computes it and publishes it, `forward` sends it over the routes
already installed. `remove()` takes all of the policy's traffic back out. A demand with several
sources is one demand, not several: how its traffic enters at the sources is its **source
allocation**, and with `ShortestPaths` `EQUAL` or `WEIGHTED`, where the routes of its sources meet,
their traffic is added and split once.
Several policies of one ledger can act in one operation, `apply_actions`: one after another, or
all at the same time with their losses shared ([several demands at once](#several-demands-at-once)).

The subsections follow a demand's life: its configuration, its roster and source allocation, the
routing models (T1 and T2 in detail), the programs they produce, the losses where traffic meets,
the calls and what they report, several demands at once, and last the
[configuration and call errors](#configuration-and-call-errors).

### Configuration

    FlowPolicyConfig(*, routing, admission=Admission.LOSSLESS)
    ShortestPaths(*, allocation=RouteAllocation.EQUAL, source_allocation=None, weight_basis=None,
                  edge_weights=None, tie_break=RouteTieBreak.FEWER_HOPS, multipath=True,
                  parallel_edges=True)
    DistributedTE(tunnel_count, *, cost_ceiling=None, min_unreserved_bandwidth=0,
                  tie_break=RouteTieBreak.FEWER_HOPS, parallel_edges=True)          # T1
    CentralizedTE(*, source_allocation, source_weights=None, support=None, tunnel_count=None,
                  prefer_headroom=False)                                          # T2 with tunnel_count
    Pinned(program, *, black_holes=None)
    TunnelCount(per_head, *, up_to=False)

`FlowPolicyConfig` holds only the routing model and the admission; every other option belongs to
the model that reads it. Every value takes keywords only after its required arguments, and every
field is read-only (a model is an immutable value; its arrays are shared by reference, so policies
built from one model read one `edge_weights` array).

| Routing | Models | Fields |
| --- | --- | --- |
| `ShortestPaths()` | IGP forwarding (as OSPF or IS-IS routers do): static shortest paths, never refreshed. `ShortestPaths()` is ECMP; `allocation=WEIGHTED` is WCMP; `allocation=OPTIMIZED` is a controller restricted to those paths | `allocation`, `source_allocation`, `weight_basis`, `edge_weights` and the route choice (`tie_break`, `multipath`, `parallel_edges`) |
| `DistributedTE(TunnelCount(N))` (T1) | distributed strict-bandwidth tunnels, as RSVP-TE head ends set them up: every head runs CSPF for its own tunnels, no controller ([T1](#t1-tunnels-set-up-by-their-heads)) | `tunnel_count`; `cost_ceiling`: `None`, `AbsoluteCost(c)` or `CostStretch(f)`; `min_unreserved_bandwidth`: a search minimum in units; `tie_break`, `parallel_edges` |
| `CentralizedTE(source_allocation=...)` | an ideal controller that sees the whole network, or only the edges of `support`; with `tunnel_count`, T2: the controller lays out equal tunnels for every head ([T2](#t2-the-controllers-equal-tunnels)) | `source_allocation` (required); `source_weights`: one integer per graph node, required exactly for `WEIGHTED` without tunnels; `support`: a `PredecessorDAG` of the graph, or `None` for every eligible edge; `tunnel_count`: a `TunnelCount(N)` or `TunnelCount(L, up_to=True)`, or `None`; `prefer_headroom`: T2's headroom objective after route cost |
| `Pinned(program)` | forwarding over a given program, computing nothing | `program`: a `ForwardingProgram`; `black_holes`: an `Eligibility` whose excluded nodes and edges are black holes, or `None` |

**Two stages of validation.** A constructor converts its arguments. A wrong type, a missing
required argument, `None` where none is allowed and an unknown keyword (an option of another model
included, such as `allocation` on `FlowPolicyConfig`) raise `TypeError` before any policy exists.
Every other check of a well-typed configuration raises `ValueError` when a policy or an evaluator
uses it, before any work, a zero offer included. The texts and their order are in
[configuration and call errors](#configuration-and-call-errors).

**The model matrix.** The options belong to each row's model; any other well-typed combination
raises `ValueError`:

| Model | route allocation | source allocation: default / allowed | `admission` | route choice | weights |
| --- | --- | --- | --- | --- | --- |
| `ShortestPaths()` | `EQUAL` (ECMP, the default) | `EQUAL` / `EQUAL` | `LOSSLESS`, `LOSSY` | `tie_break`, `multipath`, `parallel_edges` | none |
| `ShortestPaths(allocation=WEIGHTED)` | `WEIGHTED` (WCMP) | `EQUAL` / `EQUAL` | `LOSSLESS`, `LOSSY` | the same | `weight_basis`: `CAPACITY` (default), `EXPLICIT`, `RESIDUAL_SNAPSHOT` |
| `ShortestPaths(allocation=OPTIMIZED)` | `OPTIMIZED` | `OPTIMIZED` / `EQUAL`, `OPTIMIZED` | `LOSSLESS` | the same | none |
| `CentralizedTE(source_allocation=...)`, optionally `support` | optimal | required / EQUAL, OPTIMIZED, WEIGHTED | `LOSSLESS` | none | `source_weights` exactly with WEIGHTED |
| T1: `DistributedTE(TunnelCount(N))` | tunnels | `EQUAL`, fixed | `LOSSLESS`, `LOSSY` | `tie_break`, `parallel_edges` | none |
| T2: `CentralizedTE(source_allocation=..., tunnel_count=...)`, optionally `support` | tunnels | required; `EQUAL`; `OPTIMIZED` only with one source; never `WEIGHTED` | `LOSSLESS`, `LOSSY` | none (`prefer_headroom`) | none |
| `Pinned(program)` | the program's | the program's | `LOSSLESS`, `LOSSY` | none | none |

`ShortestPaths`' `source_allocation=None` resolves to the allocation's default and its
`weight_basis=None` to `CAPACITY` for `WEIGHTED`; a policy's `config` is the configuration it was
built with, these two fields resolved (the model value itself never changes).

`FlowPolicy(ledger, sources, destinations, config, eligibility=None, cache=None)` binds the policy
to its ledger. `sources` and `destinations` are one node id or a sequence; the two sets must not
overlap. Construction validates everything before the policy binds.

### The source roster

The roster is the set of sources that take part: they receive shares of the traffic. It is fixed
when the policy is built.

- For `ShortestPaths`, T1, T2 and `CentralizedTE` it is the allowed sources that have a static walk
  to an allowed destination over allowed edges, whatever their capacity, in ascending node id. A
  source whose only routes cross an edge of capacity 0 is in the roster, receives its share and
  cannot place it; a failure is expressed with an Eligibility, which removes the source.
  `CentralizedTE` with a support counts only walks inside the support.
- For `Pinned` it is the program's heads that still reach a destination under the policy's
  Eligibility ([replay](#replay-under-a-failure)).
- The roster decides which sources receive shares. For the computing models it also decides which
  records exist (T1 and T2: its heads' entries); `Pinned` keeps one record per configured program
  entry even when no head survives the replay (such a record stays at zero).
- `by_source` reports every input source, a source outside the roster with a zero row.
- A positive offer with an empty roster ends `NO_ROUTE`: a computing model creates no record, and
  `Pinned` creates or keeps its per-entry records at zero.

### The source allocation

The source allocation says how a demand's total enters at its roster sources, separately from how
it is routed:

- **`EQUAL`**: a total `x` enters as the round robin `RR(x, k, i) = x // k + (1 if i < x % k else
  0)` at roster member `i` (`k` members, ascending node id): the first members get one extra
  unit. A source whose share cannot go anywhere keeps it and loses it; the share is never moved
  to another source. Lossless therefore admits only totals that give such a source nothing, and
  Lossy loses its share.
- **`OPTIMIZED`** (controllers only): the controller chooses how much enters at each source (one
  virtual source joined to every roster source), whatever maximises the delivery and then
  minimises its cost.
- **`WEIGHTED`** (`CentralizedTE` without `tunnel_count` only): a total `x` enters as the
  highest-averages shares `H(x, w)`, `w` being the declared `source_weights` of the roster sources
  in node id order. `source_weights` is one integer per graph node in `[0, 2^62)`, copied once and
  exposed as an immutable tuple; every entry is validated, unused and masked nodes included, but
  only the roster's weights are read. A zero weight receives nothing; when every roster weight is
  0 nothing enters, with no equal fallback ([termination](#termination)). The other source
  allocations refuse the vector, and `ShortestPaths` and every tunnel configuration refuse
  `WEIGHTED`.

`DistributedTE` always allocates equally and `Pinned` takes its program's injection; neither
accepts the keyword.

Source 0 below reaches node 3 over one link of 2, and source 1 over two links of 10 (costs 1 and
2). Equal shares stop at 4: source 0 can take only 2, so source 1 gets 2 as well, at a cost of
2 x 1 + 2 x 3 = 8. The optimized controller sends 2 from source 0 and 10 from source 1, at a cost
of 2 + 30 = 32. Weights 1, 1 and 2 over three private links of 4 give `H(8, [1, 1, 2]) =
[2, 2, 4]`; a ninth unit would go to the third source, beyond its 4:

```python
S = ngc.SourceAllocation
g = graph(4, [(0, 3, 2, 1), (1, 2, 10, 1), (2, 3, 10, 2)])
for sa in (S.EQUAL, S.OPTIMIZED):
    r = ngc.FlowPolicy(ngc.FlowLedger(g), [0, 1], 3, config(ngc.CentralizedTE(source_allocation=sa))).place(12)
    print(sa.name, r.call_delivered, r.total_cost, [(row.source, row.delivered) for row in r.by_source])
g = graph(4, [(0, 3, 4, 1), (1, 3, 4, 1), (2, 3, 4, 1)])
weighted = ngc.CentralizedTE(source_allocation=S.WEIGHTED, source_weights=[1, 1, 2, 0])
r = ngc.FlowPolicy(ngc.FlowLedger(g), [0, 1, 2], 3, config(weighted)).place(float("inf"))
print("WEIGHTED", r.call_delivered, [(row.source, row.delivered) for row in r.by_source])
```

```text
EQUAL 4 8 [(0, 2), (1, 2)]
OPTIMIZED 12 32 [(0, 2), (1, 10)]
WEIGHTED 8 [(0, 2), (1, 2), (2, 4)]
```

### Routing models

**When a model computes.** A positive `place` prepares the program of an epoch that has none yet.
Every `evaluate` and every `replace`, a zero offer included, computes a whole offer afresh; only
`replace` installs what it computed. `place(0)` prepares nothing.

- **`ShortestPaths` with `EQUAL` or `WEIGHTED`** forwards like IGP routers. Its program holds the
  routes and splits used for the whole epoch, never refreshed within it:
  - each roster source's route is its static shortest-path DAG toward the destination set (with
    `multipath=False`, its canonical path), read from the
    [shortest-path cache](#9-the-shortest-path-cache);
  - the route choice is applied on the policy's service view when the program is prepared (the
    `HIGHER_AVAILABLE_CAPACITY` and `HIGHER_HEADROOM` scores read it), and `WEIGHTED` freezes its
    weights from its weight basis there ([section 5](#5-route-allocations-and-admission));
  - the program is the union of these selected routes, built exactly as
    [`ForwardingProgram.from_dag`](#supplied-dags-from_dag) builds one: at a router several
    sources pass, their traffic for the demand is added (with the router's own injection when it
    is a source) and split once over its members, by round robin (`EQUAL`) or by highest averages
    over the weights (`WEIGHTED`). Every source that passes a router has the same continuation
    there, so the union lends no source a link outside its own route;
  - it is the policy's `program`: one split context pooled over the roster, one split route per
    roster source (entry identifier = roster position, weight 1), equal sources; it lists no
    tunnel.

  Each epoch prepares on that epoch's view. With `FEWER_HOPS` or `HIGHER_CAPACITY` and `EQUAL`,
  `CAPACITY` or `EXPLICIT` weights nothing a preparation reads changes while the policy lives, so
  every epoch gets the program a fresh preparation gives; `HIGHER_AVAILABLE_CAPACITY`,
  `HIGHER_HEADROOM` and `RESIDUAL_SNAPSHOT` read the view of each preparation.
- **`ShortestPaths` with `OPTIMIZED`** is a controller restricted to shortest paths: it solves like
  `CentralizedTE` (below), but only over the union of the roster sources' routes prepared as
  above, edges of capacity 0 included (so a zero-capacity shortest route is never bypassed).
- **`DistributedTE` (T1)** sets up its tunnels for a design demand on every `replace` and
  `evaluate` and forwards them as its program; it has no incremental `place`
  ([T1](#t1-tunnels-set-up-by-their-heads)).
- **`CentralizedTE(tunnel_count=...)` (T2)** lays out its tunnels for a design demand on every
  `replace` and `evaluate` and forwards them as its program; it has no incremental `place`
  ([T2](#t2-the-controllers-equal-tunnels)). The controller rules below are `CentralizedTE`'s
  without a tunnel count.
- **`CentralizedTE`** solves the minimum-cost flow from its roster to its whole destination set
  over every eligible edge, on its service view: it may move its own traffic, never another
  owner's. The result is the largest delivery up to the target, and among those the cheapest. The
  target depends on the call:
  - a positive `place` re-solves the amount the policy already delivers plus this call's offer
    (`delivered()` plus the offer: what an earlier call could not place is not offered again). With
    the rest of the network unchanged, `place(a)` followed by `place(b)` therefore ends in the same
    allocation as `place(a + b)`;
  - `evaluate` and `replace` solve their whole offer afresh, a zero offer included;
  - `forward` solves nothing: it sends its offer over the installed paths by their weights.
  - For the same target total and service view, the final solve from zero fixes the allocation;
    the search for the largest total that fits does not affect it.
  - With `EQUAL` source allocation the virtual source's arc to roster source `i` carries exactly
    `RR(T, k, i)`: the delivered total `T` is the largest whose equal shares all fit (a total fits
    only when the solve reaches it). With `WEIGHTED` the shares are `H(T, w)`, and every source
    must send exactly its share. Each share grows with `T`, so the totals that fit are 0 to a
    largest one.
  - With `EQUAL` or `WEIGHTED` sources the result has no `min_cut` and no `marginal_cost_curve`:
    a flow with prescribed shares is not the unrestricted maximum a cut certifies, and the
    marginal costs of an equal flow need not increase.
  - A positive call on a network beyond the flow core's size limits raises `ValueError`
    ([section 8](#limits-that-refuse-a-call)).
- **`CentralizedTE(support=dag)`** is the same controller restricted to the edges of `dag`,
  trimmed to the walks from the sources to the allowed destinations inside it (an edge that leads
  nowhere inside the support is never used; a mask removes its edges as usual). It reports
  `total_cost` but neither `marginal_cost_curve` nor `min_cut`, and `evaluate_edge_failures`
  solves each deletion afresh.
- **The controllers' programs.** After every solve (a positive `place`, an `evaluate`, a `replace`)
  a controller decomposes its solved flow into paths, each roster source's share in a canonical
  order (for each source in node id order, walk the smallest positive EdgeId to the first
  destination, take the bottleneck or the remaining supply, subtract, repeat), one tunnel per path
  with its solved amount as weight (identifiers in decomposition order), also with a single source.
  With the `OPTIMIZED` source allocation the program keeps the solved supplies as **captured source
  weights** (`source_weights`: every roster source, a zero supply included). With `WEIGHTED` it
  instead keeps each roster source's declared weight, zeros included; with `EQUAL` it injects by
  round robin. Forwarded on the same service view, the program reproduces the solve exactly.
- **Records of the controllers** (`CentralizedTE`, `ShortestPaths` with `OPTIMIZED`): one record
  per roster source, holding that source's share of the solved flow (its edges and its delivery),
  extracted in a canonical order. A `CentralizedTE` record's route is the set of edges that carry
  the share, rebuilt on every positive `place` under the same FlowId (a `replace` or a `forward`
  publishes new records with new ids); a source with no share has an empty route and
  `min_route_cost` `None`. A `ShortestPaths` `OPTIMIZED` record's route is its
  source's prepared route, and the share stays inside it.
- **`Pinned`** forwards its program as given, computing nothing: injection, entries, weights and
  split tables are the program's, replayed under the policy's Eligibility and with its black
  holes ([forwarding programs](#forwarding-programs)). It keeps one record per program entry
  (`entry_id` is the entry's identifier), created when the epoch's program is first forwarded (a
  positive `place`, or any `replace`), a record of weight 0 or of a `FAILED` tunnel included (that
  one with an empty route and `min_route_cost` `None`). A positive offer with an empty replay
  roster ends `NO_ROUTE`.

### T1: tunnels set up by their heads

**In plain words.** `DistributedTE(TunnelCount(N))` models distributed RSVP-TE: there is no
controller, and each head (a roster source) sets up its own N tunnels. A design demand D is split
equally over the heads and each head's share equally over its N tunnels; each tunnel asks for
that bandwidth. One tunnel at a time, in rounds over the heads, the head keeps the links with
enough unreserved bandwidth and takes the shortest path: the tunnel comes up and holds exactly its
request on every link of its path, or it stays down and holds nothing. Real traffic then enters
equally at the heads and splits equally over each head's tunnels that are up; a head with no
tunnel up loses its share, and no one moves it elsewhere. Holds are the ledger's second
projection ([section 6](#6-the-ledger)): they decide which tunnels come up, never what traffic
gets through.

**Exact rules.**

- **Sizing.** `H` roster heads in ascending node id, `N = tunnel_count.per_head` (an "up to" count
  is refused). Head `i`'s design share is `RR(D, H, i)` (`plan.design_by_head`), its tunnel `j`'s
  request `RR(share, N, j)`: the tunnels of one head differ by at most one unit. Tunnels are
  numbered in round order (round `j = 0..N-1`, the heads ascending within a round): entry
  identifier `j * H + i`, local index `j`; this is also the setup order. An empty roster gives an
  empty program.
- **Small shares and zero requests.** When a head's share is below N units, some requests are 0.
  All N tunnels are still set up: a tunnel with request 0 comes up holding nothing, like an
  RSVP-TE tunnel with no bandwidth; its search has no bandwidth filter (unless
  `min_unreserved_bandwidth` is positive), so it may cross an edge of capacity 0, where the
  traffic later sent on it is lost on that edge. A design of 0 is the extreme case: every request
  is 0 and an offer is forwarded over the tunnels.
- **The setup view U** of an edge is the ledger's unreserved capacity with this policy's own holds
  released (never above the capacity). Every attempt runs, also after failures, on U minus the
  requests of the tunnels already up in this setup.
- **One attempt (CSPF)** with request `b` keeps the allowed edges with `U >= max(b,
  min_unreserved_bandwidth)` (a request of 0 with a minimum of 0 keeps every allowed edge, one of
  capacity 0 included), then with `parallel_edges=False` one member per pair by the member rule
  (the filter first, then the member rule, at every attempt), then takes the canonical single path
  of [section 4](#the-canonical-single-path) on static costs toward the allowed destinations: the
  lowest cost; under a `HIGHER_*` rule the larger bottleneck score; fewer hops; the smaller EdgeId
  sequence. The scores read the setup view: `HIGHER_AVAILABLE_CAPACITY` reads U,
  `HIGHER_HEADROOM` U divided by the capacity (exactly; 0 on an edge of capacity 0),
  `HIGHER_CAPACITY` the capacity.
- **The ceiling.** `AbsoluteCost(c)` admits a path of cost at most `c`. `CostStretch(n/d)` admits
  a cost `c` exactly when `c * d <= n * c0`, `c0` the smallest static shortest cost over all
  roster heads under the policy's Eligibility, edges of capacity 0 included, before any bandwidth
  filter: one baseline for the whole pool, not per head (a head whose own shortest cost is 3 can
  fail `CostStretch(2)` because another head has a zero-capacity shortest route of cost 1), fixed
  when the policy is built.
- **Statuses.** `UP`: the tunnel holds `b` on every edge of its path (its record's
  `edge_reservations`) and has weight 1. `INSUFFICIENT_BANDWIDTH`: no path passes the filter.
  `COST_CEILING_EXCEEDED`: the cheapest filtered path exceeds the ceiling. A down tunnel keeps its
  request, has no destination, path or cost (`None`), has weight 0 and holds nothing. There is no
  partial hold, retry, growth or count reduction; `min_unreserved_bandwidth` filters searches and
  is never held. A tunnel count is a number of setup attempts per head, not a number of tunnels
  guaranteed up, and two tunnels may take the same links.
- **Forwarding** is the program's: a total `x` enters as `RR(x, H, i)` at head `i`, then round
  robin over that head's up tunnels in local order. A head with no up tunnel keeps its share and
  loses it (Lossy: a source drop; Lossless admits only totals that give it nothing,
  [equal source allocation](#equal-source-allocation-and-a-source-that-cannot-pass-anything)); the
  demand is never renormalised over the heads that have tunnels. Changing only the offer never
  changes the setup. When every tunnel is up and U is at most the residual on every edge, D itself
  is carried losslessly. Holds do not police: a tunnel carries what forwarding gives it within the
  residual, more than its hold when the edges have room, and its traffic meets every other
  demand's on the residual. T1 reports no `total_cost`, `min_cut` or `marginal_cost_curve`.
- **Operations.** No `place` (`DistributedTE does not support incremental placement`, `place(0)`
  too). `replace(offer, design_demand=D)` releases the policy's own holds and traffic, sets up
  again and forwards; `evaluate` computes the same and publishes nothing (no record, no hold);
  `forward(offer)` forwards the installed program and keeps its tunnels and holds (the new records
  carry the same `edge_reservations`). `design_demand=None` means `D = offer`; an infinite offer
  needs a finite design (`DistributedTE requires a finite design demand`). To study a change of
  the offer alone, keep D explicit: `replace(offer + 1)` without it designs again for the new
  offer. A full zero offer such as `replace(0, design_demand=D)` still sets up and books the
  holds. `remove()` releases both projections. T1 is not monotone in D: more design demand can
  bring more tunnels up or fewer.
- **Records and rows.** One record per program entry in identifier order (`entry_id`), a down
  tunnel's included (zero, an empty route, `min_route_cost` `None`). `by_tunnel` shows a down
  tunnel's setup status and `UP`, or `BLACKHOLED` in stale forwarding; a `Pinned` replay of a T1
  program shows `FAILED` for an up tunnel its mask touches, keeps a down tunnel down and holds
  nothing.
- **The plan.** `result.plan` is a `TunnelPlan` of the program after the call:
  `tunnels_per_head` (N), `design_demand` (D), `design_by_head` (`(head, RR(D, H, i))` for every
  roster head) and `total_reserved_bandwidth` (the up tunnels' requests, each counted once). The
  plan describes the setup that installed the tunnels; `forward` keeps it unchanged.
  T2 also returns a `TunnelPlan`, with the fields described below; all other routings, including
  `Pinned` over a T1 program, return `None`.

A T1 call: head 0's only edge has 2 units, head 1's has 10, two tunnels per head:

```python
g = graph(3, [(0, 2, 2, 1), (1, 2, 10, 1)])
ledger = ngc.FlowLedger(g)
lossy_t1 = ngc.FlowPolicyConfig(routing=ngc.DistributedTE(ngc.TunnelCount(2)),
                                admission=ngc.Admission.LOSSY)
pol = ngc.FlowPolicy(ledger, [0, 1], 2, lossy_t1)
r = pol.replace(12, design_demand=12)
print([t.status.name for t in r.program.tunnels])
print([(t.entry_id, t.head, t.requested_bandwidth, t.reserved_bandwidth, t.edges) for t in r.program.tunnels])
print(r.call_delivered, r.call_source_drops, r.plan.total_reserved_bandwidth, r.plan.design_by_head)
print(ledger.reserved().tolist(), ledger.unreserved().tolist(), ledger.residual().tolist())
r = pol.forward(8)
print(r.call_delivered, r.call_source_drops, r.reservations_changed, ledger.reserved().tolist())
```

```text
['INSUFFICIENT_BANDWIDTH', 'UP', 'INSUFFICIENT_BANDWIDTH', 'UP']
[(0, 0, 3, 0, None), (1, 1, 3, 3, (1,)), (2, 0, 3, 0, None), (3, 1, 3, 3, (1,))]
6 6 6 ((0, 6), (1, 6))
[0, 6] [2, 4] [2, 4]
4 4 False [0, 6]
```

Each head's share is 6, so every tunnel asks for 3. Head 0's edge has 2 units unreserved: both of
its tunnels are down. Head 1's two tunnels come up on edge 1 and hold 6 of its 10. The offer of 12
enters as 6 and 6: head 0 loses its 6, head 1's tunnels carry 3 and 3. `forward(8)` keeps the
tunnels and the holds; head 0 still loses its share, now 4.

The cost ceiling's stretch baseline is one for the whole pool of heads. Below, head 0 reaches node
2 only over a link of capacity 0 and cost 1, head 1 over a link of 8 and cost 3; each head's share
of the design 4 is 2. Head 0's tunnel finds no link with 2 units unreserved. The baseline is the
smallest static shortest cost over all heads, 1 (head 0's zero-capacity route), so
`CostStretch(2)` admits a cost of at most 2 and refuses head 1's path of cost 3; an absolute
ceiling of 3 admits it:

```python
g = graph(3, [(0, 2, 0, 1), (1, 2, 8, 3)])
for ceiling in (ngc.CostStretch(2), ngc.AbsoluteCost(3)):
    t1 = ngc.FlowPolicyConfig(routing=ngc.DistributedTE(ngc.TunnelCount(1), cost_ceiling=ceiling),
                              admission=ngc.Admission.LOSSY)
    r = ngc.FlowPolicy(ngc.FlowLedger(g), [0, 1], 2, t1).replace(4, design_demand=4)
    print(repr(ceiling), [t.status.name for t in r.program.tunnels], r.call_delivered)
```

```text
CostStretch(Fraction(2, 1)) ['INSUFFICIENT_BANDWIDTH', 'COST_CEILING_EXCEEDED'] 0
AbsoluteCost(3) ['INSUFFICIENT_BANDWIDTH', 'UP'] 2
```

### T2: the controller's equal tunnels

**In plain words.** `CentralizedTE(source_allocation=..., tunnel_count=TunnelCount(N))` models a
central controller (a PCE or an SDN controller) that programs, for every head, the same number of
single-path tunnels, all planned with one common whole-unit bandwidth. The controller chooses, in
this order:

1. the largest common bandwidth;
2. for a finite design demand, no larger than the design needs;
3. the cheapest routes;
4. with `prefer_headroom=True`, the routes that leave the most spare capacity;
5. one fixed unit-by-unit procedure for any remaining tie.

Real traffic then splits equally over the heads and over each head's tunnels. A tunnel carries
what forwarding gives it wherever the links have room, more than its planned bandwidth included.
A planned bandwidth is metadata: the ledger holds nothing for it. `TunnelCount(L, up_to=True)`
gives every head one common count `k <= L`: the smallest whose program carries the live heads'
design shares; if none does, the largest live-share score, then the smallest count. It does not
promise whole-demand lossless delivery. A head that cannot reach a destination over free capacity
is dead: it keeps its share and loses it. A weak but live head limits the common size of every
tunnel, by design (heads that must not limit each other are separate demands). A `support` DAG
restricts every step to its edges.

The sizing on one head with four tunnels and two branches, of capacity 4 and 10:

```python
g = graph(4, [(0, 1, 4, 1), (0, 2, 10, 1), (1, 3, 4, 1), (2, 3, 10, 1)])
four = ngc.FlowPolicyConfig(routing=ngc.CentralizedTE(source_allocation=ngc.SourceAllocation.EQUAL,
                                                      tunnel_count=ngc.TunnelCount(4)),
                            admission=ngc.Admission.LOSSY)
r = ngc.FlowPolicy(ngc.FlowLedger(g), 0, 3, four).replace(16)
p = r.plan
print(p.max_fractional_bandwidth_per_tunnel, p.max_bandwidth_per_tunnel, p.total_planned_bandwidth,
      p.max_lossless_offer)
print([t.edges for t in r.program.tunnels], [row.assigned for row in r.by_tunnel], r.call_delivered,
      r.call_edge_drops)
```

```text
10/3 3 12 14
[(0, 2), (1, 3), (1, 3), (1, 3)] [4, 4, 4, 4] 14 [(1, 2)]
```

With `j` of the 4 tunnels on the narrow branch, the common size is at most `4 / j` there and
`10 / (4 - j)` on the wide branch: 5/2, 10/3, 2, 4/3 and 1 for `j` = 0 to 4. The best is 10/3,
with one tunnel on the narrow branch and three on the wide one, and the largest whole-unit size
is `U = 3`; the plan is 4 x 3 = 12. The lossless bound `F` is 14: `RR(14, 4) = (4, 4, 3, 3)` gives
the narrow tunnel 4 and the three wide ones 10 together, while `RR(15, 4)` would put 11 on the
wide branch. The offer 16 gives every tunnel 4, more than its plan of 3: the narrow branch carries
its 4, and the wide branch carries 10 of its 12 and drops 2.

**Exact rules.**

- **Heads.** The roster (`H` heads in ascending node id) is topological, as for T1. A head is
  **live** if it reaches an allowed destination over allowed links with free capacity `V >= 1`,
  `V` the residual with this policy's own traffic released (`plan.heads_with_capacity`); the
  others are **dead**: they keep their roster position and their share, get no tunnel and lose
  the share. `a` live heads, `k` tunnels each (`N`, or the chosen count), `K = a * k` tunnels.
- **On a support.** `CentralizedTE(source_allocation=..., support=dag, tunnel_count=...)` is T2
  with the allowed links restricted to the edges of `dag`, intersected with the policy's
  Eligibility. Every step (roster, live heads, sizing, count choice, both objectives, the
  decomposition, the program's replay basis) uses this one allowed set, and the amounts, plans,
  routes and flags equal those of the same configuration without a support and with the DAG's
  edge set composed into the Eligibility. Parallel members and EdgeIds are kept; the DAG is not
  read as separate routes.
- **Sizing.** Tunnels of a common size `s > 0` exist exactly when `K` path units flow from a
  virtual root (an arc of capacity `k` to every live head) to the destinations with link
  capacities `min(K, floor(V / s))` (paths stop at the first destination).
  `max_bandwidth_per_tunnel` U is the largest integer size that passes (0 when none does);
  `max_fractional_bandwidth_per_tunnel` is the exact fluid optimum `s*` (a `Fraction` `V[e] / j`
  with `j <= K`; `U = floor(s*)`).
- **Working size.** For a finite design D, `planned_bandwidth_per_tunnel` u is `min(U, max(1, max
  over live heads of ceil(RR(D, H, i) / k)))`, the shares taken over the whole roster, dead heads
  included; an unlimited design uses `u = U`. Every tunnel is planned at u;
  `total_planned_bandwidth` is `K * u`. Once U >= 1, every live head keeps all k tunnels even when
  its share is below k units, including a design of 0 (the tunnels are planned at 1 and carry what
  forwarding gives them).
- **Routes.** At u, the count vectors `x` (the tunnels on each link) that route `k` tunnels per
  live head with `x[e] <= min(K, floor(V[e] / u))` are compared in three steps:
  1. the smallest route-cost sum `J = sum of cost[e] * x[e]` (`total_tunnel_cost`; idle tunnels
     count);
  2. with `prefer_headroom=True`, the smallest `headroom_penalty` `Phi`: each link with static
     capacity `B > 0` charges its `j`-th tunnel `floor(2^32 * (B - V + j * u) / B)`, the fraction
     of the link then occupied, rounded down to 2^-32;
  3. the canonical augmentation: `K` single units from zero, each on the residual path from the
     root to the sink of the smallest `(J, Phi)` price, then the fewest arcs, then the smallest
     sequence of arc keys (`(0, EdgeId, direction)` for links, `(1, head, direction)` for root
     arcs, `(2, destination, direction)` for sink arcs; forward 0, reverse 1); a reverse arc
     undoes an earlier unit. Without `prefer_headroom` the second price is 0.

  `J` and `Phi` are exact Python integers and can exceed 64 bits.
- **The program.** `x` is cut, live head by live head in ascending node id, into unit walks along
  the smallest positive EdgeId to a destination; each head's `k` walks are sorted by EdgeId
  sequence (duplicates are separate tunnels). The tunnels are in round order: local index `j`,
  then roster position `i`, with the entry identifier `j * H + i` (`H` counts dead heads too, so
  a dead head leaves holes). Every tunnel is `UP` with weight 1, `planned_bandwidth = u` and
  `requested_bandwidth = reserved_bandwidth = 0`; a plan books no hold.
- **Up to L.** For each `k` in `1..L` the program at `U(k)` scores `G(min(D, F'(k)))`: `G(z)` is
  the part of a round-robin total `z` that the live heads receive, `F'(k)` the largest total whose
  live shares fit on that program (dead shares ignored); an unlimited design scores `G(F'(k))`.
  The smallest `k` whose score is the live part of D wins; if none, the largest score, then the
  smallest `k`; the chosen `k` is then laid out at its own u (`plan.tunnels_per_head`). Delivery
  is not monotone in `k`, so only two proved bounds skip a `k`: with `M = H * k` roster positions
  the score lies between `G(min(D, M * U(k)))` and `G(min(D, M * (U(k) + 1) - 1))`, and `U(k)`
  never grows with `k` (once it is 0, every later count scores 0). "Up to L" finds the smallest
  common count that carries the live heads' design shares (or the best score); it does not
  promise lossless delivery of the whole demand and does not optimise an arbitrary offer.
- **Count outcomes.** An empty roster gives an empty program (`NO_ROUTE` for a positive offer).
  No live head, or live heads with `U = 0` (`insufficient_resolution` True,
  `max_bandwidth_per_tunnel` 0, `max_fractional_bandwidth_per_tunnel` positive), give a program
  **without tunnels that keeps the roster**: `policy.program` is set, `forward` works and every
  head's share is lost at its head (Lossy `COMPLETE`, Lossless `EXHAUSTED` for a positive offer).
  Otherwise every live head has `k` tunnels. No partial program is made for one head; check the
  plan before assuming `len(tunnels) == N * H`.
- **Forwarding** is the program's: a total `x` enters as `RR(x, H, i)` over the whole roster,
  then round robin over the head's `k` tunnels in local order; a dead head's share is lost.
  `max_lossless_offer` F is the program's lossless maximum on the view it was planned on (a dead
  head's share blocks every total that gives it a unit: F is then 0 or a few remainder units).
  With `u = U` and every head live, `0 <= F - K * U < K`.
- **Configuration.** `source_allocation` is stated: `EQUAL` with any number of heads; `OPTIMIZED`
  only with one input source, the same calculation (`CentralizedTE with tunnel_count and several
  sources requires source_allocation=EQUAL`, counted before any Eligibility). `prefer_headroom`
  (`False` by default) switches the headroom objective on; T2 has no route choice of its own:
  it always keeps every parallel link, uses single-path tunnels and the canonical final tie.
  Admission `LOSSLESS` or `LOSSY`. No path-cost ceiling. The tunnel count and the tunnel limit as
  for T1 (`L` counts for "up to").
- **Operations.** No `place` (`CentralizedTE with tunnel_count does not support incremental
  placement`, `place(0)` too). `replace(offer, design_demand=D)` releases the policy's own traffic,
  lays the tunnels out again for D and forwards; `design_demand=None` means `D = offer`, and
  `float('inf')` (or an infinite Lossless offer) is an unlimited design. `evaluate` computes the
  same and publishes nothing; `forward(offer)` forwards the installed program and keeps its plan
  unchanged: its bounds and headroom values still describe the planning view, not the new offer,
  load or black holes. A planned size change over the same routes sets no flag. T2 reports no
  `total_cost`, `min_cut` or `marginal_cost_curve`.
- **Records and rows.** One record per tunnel in identifier order (a dead head has no entry and no
  record).
- **The plan.** `tunnels_per_head`, `design_demand`, `design_by_head` (empty for an unlimited
  design) and `total_reserved_bandwidth` (0), then T2's own fields: `heads_with_capacity`,
  `insufficient_resolution`, `max_fractional_bandwidth_per_tunnel`, `max_bandwidth_per_tunnel`,
  `planned_bandwidth_per_tunnel`, `total_planned_bandwidth`, `total_tunnel_cost`,
  `headroom_penalty` (`None` without `prefer_headroom`) and `max_lossless_offer`. With no live head
  every amount is 0 (`tunnels_per_head` is N, or 1 for "up to"). A plan is metadata of the
  computation that installed the tunnels, neither the installed program nor its delivery.

Two heads behind private links of 3 and 8 to node 2, two tunnels each, design 6:

```python
g = graph(3, [(0, 2, 3, 1), (1, 2, 8, 1)])
two = ngc.FlowPolicyConfig(routing=ngc.CentralizedTE(source_allocation=ngc.SourceAllocation.EQUAL,
                                                     tunnel_count=ngc.TunnelCount(2)))
r = ngc.FlowPolicy(ngc.FlowLedger(g), [0, 1], 2, two).replace(6)
p = r.plan
print(p.max_fractional_bandwidth_per_tunnel, p.max_bandwidth_per_tunnel, p.planned_bandwidth_per_tunnel,
      p.max_lossless_offer)
print([(t.entry_id, t.head, t.local_index) for t in r.program.tunnels], r.call_delivered)
```

```text
3/2 1 1 6
[(0, 0, 0), (1, 1, 0), (2, 0, 1), (3, 1, 1)] 6
```

Head 0's two tunnels share its link of 3, so the common size is at most 3/2 (head 1 alone would
allow 4): `s* = 3/2` and `U = 1`. The design 6 gives each head 3, so `u = min(1, ceil(3/2)) = 1`.
The tunnels are numbered in rounds: identifiers 0 and 1 are the two heads' first tunnels, 2 and 3
their second ones. `F = 6`: a total of 6 gives head 0 the 3 its link carries; 7 would give it 4.

A T2 call at a finer resolution: two branches from node 0 to node 3 of 2 and 10 (at q = 2^-12: 8192
and 40960 units), 16 tunnels, `prefer_headroom=True`, then the same program replayed with the 2-unit
branch's second link down:

```python
q = 2.0 ** -12
g = ngc.Graph(4, edge_sources=[0, 0, 1, 2], edge_destinations=[1, 2, 3, 3], edge_capacity_units=ngc.to_units([2, 10, 2, 10], q),
              edge_costs=[1, 1, 1, 1], resolution=q)
te16 = ngc.FlowPolicyConfig(routing=ngc.CentralizedTE(source_allocation=ngc.SourceAllocation.EQUAL,
                                                     tunnel_count=ngc.TunnelCount(16), prefer_headroom=True),
                            admission=ngc.Admission.LOSSY)
pol = ngc.FlowPolicy(ngc.FlowLedger(g), 0, 3, te16)
r = pol.replace(49152)
p = r.plan
print(p.max_bandwidth_per_tunnel, p.planned_bandwidth_per_tunnel, p.max_fractional_bandwidth_per_tunnel, p.max_lossless_offer, r.call_delivered)
print(p.total_planned_bandwidth, p.total_tunnel_cost, p.headroom_penalty, p.heads_with_capacity, p.insufficient_resolution)
print([t.edges for t in r.program.tunnels][:3], [t.planned_bandwidth for t in r.program.tunnels][:3])
fail = g.eligibility(excluded_edges=[2])
replay = ngc.FlowPolicy(ngc.FlowLedger(g), 0, 3,
                        ngc.FlowPolicyConfig(routing=ngc.Pinned(pol.program), admission=ngc.Admission.LOSSY),
                        eligibility=fail)
s = replay.replace(49152)
print([t.status.name for t in s.by_tunnel][:3], s.call_delivered)
```

```text
2925 2925 20480/7 46812 47104
46800 32 73610035200 (0,) False
[(0, 2), (0, 2), (1, 3)] [2925, 2925, 2925]
['FAILED', 'FAILED', 'UP'] 40960
```

With 2 tunnels on the 8192-unit branch and 14 on the 40960-unit one, the wide branch limits the
common size to 40960/14 = 20480/7; the narrow branch would allow 4096. The whole-unit size is 2925,
and the plan is 16 x 2925 = 46800. The lossless bound is 46812 = 16 x 2925 + 12: twelve tunnels
take one unit more, and the wide branch still carries exactly 40960. The offer 49152 gives every
tunnel 3072, more than its plan: the narrow branch carries its 6144, and the wide branch is
offered 43008, carries 40960 and drops 2048, so 6144 + 40960 = 47104 arrive. In the replay the two
tunnels on the failed branch are `FAILED` and the head splits its whole share over the other 14.

A dead head, a head too weak for one unit per tunnel, and "up to":

```python
lossy_t2 = lambda tunnel_count: ngc.FlowPolicyConfig(
    routing=ngc.CentralizedTE(source_allocation=ngc.SourceAllocation.EQUAL, tunnel_count=tunnel_count),
    admission=ngc.Admission.LOSSY)
dead = graph(3, [(0, 2, 0, 1), (1, 2, 10, 1)])
r = ngc.FlowPolicy(ngc.FlowLedger(dead), [0, 1], 2, lossy_t2(ngc.TunnelCount(2))).replace(10)
print(r.plan.heads_with_capacity, [t.entry_id for t in r.program.tunnels], r.plan.planned_bandwidth_per_tunnel, r.call_delivered, r.call_source_drops)
tiny = graph(3, [(0, 2, 1, 1), (1, 2, 100, 1)])
r = ngc.FlowPolicy(ngc.FlowLedger(tiny), [0, 1], 2, lossy_t2(ngc.TunnelCount(16))).replace(200)
print(r.plan.insufficient_resolution, r.plan.max_fractional_bandwidth_per_tunnel, r.program.tunnels, r.call_source_drops)
r = ngc.FlowPolicy(ngc.FlowLedger(tiny), [0, 1], 2, lossy_t2(ngc.TunnelCount(3, up_to=True))).replace(200)
print(r.plan.tunnels_per_head, [t.edges for t in r.program.tunnels], r.call_delivered)
```

```text
(1,) [1, 3] 3 5 5
True 1/16 () 200
1 [(0,), (1,)] 101
```

Head 0's only link has no free capacity: it is dead, keeps its share of 5 and loses it, and its
positions 0 and 2 are holes. Behind links of 1 and 100, 16 tunnels per head fit no whole unit
(`s*` = 1/16): the program has no tunnels and the whole offer is lost at the heads. Up to 3, one
tunnel per head fits a unit; it carries 1 of head 0's 100 and all of head 1's.

"Up to" compares the counts by what their programs carry of the design. Head 0 below reaches node
2 over one link of 6, head 1 over two parallel links of 3, and the design is 12. With one tunnel
per head, head 1's tunnel takes one link of 3, so the equal shares fit up to 7 (4 and 3), below the
design. With two, head 0's tunnels share its 6 and head 1's take one link each, so 12 fits, and
"up to 2" chooses 2:

```python
g = graph(3, [(0, 2, 6, 1), (1, 2, 3, 1), (1, 2, 3, 1)])
for limit in (1, 2):
    up_to = ngc.FlowPolicyConfig(routing=ngc.CentralizedTE(source_allocation=ngc.SourceAllocation.EQUAL,
                                                           tunnel_count=ngc.TunnelCount(limit, up_to=True)))
    r = ngc.FlowPolicy(ngc.FlowLedger(g), [0, 1], 2, up_to).replace(12)
    print(limit, r.plan.tunnels_per_head, r.call_delivered)
```

```text
1 1 7
2 2 12
```

The headroom objective on two equal-cost parallel links of 10, two tunnels for a design of 2:

```python
par = graph(2, [(0, 1, 10, 1), (0, 1, 10, 1)])
for on in (False, True):
    cfg = ngc.FlowPolicyConfig(routing=ngc.CentralizedTE(source_allocation=ngc.SourceAllocation.EQUAL,
                                                         tunnel_count=ngc.TunnelCount(2), prefer_headroom=on),
                               admission=ngc.Admission.LOSSY)
    r = ngc.FlowPolicy(ngc.FlowLedger(par), 0, 1, cfg).replace(12, design_demand=2)
    print(on, [t.edges for t in r.program.tunnels], r.plan.headroom_penalty, r.plan.max_lossless_offer, r.call_delivered)
```

```text
False [(0,), (0,)] None 10 10
True [(0,), (1,)] 858993458 20 12
```

Both plans size the tunnels at 1 unit (U = 10). Without the objective the cost tie goes to the
smallest keys: both tunnels on link 0, which then bounds forwarding at 10. With it, a second slot
on link 0 would pay `floor(2^32 * 2 / 10)` against `floor(2^32 / 10)` on link 1, so the tunnels
spread and the program forwards 20.

T2 on a support: three parallel links from node 0 to node 1 of (capacity, cost) (6, 1), (6, 1) and
(100, 2); the support holds the first two. The cheap pair carries one tunnel each and the third
link is never used; with the second link masked, both tunnels share the first:

```python
sup = graph(2, [(0, 1, 6, 1), (0, 1, 6, 1), (0, 1, 100, 2)])
support = ngc.PredecessorDAG.from_csr(sup, [0, 0, 2], [0, 0], [0, 1])
on_support = ngc.CentralizedTE(source_allocation=ngc.SourceAllocation.EQUAL, support=support,
                               tunnel_count=ngc.TunnelCount(2))
for mask in (None, sup.eligibility(excluded_edges=[1])):
    pol = ngc.FlowPolicy(ngc.FlowLedger(sup), 0, 1, ngc.FlowPolicyConfig(routing=on_support), mask)
    r = pol.replace(float("inf"))
    print(r.plan.max_bandwidth_per_tunnel, [t.edges for t in r.program.tunnels], r.plan.max_lossless_offer)
```

```text
6 [(0,), (1,)] 12
3 [(0,), (0,)] 6
```

The plan describes the view it was computed on. Below, T2 plans its two tunnels on one link of 10
for an unlimited design while it carries nothing yet; another demand then takes 8 of the link.
`forward(10)` delivers what the link still has, 2, and returns the installed plan unchanged:

```python
g = graph(2, [(0, 1, 10, 1)])
ledger = ngc.FlowLedger(g)
t2 = ngc.FlowPolicy(ledger, 0, 1, ngc.FlowPolicyConfig(routing=ngc.CentralizedTE(
    source_allocation=ngc.SourceAllocation.EQUAL, tunnel_count=ngc.TunnelCount(2))))
print(t2.replace(0, design_demand=float("inf")).plan.max_lossless_offer)
ngc.FlowPolicy(ledger, 0, 1, config(ngc.ShortestPaths())).place(8)
r = t2.forward(10)
print(r.call_delivered, r.plan.max_lossless_offer)
```

```text
10
2 10
```

### Several sources and destinations

- **A route ends at the first destination it reaches.** A walk that meets a destination stops
  there and never passes through it toward another destination. Passing through another source of
  the same demand is allowed and adds no traffic.
- **Destinations** (`ShortestPaths`): a route goes to the nearest destinations, and one record may
  deliver to several of them at that cost (a T1 tunnel is one path, to the allowed destination its
  CSPF reaches first). A record's `delivered` is its total; the delivery to one destination is the
  sum of the record's `edge_flows` entering it.
- **`ShortestPaths` split routing pools its sources.** Each record (one per roster source) holds the
  part of the program's traffic that is that source's. At every router the members' amounts are
  given back to the sources member by member, in member order, by highest averages over what each
  source still holds there; on an edge that drops traffic (Lossy), what the edge carries is
  divided by highest averages over what each source brought to it
  ([splitting whole units](#splitting-whole-units)). So the records sum exactly to the totals, and
  a source's record uses only edges of its own route.
- **The controllers** (`CentralizedTE`, `ShortestPaths` with `OPTIMIZED`) optimise the whole
  roster and destination set together; `CentralizedTE` may use a costlier destination to carry
  more.
- **T1 and T2 heads** share the design and the offer by round robin over the roster in ascending
  node id; a head's share goes only to its own tunnels, and what one head cannot carry is lost at
  that head, never moved to another.
- **`Pinned`** follows its program: its injection over the replay roster, then each head's entries.

Two sources with one unit each meet at node 2, whose two members lead to the destination. Pooled,
node 2 splits the two units once, one per member, and gives one member to each source; split
separately, each source would send its unit to the first member:

```python
# Edges 0: 0 -> 2, 1: 1 -> 2, then the two members 2: 2 -> 3 and 3: 2 -> 3.
g = graph(4, [(0, 2, 5, 1), (1, 2, 5, 1), (2, 3, 1, 1), (2, 3, 1, 1)])
pol = ngc.FlowPolicy(ngc.FlowLedger(g), [0, 1], 3, config(ngc.ShortestPaths()))
r = pol.place(2)
print(r.call_delivered, r.termination.name, [(row.source, row.delivered) for row in r.by_source])
print([(rec.sources, rec.edge_flows) for rec in pol.records()])
```

```text
2 COMPLETE [(0, 1), (1, 1)]
[([0], [(0, 1), (2, 1)]), ([1], [(1, 1), (3, 1)])]
```

### Equal source allocation and a source that cannot pass anything

A roster source whose every route crosses an edge with no free capacity stays in the roster: it
receives its share and loses it. Lossless admits only totals that give it nothing, which is the
number of roster sources before it, because the first remainder units go to the lowest node ids: 0
when it is the first source, 1 when it is the second (`RR(1, 2) = (1, 0)` gives it nothing,
`RR(2, 2) = (1, 1)` gives it a unit). This holds in every model. Lossy loses its share where its program
puts it: on the edge that limits it (`edge_drops`), or at the source when it has no member to send
it on (`source_drops`). An Eligibility that excludes the source or the edge removes it from the
roster instead.

```python
# Source 0 reaches node 2 only over the 0-unit edge 0; source 1 over the 10-unit edge 1.
g = graph(3, [(0, 2, 0, 1), (1, 2, 10, 1)])
equal = ngc.FlowPolicyConfig(routing=ngc.CentralizedTE(source_allocation=ngc.SourceAllocation.EQUAL))
optimized = ngc.FlowPolicyConfig(routing=ngc.CentralizedTE(source_allocation=ngc.SourceAllocation.OPTIMIZED))
for name, cfg in (("ShortestPaths", config(ngc.ShortestPaths())),
                  ("CentralizedTE EQUAL", equal),
                  ("CentralizedTE OPTIMIZED", optimized)):
    r = ngc.FlowPolicy(ngc.FlowLedger(g), [0, 1], 2, cfg).place(10)
    print(name, r.call_delivered, r.termination.name)
lossy = config(ngc.ShortestPaths(), admission=ngc.Admission.LOSSY)
r = ngc.FlowPolicy(ngc.FlowLedger(g), [0, 1], 2, lossy).place(10)
print("LOSSY", r.call_delivered, r.call_edge_drops, [(row.source, row.total_edge_drops) for row in r.by_source])
r = ngc.FlowPolicy(ngc.FlowLedger(g), [0, 1], 2, config(ngc.ShortestPaths()),
                   g.eligibility(excluded_nodes=[0])).place(10)
print("source 0 excluded", r.call_delivered)
g = graph(3, [(0, 2, 10, 1), (1, 2, 0, 1)])   # now source 1 is the one behind the 0-unit edge
print("source 1 blocked", ngc.FlowPolicy(ngc.FlowLedger(g), [0, 1], 2, config(ngc.ShortestPaths())).place(10).call_delivered)
```

```text
ShortestPaths 0 EXHAUSTED
CentralizedTE EQUAL 0 EXHAUSTED
CentralizedTE OPTIMIZED 10 COMPLETE
LOSSY 5 [(0, 5)] [(0, 5), (1, 0)]
source 0 excluded 10
source 1 blocked 1
```

### Forwarding programs

**In plain words.** A program is what a routing decided, written down: which routes each source
uses and in which integer proportions. Every computing model produces one (`policy.program`,
`result.program`). `Pinned(program)` forwards a program without computing anything: to replay a
solution on another ledger, to study a failure (the routers drop the failed routes and split over
the survivors) or to study the moment before anyone reacts (the failed elements swallow what
reaches them). `ForwardingProgram.from_routes` builds a program from routes you give, each its own
entry; `ForwardingProgram.from_dag` builds one shared hop-by-hop program on a DAG you give, its
sources pooled exactly as `ShortestPaths` pools them.

**Exact rules.**

- **The value.** A `ForwardingProgram` is immutable and shared by reference: `sources` and
  `destinations` (the normalised input sets), `tunnels` (the tunnel entries in ascending
  `entry_id`; split routes are not listed) and `source_weights` (`None` for equal sources). A
  `Tunnel` has `entry_id`, `head`, `local_index` (its position among its head's entries of both
  kinds), `destination`, `edges` (in path order), `path_cost` (these three `None` for a tunnel
  T1's setup left down), `status` (`UP`, or a down tunnel's `INSUFFICIENT_BANDWIDTH` or
  `COST_CEILING_EXCEEDED`), `requested_bandwidth` and `reserved_bandwidth` (T1's request and what
  it holds on every edge of its path: the request when up, else 0; 0 for the other producers),
  `planned_bandwidth` (T2's planned bandwidth u, never held; 0 for the other producers) and
  `weight`. It has no public constructor.
- **What a program shows.** A program can be replayed without exporting its tables. The public
  value does not expose per-router effective weights, the split tables or a reachable-source list
  (its `sources` are the normalised input sources); a split program may have no tunnels while it
  carries traffic. Records and `OutputOptions` give the carried and dropped amounts per edge.
- **Injection.** A total `x` enters at the roster (the heads with an entry; for T1 and T2 the
  policy's roster, a head without an up tunnel or without any tunnel included) either by round
  robin (**equal sources**) or by highest averages over the captured source weights; each head
  then splits its share by highest averages over its entries' weights (equal weights give round
  robin). A head whose entries all have weight 0 loses its share at the source. Identical routes
  are never merged: each entry takes its own share. `ShortestPaths` split routing, `from_dag`, T1,
  T2, equal controllers and supplied routes use equal sources; optimized controllers use their
  solved supplies.
- **Entry identifiers.** Every entry has an integer `entry_id` that identifies it within one
  program (not across recomputations): T1 and T2 tunnels are numbered `j * H + i` (round `j`,
  roster position `i`, `H` heads; T2's dead heads leave holes); controller paths are numbered in
  decomposition order; a `ShortestPaths` split program and a `from_dag` program number their routes
  by roster position; supplied routes keep their position in the `routes` list (a split route uses
  its position too, so tunnel identifiers can skip it). Records (`FlowRecord.entry_id`) and
  `by_tunnel` rows (`TunnelFlow.entry_id`) use the same identifier.
- **Split contexts.** A split route given to `from_routes` is its own context: its routers split
  only that route's share, even where another split route reaches the same router with the same
  members. A `ShortestPaths` `EQUAL` or `WEIGHTED` program and a `from_dag` program pool all
  their sources in one context. A controller's program (`ShortestPaths` `OPTIMIZED` included) is
  weighted path tunnels: its sources were solved jointly, and forwarding them splits no pooled
  router.
- **`ForwardingProgram.from_routes(graph, sources, destinations, routes, *, allocation=None,
  weight_basis=None, edge_weights=None, residual=None)`** with `RouteDefinition(head, route, *,
  weight=1)` (`route` a `Path` with one edge per hop, a tunnel, or a `PredecessorDAG`, a split
  route). Order and duplicates are kept; the roster is the sources with a route. `allocation` is
  the interior split of the `PredecessorDAG` routes, `EQUAL` (`None`) or `WEIGHTED`, whose weights
  `weight_basis`, `edge_weights` and `residual` choose as in
  [section 5](#5-route-allocations-and-admission) (`residual` is the free capacity
  `RESIDUAL_SNAPSHOT` reads; it has no default). A split route keeps the edges of its DAG on the
  walks from its head to a destination. A route of weight 0 is validated in full and carries
  nothing. The first failure, in this order, raises `ValueError`:
  1. the endpoint sets;
  2. `routes must not be empty`;
  3. when no route is a `PredecessorDAG`, any option given: `<allocation, weight_basis,
     edge_weights or residual> applies only to PredecessorDAG routes`;
  4. `OPTIMIZED`: `from_routes allocation must be EQUAL or WEIGHTED`;
  5. the weights' combination: `weight_basis requires allocation=WEIGHTED`,
     `weight_basis=EXPLICIT requires edge_weights`, `edge_weights require weight_basis=EXPLICIT`,
     `weight_basis=RESIDUAL_SNAPSHOT requires residual`, `residual requires
     weight_basis=RESIDUAL_SNAPSHOT`;
  6. `edge_weights`' length and entries, then `residual`'s (the ledger's baseline messages);
  7. each route in list order: `routes[<i>]: route belongs to another graph`, `routes[<i>]: head
     <h> is not a source`; for a `Path`, `routes[<i>]: a tunnel has one edge per hop`,
     `routes[<i>]: the path starts at <v>, not at head <h>`, `routes[<i>]: the path ends at <v>,
     not at a destination`, `routes[<i>]: the path meets destination <v> before its end`; for a
     `PredecessorDAG`, `routes[<i>]: no walk from head <h> to a destination`; then
     `routes[<i>]: weight must be an integer in [0, 2^62)`.
- **The lossless bound** of a program is the one of
  [section 5](#5-route-allocations-and-admission): each entry's share forwarded without drops,
  the loads added per edge, the largest total that leaves no share without a route and puts no
  edge above its free capacity.

#### Supplied DAGs: from_dag

`ForwardingProgram.from_dag(graph, sources, destinations, dag, *, allocation=None,
weight_basis=None, edge_weights=None, residual=None)` builds one shared hop-by-hop program on a
`PredecessorDAG` of the graph, trimmed to the walks from each source to its first destinations.
It is the mechanism `ShortestPaths` `EQUAL`/`WEIGHTED` uses for its own routes (one builder):
sources are allocated equally and pooled at every shared router with exactly the same splitting
and attribution, in one split context. `allocation` is `EQUAL` (`None`) or `WEIGHTED`; the weights
and the snapshot follow [section 5](#5-route-allocations-and-admission), and `RESIDUAL_SNAPSHOT`
needs an explicit `residual`. The roster is the sources that reach a destination in the DAG,
ignoring capacity and weights; an empty roster is allowed (an empty program). There is one split
route per roster source, identified by roster position, and no tunnel; a `Pinned` policy keeps one
record per route, attributed to its source. A replay rebuilds the roster on the DAG and keeps the
pooled context; stale forwarding keeps the tables and the roster. The first failure, in this
order, raises `ValueError`: the endpoint sets; `dag belongs to another graph`; the allocation
(`from_dag allocation must be EQUAL or WEIGHTED`) and the weights' combination as for
`from_routes`; `edge_weights`' and `residual`'s lengths and entries. A `dag` of another type is a
`TypeError` (`dag must be a PredecessorDAG (got <type>)`).

Two sources whose routes meet at node 2, which has two links of 1 to the destination. `from_dag`
pools the two arrivals and splits them once, one per link; two separate split routes of
`from_routes` each send their unit to the first link:

```python
g = graph(4, [(0, 2, 10, 1), (1, 2, 10, 1), (2, 3, 1, 1), (2, 3, 1, 1)])
dag = ngc.PredecessorDAG.from_csr(g, [0, 0, 0, 2, 4], [0, 1, 2, 2], [0, 1, 2, 3])
pooled = ngc.ForwardingProgram.from_dag(g, [0, 1], 3, dag)
separate = ngc.ForwardingProgram.from_routes(g, [0, 1], 3, [ngc.RouteDefinition(0, dag), ngc.RouteDefinition(1, dag)])
for name, prog in (("from_dag", pooled), ("from_routes", separate)):
    cfg = ngc.FlowPolicyConfig(routing=ngc.Pinned(prog))
    print(name, prog.tunnels, ngc.FlowPolicy(ngc.FlowLedger(g), [0, 1], 3, cfg).place(float("inf")).call_delivered)
sp = ngc.FlowPolicy(ngc.FlowLedger(g), [0, 1], 3, config(ngc.ShortestPaths())).place(float("inf"))
print("ShortestPaths", sp.call_delivered)
```

```text
from_dag () 2
from_routes () 1
ShortestPaths 2
```

#### Replay under a failure

A `Pinned` policy's Eligibility is the failure being studied (local repair):

- A tunnel that touches an excluded node or edge is `FAILED` and gets nothing; replay never
  computes a replacement path. A tunnel T1's setup left down stays down with its status: replay
  never brings one up.
- A split route loses its excluded members and the branches left without a destination; members
  of capacity 0 stay (`EQUAL` still sends them their share).
- The roster is rebuilt from the program's **basis**, together with the replaying policy's
  exclusions, ignoring capacities and weights:

  | Producer of the program | Basis |
  | --- | --- |
  | unrestricted `CentralizedTE`; T1; T2 without a support | the topology, under the producing policy's exclusions |
  | `ShortestPaths`; `from_dag`; `CentralizedTE` with a support (with or without a tunnel count); `from_routes` | the original routes or support |

- A head that no longer reaches a destination leaves the pool: with equal sources its share moves
  to the others. A head that still reaches one but has no surviving entry keeps its share and
  loses it.
- Captured source weights keep the remaining heads' saved weights. When all of them are 0, a
  positive offer is lost at the sources, and the call ends `NO_ACTIVE_WEIGHT` when the program has
  an entry with a path. A program without one (the program of a controller that sent nothing, such
  as an all-zero `WEIGHTED` one) ends `COMPLETE` under `LOSSY`, the offer counted in the demand's
  `source_drops` while every per-source row stays zero; under `LOSSLESS` it admits nothing and ends
  `EXHAUSTED`.
- An edge whose capacity is now 0 is not a failure: its traffic is lost on that edge, the tunnel
  stays `UP`. Weights stay as the program froze them.

#### Stale forwarding: black holes

`Pinned(program, black_holes=e)` and `forward(offer, black_holes=e)` forward the program unchanged
while the nodes and edges that `e` excludes swallow what reaches them (`None`, or an all-allowed
value, means none):

- **The program is unchanged.** Nothing is pruned, re-split or reweighted: every router keeps its
  members and weights and applies the program's whole-unit rule to the traffic that actually
  arrives.
- **Where traffic disappears.** A black-holed edge carries nothing: all that arrives is dropped
  there. A black-holed node acts the same on every edge entering it. A black-holed source loses its
  own share as a source drop.
- **`LOSSY`** drops what reaches a black hole where it arrives. That traffic takes no part in
  sharing an edge ([the loss rule](#losses-on-shared-edges)) and makes no edge wait, and the
  traffic behind it can carry nothing.
- **`LOSSLESS`** admits only totals that put nothing on a black hole (a black-holed edge serves 0).
- **Statuses.** A tunnel that is up but whose path reaches a black hole is `BLACKHOLED` in
  `by_tunnel`, its traffic dropped at the first black hole on its path. A down tunnel keeps its
  setup status.
- **With a mask.** The policy's Eligibility applies first: an element both excluded and
  black-holed is excluded.
- A program forwarded by `Pinned` books no hold, whatever its producer.

```python
# Two parallel edges of 4 and 8: the duplicate tunnel takes its own share.
g = graph(2, [(0, 1, 4, 1), (0, 1, 8, 1)])
def tunnel(g, head, edges, weight=1):
    return ngc.RouteDefinition(head, ngc.Path(g, [[e] for e in edges]), weight=weight)
def pinned(program, admission=ngc.Admission.LOSSLESS, black_holes=None):
    return ngc.FlowPolicyConfig(routing=ngc.Pinned(program, black_holes=black_holes),
                                admission=admission)
for routes in ([tunnel(g, 0, [0]), tunnel(g, 0, [1]), tunnel(g, 0, [1])],
               [tunnel(g, 0, [0]), tunnel(g, 0, [1])]):
    prog = ngc.ForwardingProgram.from_routes(g, 0, 1, routes)
    r = ngc.FlowPolicy(ngc.FlowLedger(g), 0, 1, pinned(prog)).place(float("inf"))
    print(r.call_delivered, [(t.entry_id, t.assigned) for t in r.by_tunnel])

# Two tunnels 0 -> 1 -> 3 and 0 -> 2 -> 3, every edge 10: node 2 fails.
g = graph(4, [(0, 1, 10, 1), (0, 2, 10, 1), (1, 3, 10, 1), (2, 3, 10, 1)])
prog = ngc.ForwardingProgram.from_routes(g, 0, 3, [tunnel(g, 0, [0, 2]), tunnel(g, 0, [1, 3])])
down = g.eligibility(excluded_nodes=[2])
lossy = pinned(prog, ngc.Admission.LOSSY)
r = ngc.FlowPolicy(ngc.FlowLedger(g), 0, 3, lossy, down).replace(10)
print("repair", r.call_delivered, [(t.status.name, t.assigned) for t in r.by_tunnel])
r = ngc.FlowPolicy(ngc.FlowLedger(g), 0, 3, pinned(prog, ngc.Admission.LOSSY, down)).replace(10)
print("stale", r.call_delivered, r.call_edge_drops, [(t.status.name, t.assigned) for t in r.by_tunnel])
print("stale lossless", ngc.FlowPolicy(ngc.FlowLedger(g), 0, 3, pinned(prog, black_holes=down))
      .place(float("inf")).call_delivered)
```

```text
12 [(0, 4), (1, 4), (2, 4)]
8 [(0, 4), (1, 4)]
repair 10 [('UP', 10), ('FAILED', 0)]
stale 5 [(1, 5)] [('UP', 5), ('BLACKHOLED', 5)]
stale lossless 1
```

The repaired head sends all 10 over its surviving tunnel. Stale, node 0 still splits 5 and 5 and
the second tunnel's 5 are dropped on edge 1, the edge entering node 2; Lossless admits 1 unit (2
units would put one on the black hole).

A program of a controller replays its solve; under a mask its failed tunnel's share moves to the
head's surviving tunnel:

```python
g = graph(4, [(0, 1, 2, 1), (0, 2, 10, 1), (1, 3, 2, 1), (2, 3, 10, 1)])
solved = ngc.FlowPolicy(ngc.FlowLedger(g), 0, 3, ngc.FlowPolicyConfig(routing=ngc.CentralizedTE(source_allocation=ngc.SourceAllocation.OPTIMIZED)))
prog = solved.replace(12).program
print([(t.entry_id, list(t.edges), t.weight) for t in prog.tunnels], prog.source_weights)
r = ngc.FlowPolicy(ngc.FlowLedger(g), 0, 3, pinned(prog, ngc.Admission.LOSSY),
                   g.eligibility(excluded_edges=[2])).replace(12)
print(r.call_delivered, r.call_edge_drops, [(t.status.name, t.assigned, t.delivered) for t in r.by_tunnel])
```

```text
[(0, [0, 2], 2), (1, [1, 3], 10)] ((0, 12),)
10 [(1, 2)] [('FAILED', 0, 0), ('UP', 12, 10)]
```

The first tunnel is `FAILED`; the head's 12 go to its other tunnel by its weight 10 alone, which
carries 10 and drops 2 on edge 1.

#### Failure studies

A failure has three states, all built from these primitives: **stale forwarding** (no router has
reacted yet: the black holes above), **local repair** (the routers drop the failed routes and split
over the survivors: the mask of a `Pinned` policy) and **converged routing** (a new computation on
the failed network: a new policy with that mask). Save the baseline program before the failure.
`evaluate_edge_failures` recomputes or repairs under each deletion
([section 10](#10-evaluators)); it never turns a deletion into a black hole, so a stale study of
every link is a loop over `Pinned(base, black_holes=...)`. A failure never changes the graph or
another policy's mask; on a shared ledger, the old owner is removed or replaced before its
replacement is modelled, since a new policy is a different demand that does not release the old
one's traffic or holds.

### Losses on shared edges

**In plain words.** Where `LOSSY` traffic meets on an edge that cannot carry it all, the edge
carries as much as it can and shares that amount out by highest averages: first between the
demands, by what each brings, then within a demand between its sources, then within a source
between its entries. Each part goes on with its share; the rest is dropped on that edge. Sharing
demand by demand first means that a demand gains nothing by spreading its traffic over more
tunnels.

Two demands of 8 and 4 meet on a link of 9: the link carries `T = min(9, 12) = 9`, and
`H(9, [8, 4]) = (6, 3)`, so they lose 2 and 1. Offers of 5, 1 and 1 on a link of 3 give
`H(3, [5, 1, 1]) = (3, 0, 0)`: the first offer's quotients 5, 5/2 and 5/3 all exceed the others'
1. On a link of 1, demand A with two tunnels of 1 unit and demand B with one tunnel of 2 tie at
the demand level, 2 against 2. Ties go to the demand created first, A, whatever the order of the
list, and within A to its first tunnel:

```python
def at_once(capacity, offers):
    """One link 0 -> 1; each offer is its own ECMP demand, all forwarded at the same time."""
    g = graph(2, [(0, 1, capacity, 1)])
    ledger = ngc.FlowLedger(g)
    demands = [ngc.FlowPolicy(ledger, 0, 1, config(ngc.ShortestPaths(), admission=ngc.Admission.LOSSY))
               for _ in offers]
    results = ngc.apply_actions(ledger, [ngc.Replace(d, x) for d, x in zip(demands, offers)],
                                mode=ngc.BatchMode.SIMULTANEOUS)
    return [r.call_delivered for r in results]

print(at_once(9, [8, 4]), at_once(3, [5, 1, 1]))
g = graph(2, [(0, 1, 1, 1)])
ledger = ngc.FlowLedger(g)
a_program = ngc.ForwardingProgram.from_routes(g, 0, 1, [tunnel(g, 0, [0]), tunnel(g, 0, [0])])
b_program = ngc.ForwardingProgram.from_routes(g, 0, 1, [tunnel(g, 0, [0])])
a = ngc.FlowPolicy(ledger, 0, 1, pinned(a_program, ngc.Admission.LOSSY))   # created first
b = ngc.FlowPolicy(ledger, 0, 1, pinned(b_program, ngc.Admission.LOSSY))
rb, ra = ngc.apply_actions(ledger, [ngc.Replace(b, 2), ngc.Replace(a, 2)], mode=ngc.BatchMode.SIMULTANEOUS)
print(ra.call_delivered, [t.delivered for t in ra.by_tunnel], rb.call_delivered)
```

```text
[6, 3] [3, 0, 0]
1 [1, 0] 0
```

**Exact rules.**

- **The local rule.** At an edge `e` that receives more than its free capacity `V[e]`, with
  arrivals `a[d, s, t]` by demand `d`, source `s` and entry `t`, the edge carries
  `T = min(V[e], sum of a)`. Each demand gets its share of `H(T, arrivals by demand)`, the demands
  taken in the order their policies were created in the ledger, never by list position (one
  program is one demand; in a `SIMULTANEOUS` list each policy is a demand). Within a demand its
  sources share by `H` over their arrivals, and within a source its entries the same way. Each
  part carries only its share onward; the rest is dropped on this edge. A pooled `ShortestPaths`
  demand has no entry level: it splits its total once at each router, and each member's amount
  is shared among the sources by highest averages over what each still brings.
- **The whole forwarding state** is evaluated at once:
  1. Every demand is forwarded with its fixed injections and splits, sharing no edge and dropping
     its traffic only at its own black holes, and the loads are added. Traffic that carries
     nothing in this forward is **dead**: an entry, a split member or a pooled source whose
     whole-unit share is 0, or traffic that can only come through its own demand's black hole.
     Dead traffic takes no part in what follows: it loads no edge, makes no edge wait and joins no
     series of edges below. So adding a demand that can carry nothing never changes another
     demand's result, error or round count.
  2. **Edges in series.** Two consecutive edges of equal free capacity that all the live traffic
     crosses together (nothing enters, ends, merges, branches or is re-attributed between them)
     act as one: in a chain of such edges only the first keeps its service, and the others pass
     traffic through and report no drops (the records keep every physical edge). So splitting an
     edge into two in series with the same capacity changes nothing (a black hole on it becomes
     one on its first part, or on both).
  3. An edge whose load fits its free capacity can never drop: it passes traffic through and
     waits for nothing. Every other shared edge waits for all the traffic that can reach it, and a
     router waits for what it receives.
  4. The edges are evaluated once each, in dependency order, by the local rule; the order of
     independent edges does not change the result.
  5. Edges that wait on each other in a cycle are solved as a group (below).

  Traffic that reaches its own demand's black hole is dropped there and takes no part in the
  sharing. This is a steady proportional-sharing convention, not a packet, buffer or
  retransmission simulation.
- **Cycles.** When edges wait on each other in a cycle (one tunnel meets edge X before Y, another
  Y before X, and both can drop; the tunnels may belong to different demands), only that group is
  recalculated, round by round. Round 1 starts from nothing carried. Every round re-injects the
  full offers, recomputes all the group's edges at once from what the previous round carried, and
  never offers dropped traffic again. Each round's carried amounts are compared with the earlier
  rounds', and the first case that applies ends the group:

  | After a round | Outcome |
  | --- | --- |
  | it carries exactly what the round before carried | the group is solved |
  | it carries what the round two before carried (the rounds alternate between two states) | `ValueError("lossy forwarding iteration entered a period-2 cycle")` |
  | the group has used 64 rounds per waiting edge in it | `ValueError("lossy forwarding iteration did not reach a fixed point within the round limit")` |

  A solved answer is checked independently against the rule before anything is published; a
  failed check, an implementation fault, raises `RuntimeError("lossy forwarding fixed-point
  verification failed")`. Each error publishes nothing and does not claim that no answer exists.

```python
# Two tunnels of one head share the second edge of 2 on tunnel 0 only.
g = graph(4, [(0, 1, 10, 1), (0, 2, 10, 1), (1, 3, 2, 1), (2, 3, 10, 1)])
prog = ngc.ForwardingProgram.from_routes(g, 0, 3, [tunnel(g, 0, [0, 2]), tunnel(g, 0, [1, 3])])
r = ngc.FlowPolicy(ngc.FlowLedger(g), 0, 3, pinned(prog, ngc.Admission.LOSSY)).place(12)
print(r.call_delivered, r.call_edge_drops, [(t.assigned, t.delivered, t.total_edge_drops) for t in r.by_tunnel])

# One edge of 9 shared by two tunnels of weights 2 and 1: the offer splits 8 and 4, the edge 6 and 3.
g = graph(2, [(0, 1, 9, 1)])
prog = ngc.ForwardingProgram.from_routes(g, 0, 1, [tunnel(g, 0, [0], 2), tunnel(g, 0, [0], 1)])
r = ngc.FlowPolicy(ngc.FlowLedger(g), 0, 1, pinned(prog, ngc.Admission.LOSSY)).place(12)
print(r.call_delivered, [(t.assigned, t.delivered, t.total_edge_drops) for t in r.by_tunnel])

# Two tunnels cross edges X (e0, 1 unit) and Y (e3, 2 units) in opposite orders, weights 1 and 2.
g = graph(6, [(0, 1, 1, 1), (1, 2, 1000, 1), (1, 5, 1000, 1), (2, 3, 2, 1), (3, 0, 1000, 1),
              (3, 5, 1000, 1), (4, 0, 1000, 1), (4, 2, 1000, 1)])
prog = ngc.ForwardingProgram.from_routes(
    g, 4, 5, [tunnel(g, 4, [6, 0, 1, 3, 5], 1), tunnel(g, 4, [7, 3, 4, 0, 2], 2)])
try:
    ngc.FlowPolicy(ngc.FlowLedger(g), 4, 5, pinned(prog, ngc.Admission.LOSSY)).place(3)
except ValueError as error:
    print(error)
```

```text
8 [(2, 4)] [(6, 2, 4), (6, 6, 0)]
9 [(8, 6, 2), (4, 3, 1)]
lossy forwarding iteration entered a period-2 cycle
```

In the first program the second edge of 2 on tunnel 0 drops 4 of its 6. In the cycle, tunnel 0
sends 1 unit over X (edge 0, capacity 1) and then Y (edge 3, capacity 2), and tunnel 1 sends 2
units over Y and then X. What tunnel 0 keeps on X decides what reaches Y, and what tunnel 1 keeps
on Y decides what reaches X. Shares are by highest averages over the arrivals, ties to tunnel 0:

| Round | At X: arrivals of tunnels 0, 1 | X carries | At Y: arrivals of tunnels 0, 1 | Y carries |
| --- | --- | --- | --- | --- |
| 1 | 1, 0 | 1, 0 | 0, 2 | 0, 2 |
| 2 | 1, 2 | 0, 1 | 1, 2 | 1, 1 |
| 3 | 1, 1 | 1, 0 | 0, 2 | 0, 2 |

Round 3 repeats round 1, so the call raises the period-2 error and publishes nothing.

### The placement call

`place(offer, outputs=None)` takes an offer in units (an `int`) or `float('inf')`.

1. **Validation** happens before anything changes. Within an epoch, a finite offer after an
   infinite one is rejected (`finite placement after infinite placement requires a new epoch`, the
   cumulative request would be undefined), `place(0)` included.
2. **An offer of 0** is `COMPLETE` and does nothing: it creates no flow and prepares nothing.
   Before any program exists (a new policy, or after `remove()`) its result is all zeros, every
   flag false; with a program installed it reports the demand's current state with every flag
   false. It is not a removal: `remove()` is.
3. **`ShortestPaths` split routing and `Pinned`.** `ShortestPaths` prepares the program at the
   epoch's first positive call (one record per roster source, also at zero delivery); `Pinned`
   creates its records, one per entry, there. The two admissions count the new offer differently:
   - **Lossless** re-solves the whole target: it admits the largest total `a` between what it
     already delivers and that plus the offer that fits on its service view
     ([the lossless bound](#5-route-allocations-and-admission)), forwarded from zero. A unit
     refused by an earlier call is not offered again: the next call targets what is delivered
     plus its own offer.
   - **Lossy** sends only the new units: the sources inject this call's increase of the
     injection, `RR(offered + v) - RR(offered)`, and it is forwarded hop by hop on the free
     capacity; each edge carries the smaller of what reaches it and its residual and drops the
     rest, and old drops are never offered again. With captured source weights the increase is
     `H(offered + v, w) - H(offered, w)`, so after any number of calls the sources have received
     exactly `H(total offered, w)`.

   An infinite Lossy offer is rejected for `Pinned` (`LOSSY forwarding requires a finite offer`);
   for `ShortestPaths` it is the first-hop experiment of
   [section 5](#5-route-allocations-and-admission).
4. **The controllers** re-solve as described above.
5. **T1 and T2** have no `place`: every call raises `DistributedTE does not support incremental
   placement` or `CentralizedTE with tunnel_count does not support incremental placement`, before
   anything else (an offer of 0 too); their calls are `replace`, `evaluate` and `forward`.
6. **Commit.** All changes to the ledger and the policy are committed together at the end, or
   none are ([section 8](#8-failures-threads-and-lifetimes)).

What the call guarantees:

- `call_delivered` is the increase in delivery made by this call, and `0 <= call_delivered <=
  offer`. For Lossy it is this call's delivery. `call_undelivered = offer - call_delivered`
  (`float('inf')` for an infinite offer; for Lossy it includes what was lost); it is never a
  backlog to retry.
- A `place` call never reduces `delivered()`, the sum of `delivered` over the policy's records.
- The first positive offer creates the policy's records even when none of them can deliver
  anything yet. Such a zero-delivery record keeps its FlowId for the rest of the epoch and can fill
  on a later call.
- **Additivity.** With the rest of the network unchanged between the calls, `place(a); place(b)`
  equals `place(a + b)` exactly for Lossless `ShortestPaths` split routing, the controllers (value,
  allocation and cost curve) and Lossless `Pinned`. Lossy calls are not additive: each forwards its
  own injection on the free capacity the earlier calls left, and nothing is offered again.

Two sources behind private links of 1 and 3, and calls of 3 and then 1. Lossless: the first call
tries 3, split (2, 1), but source 0's link carries only 1, so it admits 2, split (1, 1). The
second call targets what is delivered plus its offer, 2 + 1 = 3 again, which still does not fit,
and delivers 0. Lossy: the first call sends (2, 1) and drops 1 on source 0's link; the second sends
only the new units, `RR(4, 2) - RR(3, 2) = (2, 2) - (2, 1) = (0, 1)`, and delivers that 1:

```python
g = graph(3, [(0, 2, 1, 1), (1, 2, 3, 1)])
for admission in (ngc.Admission.LOSSLESS, ngc.Admission.LOSSY):
    pol = ngc.FlowPolicy(ngc.FlowLedger(g), [0, 1], 2, config(ngc.ShortestPaths(), admission=admission))
    calls = [pol.place(v) for v in (3, 1)]
    print(admission.name, [(r.call_delivered, r.call_edge_drops) for r in calls], calls[-1].delivered)
```

```text
LOSSLESS [(2, []), (0, [])] 2
LOSSY [(2, [(0, 1)]), (1, [])] 3
```

### Epochs

`remove()` takes every flow of the policy out of the ledger (allocations reverted, records gone,
ids retired), clears the installed program and starts a new epoch: the cumulative request, the
round-robin position of the source split and the finite-or-infinite state are reset. The next
positive offer prepares again and creates new records with new ids; ids are never reused. After
`remove()` the policy's `program` reads `None` until the next preparation. A repeated `remove()`
does nothing. A policy destroyed without `remove()` leaves its traffic and holds in the ledger as
orphans until `FlowLedger.reset()`; nothing is redistributed.
`replace` and `forward` also start a new epoch (below).

```python
g = graph(4, [(0, 1, 10, 1), (1, 3, 10, 1), (0, 2, 2, 1), (2, 3, 2, 1)])
ledger = ngc.FlowLedger(g)
pol = ngc.FlowPolicy(ledger, [0, 1], 3, config(ngc.ShortestPaths()))
for v in (1, 1, 0):
    r = pol.place(v)
    print(v, r.call_delivered, r.routes_changed, [(x.flow_id[1], x.delivered) for x in pol.records()])
pol.remove()
r = pol.place(1)
print(1, r.call_delivered, r.routes_changed, [(x.flow_id[1], x.delivered) for x in pol.records()])
```

```text
1 1 True [(0, 1), (1, 0)]
1 1 False [(0, 1), (1, 1)]
0 0 False [(0, 1), (1, 1)]
1 1 True [(2, 1), (3, 0)]
```

The first unit goes to source 0 and the second to source 1 (round robin of the cumulative
offer); the first call creates both records, source 1's at zero delivery. After `remove()` the
split starts again at source 0, under new ids.

### Operations: evaluate, replace and forward

Besides `place`, a policy has three full-offer calls. Each starts from the policy's **released
view**: its service view, the ledger's residual with the policy's own traffic given back (T1's
setup uses its setup view, with its own holds given back).

| Call | What it does | Routings |
| --- | --- | --- |
| `evaluate(offer, *, design_demand=None, outputs=None)` | computes and forwards `offer` from the released view as a fresh epoch would (T1: sets up on its setup view), and reports it, the computed program included; publishes nothing: no record, no hold, and the policy's program and epoch stay as they were | all |
| `replace(offer, *, design_demand=None, outputs=None)` | the same computation, then published at once: the old records and holds are retired, new ones created, the new program installed, a new epoch begun with `requested = offer`; delivery may go down | all |
| `forward(offer, *, black_holes=None, outputs=None)` | forwards a full new offer over the installed program with no route computation and no setup (weights stay frozen, a controller does not re-solve, T1 keeps its tunnels and holds, T2 its tunnels and plan), on the released view, then publishes like `replace`; `black_holes` forwards it stale for this call | every computing model with an installed program; `Pinned` uses `replace` |

- `call_delivered` of a full-offer call is its whole new delivery and `delivery_delta` the new delivery
  minus the old (it may be negative). A later `place` continues the epoch a `replace` or a
  `forward` began (its records keep their ids).
- With other owners fixed, identical full inputs give identical results, except fresh record ids.
- **`design_demand`** belongs to T1 and T2 alone: an integer in `[0, 2^62)` or `float('inf')`;
  `None` means equal to the offer. T1 needs it finite, also for an infinite offer; T2 takes an
  unlimited one (its largest common size). It never accumulates. The texts of these checks and of
  `forward`'s are in [configuration and call errors](#configuration-and-call-errors).
- A program reused elsewhere, or read again under a different mask, goes through a new `Pinned`
  policy; a program forwarded by `Pinned` never books holds, whatever its producer (read a
  result's `edge_reservations` or the ledger for actual holds).

```python
g = graph(2, [(0, 1, 10, 1), (0, 1, 10, 1)])
ledger = ngc.FlowLedger(g)
pol = ngc.FlowPolicy(ledger, 0, 1, config(ngc.ShortestPaths()))
r = pol.replace(6)
print(r.call_delivered, [x.flow_id[1] for x in pol.records()])
r = pol.evaluate(12)
print(r.call_delivered, r.delivery_delta, ledger.carried().tolist())
r = pol.replace(4)
print(r.call_delivered, r.delivery_delta, [x.flow_id[1] for x in pol.records()])
r = pol.forward(7)
print(r.call_delivered, r.routes_changed, r.weights_changed, ledger.carried().tolist())
```

```text
6 [0]
12 6 [3, 3]
4 -2 [1]
7 False False [4, 3]
```

### Termination

The call's `termination` (of `place` and of the full-offer calls alike) is the first row that
applies (an invalid call raises instead):

| Termination | When |
| --- | --- |
| `COMPLETE` | the offer is 0 |
| `NO_ROUTE` | the offer is positive and the roster is empty (`Pinned`: the replay roster) |
| `NO_ACTIVE_WEIGHT` | the program has an entry with a path (an up tunnel, a controller path or a split route), but nothing can be sent anywhere: every `WEIGHTED` member has weight 0, or no head with a positive injection weight has a surviving entry of positive weight (routes of weight 0, every up tunnel failed). A program whose tunnels are all down (T1) or that has no tunnels (T2 with every head dead, or U = 0) is not this case: it goes on to the rows below, every share lost at its head |
| `NO_ACTIVE_WEIGHT` | a live `WEIGHTED` controller whose roster weights are all 0: nothing enters, nothing is delivered and a finite offer is wholly unadmitted; its program has captured weights that are all 0 and no tunnel |
| `COMPLETE` | Lossy with a finite offer: the whole offer was forwarded, and delivery and losses account for all of it (it does not mean nothing was lost) |
| `COMPLETE` | Lossless: the whole offer was admitted |
| `EXHAUSTED` | otherwise: Lossless below the offer, any infinite offer (the first-hop experiment included), a controller whose maximum is 0 over existing routes |

`call_undelivered` is an amount, never a status: a Lossy call can be `COMPLETE` with a positive
`call_undelivered`. A `Pinned` replay of the all-zero `WEIGHTED` controller's program follows the
rows of a fixed program instead ([replay](#replay-under-a-failure)): Lossy `COMPLETE` with the
whole offer lost at the sources, Lossless `EXHAUSTED`. That is the one case in which a live call
and its replay end differently.

```python
def show(r):
    print(r.call_delivered, r.termination.name, [(row.source, row.delivered) for row in r.by_source])

# Every tunnel down: five edges, the design 14 asks 7 per tunnel and no path has 7.
g = graph(4, [(0, 1, 6, 1), (0, 2, 6, 3), (1, 2, 4, 1), (1, 3, 6, 3), (2, 3, 6, 1)])
for admission in (ngc.Admission.LOSSY, ngc.Admission.LOSSLESS):
    all_down = ngc.FlowPolicyConfig(routing=ngc.DistributedTE(ngc.TunnelCount(2)), admission=admission)
    r = ngc.FlowPolicy(ngc.FlowLedger(g), 0, 3, all_down).replace(14, design_demand=14)
    print(admission.name, [t.status.name for t in r.by_tunnel], r.call_source_drops, end=" ")
    show(r)

# One edge for source 1 to node 4.
g = graph(5, [(0, 2, 1, 1), (0, 3, 1, 1), (2, 4, 1, 1), (3, 4, 1, 1), (1, 4, 1, 1)])

# Weights from the residual where another demand holds the only edge: no positive weight.
ledger = ngc.FlowLedger(g)
ngc.FlowPolicy(ledger, 1, 4, config(ngc.ShortestPaths())).place(1)
weighted = ngc.FlowPolicyConfig(routing=ngc.ShortestPaths(allocation=ngc.RouteAllocation.WEIGHTED,
                                                          weight_basis=ngc.WeightBasis.RESIDUAL_SNAPSHOT))
show(ngc.FlowPolicy(ledger, 1, 4, weighted).place(1))
```

```text
LOSSY ['INSUFFICIENT_BANDWIDTH', 'INSUFFICIENT_BANDWIDTH'] 14 0 COMPLETE [(0, 0)]
LOSSLESS ['INSUFFICIENT_BANDWIDTH', 'INSUFFICIENT_BANDWIDTH'] 0 0 EXHAUSTED [(0, 0)]
0 NO_ACTIVE_WEIGHT [(1, 0)]
```

Both of T1's tunnels are down, so the program has no entry with a path: Lossy forwards the whole
offer and loses it at the head (`COMPLETE`), Lossless admits nothing (`EXHAUSTED`). The Weighted
route has an entry with a path whose only member has weight 0: `NO_ACTIVE_WEIGHT`.

### Flags

`routes_changed`, `weights_changed` and `reservations_changed` start false in every call and
record only changes that were committed; a call that raises publishes nothing (`evaluate` reports
what `replace` would set). They compare the installed program before and after the call (no
program counts as an empty one): `routes_changed` the ordered heads, each entry's kind and route
(a tunnel's path, a split route's members), duplicates included; `weights_changed` the injection
(equal or the captured weights, by value) and every entry and split-table weight, an entry absent
from one program counting as weight 0; `reservations_changed` every tunnel's `requested_bandwidth`
and `reserved_bandwidth` by identifier, an absent tunnel counting as 0 (only T1 sets it: T2's planned
bandwidth is no reservation, and a T2 layout over the same routes at another planned size sets
no flag). Hence:

| Committed change | `routes_changed` | `weights_changed` |
| --- | --- | --- |
| the epoch's first nonempty program (`ShortestPaths` and the controllers' first positive call, a `replace`, `replace(0)` included) | true | true when some entry has a positive weight (a T1 program whose tunnels are all down: false) |
| a T1 setup whose paths or up and down tunnels differ from the installed one's (requests changed too: `reservations_changed`) | true | true when a tunnel came up or went down |
| an empty program over no program or an empty one (a call with an empty roster, `replace(0)` included) | false | false |
| a controller's re-solve that moves a source's share to other edges, also with `call_delivered == 0` | true | true unless the entries keep their number and amounts |
| a controller's re-solve that changes only a path's amount (`CentralizedTE` `place(5)` then `place(5)` on one path) | false | true |
| an offer-only call over an unchanged program: `forward`, `Pinned`, `ShortestPaths` increments, an identical `replace` (new ids) | false | false |
| nothing of the above: `place(0)` (an offer of 0, which prepares and installs nothing, unlike a full zero offer such as `replace(0)`, which follows the rows above), rejected candidates only | unchanged | unchanged |

`reservations_changed` is true exactly when some tunnel's requested or reserved amount changed by
identifier: a T1 `replace` (or `evaluate`) whose design or setup differs, even over unchanged
paths; it is never set by `forward`, by `Pinned` or by the other routings. The flags describe
routes, weights and holds only: delivered and dropped amounts can change with all three false.

```python
# A first offer on a route whose only unit another demand holds.
g = graph(5, [(0, 2, 1, 1), (0, 3, 1, 1), (2, 4, 1, 1), (3, 4, 1, 1), (1, 4, 1, 1)])
ledger = ngc.FlowLedger(g)
ngc.FlowPolicy(ledger, 1, 4, config(ngc.ShortestPaths())).place(1)
pol = ngc.FlowPolicy(ledger, 1, 4, config(ngc.ShortestPaths()))
r = pol.place(1)
print(r.call_delivered, r.termination.name, r.routes_changed, r.weights_changed,
      [x.delivered for x in pol.records()])
```

```text
0 EXHAUSTED True True [0]
```

### The result

`place`, `evaluate`, `replace` and `forward` (and every evaluator) return a `PlacementResult`. Its
fields describe either this call or the demand's state after it (totals over the epoch); a
per-call loss field carries the `call_` prefix, a state field does not, so the two are never
confused.

- **This call**: `call_delivered`, `call_undelivered`, `termination`, `routes_changed`,
  `weights_changed` and `reservations_changed` as above; `delivery_delta`, the change of the
  delivery (for `place` it equals `call_delivered`); `call_edge_drops`, this call's drops per edge
  as `[(EdgeId, units)]`; `call_source_drops`, the units this call lost at sources.
- **The demand's state**: `program` (the installed `ForwardingProgram`, `evaluate`'s computed one;
  `None` before the epoch's first preparation and after `remove()`), `plan` (T1 and T2: the
  `TunnelPlan` of that program, [T1](#t1-tunnels-set-up-by-their-heads),
  [T2](#t2-the-controllers-equal-tunnels); `None` for every other routing, `Pinned` included),
  `requested` (the sum of the offers, `float('inf')` after an infinite one), `offered`,
  `admitted`, `delivered`, `unadmitted`, `source_drops`, `edge_drops` (per edge, as
  `[(EdgeId, units)]`), `by_source`, `by_tunnel` and `total_cost`. Every amount is an exact Python
  integer and may exceed 64 bits.
- **Accounting.** Lossless: `offered` is `None`, `source_drops` 0, no edge drops, `admitted =
  delivered`, and `requested = unadmitted + delivered` for a finite offer (`unadmitted` is `None`
  after an infinite one). Lossy: `offered` is the traffic actually sent, `admitted = offered`,
  `unadmitted = 0`, and `offered = source_drops + sum(edge_drops) + delivered`; `offered =
  requested` except in the first-hop experiment, where `requested` is infinite and `offered` the
  finite injection. `delivered + sum(edge_drops)` is what entered the network. `edge_flows` are
  carried loads per edge: summing them counts every hop, not the delivery.
- **`by_source`** is a tuple of `SourceFlow` rows, one per input source in ascending node id, a
  source outside the roster included (a zero row): `source`, `offered` (`None` when Lossless),
  `admitted`, `delivered`, `source_drops` and `total_edge_drops` (its share of `edge_drops`). The
  rows sum to the totals, with two exceptions, both Lossy with a positive offer: an empty roster,
  and a forwarding whose remaining captured source weights are all 0 (the program of a controller
  that sent nothing, such as an all-zero `WEIGHTED` one, or a mask that left only zero-weight
  heads). In both the whole offer is the demand's `source_drops`, lost before entering, while every
  row stays zero; read the demand's totals from the result's own fields, not from the rows.
- **`by_tunnel`** is a tuple of `TunnelFlow` rows, one per program tunnel in ascending `entry_id`
  (empty without tunnels): `entry_id`, `head`, `status` (`UP`; `INSUFFICIENT_BANDWIDTH` or
  `COST_CEILING_EXCEEDED` for a tunnel T1's setup left down; `FAILED` or `BLACKHOLED` in this
  forwarding), `assigned` (its injected share), `delivered` and `total_edge_drops`. It is another
  view of the same traffic as `by_source`, never added to it; a head's source drops are not
  attributed to one of its down tunnels.
- **`total_cost`**: the controllers only, the exact `sum(cost[e] * flow[e])` of the policy's
  allocation; `None` otherwise (T1, T2 and `Pinned` included).
- **`marginal_cost_curve` and `min_cut`** exist only for unrestricted `CentralizedTE` with
  `OPTIMIZED` source allocation. For every other model, for `EQUAL` or `WEIGHTED` sources and with
  a supplied support both are `None`: the marginal costs of a flow with prescribed shares need not
  increase, and such a flow is not the unrestricted maximum a cut certifies. For that controller
  the call decides what is reported:

  | Call | `marginal_cost_curve` | `min_cut` (when requested) |
  | --- | --- | --- |
  | a positive `place`, an `evaluate` or a `replace`, a zero offer included | the curve of the call's solve: of the policy's cumulative flow for `place`, of the whole offer for `evaluate` and `replace` | the cut, when the solve reached the maximum on its view; `None` when the offer stopped it below |
  | `place(0)` | the stored curve unchanged: `None` before the epoch's first solve or forward; after a `forward`, `None` for a nonempty roster and `[]` for an empty one; otherwise the last solve's | `None` |
  | `forward` with a nonempty roster | `None`: nothing is solved | `None` |
  | any of these but `place(0)`, with an empty roster | `[]` (the flow is 0 and maximal) | `[]` |

- **Histograms for the other models.** The sum of `(min_route_cost, delivered)` over the records is
  a histogram of delivery by each record's cheapest walk cost: the cost every unit travels at for
  `ShortestPaths` split routing, T1 and T2, whose routes have one cost, but not for a split route
  whose walks differ in cost.
- **Optional outputs**, `None` unless requested through `OutputOptions(edge_flows=False,
  residual=False, saturated_boundary_edges=False, edge_reservations=False, min_cut=False)`:
  `edge_flows` is the policy's own allocation per edge after the call (the ledger's total is
  `carried()`); `residual` and `saturated_boundary_edges` are the ledger's after the call;
  `edge_reservations` is the policy's own holds per edge after the call (zeros for a routing
  without holds; the ledger's total is `reserved()`); `min_cut` is described below.
- **`saturated_boundary_edges` and `min_cut` are different things.** `saturated_boundary_edges` is
  the saturation report of [section 6](#6-the-ledger), available in every model. `min_cut`
  (requested by `OutputOptions(min_cut=True)`, reported as the table above says) is the edges of a
  minimum cut of the policy's service view, the one on the sources' side: the allowed edges that
  leave the nodes the sources can still reach once the flow is maximal (counting the option to undo
  flow on an edge), in ascending EdgeId. An edge with nothing free in the view (capacity 0, or held
  by other policies) is included: it adds 0 to the cut's value and is listed because every allowed
  edge leaving that set is listed. An edge the policy's Eligibility excludes is never listed.

```python
g = graph(4, [(0, 1, 10, 1), (1, 3, 10, 1), (0, 2, 5, 2), (2, 3, 5, 2)])
ledger = ngc.FlowLedger(g)
pol = ngc.FlowPolicy(ledger, 0, 3, ngc.FlowPolicyConfig(routing=ngc.CentralizedTE(source_allocation=ngc.SourceAllocation.OPTIMIZED)))
want = ngc.OutputOptions(saturated_boundary_edges=True, min_cut=True)
r = pol.place(12, want)
print(r.call_delivered, r.marginal_cost_curve, r.min_cut, r.saturated_boundary_edges, r.termination.name)
r = pol.place(float("inf"), want)
print(r.call_delivered, r.marginal_cost_curve, r.min_cut, r.saturated_boundary_edges, r.termination.name)
print(r.requested, r.admitted, r.unadmitted, r.total_cost)
r = pol.forward(5, outputs=want)
print(r.call_delivered, r.total_cost, r.marginal_cost_curve, r.min_cut)
```

```text
12 [(2, 10), (4, 2)] None [0] COMPLETE
3 [(2, 10), (4, 5)] [0, 1] [0, 1] EXHAUSTED
inf 15 None 40
5 12 None None
```

After 12 units the cheap path (cost 2) is full and the solve stopped below the maximum, so there is
no cut; `saturated_boundary_edges` already reports the full edge 0. The second call adds 3 units at
cost 4 and reaches the maximum of 15, at a total cost of 10 x 2 + 5 x 4 = 40. `forward(5)` sends 5
over the installed paths by their weights 10 and 5, `H(5, [10, 5]) = (4, 1)`, at a cost of 4 x 2 + 1
x 4 = 12; it solves nothing, so it reports no curve and no cut.

### Several demands at once

    apply_actions(ledger, actions, *, mode=BatchMode.ORDERED, outputs=None) -> list[PlacementResult]
    Place(policy, offer)
    Replace(policy, offer, *, design_demand=None)
    Forward(policy, offer, *, black_holes=None)

`apply_actions` runs the actions of several policies of one ledger as one operation and returns
one `PlacementResult` per action, in input order, each with the fields of that policy's scalar
call (`place`, `replace` or `forward`). The whole list is validated first and published once: if
anything fails, the ledger and every policy are as they were. A policy appears at most once, and
every policy belongs to `ledger`. An action is an immutable value that keeps its policy (and its
black holes) alive; its offer and design demand are kept as given (an int or a float: another type
raises `TypeError` when the action is built), and the call checks their domain.

- **`BatchMode.ORDERED`** runs each action in list order exactly as its scalar call would, on the
  state the earlier actions left (their traffic and their holds); the policies not reached yet keep
  their old state. Normal outcomes (exhaustion, source drops, a down tunnel) are kept. Earlier
  actions are served first: on one edge of 10, an ECMP replace of 10 and then a T1 tunnel's replace
  of 10 deliver 10 and 0, though the tunnel holds its 10, because a hold reserves setup bandwidth,
  not service. New records take FlowIds in list order.
- **`BatchMode.SIMULTANEOUS`** forwards every action at the same time on fixed routes, the losses
  shared across the demands by [the loss rule](#losses-on-shared-edges). It takes `Replace` of
  `ShortestPaths` `EQUAL` or `WEIGHTED` (its routes prepared on the list's view) and of `Pinned`,
  and `Forward` of any policy with an installed program; every participant is LOSSY with a finite
  offer. Controllers compute nothing here and T1 sets nothing up: prepare them first (scalar calls
  or an ORDERED list), then forward them; a T1 owner keeps its holds. The view is the residual
  plus every participant's old traffic, released together: owners outside the list and orphans
  stay as fixed load, and holds are never added. Repeating a list with the same outside load gives
  the same state, and permuting it changes nothing. Each participant's black holes are its own:
  what reaches one is dropped there and takes no part in sharing that edge, while the others share
  it as usual. Removing a participant later frees only its own traffic; the others do not
  rebalance until they are forwarded again.
- **Outputs.** `outputs` applies to every entry. `residual` and `saturated_boundary_edges` read the
  ledger after that action (ORDERED) or after the whole list (SIMULTANEOUS: one residual for every
  entry, each entry's boundary from its own sources and Eligibility).
- **Validation, in this order:** the arguments' types (`TypeError`: `ledger must be a FlowLedger`,
  `actions must be a sequence of Place, Replace or Forward (got <type>)`, `actions[<i>] must be a
  Place, Replace or Forward (got <type>)`, `mode must be a BatchMode (got <type>)`, `outputs must be
  an OutputOptions (got <type>)`); the mode's value (`mode: <v> is not a value of BatchMode`); an
  empty list checks ledger usability under the gate before returning `[]`. Each action's argument
  checks run in list order: `list policy belongs to another ledger`, `duplicate policy in list`, the
  SIMULTANEOUS rules (`SIMULTANEOUS does not accept Place`, `SIMULTANEOUS Replace requires
  ShortestPaths EQUAL or WEIGHTED, or Pinned`, `SIMULTANEOUS requires admission=LOSSY`,
  `SIMULTANEOUS requires finite offers`), then the scalar call's argument checks. After all actions
  pass, check ledger usability and every action's mutable state under the gate, in list order. Only
  then compute anything. An error of the computation (a loss cycle that does not settle, a solver or
  forwarding size limit, an allocation failure) publishes nothing.

```python
# Two ECMP demands of 10 on one edge of 10, at the same time and one after the other.
g = graph(2, [(0, 1, 10, 1)])
ecmp = ngc.FlowPolicyConfig(routing=ngc.ShortestPaths(), admission=ngc.Admission.LOSSY)
for mode in (ngc.BatchMode.SIMULTANEOUS, ngc.BatchMode.ORDERED):
    led = ngc.FlowLedger(g)
    a = ngc.FlowPolicy(led, 0, 1, ecmp)
    b = ngc.FlowPolicy(led, 0, 1, ecmp)
    res = ngc.apply_actions(led, [ngc.Replace(a, 10), ngc.Replace(b, 10)], mode=mode)
    print(mode.name, [r.call_delivered for r in res], [r.call_edge_drops for r in res])

# Edges 0: 0 -> 2 of 5, 1: 1 -> 2 of 10, 2: 2 -> 3 of 10; a from 0 and b from 1, 5 each.
g = graph(4, [(0, 2, 5, 1), (1, 2, 10, 1), (2, 3, 10, 1)])
led = ngc.FlowLedger(g)
a, b = (ngc.FlowPolicy(led, s, 3, ecmp) for s in (0, 1))
out = ngc.OutputOptions(residual=True, saturated_boundary_edges=True)
for r in ngc.apply_actions(led, [ngc.Replace(a, 5), ngc.Replace(b, 5)],
                           mode=ngc.BatchMode.SIMULTANEOUS, outputs=out):
    print(r.call_delivered, r.residual.tolist(), r.saturated_boundary_edges)
```

```text
SIMULTANEOUS [5, 5] [[(0, 5)], [(0, 5)]]
ORDERED [10, 0] [[], [(0, 10)]]
5 [0, 5, 0] [0]
5 [0, 5, 0] [2]
```

At the same time the edge shares its 10 between the two demands; one after the other the first
takes it all. In the second list both entries report the residual after the whole list, and each
its own boundary: a's source reaches nothing beyond its full edge 0, b's reaches node 2, whose exit
is full.

### Configuration and call errors

These are the checks of a configuration and of a call, with their texts, in the order they run.
The first that fails is reported, before any work and with nothing changed, a zero offer included.
[Exception classes](#exception-classes) gives the rule for `TypeError` and `ValueError`.

**At construction** (`TypeError` unless stated):

- An unknown keyword (also every field that belongs to another model, such as `allocation` on
  `FlowPolicyConfig` or `source_allocation` on `DistributedTE`), a missing required argument
  (`routing`, `tunnel_count`, `CentralizedTE`'s `source_allocation`, `program`), `None` where none
  is allowed, and a value of the wrong type raise `TypeError` before any policy exists.
- A bool field (and a bool search option) takes only a `bool` or a `numpy.bool_`, known by its
  type: `multipath=None`, `multipath=1` and an object with `__bool__` are `TypeError`.
- An integer field of a configuration value (`TunnelCount` `per_head`, `AbsoluteCost` `cost`,
  `min_unreserved_bandwidth`, a `RouteDefinition`'s `head` and `weight`) given a value of another
  type raises `TypeError` (`TunnelCount per_head must be an integer (got bool)`); an integer
  outside 64 bits raises `ValueError` (`TunnelCount per_head = <v> is out of range [...]`).
- An `edge_weights` that is not a 1-D sequence of integers (a scalar, a string, a 2-D array, a
  bool or float array, a bool or float entry) raises `TypeError`
  (`edge_weights[0] must be an integer (got bool)`); an entry beyond 64 bits raises
  `ValueError('edge_weights[<i>] must be an integer in [0, 2^62)')`.
- `source_weights` is converted as `edge_weights` is, in three passes over the whole vector, each
  reporting its first defect. First the shape: an array of another rank, a string, bytes or a
  non-sequence, then an entry that is itself a sequence. Then the entry types:
  `source_weights[<i>] must be an integer (got <type>)` (`TypeError`). Last, an integer that does
  not fit in signed 64 bits: `source_weights[<i>] must be an integer in [0, 2^62)` (`ValueError`;
  an unsigned array is checked like any other). So `[1.5, [1]]` fails on entry 1's shape, and
  `[2**64, 1.5]` on entry 1's type. Any sequence of integers is accepted, a `range` or a tuple
  too.

**When a policy or an evaluator uses the configuration** (`ValueError`). The model matrix, in the
order checked:

1. the enumeration values (`admission`, then the model's `allocation`, `source_allocation`,
   `weight_basis`, `tie_break`), each `<field>: <v> is not a value of <Enum>` (reachable with a
   value built from an integer);
2. the route allocation and the source allocation: first `source_allocation=WEIGHTED requires
   CentralizedTE without tunnel_count`; then `ShortestPaths with EQUAL or WEIGHTED requires
   source_allocation=EQUAL`; for T2, `CentralizedTE with tunnel_count and several sources requires
   source_allocation=EQUAL` (`OPTIMIZED` with more than one input source, counted after removing
   repeats and before any Eligibility; with one source it is the same calculation as `EQUAL`);
3. the admission: `ShortestPaths OPTIMIZED requires admission=LOSSLESS`; `CentralizedTE without
   tunnel_count requires admission=LOSSLESS`;
4. the route choice: `ShortestPaths tie_break requires multipath=False or parallel_edges=False` (a
   `HIGHER_*` rule with every route and every parallel edge kept, also for `OPTIMIZED`, whose route
   the rule selects); then the **pooled single-path rule**: `ShortestPaths with multiple input
   sources requires tie_break=FEWER_HOPS when multipath=False` (more than one input source, counted
   after removing repeats and before any Eligibility: a bottleneck choice could differ per source
   at a router the sources share);
5. `prefer_headroom requires tunnel_count`;
6. the weights: `weight_basis requires allocation=WEIGHTED`; `weight_basis=EXPLICIT requires
   edge_weights`; `edge_weights require weight_basis=EXPLICIT`; for CentralizedTE,
   `source_allocation=WEIGHTED requires source_weights` and `source_weights require
   source_allocation=WEIGHTED`;
7. last, a cache given to a model that does not read one: `<DistributedTE or CentralizedTE or
   Pinned> does not use ShortestPathCache` (only `ShortestPaths` reads a `ShortestPathCache`).

Then each model's own fields: the tunnel count in `[1, 10000]` (`tunnel count must be in [1,
10000]`; `TunnelCount(per_head)` itself takes any 64-bit integer); for T1 an "up to" count
(`DistributedTE requires an exact tunnel count`), then the ceiling (`cost_ceiling AbsoluteCost
cost must be >= 0 (got <c>)`, `cost_ceiling CostStretch factor must be >= 1 (got <n>/<d>)`), then
`min_unreserved_bandwidth must be an integer in [0, 2^62)`; `edge_weights`' length and entries
(`edge_weights: expected <E> entries`, `edge_weights[<i>] must be an integer in [0, 2^62)`).
Then `source_weights`: its length against the graph's node count
(`source_weights: expected <N> entries`), then every entry in node id order
(`source_weights[<i>] must be an integer in [0, 2^62)`). Each weight is bounded; their sum may
exceed 64 bits.

**`FlowPolicy` construction**, in this order:

1. the endpoint sets (empty, out of range, the lowest node in both);
2. the graph of each object the policy is given: the `Eligibility` (`foreign Eligibility: built
   for another graph`), the program (`program belongs to another graph`), the support (`support
   belongs to another graph`, with or without a tunnel count, an empty support included), the
   `ShortestPathCache` (`foreign ShortestPathCache: bound to another graph`), the black holes
   (`black_holes belongs to another graph`);
3. the model matrix and the models' fields above;
4. the tunnel limit (T1, T2: the roster's heads times the tunnel count, for "up to" the largest
   count, is at most 2^30 - 1: `tunnel limit exceeded: <H> heads x <N> tunnels > 1073741823`; no
   tunnel is made at construction);
5. `Pinned`'s endpoints (`program endpoints do not match policy endpoints`; any number of each).

**Calls.**

- `evaluate` and `replace` check the offer (`offer must be an integer in [0, 2^62) or permitted
  +inf`), then a Lossy `+inf` (`LOSSY forwarding requires a finite offer`), then the design
  demand: the routing (`<model> does not accept design_demand`, the model `ShortestPaths`,
  `CentralizedTE` or `Pinned`), its domain (`design_demand must be an integer in [0, 2^62) or
  permitted +inf`), its finiteness for T1 (`DistributedTE requires a finite design demand`, also
  for an infinite offer without a design).
- `forward` checks, in order: `Pinned forwards with replace`; the offer; a Lossy `+inf`
  (`LOSSY forwarding requires a finite offer`, as for `evaluate`, `replace` and Lossy `Pinned`);
  the black holes' type and graph (`black_holes belongs to another graph`); then, under the ledger
  gate, `forward requires an installed program` (also after `remove()`).
- `place`'s checks are in [the placement call](#the-placement-call), and a list's in
  [several demands at once](#several-demands-at-once).

## 8. Failures, threads and lifetimes

**In plain words.** A call either does everything or changes nothing: if it fails, for any
reason, every ledger, policy and cache is as it was and you get one exception. There are two
stated exceptions to this, described below. Immutable values can be shared freely between
threads. Calls on one ledger serialize; different ledgers run in parallel. Cache leases serialize
building calls and admit frozen readers together.

**Exact rules.**

- **All or nothing.** A failure anywhere in a public call (validation, memory allocation, a
  worker thread, the construction of the Python result) leaves every ledger, policy, cache and
  internal buffer as it was before the call and raises one exception in the caller's thread.
  Worker threads are always joined first; an exception never escapes a worker. When several pairs
  or candidates of an evaluator fail, the exception of the lowest index is raised.
- **Limit 1: pybind11's own object machinery.** Constructing a bound object from Python (any of
  the library's classes), parsing the arguments of a call that passes keywords, and converting the
  result of a method into a new Python object outside the library's checked construction (for
  example `PredecessorDAG.from_csr`, `Graph.eligibility`) use pybind11's own allocation paths,
  which do not all check for failure. A memory failure at exactly that point may crash the process
  instead of raising. Library state is never left half-changed by it. The results of `spf`, `ksp`,
  `resolve_to_paths`, `ForwardingProgram.from_routes` and `from_dag`, `FlowPolicy.place`,
  `evaluate`, `replace` and `forward`, `apply_actions`, `evaluate_demand`, `evaluate_pairs` and
  `evaluate_edge_failures` are built through a checked path and stay inside the guarantee, and so
  is every read of a field or property and every `repr` of the library's values (the configuration
  values and their enumerations, `FlowPolicy.config`, `program`, `plan`, `by_source`, `by_tunnel`,
  `tunnels`, `termination`, a tunnel's `status`, a record's `route`, `SensitivityResult.baseline`
  and `entries`, among others), and `records()` and `record()` of a policy and a ledger,
  `FlowLedger.graph` and `Eligibility.intersect`.
- **Limit 2: a `ShortestPathCache` passed to the call.** A building cache is a memo: when a call fails
  after some cache lookups succeeded (a policy constructor failing after its lookup, an evaluator
  failing midway), the cache keeps those complete entries and their statistics. It never holds a
  partial entry, every value it handed out stays valid, and a kept entry never changes any later
  result. A failed lookup counts one miss. A frozen cache never changes its entries; a failed
  call changes none of its statistics, except one miss when the call raised `ShortestPathCacheMiss`.
- **Threads and the GIL.** Each ledger has one native gate protecting its projections and every
  bound policy, including reads, lists, removal and policy destruction. Calls on one ledger take
  effect one at a time; different ledgers can compute in parallel. Native work releases the GIL,
  including policy calls, searches, factories and array conversions. Immutable reads keep the GIL
  and take no gate: graph, mask, route, program, configuration and result fields, `policy.config`,
  cache `graph_uid()` and `stats()`. Each statistics field is read atomically, without promising a
  consistent snapshot of all fields.
- **Python and the native locks.** Python conversions and callbacks (an argument's `__index__`,
  `__eq__`, `__hash__`, `repr`) run with no native lock held. Waiting for a ledger gate or a cache
  lease releases the GIL, or detaches the thread state on free-threaded Python. Native checks that
  grow with the input release the GIL too: the tunnel-limit roster search of policy construction and
  of the evaluators (once per pair of `evaluate_pairs`), and `from_routes`' walk of each supplied
  `PredecessorDAG`. A call takes at most one cache lease, then the ledger gate.
- **Results and publication.** A call that publishes (a scalar call, a list, an evaluator, a
  policy's construction) builds every Python result object it returns before it changes
  anything, so a failure there (`MemoryError`) never leaves a published change behind. A call that
  publishes nothing (a read, a search, a factory) builds its result after releasing every lock,
  from a copy it owns.
- **Validation passes.** In the documented validation order, all argument and immutable checks
  precede all mutable-state checks. Thus an invalid offer or foreign black holes precede
  `forward requires an installed program`; a later invalid list action precedes an earlier
  action's epoch error. Installed programs, finite `place` calls after infinity, requested limits,
  cache mode, ledger usability and identity availability are checked under the gate or lease,
  in their documented order, action by action. No list action computes before both passes finish.
- **Free-threaded Python.** The cp314t build follows the same concurrency rules. It detaches its
  Python thread state before waiting for a gate or lease. A call still converts its arguments and
  copies every borrowed array before taking any lock; changing an argument object from another
  thread while a call converts it is the caller's race.
- **Lifetimes.** The graph outlives the ledger, the ledger outlives its policies; Python keeps
  the graph, ledger, Eligibility and cache alive while a policy needs them.
- **Immutable values.** The configuration values (`FlowPolicyConfig`, the four models,
  `TunnelCount`, `AbsoluteCost`, `CostStretch`, `OutputOptions`, `EvaluationOptions`,
  `RouteDefinition`), the search options (`SpfOptions`, `KspOptions`) and the actions of a list
  (`Place`, `Replace`, `Forward`) have read-only fields: no option object changes after it is
  built. A policy validates its configuration once and keeps that snapshot (`policy.config`); an
  evaluator validates and copies its options at call entry, and `apply_actions` its actions.
```python
g = graph(4, [(0, 1, 10, 1), (1, 3, 10, 1), (0, 2, 2, 1), (2, 3, 2, 1)])
ledger = ngc.FlowLedger(g)
pol = ngc.FlowPolicy(ledger, 0, 3, config(ngc.ShortestPaths()))
print(pol.place(float("inf")).call_delivered)
try:
    pol.place(1)          # a finite offer after +inf in the same epoch
except ValueError:
    print("ValueError; unchanged:", pol.delivered(), ledger.carried().tolist())
```

```text
5
ValueError; unchanged: 5 [3, 2, 3, 2]
```

### Exception classes

One rule for the class of an input error, at every public entry (constructors, search options,
call arguments, evaluator options): a Python type the position never takes is a `TypeError`; an
accepted type whose value is outside its domain is a `ValueError`. The text is the same either
way; a text that covers both causes (`offer must be an integer in [0, 2^62) or permitted +inf`,
`a FlowId is a tuple (ledger uid, id)`, `each pair must be (sources, destinations)`) takes the
class of its cause.

| Position | Takes | `TypeError` | `ValueError` |
| --- | --- | --- | --- |
| an integer (a count, a node or edge id, an entry of an edge array, a FlowId part, `as_float`'s amounts) | any object with `__index__`, read exactly; a NumPy integer; an `int` at or above `2^63` is an integer | a bool (Python's or `numpy.bool_`); a float, `k=1.0` and `k=float('inf')` included; a string; `None` where none is allowed; an object of another class | out of range; beyond 64 bits where the range is bounded |
| an offer, a design demand | an integer, or `float('inf')` where the call allows it | a bool, a string, `None` where one is required, another non-number | a finite float, NaN, `-inf`, an integer out of range, `+inf` where the call refuses it (an evaluator's `offer`) |
| a non-integer number: a `resolution` (`Graph`, `to_units`, `as_float`) and `to_units`' values | exactly an `int` (not a bool), a `float`, a NumPy integer or floating scalar, an `int` or `float` subclass by its value (never through its own `__int__`, `__float__` or `__array__`); for `to_units` also an array of an integer or floating dtype or a sequence of those scalars | a `Decimal`, a `Fraction`, an object known by `__float__`, `__index__` or `__array__`, a 0-d array as a resolution; a bool entry beside numbers in `to_units` (NumPy would read it as a number; an all-bool sequence keeps its dtype text) | NaN, an infinity, a negative amount, a result at or above `2^62` units, a resolution that is not a power of two in `[2^-32, 2^32]` (an integer beyond the float range included: `resolution must be a power of two in [2^-32, 2^32]`), an integer too large for NumPy beside numbers |
| a bool field or option | a `bool` or a `numpy.bool_` | anything else | - |
| a vector | a 1-D sequence or NumPy array of integers | a non-sequence, a string, a nested entry, a NumPy array of another rank or of a non-integer dtype, an entry of a wrong type | a wrong length, an entry out of range |
| a `fractions.Fraction` | only by `CostStretch` ([cost ceilings](#cost-ceilings-and-the-exact-stretch)) | everywhere else | - |

The other classes:

- `ValueError` also for an invalid combination, a foreign object (a value of another graph,
  ledger or policy) or a call in the wrong state.
- `TypeError` also for a field or option that does not exist (for example
  `ShortestPaths(tunnel_count=ngc.TunnelCount(2))`), a missing required argument, an argument
  that cannot be converted to its parameter's type, a bool field given anything but a `bool` or
  `numpy.bool_`, and an integer field or `edge_weights` of a configuration value given a
  non-integer.
- Every converted value passes on unchanged any exception raised by the value's own conversion
  (its `__float__` or `__index__`), whatever its class.
- `ShortestPathCacheMiss` for a miss in a frozen cache (its own class, neither a `ValueError` nor
  a `RuntimeError`); `MemoryError` for an allocation failure (everything rolled back);
  `RuntimeError` for an internal invariant failure and for an exhausted identity counter.

### Limits that refuse a call

Each of these raises before anything changes.

| Limit | What it counts | Error |
| --- | --- | --- |
| a flow solve (`CentralizedTE` routing, the `OPTIMIZED` allocation and the evaluators that use them) | at most `INT32_MAX - 2` nodes, and at most `2^30 - 1` edges, sources and destinations together. A `CentralizedTE` solve counts the graph's nodes, the eligible edges, the roster sources and the eligible destinations; a `ShortestPaths` `OPTIMIZED` solve counts the nodes of its route's edges and endpoints, those edges, its sources and its destinations | `ValueError('flow solver node limit exceeded')`, `ValueError('flow solver arc limit exceeded')` |
| a `LOSSY` forwarding | each at most `2^31 - 1`, checked in this order: `entry` (the forwarded program entries), `occurrence` (every link visit of every entry: a tunnel's path length, a split route's member links), `router` (the routers of split routes), `event` (the links that can drop plus the routers, counted on the actual traffic, each demand's stopped at its own black holes) and `dependency` (the arcs between events). A list of several demands counts its one joint network, every participant together | `ValueError('forwarding <counter> limit exceeded: <n> > 2147483647')` with the exact count `<n>` |
| cumulative counters (requests, `offered`, drops) | below `2^127` ([section 1](#the-conversion-exactly)) | `ValueError('<counter> limit exceeded')` |
| policy identifiers of a ledger | at most `2^63 - 1` | `RuntimeError('policy identity limit exceeded')` |
| FlowId counters of a ledger | at most `2^63 - 1`; a publication's new records counted together, a list's included | `RuntimeError('record identity limit exceeded')` |
| identities of graphs, ledgers and masks in the process | at most `2^63 - 1` (0, the all-allowed mask's, is never issued) | `RuntimeError('object identity limit exceeded')` |

Only graphs of billions of nodes or about a billion edges reach the flow-solve limits.

## 9. The shortest-path cache

**In plain words.** Many placements need the same static shortest-path DAG toward the same
destinations. A `ShortestPathCache` remembers those DAGs for one graph. It stores the full DAG and
applies each caller's route choice afterwards, so one entry serves every route choice. A cache is either
building (one writer adds entries) or frozen (read-only, shared by many threads, and a miss is an
error).

**Exact rules.**

- **Binding and key.** `ShortestPathCache(graph)` is bound to one graph (`graph_uid()`). The key is
  `(the root set as a sorted set, the direction, the Eligibility uid)`. The value is the full
  static shortest-path DAG: the DAG of every equal-cost walk (`multipath=True`,
  `parallel_edges=True`) with no residual, edges of capacity 0 included. Values are
  immutable and stay valid after `clear()` and after the cache is destroyed.
- **The route choice is never cached.** The single path and the member rule are applied when a
  program is prepared, on the cached DAG, reading the scores of `HIGHER_AVAILABLE_CAPACITY` and
  `HIGHER_HEADROOM` from the residual at that moment and `HIGHER_CAPACITY` from static capacities.
  Because the cached DAG keeps every equal-cost parallel member, a member that a static choice
  would have dropped is still available when its residual makes it the right one.
- **Who reads it.** `ShortestPaths` policy construction reads the cache (keyed by its destination
  set, reverse direction, and Eligibility uid) and retains the support DAG. Later scalar calls and
  lists use that retained support without a cache lease. The evaluators read their supplied cache
  when they run a `ShortestPaths` configuration. The other models do not read it and reject a
  supplied cache (`ValueError`, for example `CentralizedTE does not use ShortestPathCache`), as the
  last check of the model matrix.

| Operation | Building | Frozen |
| --- | --- | --- |
| lookup hit | returns the value | returns the value; safe from many threads |
| lookup miss | searches, inserts, returns | raises `ShortestPathCacheMiss` |
| `freeze()` | becomes frozen | no effect |
| `clear()` | empties the cache; hit and miss counters persist | `ValueError` |
| `thaw()` | no effect | waits for active calls, then becomes building |
| destruction | callers keep it alive during calls | callers keep it alive during calls; returned values stay valid |

- **Leases.** Each building call owns one exclusive lease for its whole computation; frozen
  calls own shared leases. `freeze()`, `thaw()` and `clear()` take exclusive admission and wait
  for active calls. Admission is FIFO: a queued lifecycle call cannot be overtaken by later
  builders or readers. Evaluator workers borrow their outer call's lease and never acquire it
  again. The cache lease precedes any ledger gate.
- **Statistics.** `stats()` returns `hits`, `misses`, `entries`, `bytes` (the entries' retained
  arrays, without container overhead) and `live_values` (values alive anywhere, cache entries
  included). Each field is atomic; together they are not one consistent snapshot. An invalid
  query (foreign graph, bad root) counts nothing. A frozen cache's statistics follow the public
  call: its hits count only when the call succeeds, its Python result included; a failed call
  changes no statistic, except one miss for a call that raises `ShortestPathCacheMiss`.
- **Workers.** A supplied building cache is accepted only with `workers == 1` (`ValueError`
  otherwise). A frozen cache may be shared by every worker. Without a supplied cache, each
  worker uses its own private building cache.
- **Uids and frozen caches.** A per-pair Eligibility that the evaluator must intersect with a
  common one gets a new uid when the result differs from both, and that key cannot be in a
  frozen cache. To reuse prefilled entries, pass already-combined per-pair values with
  `eligibility=None`, or a common value that contains them (so the intersection returns the
  prefilled operand).
- **Scope.** Use one cache per series of analyses on one graph, and drop it at the end of the
  series.

```python
g = graph(5, [(0, 1, 10, 1), (1, 3, 10, 1), (0, 2, 2, 1), (2, 3, 2, 1), (3, 4, 9, 1)])
cache = ngc.ShortestPathCache(g)
opts = ngc.EvaluationOptions(config(ngc.ShortestPaths()), cache=cache)
rs = ngc.evaluate_pairs(g, [(0, 3), (1, 3)], opts)          # building: workers defaults to 1
print([r.call_delivered for r in rs], cache.stats().hits, cache.stats().misses)
cache.freeze()
rs = ngc.evaluate_pairs(g, [(0, 3), (1, 3)], opts, workers=2)  # frozen: shared by 2 workers
print([r.call_delivered for r in rs], cache.stats().hits, cache.stats().misses)
try:
    ngc.evaluate_demand(g, 0, 4, opts)                # destination 4: not in the cache
except ngc.ShortestPathCacheMiss:
    print("ShortestPathCacheMiss; misses =", cache.stats().misses)
```

```text
[5, 10] 1 1
[5, 10] 3 1
ShortestPathCacheMiss; misses = 2
```

Both pairs go to destination 3, so they read one key (destination 3, reverse search, the
all-allowed Eligibility). Building, the first pair misses and fills it and the second hits: 1 hit,
1 miss. Frozen, both pairs hit: 3 hits. Destination 4 has no entry, so the frozen cache raises and
counts one more miss.

## 10. Evaluators

**In plain words.** The evaluators answer questions about a model without touching your ledgers.
`evaluate_demand` asks how much one model carries between two node sets on a fresh, empty ledger,
either as much as possible or at a given offer. `evaluate_pairs` asks that for many independent
pairs. `evaluate_edge_failures` asks, for each link, how much the answer drops when that link is
removed, at the same offer. A fresh ledger means these answers can differ from placing into a
ledger that other demands already occupy: to ask about a loaded network, use `FlowPolicy.evaluate`
on a ledger with that load.

**Exact rules.**

- **`EvaluationOptions(config, *, eligibility=None, cache=None, outputs=None, design_demand=None,
  offer=None)`** is the one options type of the three evaluators: the configuration, the static
  Eligibility (`None`: everything allowed), the `ShortestPathCache` (`None`: a private cache per
  worker), the optional outputs of each result (`evaluate_edge_failures`: of its baseline), the
  tunnel models' design (`design_demand`) and the `offer`. The design and the offer are kept as
  given (an `int` or a float; another type, a `bool` or a string, raises `TypeError` at
  construction); the evaluator checks them.
- **A finite offer for every model.** `offer`, when given, is a finite integer in `[0, 2^62)` for
  every model under each admission it supports, in all three evaluators: the offer of every pair,
  of the baseline and of every failed network. Each evaluation computes the same amounts as a
  fresh policy's `evaluate(offer, design_demand=...)`: Lossless delivers up to the offer, Lossy
  accounts for exactly that offer. An absent offer is the uncapped Lossless question or, only for
  `ShortestPaths` `EQUAL`/`WEIGHTED` with Lossy, the first-hop experiment.

  | Model | Lossless | Lossy |
  | --- | --- | --- |
  | `ShortestPaths` `EQUAL` / `WEIGHTED` | `offer` `None`: the uncapped bound; finite: delivery up to it | `offer` finite: that offer; `None`: the first-hop experiment |
  | `ShortestPaths` `OPTIMIZED`, `CentralizedTE` without a tunnel count | `offer` `None`: the uncapped optimum; finite: the optimum up to it | refused by the matrix |
  | T1 | `design_demand` required and finite; `offer` `None`: the lossless bound of those tunnels; finite: up to it | `design_demand` required and finite; `offer` required |
  | T2 | `design_demand` `None`: the offer if it is finite, else unlimited; `offer` `None`: the bound at that design | `offer` finite; `None`: the design if finite, else an error |
  | `Pinned` | `offer` `None`: the fixed-program bound; finite: up to it | `offer` required |

- **Validation order** (before any worker starts): the policy's construction checks (the
  configuration, [section 7](#configuration)); then an explicitly given offer (`offer must be an
  integer in [0, 2^62)`, an explicit `float('inf')` included); then the design: only T1 and T2
  take one (`<model> does not accept design_demand`, checked before its range), in `[0, 2^62)` or
  `float('inf')` (`design_demand must be an integer in [0, 2^62) or permitted +inf`), T1's finite
  (`DistributedTE requires a finite design demand`); then T1's missing design (`DistributedTE
  evaluation requires design_demand`); then the missing Lossy offer (`DistributedTE LOSSY
  evaluation requires a finite offer`, `Pinned LOSSY evaluation requires a finite offer`, and
  `CentralizedTE LOSSY evaluation requires a finite offer` for T2 when no finite design gives one).
  The model in these texts is the class name, T2 included.
- **`evaluate_demand(graph, sources, destinations, options)`** evaluates a `FlowPolicy` with
  `options.config` on a fresh ledger of the graph and returns its `PlacementResult`; it commits
  nothing anywhere. Without an offer (the uncapped question) the configuration
  `CentralizedTE(source_allocation=OPTIMIZED)` gives the exact maximum flow: its result carries the
  marginal cost curve and, with `OutputOptions(min_cut=True)`, the minimum cut. A finite offer
  below the maximum stops the solve there, with no cut ([the result](#the-result)). With `EQUAL`
  source allocation it is the largest total whose equal shares all fit.
- **`evaluate_pairs(graph, pairs, options, per_pair_eligibility=None, workers=None)`** runs
  `evaluate_demand` for each `(sources, destinations)` pair and returns the results in pair order.
  Nothing is shared between pairs (no reservations, no capacity): the results are independent
  capacities, not allocations that coexist; demands that share capacity are policies on one ledger.
  Pair `i` uses `options.eligibility` intersected with `per_pair_eligibility[i]` (entries may be
  `None`; a given list has one entry per pair). All checks finish before evaluation workers start:
  first the arguments (the worker count, the graph of `options.eligibility` and of the cache,
  `Pinned`, which is rejected because a program fixes its endpoints (`evaluate_pairs does not accept
  Pinned`), the configuration, the offer and design rules above, then the length of
  `per_pair_eligibility`), then each pair completely (its endpoint sets, their overlap, its
  Eligibility's graph, the pooled single-path rule of [section 7](#configuration), T2's
  several-sources rule and, for the tunnel models, the tunnel limit for its sources) in pair order,
  so the error of the lowest failing pair is the one raised. A supplied building cache's `workers ==
  1` rule reads the cache's mode and comes last, under the cache's lease. An empty batch returns
  `[]`.
- **`evaluate_edge_failures(graph, sources, destinations, options, workers=None)`** returns a
  `SensitivityResult`: `baseline` (the result of `evaluate_demand`, with the requested outputs)
  and `entries`, one `SensitivityEntry(edge, delivery_reduction)` per candidate edge, in ascending
  EdgeId.
  - **Candidates.** The edges that get an entry depend on the configuration; an edge without an
    entry has reduction 0:

    | Configuration | Edges with an entry |
    | --- | --- |
    | unrestricted `CentralizedTE` with `OPTIMIZED` sources | the edges that carry flow in the baseline: removing any other edge keeps the baseline feasible, so the delivery stays, at any offer |
    | every other configuration (`EQUAL` or `WEIGHTED` sources, a support, `ShortestPaths`, T1, T2, `Pinned`) | every allowed edge, whatever its capacity: an edge of capacity 0 can carry a route and keep a source in the roster, so deleting it can change the result either way. `ShortestPaths` gives an edge outside the static shortest-path support reachable from the roster sources reduction 0 without recomputing (no route can use it) |

  - **`delivery_reduction`** is `baseline.delivered` minus the delivery of a fresh
    `evaluate_demand` with that edge also excluded, at the same offer, design and other options,
    in units, signed and never rounded or thresholded. It is not the traffic the deleted edge
    carried. A negative value means the model carries more without the edge (for example, equal
    splitting no longer sends a share down a thin path, or a source behind it no longer holds the
    equal allocation back). The test suite compares every reduction of the optimized controller
    with such a fresh evaluation, on small networks whose expected values are written out in the
    tests.
  - **Each deletion.** T1 sets up again for the same design under each deletion; T2 redoes its
    sizing, count choice and objectives, on the same supplied support when there is one.
    `CentralizedTE` with a support solves each deletion afresh. `Pinned` replays the same program
    under each deletion (a repaired fixed-route question, not a new optimum and not stale
    forwarding): deleting an edge of a tunnel fails that tunnel.
  - **Cache.** Each deletion uses a new Eligibility uid, so a `ShortestPaths` deletion that has to
    be placed raises `ShortestPathCacheMiss` through a frozen cache. The other models reject a
    supplied cache. Private caches are cleared between candidates.
- **Workers.** An explicit `workers` (an integer of at least 1) is always honoured; the
  evaluator uses at most one thread per pair or candidate. With `workers=None` the default is 1
  when a cache is supplied, and otherwise comes from the environment:

  | Variable | Default worker count of | Accepted values |
  | --- | --- | --- |
  | `NGRAPH_CORE_BATCH_THREADS` | `evaluate_pairs` | an integer in `[1, 4096]`; unset or anything else means 1 |
  | `NGRAPH_CORE_SENSITIVITY_THREADS` | `evaluate_edge_failures` | the same |

  The variables only supply defaults; they never override an explicit `workers`. The C++
  functions take the worker count as an explicit argument.
- **Failures** follow [section 8](#8-failures-threads-and-lifetimes): a failed evaluator call raises
  (also a solve beyond the flow core's size limits of [section 8](#limits-that-refuse-a-call), with
  `CentralizedTE` or `OPTIMIZED`), and a failed deletion fails the whole call (the exception of the
  lowest failing candidate; no entry is returned); a supplied building cache keeps the entries of
  the lookups that succeeded before the failure; a supplied frozen cache's statistics are unchanged,
  except one miss when the call raises `ShortestPathCacheMiss`.

```python
g = graph(4, [(0, 1, 10, 1), (1, 3, 10, 1), (0, 2, 2, 1), (2, 3, 2, 1)])
for name, cfg in (("ShortestPaths", config(ngc.ShortestPaths())),
                  ("CentralizedTE", ngc.FlowPolicyConfig(routing=ngc.CentralizedTE(source_allocation=ngc.SourceAllocation.OPTIMIZED)))):
    s = ngc.evaluate_edge_failures(g, 0, 3, ngc.EvaluationOptions(cfg))
    print(name, s.baseline.call_delivered, [(e.edge, e.delivery_reduction) for e in s.entries])
```

```text
ShortestPaths 5 [(0, 3), (1, -5), (2, 3), (3, -5)]
CentralizedTE 12 [(0, 10), (1, 2), (2, 10), (3, 2)]
```

Under equal splitting the thin path (edges 1 and 3) holds the whole diamond back to 5 units;
without it the wide path alone carries 10, a reduction of -5.

A failure study at a fixed offer, against the [first-hop
experiment](#5-route-allocations-and-admission): two parallel links, ECMP, Lossy. On links of 10 and
10, an offer of 4 is delivered by every network; without an offer the experiment injects each
remaining link's capacity, 20 before a deletion and 10 after, a different question. On links of 4
and 10 without an offer, the experiment delivers 14 and the deletions leave 10 or 4: reductions 4
and 10. With an offer of 12 each link is offered 6, the 4-unit link drops 2 and the baseline
delivers 10; the deletions deliver 10 and 4: reductions 0 and 6:

```python
lossy = config(ngc.ShortestPaths(), admission=ngc.Admission.LOSSY)
for capacities, offers in (((10, 10), (4, None)), ((4, 10), (None, 12))):
    g = graph(2, [(0, 1, c, 1) for c in capacities])
    for offer in offers:
        s = ngc.evaluate_edge_failures(g, 0, 1, ngc.EvaluationOptions(lossy, offer=offer))
        print(capacities, offer, s.baseline.call_delivered, [e.delivery_reduction for e in s.entries])
try:
    ngc.evaluate_demand(g, 0, 1, ngc.EvaluationOptions(lossy, offer=float("inf")))
except ValueError as error:
    print(error)
```

```text
(10, 10) 4 4 [0, 0]
(10, 10) None 20 [10, 10]
(4, 10) None 14 [4, 10]
(4, 10) 12 10 [0, 6]
offer must be an integer in [0, 2^62)
```

### Limits

- **One route allocation per program.** A program applies one route allocation at every router.
  Weighted ingress above routers that should split equally is expressible exactly when every
  member of those downstream equal splits has positive static capacity: give each member of such a
  split the same positive `EXPLICIT` weight, since `H(x, [w] * k) = RR(x, k)` for any positive `w`.
  A downstream equal split with a zero-capacity member differs (`WEIGHTED` prunes that member,
  `EQUAL` keeps it) and is not expressible.
- **Certificates and exports.** The public optimality certificates are `min_cut` and
  `marginal_cost_curve` of unrestricted optimized `CentralizedTE`. A program has no serialized
  form for moving it between processes, and its per-router tables are not exported.
- **Not modelled:** traffic classes and priority queues, RSVP setup and holding priorities and
  preemption, bandwidth policing, automatically disjoint or SRLG-diverse tunnels, a stretch
  baseline per head inside one pooled T1 demand, growth of a tunnel count over time, and a joint
  optimum or joint bound over several demands. Per-pair offers, designs or configurations are
  separate `evaluate_demand` calls, and a study of the three failure states is a loop of the
  caller's.

## 11. C++ names

The C++ library (namespace `netgraph::core`) has the same semantics. Its public headers are
self-contained:

| Header | Declares |
| --- | --- |
| `netgraph/core/types.hpp` | `NodeId`, `EdgeId` (`int32`), the cost and amount aliases (`int64`), `kDomainLimit` (`2^62`), `kUnreachable` |
| `netgraph/core/graph.hpp` | `Graph`, `to_units`, `as_float` |
| `netgraph/core/eligibility.hpp` | `Eligibility` (`allows_node`, `allows_edge`) |
| `netgraph/core/route_structures.hpp` | `PredecessorDAG`, `Path` |
| `netgraph/core/options.hpp` | `Direction`, `AbsoluteCost`, `CostStretch`, `RouteTieBreak`, `RouteAllocation`, `WeightBasis`, `Admission`, `SourceAllocation` |
| `netgraph/core/forwarding_program.hpp` | `ForwardingProgram` (`from_routes`, `from_dag`), `Tunnel`, `TunnelStatus`, `RouteDefinition`, `RouteOptions` (the options of `from_routes` and `from_dag`) |
| `netgraph/core/search.hpp` | `spf`, `ksp`, `resolve_to_paths`, `SpfOptions`, `SpfResult`, `KspOptions`, `SpfQueue` |
| `netgraph/core/flow_record.hpp` | `FlowId`, `Counter128` (the 128-bit cumulative counter), `FlowRecord` |
| `netgraph/core/flow_ledger.hpp` | `FlowLedger` |
| `netgraph/core/flow_policy.hpp` | `FlowPolicy`, its configuration (`FlowPolicyConfig`, `ShortestPaths`, `DistributedTE`, `CentralizedTE`, `Pinned`, `TunnelCount`, `OutputOptions`) and result types (`PlacementResult`, `SourceFlow`, `TunnelFlow`, `TunnelPlan`), `Termination`; the lists: `apply_actions`, `BatchMode`, `Place`, `Replace`, `Forward`, `Action` |
| `netgraph/core/shortest_path_cache.hpp` | `ShortestPathCache`, `ShortestPathCacheMiss` |
| `netgraph/core/evaluators.hpp` | `evaluate_demand`, `evaluate_pairs`, `evaluate_edge_failures`, `EvaluationOptions`, `SensitivityEntry`, `SensitivityResult` |

Spelling differences: Python enum members and constants are upper case
(`RouteTieBreak.FEWER_HOPS`, `UNREACHABLE`), C++ enumerators are the same words in PascalCase
(`RouteTieBreak::FewerHops`, `RouteAllocation::Optimized`) and C++ constants carry a `k` prefix
(`kUnreachable`). Where Python returns `None` for an unreachable cost (`target_cost`,
`min_route_cost`), C++ returns `kUnreachable`; an infinite offer is an empty `std::optional`, and
its `call_undelivered` and `requested` are empty too. An optional field that Python reports as
`None` (`ShortestPaths`' `source_allocation` and `weight_basis` before resolution, `offered`,
`unadmitted`, `total_cost`, `marginal_cost_curve`, `plan` and T2's fields of a `TunnelPlan`,
`CentralizedTE`'s `tunnel_count`, `FlowRecord.entry_id`, a down tunnel's `destination`, `edges`
and `path_cost`, the optional outputs) is an empty `std::optional` in C++;
`max_fractional_bandwidth_per_tunnel` is a reduced `(numerator, denominator)` pair,
`total_tunnel_cost` and `headroom_penalty` are `Counter128`; `CostStretch` is constructed from a
numerator and a denominator (reduced on construction, `numerator()` and `denominator()`); a
`design_demand` or an evaluator's `offer` is a `std::optional<std::optional<Amount>>`, absent for
`None` and an empty inner value for `float('inf')`; and a `None` program, support or black holes
a null `std::shared_ptr` (`edge_weights` and `source_weights` are shared arrays, `std::shared_ptr<const
std::vector<Amount>>`). `CentralizedTE` has no default constructor: it is built from its
`SourceAllocation`. The epoch's amounts are `Counter128` values. A Python `ValueError` is a
`std::invalid_argument` in C++. `spf` takes a trailing C++-only `SpfQueue` argument that chooses
the search's internal queue; the results never depend on it. `apply_actions` takes a
`std::span<const Action>`, where `Action` is a `std::variant<Place, Replace, Forward>` whose
values borrow their policy (a `std::reference_wrapper<FlowPolicy>`) and black holes (a `const
Eligibility*`, null for none), which must outlive the call.
