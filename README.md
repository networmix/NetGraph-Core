# NetGraph-Core

C++ graph engine for network flow analysis, traffic engineering simulation, and capacity planning.

> **Unreleased, breaking.** This version replaces the 0.11 API (see the Unreleased section of
> [CHANGELOG.md](CHANGELOG.md)): the graph on integer units, the minimum-cost flow core, the route
> structures (`PredecessorDAG`, `Path`), the integer ledger `FlowLedger` with the route allocations, the
> policy engine `FlowPolicy` with its four routing models and its lists of several policies' actions
> (`apply_actions`), forwarding programs, the searches (`spf`,
> `ksp`, `resolve_to_paths`), the `ShortestPathCache` in building and frozen modes, and the evaluators
> `evaluate_demand`, `evaluate_pairs` and `evaluate_edge_failures`.
> The previous API (0.11 and earlier) is in git history and on PyPI.

## Overview

- **Exact integer amounts**: a graph carries an explicit resolution `q` (a power of two in
  `[2^-32, 2^32]`); capacities and flows are `int64` units of `q`. Physical values enter only
  through `to_units` (round half up) and leave through `as_float`.
- **Positive costs**: edge costs are integers `>= 1`, so the shortest-path DAG (every edge on
  some shortest path) has no cycle and does not depend on search order or tie rules.
- **Zero capacity is not absence**: an edge of capacity 0 stays in routing, as a configured but
  empty link does: searches, routes and rosters use it, an equal split still sends it a share
  (which is lost), a weighted split gives it weight 0. Only an explicit residual or an
  `Eligibility` removes an edge from a search.
- **Stable edge ids**: edges are numbered by sorting on `(source, cost, destination)`, ties kept in input
  order; that number is the public `EdgeId`, and `input_positions()` maps it back to the input position.
- **One flow core**: an exact integer minimum-cost flow solver (the most traffic, sent as cheaply
  as possible). It returns the flow value, the marginal cost curve (the extra minimum total cost
  per additional delivered unit, as `(cost, units)` segments), the flow on every edge, a minimum
  cut (the full edges that separate the sources from the destinations) and the node potentials
  (per-node prices that prove the flow is cheapest). It works one cost level at a time: shortest
  paths with potentials find the next cheapest cost, then Dinic phases fill every path of that
  cost.

## Core Features

The exact behaviour of every call (units and rounding, edge ids, which route a search returns,
how the route allocations split traffic, outcomes and flags, failures, the cache and the evaluators) is
specified in [docs/semantics.md](docs/semantics.md).

### Graph

- **`Graph(num_nodes, edge_sources, edge_destinations, edge_capacity_units, edge_costs,
  resolution)`**: immutable directed multigraph that stores each node's out-edges and in-edges as
  arrays; parallel edges and self-loops are allowed. Validation reports the first defect in a
  fixed order (the resolution, the unit values, the shapes, zero or negative costs, the remaining
  bounds). Accessors (`edge_capacity_units()`, `edge_costs()`, `edge_sources()`,
  `edge_destinations()`, `input_positions()`, `outgoing_offsets()`, `incoming_offsets()`,
  `incoming_edges()`) return read-only NumPy views that keep the graph alive.
- **`Eligibility`**: which nodes and edges of one graph may be used (immutable), built by
  `graph.eligibility(excluded_nodes, excluded_edges)`; `Eligibility.intersect(a, b)` returns an
  operand when the result equals it, otherwise a new value.

### Route structures

- **`PredecessorDAG`**: the route structure a search returns, a predecessor DAG (directed, no
  cycle): for each node, the edges that reach it on the route and the node each comes from.
  `parents()` and `via_edges()` list these entries, in ascending EdgeId within a node;
  `parent_offsets()` (built on each call) says where each node's entries start.
  `PredecessorDAG.from_csr(graph, parent_offsets, parents, via_edges)` validates the shape, the
  ranges and endpoints, duplicates, acyclicity and the canonical order before accepting a DAG.
- **`Path(graph, hops)`**: a route given hop by hop, where a hop may hold several parallel edges
  between the same two nodes; `cost()` is the sum over hops of the minimum member cost,
  `source()` and `destination()` its ends, and `to_dag()` gives the equivalent `PredecessorDAG`.
- **Route choice**: one tie rule, `RouteTieBreak` (`FEWER_HOPS`, the default, `HIGHER_CAPACITY`,
  `HIGHER_AVAILABLE_CAPACITY`, `HIGHER_HEADROOM`), with `multipath` (every equal-cost route, or the
  one canonical path) and `parallel_edges` (every parallel link, or one per node pair). Cost always
  comes first; the tie rule only orders equal-cost choices.

### Ledger and route allocations

- **`FlowLedger(graph, service_baseline=None, *, reservation_baseline=None)`**: the integer ledger
  of one graph, in two projections that are never added: the carried traffic (`carried()` is the
  exact sum of the records' edge flows and `residual()` the service baseline minus it) and the
  strict holds of T1 tunnels (`reserved()` is the sum of the records' holds and `unreserved()` the
  reservation baseline minus it). Holds never police traffic. The records (`FlowRecord`: flow id,
  `policy_id`, endpoints, admission, `entry_id`, `route`, `min_route_cost`, offered, delivered,
  `edge_flows`, `edge_drops` and `edge_reservations`) are the flows placed on it;
  `saturated_boundary_edges(sources)` lists the full edges leaving the set of nodes the sources
  can still reach. `reset()` and `copy()` are rejected while a policy is bound to the ledger.
  Policies place, hold and remove; `reset(service_baseline=None, *, reservation_baseline=None)`
  replaces an unbound ledger's baselines and clears its records.
- **Route allocations**, the rules that spread one flow over its route (`RouteAllocation`: `OPTIMIZED`,
  `WEIGHTED`, `EQUAL`; `Admission`: `LOSSLESS`, `LOSSY`): `OPTIMIZED` re-solves the flow with the
  flow core on its route, over the free capacity plus its own allocation; `WEIGHTED` splits at
  every node in proportion to integer weights frozen when the route is prepared, from a
  `WeightBasis` (`CAPACITY`, `EXPLICIT`, `RESIDUAL_SNAPSHOT`) (highest averages (D'Hondt): assign
  each next unit to the largest weight / (assigned + 1), ties to the lower EdgeId); `EQUAL` splits equally (round robin); `LOSSY` forwards hop
  by hop and records what it drops.

### Policies

- **`FlowPolicy(ledger, sources, destinations, config, eligibility=None, cache=None)`**: the
  placement engine for one demand, bound to its ledger. `FlowPolicyConfig(*, routing,
  admission=LOSSLESS)` holds the routing model and the admission; every other option belongs to
  the model that reads it. A model's source allocation is `EQUAL` (every source sends the same
  round-robin share) or `OPTIMIZED` (a controller chooses how much enters where). The routing is
  one of:
  - `ShortestPaths(*, allocation=EQUAL, source_allocation=None, weight_basis=None,
    edge_weights=None, tie_break=FEWER_HOPS, multipath=True, parallel_edges=True)`: static
    shortest paths, as routers forward with an IGP such as OSPF. `ShortestPaths()` is ECMP;
    `allocation=WEIGHTED` is WCMP with its `weight_basis` and `edge_weights`; the sources of one
    demand are pooled (where their routes meet, a router splits their traffic once);
    `allocation=OPTIMIZED` is a controller restricted to those paths;
  - `DistributedTE(TunnelCount(N), *, cost_ceiling=None, min_unreserved_bandwidth=0,
    tie_break=FEWER_HOPS, parallel_edges=True)` (T1): every head sets up its own N
    strict-bandwidth tunnels, as RSVP-TE head ends do, with no controller: a design volume is split
    equally over the heads and their tunnels, each tunnel runs CSPF (the shortest path among the
    links with enough unreserved bandwidth, within an `AbsoluteCost` or `CostStretch` cost ceiling)
    and either comes up holding its request on every link or stays down
    (`INSUFFICIENT_BANDWIDTH`, `COST_CEILING_EXCEEDED`); traffic then enters equally at the heads
    and splits equally over their up tunnels;
  - `CentralizedTE(*, source_allocation, support=None, tunnel_count=None, prefer_headroom=False)`:
    without a tunnel count, an ideal controller: the minimum-cost flow from all the policy's
    sources to its destinations, on the free capacity plus its own allocation, optionally only
    over the edges of a support DAG. With `tunnel_count=TunnelCount(N)` (T2), a central controller
    (as a PCE) lays out, for every head that reaches a destination over free capacity, N
    single-path tunnels with one common whole-unit planned bandwidth: as large as possible, then
    (for a finite design) no larger than needed, then the cheapest routes, then with
    `prefer_headroom=True` the routes that leave the most spare capacity, then one fixed tie;
    `TunnelCount(L, up_to=True)` chooses the smallest common count up to L whose program carries
    the live heads' design shares (if none does, the largest live-share score, then the smallest
    count; it does not promise whole-demand lossless delivery); a `support` restricts the tunnels to its edges. The planned bandwidth is never held.
    The source allocation is always stated. NetGraph's TE presets map onto it: `TE_ECMP_16_LSP` is
    `CentralizedTE(source_allocation=EQUAL, tunnel_count=TunnelCount(16))` and
    `TE_ECMP_UP_TO_256_LSP` the same with `TunnelCount(256, up_to=True)`, both with
    `prefer_headroom=True`;
  - `Pinned(program, *, black_holes=None)`: forwarding over a given `ForwardingProgram`, computing
    nothing; the policy's Eligibility fails the tunnels it touches (their heads split over the
    survivors) and `black_holes` makes failed elements swallow what reaches them (stale
    forwarding).
- **`ForwardingProgram`**: the routes and integer weights of one demand. Every computing model
  produces one (`policy.program`, `result.program`); `ForwardingProgram.from_routes(graph,
  sources, destinations, routes)` builds one from `RouteDefinition(head, route, weight=1)` entries
  (a `Path` is a tunnel, a `PredecessorDAG` a split route of its own), and
  `ForwardingProgram.from_dag(graph, sources, destinations, dag)` builds one shared hop-by-hop
  program on a DAG whose sources are pooled exactly as `ShortestPaths` pools them.
- **`place(offer, outputs=None)`** adds `offer` to the policy's request as one all-or-nothing step
  (`ShortestPaths` forwards over the program it prepares once per epoch; the controllers re-solve
  what they already deliver plus this request; T1 and T2 have no `place`) and returns a
  `PlacementResult`: this call's `call_delivered`, `call_undelivered`, `delivery_delta`,
  `call_edge_drops`, `call_source_drops`, `termination` and flags; the demand's `requested`,
  `offered`, `admitted`, `delivered`, `unadmitted`, `source_drops`, `edge_drops`, per-source rows
  `by_source`, `total_cost` (controllers), `marginal_cost_curve` (optimized `CentralizedTE`), the
  `program`, the tunnel models' `plan`, and per-tunnel rows `by_tunnel`; and, when requested
  through `OutputOptions`, `edge_flows`, `residual`, `saturated_boundary_edges`,
  `edge_reservations` and `min_cut` (optimized `CentralizedTE`). `evaluate(offer, *,
  design_demand=None)` computes a whole new offer without publishing it, `replace(offer, *,
  design_demand=None)` publishes it (T1 sets its tunnels up and T2 lays them out for the design,
  `None` meaning the offer), `forward(offer, *, black_holes=None)` sends it over the installed
  program (T1 keeps its tunnels and holds, T2 its tunnels and plan); `remove()` reverts the
  policy's flows and holds; `records()`, `record(flow_id)`, `delivered()`, `config` and `program`
  read its state.
- **`apply_actions(ledger, actions, *, mode=BatchMode.ORDERED, outputs=None)`**: the actions of
  several policies of one ledger (`Place(policy, offer)`, `Replace(policy, offer, *,
  design_demand=None)`, `Forward(policy, offer, *, black_holes=None)`, immutable values) as one
  transaction, validated whole and published once, one `PlacementResult` per action in input
  order. `BatchMode.ORDERED` runs each as its scalar call on the state the earlier ones left
  (earlier actions are served first); `BatchMode.SIMULTANEOUS` forwards LOSSY `Replace` actions of
  `ShortestPaths` `EQUAL` or `WEIGHTED` and `Pinned` and `Forward` actions of installed programs at
  the same time on one view (the residual plus every participant's old traffic), a full edge
  shared by highest averages over the demands, then their heads, then their entries
  ([docs/semantics.md](docs/semantics.md#several-demands-at-once)).
- **Failure guarantee**: every public call either completes or raises one exception and leaves
  every ledger, policy, cache and scratch buffer as it was, whatever fails (validation, memory, a
  worker thread, building the Python result). Two limits: pybind11 does not check every
  allocation when it constructs a bound object from Python, parses the arguments of a call that
  passes keywords, or converts a returned value into a new Python object outside the library's
  checked result construction, so running out of memory
  exactly there may crash the process (Core state stays unchanged); and a `ShortestPathCache`
  passed to the call: a building cache keeps the complete entries and statistics of the lookups
  that succeeded before the failure, and a frozen cache keeps its entries and statistics
  unchanged, except that a call raising `ShortestPathCacheMiss` counts one miss
  ([docs/semantics.md](docs/semantics.md#8-failures-threads-and-lifetimes)).

### Search

`spf(graph, roots, SpfOptions(direction, *, targets=None, tie_break=FEWER_HOPS, multipath=True,
parallel_edges=True, residual=None, min_residual=1, eligibility=None))` searches from a root set,
forward or reverse. Edges of capacity 0 are searched like any other; only a supplied residual
filters by capacity.
- **Without targets:** it returns `costs` (`UNREACHABLE` for nodes not reached) and the full
  DAG of tight edges, the edges that lie on some shortest path.
- **With targets:** it stops at the nearest targets and returns `target_cost` and the DAG of
  the walks between the roots and the nearest targets (for a reverse search, from a target to
  a root).
- **With `multipath=False`:** it returns the canonical single path: the walk of lowest cost;
  among those, under a `HIGHER_*` tie rule, the one of highest bottleneck (the smallest score
  along the walk); then the fewest hops; then the lexicographically smallest sequence of EdgeIds
  ([docs/semantics.md](docs/semantics.md#the-canonical-single-path)).
- **With `parallel_edges=False`:** each node pair is first reduced to one parallel member (lowest
  cost, then under a `HIGHER_*` tie rule the highest score, then the lowest EdgeId).
- A `HIGHER_*` tie rule with `multipath=True` and `parallel_edges=True` has nothing to choose and
  is rejected (`spf tie_break requires multipath=False or parallel_edges=False`).

`ksp(graph, source, destination, KspOptions(k, *, tie_break=FEWER_HOPS, parallel_edges=True,
cost_ceiling=None, residual=None, eligibility=None))` returns the K shortest simple paths as
`(cost, Path)` pairs, in the order of [docs/semantics.md](docs/semantics.md#k-shortest-paths) and
up to an inclusive cost ceiling (`AbsoluteCost(c)`, or `CostStretch(f)` with `f` an exact
fraction: a float means the decimal it prints, so `CostStretch(1.2)` is exactly 6/5).
`resolve_to_paths(graph, dag, source, destination, split_parallel_edges=False, max_paths=None)`
lists a DAG's walks as `Path` values.

### ShortestPathCache

`ShortestPathCache(graph)` keeps the static shortest-path DAGs (from costs and the Eligibility,
not traffic) that policies build, keyed by the root set, the direction and the Eligibility. Values
are immutable and shared, and `stats()` reports hits, misses, entries, retained bytes and live
values.
- **Building mode (new caches):** one writer; a miss searches and inserts. An evaluator given a
  building cache runs on one worker (an explicit `workers` above 1 is rejected).
- **Frozen mode (`freeze()`):** read-only and shared by every worker of an evaluator (pass
  `workers` explicitly; with a cache the default is 1); a miss raises `ShortestPathCacheMiss`.
  `thaw()` returns to building mode once every reader has finished.

### Evaluators

All three take one options type, `EvaluationOptions(config, *, eligibility=None, cache=None,
outputs=None, design_demand=None, offer=None)`. `evaluate_demand(graph, sources, destinations,
options)` evaluates a fresh policy on a fresh ledger at full capacity: with no `offer`, the
uncapped question (Lossless) or, for `ShortestPaths` `EQUAL`/`WEIGHTED` Lossy, the first-hop
experiment; with a finite `offer`, that offer, for every model (T1 always names a finite
`design_demand`; T2 takes the offer as its design when none is given).
`CentralizedTE(source_allocation=OPTIMIZED)` gives the exact maximum flow with its marginal cost
curve, and its minimum cut with `outputs=OutputOptions(min_cut=True)`. `evaluate_pairs(graph,
pairs, options, per_pair_eligibility=None, workers=None)` evaluates independent pairs in pair order
(the lowest failing index raises; Pinned, whose program fixes its endpoints, is rejected).
`evaluate_edge_failures(graph, sources, destinations, options, workers=None)` reports, for every
candidate edge, the signed `delivery_reduction` when the edge is removed, at the same offer and
design. With `workers=None` the worker count is 1 when a cache is supplied, otherwise it comes
from `NGRAPH_CORE_BATCH_THREADS` and `NGRAPH_CORE_SENSITIVITY_THREADS`.

### Flow core

One native integer minimum-cost flow core serves `CentralizedTE` routing, the `OPTIMIZED` allocation,
the evaluators and T2's sizing (there with every cost taken as 0: a maximum flow): the largest flow
up to a cap from a source set to a destination set over the eligible edges, sent at the lowest
total cost, starting from zero flow. For `evaluate_edge_failures` of unrestricted `CentralizedTE`
with `OPTIMIZED` source allocation it repairs the baseline's solution for each removed edge instead
of solving again, when that baseline is a maximum flow (no offer, or an offer above it); at an
offer that binds each deletion is solved at the offer. Its private Python entry
`_netgraph_core._solve_min_cost_flow(graph, view, sources, destinations, eligibility, cap)` returns
the payload `{flow_value, marginal_cost_curve, edge_flows, min_cut, potentials, counters}` for the flow
checks of the test suites.

## Quick Start

```python
import netgraph_core as ngc

q = 2.0**-12
# Two parallel paths from node 0 to node 3 (physical capacities 10 and 5).
caps = ngc.to_units([10.0, 10.0, 5.0, 5.0], q)
graph = ngc.Graph(4, [0, 1, 0, 2], [1, 3, 2, 3], caps, [1, 1, 2, 2], q)

config = ngc.FlowPolicyConfig(routing=ngc.CentralizedTE(source_allocation=ngc.SourceAllocation.OPTIMIZED))
result = ngc.evaluate_demand(graph, 0, 3, ngc.EvaluationOptions(config))
print(ngc.as_float(result.call_delivered, q))   # 15.0
print(result.marginal_cost_curve)         # [(2, 40960), (4, 20480)]
print(result.total_cost)                # 163840

# T1: two tunnels sized for 12; the second finds no path with 6 units unreserved.
t1 = ngc.FlowPolicyConfig(routing=ngc.DistributedTE(ngc.TunnelCount(2)), admission=ngc.Admission.LOSSY)
policy = ngc.FlowPolicy(ngc.FlowLedger(graph), 0, 3, t1)
r = policy.replace(ngc.to_units(12.0, q), design_demand=ngc.to_units(12.0, q))
print([t.status.name for t in r.program.tunnels])  # ['UP', 'INSUFFICIENT_BANDWIDTH']
print(ngc.as_float(r.call_delivered, q))                    # 10.0
print(ngc.as_float(r.plan.total_reserved_bandwidth, q))       # 6.0

# T2: 16 tunnels of one common size, 11 on the 10-unit path and 5 on the 5-unit one.
te16 = ngc.FlowPolicyConfig(routing=ngc.CentralizedTE(source_allocation=ngc.SourceAllocation.EQUAL,
                                                     tunnel_count=ngc.TunnelCount(16), prefer_headroom=True))
r = ngc.FlowPolicy(ngc.FlowLedger(graph), 0, 3, te16).replace(float("inf"))
print(r.plan.max_fractional_bandwidth_per_tunnel, r.plan.max_bandwidth_per_tunnel, r.call_delivered)  # 40960/11 3723 59575
```

## Installation

```bash
pip install netgraph-core
```

Or from source (the shipped configuration; `pip install -e .` is the editable development build,
which carries the private test hooks):

```bash
pip install .
```

The package is typed: it ships `py.typed` and a type stub generated from the compiled module.

### Build Optimizations

Default builds include LTO and loop unrolling. For local development:

```bash
make install-native   # CPU-specific optimizations (not portable)
```

## Repository Structure

```
src/                    # C++ implementation
include/netgraph/core/  # Public C++ headers
bindings/python/        # pybind11 bindings
python/netgraph_core/   # Python package
tests/cpp/              # C++ tests (googletest)
tests/py/               # Python tests (pytest)
dev/                    # Developer scripts (not shipped): checks, coverage,
dev/perf/               #   and the A/B/A timing tool
```

## Development

```bash
make dev           # Setup: venv, dependencies, pre-commit hooks
make check         # Rebuild the editable extension, then everything check-ci runs
make check-ci      # Lint, types, C++ tests (hooks on and off), Python tests (CI entrypoint)
make test          # Python tests
make cpp-test      # C++ tests (private hooks on)
make cpp-test-off  # C++ tests without the hooks (the shipped configuration)
make sanitize-test # C++ tests under ASan/UBSan
make tsan-test     # C++ tests under ThreadSanitizer
make cov           # Combined coverage report (C++ + Python; its own venv and builds)
make stubs-check   # The generated type stub: drift, package contents, typed clients
make format        # Apply ruff formatting (make hooks runs every pre-commit fixer)
```

`make check` and `make check-ci` only verify; they never edit a file.

Performance claims are timed with the A/B/A tool (old build, new build, then the old build again
on one machine, so drift shows), described in [dev/perf/README.md](dev/perf/README.md)
(`make perf-tool-test` tests the tool itself; CI runs it on one macOS job).

The editable development build (`make install`, `make dev`) carries the private test hooks
(CMake option `NETGRAPH_CORE_PRIVATE_HOOKS`; a wheel build does not); the test suites skip the
tests that need them when they are absent.

## Environment Variables

| Variable | Effect |
| --- | --- |
| `NGRAPH_CORE_BATCH_THREADS` | Default worker count of `evaluate_pairs` (an integer in [1, 4096]; 1 when unset or invalid; not read when a cache is supplied, which defaults to 1; an explicit `workers` wins). |
| `NGRAPH_CORE_SENSITIVITY_THREADS` | Default worker count of `evaluate_edge_failures` (as above). |

## Requirements

- **C++:** C++20 compiler (GCC 10+, Clang 12+, MSVC 2019+)
- **Python:** 3.11+
- **Build:** CMake 3.23+, scikit-build-core
- **Dependencies:** pybind11, NumPy

## License

[MIT License](LICENSE)
