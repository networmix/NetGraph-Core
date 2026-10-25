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

## How it fits together

A demand (traffic from a set of sources to a set of destinations) goes through six things: a
`Graph` on integer units, a `FlowLedger` of what its edges carry and what tunnels hold, a
`FlowPolicy` that routes the demand with one routing model and one admission, the
`ForwardingProgram` the model computes, the forwarding of the offer over that program, and the
`PlacementResult` that reports what was delivered and lost. [docs/semantics.md](docs/semantics.md)
explains each with worked examples and gives the exact rules.

## Overview

- **Exact integer amounts**: a graph carries an explicit resolution `q` (a power of two in
  `[2^-32, 2^32]`); capacities and flows are `int64` units of `q`. Physical values enter only
  through `to_units` (the binary64 value rounded half up to whole units, exactly) and leave through
  `as_float`.
- **Positive costs**: edge costs are integers `>= 1`, so the shortest-path DAG (every edge on
  some shortest path) has no cycle; for a fixed set of usable edges it does not depend on search
  order or tie rules.
- **Zero capacity is not absence**: an edge of capacity 0 stays in routing, as a configured but
  empty link does: searches, routes and rosters use it, an equal split still sends it a share
  (which is lost), a weighted split gives it weight 0. Zero capacity alone never excludes an edge:
  an `Eligibility` excludes it, the route choice can keep another of its parallel members, and an
  explicit residual adds a capacity threshold.
- **Edge ids by sort**: edges are numbered by sorting on `(source, cost, destination)`, ties kept
  in input order; that number is the public `EdgeId`, and `input_positions()` maps it back to the
  input position.
- **One flow core**: an exact integer minimum-cost flow solver (the most traffic, sent as cheaply
  as possible). A solve of the ideal controller reports its marginal cost curve (the extra minimum
  total cost per additional delivered unit, as `(cost, units)` segments) and, on request, a
  minimum cut when it reaches the maximum ([the result](docs/semantics.md#the-result)).

## Core Features

The exact behaviour of every call (units and rounding, edge ids, which route a search returns,
how the route allocations split traffic, outcomes and flags, failures, the cache and the evaluators) is
specified in [docs/semantics.md](docs/semantics.md).

### Graph

- **`Graph(num_nodes, edge_sources, edge_destinations, edge_capacity_units, edge_costs,
  resolution)`**: an immutable directed multigraph (parallel edges and self-loops allowed).
  Validation reports the first defect in a fixed order. Accessors (`edge_capacity_units()`,
  `edge_costs()`, `edge_sources()`, `edge_destinations()`, `input_positions()`,
  `outgoing_offsets()`, `incoming_offsets()`, `incoming_edges()`) return read-only NumPy views
  that keep the graph alive ([the graph](docs/semantics.md#2-the-graph)).
- **`Eligibility`**: which nodes and edges of one graph may be used (immutable), built by
  `graph.eligibility(excluded_nodes, excluded_edges)`; `Eligibility.intersect(a, b)` combines two.

### Route structures

- **`PredecessorDAG`**: the route structure a search returns, a predecessor DAG (directed, no
  cycle): for each node, the edges that reach it on the route and the node each comes from;
  `PredecessorDAG.from_csr(graph, parent_offsets, parents, via_edges)` builds and validates one.
- **`Path(graph, hops)`**: a route given hop by hop, where a hop may hold several parallel edges
  between the same two nodes; `cost()` is the sum over hops of the cheapest member's cost.
- **Route choice**: one tie rule, `RouteTieBreak` (`FEWER_HOPS`, the default, `HIGHER_CAPACITY`,
  `HIGHER_AVAILABLE_CAPACITY`, `HIGHER_HEADROOM`), with `multipath` (every equal-cost route, or the
  one canonical path) and `parallel_edges` (every parallel link, or one per node pair). Cost always
  comes first; the tie rule only orders equal-cost choices
  ([the route choice](docs/semantics.md#the-route-choice)).

### Ledger and route allocations

- **`FlowLedger(graph, service_baseline=None, *, reservation_baseline=None)`**: the integer ledger
  of one graph, in two projections that are never added: the carried traffic (`carried()`,
  `residual()`) and the strict holds of T1 tunnels (`reserved()`, `unreserved()`). Holds never
  police traffic. Its `FlowRecord`s are the routed amounts of its policies
  ([the ledger](docs/semantics.md#6-the-ledger)).
- **Route allocations** (`RouteAllocation`) say how a demand's traffic is split at every router of
  its routes: `EQUAL` (ECMP, round robin), `WEIGHTED` (WCMP, highest averages over integer weights
  from a `WeightBasis`) or `OPTIMIZED` (a controller restricted to those routes). **Admission**
  (`Admission`) says what happens to traffic that does not fit: `LOSSLESS` admits only what the
  split carries without loss, `LOSSY` sends everything and records where it is lost
  ([route allocations and admission](docs/semantics.md#5-route-allocations-and-admission)).

### Policies

- **`FlowPolicy(ledger, sources, destinations, config, eligibility=None, cache=None)`**: one
  demand, bound to its ledger. `FlowPolicyConfig(*, routing, admission=LOSSLESS)` holds the routing
  model, which holds every option it reads. The models:
  - `ShortestPaths(...)`: static shortest paths, as routers forward with an IGP such as OSPF.
    `ShortestPaths()` is ECMP and `allocation=WEIGHTED` is WCMP; with either, the sources of one
    demand are pooled: where their routes meet, a router splits their traffic once.
    `allocation=OPTIMIZED` is a controller restricted to those paths, solving the sources jointly.
  - `DistributedTE(TunnelCount(N), ...)` (T1): every head sets up its own N strict-bandwidth
    tunnels by CSPF, as RSVP-TE head ends do, with no controller
    ([T1](docs/semantics.md#t1-tunnels-set-up-by-their-heads)).
  - `CentralizedTE(*, source_allocation, ...)`: an ideal controller, the minimum-cost flow of the
    demand over every eligible edge or a support DAG. With `tunnel_count=TunnelCount(N)` (T2) it
    lays out N single-path tunnels per head with one common whole-unit planned bandwidth
    ([T2](docs/semantics.md#t2-the-controllers-equal-tunnels)).
  - `Pinned(program, *, black_holes=None)`: forwarding over a given `ForwardingProgram`, computing
    nothing: repaired under the policy's Eligibility, or stale with black holes
    ([forwarding programs](docs/semantics.md#forwarding-programs)).
- **Source allocation** (`SourceAllocation`): `EQUAL` (every source sends the same round-robin
  share), `OPTIMIZED` (a controller chooses how much enters where) or `WEIGHTED` (`CentralizedTE`
  without tunnels: highest-averages shares of declared source weights).
- **`ForwardingProgram`**: the routes and integer weights of one demand. Every computing model
  produces one (`policy.program`, `result.program`); `ForwardingProgram.from_routes` builds one
  from given routes, and `ForwardingProgram.from_dag` one shared hop-by-hop program on a DAG.
- **Calls**: `place(offer)` adds to the demand (T1 and T2 have none); `evaluate(offer)` computes a
  whole new offer without publishing it, `replace(offer)` publishes it, `forward(offer)` sends it
  over the installed program, and `remove()` takes the demand's traffic out. The first four return
  a `PlacementResult` (`remove()` returns `None`): this call's delivery and losses, the demand's
  totals, per-source and per-tunnel rows, the program and, for T1 and T2, their `TunnelPlan`
  ([policies](docs/semantics.md#7-policies)).
- **`apply_actions(ledger, actions, *, mode=BatchMode.ORDERED, outputs=None)`**: the `Place`,
  `Replace` and `Forward` actions of several policies of one ledger as one transaction, one after
  another (`ORDERED`) or at the same time with the losses on shared links divided between the
  demands (`SIMULTANEOUS`) ([several demands at once](docs/semantics.md#several-demands-at-once)).
- **Failure guarantee**: every public call either completes or raises one exception and leaves
  every ledger, policy and cache as it was, with two stated limits (running out of memory inside
  pybind11's own object construction, and a `ShortestPathCache` passed to the call)
  ([failures, threads and lifetimes](docs/semantics.md#8-failures-threads-and-lifetimes)).

Native computation releases the GIL. Calls on one ledger take effect one at a time; different
ledgers run in parallel. The free-threaded Python build follows the same rules.

### Search

- `spf(graph, roots, SpfOptions(direction, *, targets=None, tie_break=FEWER_HOPS, multipath=True,
  parallel_edges=True, residual=None, min_residual=1, eligibility=None))` searches from a root
  set, forward or reverse, and returns the shortest-path DAG (to every node, or to the nearest
  targets) or, with `multipath=False`, the canonical single path
  ([the canonical single path](docs/semantics.md#the-canonical-single-path)). Edges of capacity 0
  are searched like any other; only a supplied residual adds a capacity threshold.
- `ksp(graph, source, destination, KspOptions(k, *, tie_break=FEWER_HOPS, parallel_edges=True,
  cost_ceiling=None, residual=None, eligibility=None))` returns the K shortest simple paths as
  `(cost, Path)` pairs, in the order of [docs/semantics.md](docs/semantics.md#k-shortest-paths),
  up to an inclusive cost ceiling (`AbsoluteCost(c)`, or `CostStretch(f)` with `f` an exact
  fraction: a float means the decimal it prints, so `CostStretch(1.2)` is exactly 6/5).
- `resolve_to_paths(graph, dag, source, destination, split_parallel_edges=False, max_paths=None)`
  lists a DAG's walks as `Path` values.

### ShortestPathCache

`ShortestPathCache(graph)` keeps the static shortest-path DAGs (from costs and the Eligibility,
not traffic) that policies use, keyed by the root set, the direction and the Eligibility. In
building mode a miss searches and inserts; `freeze()` makes it read-only and shareable by every
worker of an evaluator, and a miss then raises `ShortestPathCacheMiss`
([the shortest-path cache](docs/semantics.md#9-the-shortest-path-cache)).

### Evaluators

All three take one options type, `EvaluationOptions(config, *, eligibility=None, cache=None,
outputs=None, design_demand=None, offer=None)`, and evaluate on a fresh ledger at full capacity:
`evaluate_demand(graph, sources, destinations, options)` one demand, `evaluate_pairs(graph, pairs,
options, per_pair_eligibility=None, workers=None)` many independent pairs, and
`evaluate_edge_failures(graph, sources, destinations, options, workers=None)` the signed
`delivery_reduction` when each edge is removed, at the same offer and design.
Without an offer, `CentralizedTE(source_allocation=OPTIMIZED)` gives the exact maximum flow with
its marginal cost curve, and its minimum cut with `outputs=OutputOptions(min_cut=True)`; a finite
offer below the maximum reports no cut ([evaluators](docs/semantics.md#10-evaluators)).

### Flow core

One native integer minimum-cost flow core serves `CentralizedTE` routing, the `OPTIMIZED`
allocation, the evaluators and T2's sizing (there with every cost taken as 0: a maximum flow): the
largest flow up to a cap from a source set to a destination set over the eligible edges, sent at
the lowest total cost.

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

# ECMP over the same graph: only the cost-2 path is shortest; 2 of the 12 are dropped on its first edge.
ecmp = ngc.FlowPolicyConfig(routing=ngc.ShortestPaths(), admission=ngc.Admission.LOSSY)
r = ngc.FlowPolicy(ngc.FlowLedger(graph), 0, 3, ecmp).place(ngc.to_units(12.0, q))
print(ngc.as_float(r.call_delivered, q), r.call_edge_drops)  # 10.0 [(0, 8192)]

# T1: two tunnels sized for 12; the second finds no path with 6 units unreserved.
t1 = ngc.FlowPolicyConfig(routing=ngc.DistributedTE(ngc.TunnelCount(2)), admission=ngc.Admission.LOSSY)
policy = ngc.FlowPolicy(ngc.FlowLedger(graph), 0, 3, t1)
r = policy.replace(ngc.to_units(12.0, q), design_demand=ngc.to_units(12.0, q))
print([t.status.name for t in r.program.tunnels])  # ['UP', 'INSUFFICIENT_BANDWIDTH']
print(ngc.as_float(r.call_delivered, q))                    # 10.0
print(ngc.as_float(r.plan.total_reserved_bandwidth, q))       # 6.0

# T2: 16 tunnels of one common size, 11 on the 10-unit path and 5 on the 5-unit one.
# The 10-unit path (40960 units) bounds the size at 40960/11, whose whole-unit floor is 3723;
# 59575 is the largest total whose round-robin split over the 16 tunnels still fits.
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
(`make perf-tool-test` tests the tool itself).

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
