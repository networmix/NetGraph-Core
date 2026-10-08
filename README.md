# NetGraph-Core

C++ graph engine for network flow analysis, traffic engineering simulation, and capacity planning.

> **Streamline in progress.** The `streamline` branch rebuilds Core. This build carries the
> graph on integer units, the minimum-cost flow core, the route structures (`PredDAG`, `Path`),
> the integer ledger `FlowLedger` with the placement kernels, the policy engine `FlowPolicy`, the
> `DagCache` (building mode) and the evaluators `capacity`, `batch_capacity` and `sensitivity`.
> The public search API (shortest paths, k shortest paths, path resolution) and the frozen cache
> are not part of it yet. The previous API (0.11 and earlier) is in git history and on PyPI.

## Overview

- **Exact integer amounts**: a graph carries an explicit resolution `q` (a power of two in
  `[2^-32, 2^32]`); capacities and flows are `int64` units of `q`. Physical values enter only
  through `to_units` (round half up) and leave through `as_float`.
- **Positive costs**: edge costs are integers `>= 1`, so the shortest-path DAG is unique.
- **Stable edge ids**: the public `EdgeId` is the rank of an edge in the order
  `(src, cost, dst, input)`; `input_index()` maps ids back to input positions.
- **One flow core**: an integer minimum-cost flow (successive shortest paths with potentials,
  Dinic phases per cost tier) returning the flow value, the marginal cost curve, per-edge
  flows, the residual minimum cut and the dual potentials.

## Core Features

### Graph

- **`Graph(num_nodes, src, dst, capacity_units, cost, resolution)`**: immutable
  directed multigraph with forward and reverse CSR adjacency; parallel edges and self-loops are
  allowed. Validation reports the first defect in a fixed order (the resolution, the unit
  values, the shapes, zero or negative costs, the remaining bounds). Accessors (`capacity_units()`, `cost()`, `edge_src()`, `edge_dst()`,
  `input_index()`, `row_offsets()`, `in_row_offsets()`, `in_adj_edge_index()`) return read-only
  NumPy views that keep the graph alive.
- **`Eligibility`**: immutable, graph-bound static eligibility built by
  `graph.eligibility(excluded_nodes, excluded_edges)`; `Eligibility.intersect(a, b)` returns an
  operand when the result equals it, otherwise a new value.

### Route structures

- **`PredDAG`**: immutable, graph-bound predecessor DAG, stored over its member nodes and
  presented in CSR form (`parent_offsets()` built on demand, `parents()`, `via_edges()`),
  entries in ascending via EdgeId within a node.
  `PredDAG.from_csr(graph, offsets, parents, via)` validates the shape, the ranges and endpoints,
  duplicates, acyclicity and the canonical order before accepting a DAG.
- **`Path(graph, hops)`**: an explicit walk whose hops may hold several parallel members;
  `cost()` is the sum over hops of the minimum member cost and `to_dag()` gives the equivalent
  `PredDAG`.
- **`PathSelection(tie_break, priority, multipath=True, multi_edge=True)`** with
  `EdgeTieBreak` and `EdgePriority`: the selection axes the searches and policies read.

### Ledger and kernels

- **`FlowLedger(graph, residual=None)`**: the integer ledger of one graph: an immutable baseline
  in units and the records (`FlowRecord`: flow id, owner, endpoints, placement, admission,
  bundle, route cost, offered, delivered, per-edge allocation and drops) of the flows placed on
  it. `allocated()` is the exact sum of the records and `residual()` the baseline minus it;
  `bottleneck(sources)` lists the saturated frontier; `reset()` and `copy()` follow the lifetime
  rules (rejected while a policy is bound). The ledger is mutated only through policies.
- **Kernels** (`FlowPlacement`: `MAX_FLOW`, `WEIGHTED`, `EQUAL_BALANCED`; `Admission`:
  `LOSSLESS`, `LOSSY`): MaxFlow is the flow core on the flow's DAG (a replacement on its own
  allocation plus the residual), Weighted splits by highest averages over weights frozen at
  creation, EqualBalanced by round robin; Lossy forwards hop by hop and records drops.

### Policies

- **`FlowPolicy(ledger, sources, destinations, config, elig=None, cache=None)`**: the one
  placement engine, bound to its ledger. `FlowPolicyConfig(routing=..., placement=...,
  admission=..., selection=...)` has four explicit axes; the routing is `CostOnly()` (static
  shortest paths), `CapacityAware(flow_count=1, max_flow_count=None, max_path_cost=None,
  min_flow_bandwidth=0, reoptimize=False, diminishing=None)` (residual-aware searches with an
  `Absolute` or `Stretch` cost ceiling, a `DiminishingReturns(window, fraction)` guard, and at
  most `MAX_FLOWS_PER_CALL` (10,000) creations per call), `Global()` (the minimum-cost flow of the policy's
  sources to its destinations on the residual plus its own allocation) or `Pinned(bundles)`
  (fixed `PredDAG`s).
- **`place(volume, outputs=None)`** runs one transaction per call (per source for CostOnly and
  CapacityAware, cumulative on Global) and returns a `PlacementResult` (`placed`, `remaining`,
  `termination`, `local_outcomes`, `routes_changed`, `weights_changed`, `drops` and the
  optional `cost_distribution` (Global) and, when requested through `OutputOptions`, `edge_flows`,
  `residual`, `bottleneck` and, on Global with `OutputOptions(bottleneck=True)`, `min_cut`).
  `remove()` reverts the policy's flows; `records()`, `record(flow_id)`, `placed_demand()` and
  `usable_bundles()` read its state.
- **Failure guarantee**: every public call is all-or-nothing up to its caller: a failure
  anywhere (validation, allocation, a worker thread, the construction of the Python result)
  leaves every ledger, policy, cache and scratch buffer as it was and raises one exception.
  Two limits: the pybind11 object machinery is outside the guarantee (constructing a bound
  object from Python, or a conversion into a new bound object or container outside the staged
  calls, allocates through pybind11's `make_new_instance`, which does not check `tp_alloc`, so a
  Python allocation failure there crashes the process, with Core state unchanged); and a
  supplied building `DagCache` is a memo that keeps the complete entries and statistics of the
  lookups that succeeded before a later failure of the call.

### Search

`spf(graph, roots, SpfOptions(direction, selection, targets=None, residual=None,
min_residual=1, elig=None))` searches from a root set, forward or reverse.
- **Without targets:** it returns `costs` (`UNREACHABLE` for nodes not reached) and the full
  DAG of tight edges, the edges that lie on some shortest path.
- **With targets:** it stops at the nearest targets and returns `target_cost` and the DAG of
  the walks between the roots and the nearest targets (for a reverse search, from a target to
  a root).
- **With `multipath=False`:** it returns the canonical single path, the one walk that minimises
  cost, then bottleneck, then hops, then edge ids.

`ksp(graph, src, dst, KspOptions(k, selection, max_path_cost=None, residual=None, elig=None))`
returns the K shortest simple paths as `(cost, Path)` pairs, in the contract's order and up to
an inclusive cost ceiling. `resolve_to_paths(graph, dag, src, dst, split_parallel_edges=False,
max_paths=None)` lists a DAG's walks as `Path` values.

### DagCache

`DagCache(graph)` keeps the static shortest-path supports that policies build, keyed by the
root set, the direction and the Eligibility. Values are immutable and shared, and `stats()`
reports hits, misses, entries, retained bytes and live values.
- **Building mode (new caches):** one writer; a miss searches and inserts.
- **Frozen mode (`freeze()`):** read-only and shared by every worker of an evaluator; a miss
  raises `DagCacheMiss`. `thaw()` returns to building mode once every reader has finished.

### Evaluators

`capacity(graph, sources, destinations, CapacityOptions(config, elig=None, cache=None,
outputs=None))` is `place(+inf)` of a policy on a fresh ledger; Global with MaxFlow gives the
exact maximum flow with its marginal cost curve, and its minimum cut with
`outputs=OutputOptions(bottleneck=True)`. `batch_capacity(graph, pairs,
opts, per_pair_elig=None, workers=None)` evaluates pairs in pair order (the lowest failing
index raises). `sensitivity(graph, sources, destinations, SensitivityOptions(...),
workers=None)` reports, for every candidate edge, the signed loss in placed units when the
edge is removed. The default worker counts come from `NGRAPH_CORE_BATCH_THREADS` and
`NGRAPH_CORE_SENSITIVITY_THREADS`.

### Flow core

One native integer minimum-cost flow core serves Global routing, the MaxFlow kernel and the
evaluators: the minimum-cost flow of the maximum value attainable under a cap from a source set
to a destination set over the eligible edges, starting from zero flow; the exact sensitivity
model re-solves each deletion from the baseline solve's final state. Its private Python entry
`_netgraph_core._place_global_units(graph, view, sources, destinations, elig, cap)` returns the
payload `{placed, cost_distribution, edge_flows, min_cut, potentials, counters}` for the
equivalence recorder and the oracle checks.

## Quick Start

```python
import netgraph_core as ngc

q = 2.0**-12
# Two parallel paths from node 0 to node 3 (physical capacities 10 and 5).
caps = ngc.to_units([10.0, 10.0, 5.0, 5.0], q)
graph = ngc.Graph(4, [0, 1, 0, 2], [1, 3, 2, 3], caps, [1, 1, 2, 2], q)

config = ngc.FlowPolicyConfig(
    routing=ngc.Global(),
    placement=ngc.FlowPlacement.MAX_FLOW,
    admission=ngc.Admission.LOSSLESS,
    selection=ngc.PathSelection(ngc.EdgeTieBreak.DETERMINISTIC, ngc.EdgePriority.CAPACITY),
)
result = ngc.capacity(graph, 0, 3, ngc.CapacityOptions(config))
print(ngc.as_float(result.placed, q))   # 15.0
print(result.cost_distribution)         # [(2, 40960), (4, 20480)]
```

## Installation

```bash
pip install netgraph-core
```

Or from source:

```bash
pip install -e .
```

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
dev/perf/               #   performance tooling and recorded A/B/A results,
dev/research/           #   research prototypes and their evidence
```

## Development

```bash
make dev        # Setup: venv, dependencies, pre-commit hooks
make check      # Run all tests and linting (auto-fix formatting)
make check-ci   # Strict checks without auto-fix (for CI)
make test       # Python tests
make cpp-test   # C++ tests only
make cov        # Combined coverage report (C++ + Python)
```

Performance work (harnesses, old/new/old protocol, recorded measurements) is
described in [dev/perf/README.md](dev/perf/README.md).

The editable development build (`make install`, `make dev`) carries the private test hooks
(CMake option `NETGRAPH_CORE_PRIVATE_HOOKS`; a wheel build does not); the test suites skip the
tests that need them when they are absent.

## Environment Variables

| Variable | Effect |
| --- | --- |
| `NGRAPH_CORE_BATCH_THREADS` | Default worker count of `batch_capacity` (an integer in [1, 4096]; 1 when unset or invalid; an explicit `workers` wins). |
| `NGRAPH_CORE_SENSITIVITY_THREADS` | Default worker count of `sensitivity` (as above). |

## Requirements

- **C++:** C++20 compiler (GCC 10+, Clang 12+, MSVC 2019+)
- **Python:** 3.11+
- **Build:** CMake 3.23+, scikit-build-core
- **Dependencies:** pybind11, NumPy

## License

[MIT License](LICENSE)
