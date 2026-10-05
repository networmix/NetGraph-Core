# NetGraph-Core

C++ graph engine for network flow analysis, traffic engineering simulation, and capacity planning.

> **Streamline in progress.** The `streamline` branch rebuilds Core in sub-PRs
> (`dev/design/2026-10-streamline/`). This build (sub-PR 1a) carries the graph on integer
> units and the minimum-cost flow core behind a private entry point; the ledger, placement
> kernels, policies, searches and evaluators return in the following sub-PRs. The previous API
> (0.12 and earlier) is in git history and on PyPI.

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

- **`StrictMultiDiGraph(num_nodes, src, dst, capacity_units, cost, resolution)`**: immutable
  directed multigraph with forward and reverse CSR adjacency; parallel edges and self-loops are
  allowed. Validation reports the first defect in a fixed order (unit values, shapes, costs,
  remaining bounds). Accessors (`capacity_units()`, `cost()`, `edge_src()`, `edge_dst()`,
  `input_index()`, `row_offsets()`, `in_row_offsets()`, `in_adj_edge_index()`) return read-only
  NumPy views that keep the graph alive.
- **`Eligibility`**: immutable, graph-bound static eligibility built by
  `graph.eligibility(excluded_nodes, excluded_edges)`.

### Flow core

`_netgraph_core._place_global_units(graph, view, sources, destinations, elig, cap)` computes
the minimum-cost flow of the maximum value attainable under `cap` (`None`: uncapped) from a
source set to a destination set over the eligible edges, on the capacities `view` (units, id
order), starting from zero flow. It returns `{placed, cost_distribution, edge_flows, min_cut,
potentials, counters}`. It is private: the recorder and the gates use it; the public evaluators
built on it arrive with the policy engine.

## Quick Start

```python
import netgraph_core as ngc
import _netgraph_core

q = 2.0**-12
# Two parallel paths from node 0 to node 3 (physical capacities 10 and 5).
caps = ngc.to_units([10.0, 10.0, 5.0, 5.0], q)
graph = ngc.StrictMultiDiGraph(4, [0, 1, 0, 2], [1, 3, 2, 3], caps, [1, 1, 2, 2], q)

payload = _netgraph_core._place_global_units(
    graph, graph.capacity_units(), [0], [3], None, None
)
print(ngc.as_float(payload["placed"], q))   # 15.0
print(payload["cost_distribution"])         # [[2, 40960], [4, 20480]]
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

## Environment Variables

| Variable | Effect |
| --- | --- |
| `NGRAPH_CORE_PROFILE=1` | Enable the C++ profiling scopes (`profiling.hpp`); the Python accessors return with `netgraph_core.diagnostics`. |

## Requirements

- **C++:** C++20 compiler (GCC 10+, Clang 12+, MSVC 2019+)
- **Python:** 3.11+
- **Build:** CMake 3.23+, scikit-build-core
- **Dependencies:** pybind11, NumPy

## License

[MIT License](LICENSE)
