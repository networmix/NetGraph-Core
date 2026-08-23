# NetGraph-Core

C++ graph engine for network flow analysis, traffic engineering simulation, and capacity planning.

## Overview

NetGraph-Core provides a specialized graph implementation for networking problems. Key design priorities:

- **Determinism**: Guaranteed reproducible edge ordering by (cost, src, dst).
- **Flow Modeling**: Native support for multi-commodity flow state, residual tracking, and ECMP/WCMP placement.
- **Performance**: Immutable CSR (Compressed Sparse Row) adjacency, and zero-copy NumPy views over mutable flow state.

## Core Features

### 1. Graph Representations

- **`StrictMultiDiGraph`**: Immutable directed multigraph using CSR adjacency. Supports parallel edges (multi-graph), essential for network topologies.
- **`FlowGraph`**: Topology overlay managing mutable flow state, per-flow edge allocations, and residual capacities.

### 2. Network Algorithms

- **Shortest Paths (SPF)**:
  - Modified Dijkstra returns a **Predecessor DAG** to capture all equal-cost paths.
  - Supports **ECMP** (Equal-Cost Multi-Path) routing.
  - Features **node/edge masking** and **residual-aware tie-breaking**.

- **K-Shortest Paths (KSP)**:
  - Yen's algorithm returning DAG-wrapped paths.
  - Configurable constraints on cost factors (e.g., paths within 1.5x of optimal).

- **Explicit Paths**:
  - `PredDAG.from_edges(graph, edges)` builds a path bundle from an operator-supplied
    edge sequence, usable anywhere a `PredDAG` is accepted.

- **Max-Flow**:
  - **Algorithm**: Iterative augmentation using Successive Shortest Path on residual graphs, pushing flow across full ECMP/WCMP DAGs at each step. For `Proportional` placement with `require_capacity=True`, a final residual completion phase augments with reverse arcs so the result is a true maximum flow; the other modes are placement models and may report less.
  - **Traffic Engineering (TE) Mode**: Routing adapts to residual capacity (progressive fill).
  - **IP Routing Mode**: Cost-only routing (ECMP/WCMP) ignoring capacity constraints.

- **Analysis**:
  - **Sensitivity Analysis**: Identifies critical edges by removing each saturated edge and measuring how much total flow is *lost*. Supports `shortest_path` mode to analyze only edges used under ECMP routing (IP/IGP networks) vs. full max-flow (SDN/TE networks).
  - **Min-Cut**: Computes minimum cuts on residual graphs.

### 3. Flow Policy Engine

Unified configuration object (`FlowPolicy`) that models diverse routing behaviors:

- **Modeling**: Unified configuration for **IP Routing** (static costs) and **Traffic Engineering** (dynamic residuals).
- **Placement Strategies**:
  - `EqualBalanced`: **ECMP** (equal splitting) - equal distribution across next-hops and parallel edges.
  - `Proportional`: **WCMP** (weighted splitting) - distribution proportional to residual capacity.
- **Lifecycle Management**: Handles demand placement, re-optimization of existing flows, and constraints (path cost, stretch factor, flow counts).
- **Static (Pinned) Paths**: `FlowPolicy.set_static_paths()` pins a demand to explicit path bundles (MPLS-style). One flow per usable bundle; a bundle with no path surviving the failure masks is *down* and carries nothing, since a pinned path does not reroute.

### 4. Python Integration

- **Zero-Copy**: `FlowState` and `FlowGraph` `*_view()` methods expose C++ buffers as read-only NumPy arrays that track mutation in place. `StrictMultiDiGraph.*_view()` returns copies instead, so the graph's immutability cannot be violated.
- **Concurrency**: Releases the Python GIL during the long-running algorithms (SPF, KSP, max-flow, placement) to enable threading. `PredDAG.resolve_to_paths` holds the GIL. `batch_max_flow` and `sensitivity_analysis` also use internal worker threads, sized by `NGRAPH_CORE_BATCH_THREADS` / `NGRAPH_CORE_SENSITIVITY_THREADS`.

## Quick Start

```python
import numpy as np
import netgraph_core as ngc

# Two parallel paths from node 0 to node 3.
graph = ngc.StrictMultiDiGraph.from_arrays(
    num_nodes=4,
    src=np.array([0, 1, 0, 2], dtype=np.int32),
    dst=np.array([1, 3, 2, 3], dtype=np.int32),
    capacity=np.array([10.0, 10.0, 5.0, 5.0], dtype=np.float64),
    cost=np.array([1, 1, 1, 1], dtype=np.int64),
    ext_edge_ids=np.arange(4, dtype=np.int64),  # your own stable edge ids
)

algs = ngc.Algorithms(ngc.Backend.cpu())
handle = algs.build_graph(graph)

total, summary = algs.max_flow(handle, 0, 3)
print(total)                            # 15.0
print(summary.min_cut.edges.tolist())   # [0, 1] -- bottleneck edge ids
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

## Environment Variables

| Variable | Effect |
| --- | --- |
| `NGRAPH_CORE_PROFILE=1` | Enable profiling of C++ hot paths (`profiling_dump()` / `profiling_reset()`). |
| `NGRAPH_CORE_BATCH_THREADS` | Worker threads for `batch_max_flow` (default: hardware concurrency). Set to `1` when calling from your own worker pool. |
| `NGRAPH_CORE_SENSITIVITY_THREADS` | Worker threads for `sensitivity_analysis` (default: hardware concurrency). |

## Requirements

- **C++:** C++20 compiler (GCC 10+, Clang 12+, MSVC 2019+)
- **Python:** 3.11+
- **Build:** CMake 3.23+, scikit-build-core
- **Dependencies:** pybind11, NumPy

## License

[MIT License](LICENSE)
