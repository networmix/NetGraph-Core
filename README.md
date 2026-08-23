# NetGraph-Core

C++ graph engine for network flow analysis, traffic engineering simulation, and capacity planning.

## Overview

NetGraph-Core provides a specialized graph implementation for networking problems. Key design priorities:

- **Determinism**: Guaranteed reproducible edge ordering by (cost, src, dst).
- **Flow Modeling**: Native support for multi-commodity flow state, residual tracking, and ECMP/WCMP placement.
- **Performance**: Immutable CSR (Compressed Sparse Row) adjacency and zero-copy NumPy views.

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

### 4. Python Integration

- **Zero-Copy**: `FlowState` and `FlowGraph` `*_view()` methods expose C++ buffers as read-only NumPy arrays that track mutation in place. `StrictMultiDiGraph.*_view()` returns copies instead, so the graph's immutability cannot be violated.
- **Concurrency**: Releases the Python GIL during the long-running algorithms (SPF, KSP, max-flow, placement) to enable threading. `PredDAG.resolve_to_paths` holds the GIL. `batch_max_flow` and `sensitivity_analysis` also use internal worker threads, sized by `NGRAPH_CORE_BATCH_THREADS` / `NGRAPH_CORE_SENSITIVITY_THREADS`.

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

## Requirements

- **C++:** C++20 compiler (GCC 10+, Clang 12+, MSVC 2019+)
- **Python:** 3.11+
- **Build:** CMake 3.23+, scikit-build-core
- **Dependencies:** pybind11, NumPy

## License

[MIT License](LICENSE)
