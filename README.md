# NetGraph-Core

C++ implementation of graph algorithms for network flow analysis and traffic engineering with Python bindings.

## Features

- **Algorithms:** Shortest paths (Dijkstra), K-shortest paths (Yen), max-flow (successive shortest paths)
- **Graph representation:** Immutable directed multigraph with CSR adjacency
- **Flow placement:** Tunable policies (proportional to capacity, equal-balanced across paths)
- **Python bindings:** NumPy integration, GIL released during computation
- **Deterministic:** Reproducible edge ordering by (cost, src, dst)

## Installation

```bash
pip install netgraph-core
```

Or from source:

```bash
pip install -e .
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
make cpp-test   # C++ tests only
make py-test    # Python tests only
make cov        # Coverage report (C++ + Python)
```

## Requirements

- **C++:** C++20 compiler (GCC 10+, Clang 12+, MSVC 2019+)
- **Python:** 3.9+
- **Build:** CMake 3.15+, scikit-build-core
- **Dependencies:** pybind11, NumPy

## License

AGPL-3.0-or-later
