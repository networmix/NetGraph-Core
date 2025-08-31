NetGraph-Core

High-performance C++ core for NetGraph with Python bindings. CPU-first, GPU-ready.

- Native C++ kernels for shortest paths, k-shortest paths, and max-flow
- Deterministic, thread-parallel, batch execution with GIL released
- Clean Python API compatible with NetGraph’s public surface

Development

- Setup: `make dev` (creates venv, installs dev deps, sets up pre-commit)
- Local checks: `make check` (auto-fix with pre-commit, run C++/Python tests, then lint)
- CI checks: `make check-ci` (strict, non-mutating: lint → C++ tests → Python tests)
- Coverage: `make cov` (prints unified summary and writes XML + single-page HTML under build/coverage/)
