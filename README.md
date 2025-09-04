NetGraph-Core

C++ core with Python bindings.

- Native C++ kernels for shortest paths, k-shortest paths, and max-flow
- Deterministic, thread-parallel, batch execution with GIL released
- Python API with a stable, well-documented surface

Development

- Setup: `make dev` (creates venv, installs dev deps, sets up pre-commit)
- Local checks: `make check` (auto-fix with pre-commit, run C++/Python tests, then lint)
- CI checks: `make check-ci` (strict, non-mutating: lint → C++ tests → Python tests)
- Coverage: `make cov` (prints unified summary and writes XML + single-page HTML under build/coverage/)
