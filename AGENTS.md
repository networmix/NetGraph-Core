# Development

- Setup: `bash .superset/workspace.sh setup`. Each worktree gets its own venv;
  setup does not replace the Git hooks shared by worktrees.
- Superset **Run**: `bash .superset/workspace.sh check` rebuilds the Python
  extension, then runs lint, types, the C++ tests with and without the private
  hooks, and the Python tests. Editable Python imports alone do not rebuild
  changed C++.
- For C++ changes, also run `make sanitize-test` (ASan/UBSan) and
  `make tsan-test` (ThreadSanitizer); errors must fail.
- For downstream changes, run NetGraph's `dev/check_core_integration.sh` with
  this worktree's absolute path. Test the chosen build, not the published wheel.

Confirm defects with a reproducer. Max-flow tests need feasibility, conservation
and optimality/cut checks with nonuniform costs, excluded nodes and edges, and
cost ceilings (Absolute, Stretch).
A binding signature or docstring change needs the regenerated type stub
(`bash dev/stubs/run.sh generate`; `make stubs-check` checks drift, package
contents and typed clients). Performance claims need correctness checks and
A/B/A on one quiet machine (`dev/perf/README.md`; `make perf-tool-test` after
changing that tool; CI runs it on one macOS job, as a library change can break
it). Local checks do not replace CI.
