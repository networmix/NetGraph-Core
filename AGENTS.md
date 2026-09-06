# Development

- Setup: `bash .superset/workspace.sh setup`. Each worktree gets its own venv;
  setup does not replace the Git hooks shared by worktrees.
- Superset **Run**: `bash .superset/workspace.sh check` rebuilds the Python
  extension, then runs lint, types, C++ and Python tests. Editable Python imports
  alone do not rebuild changed C++.
- For C++ changes, also run `make sanitize-test` (ASan/UBSan; errors must fail).
- For downstream changes, run NetGraph's `dev/check_core_integration.sh` with
  this worktree's absolute path. Test the chosen build, not the published wheel.

Confirm defects with a reproducer. Max-flow tests need feasibility, conservation
and optimality/cut checks with nonuniform costs, masks and cost-factor cases.
Keep bindings and `_docs.py` consistent. Performance claims need correctness
checks and A/B/A on one quiet machine. Local checks do not replace CI.
