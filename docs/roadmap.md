# Roadmap

This page is for maintainers. It is the product's rolling roadmap and the only place that records planned and
deferred work and its state. The user documentation describes what the product does today; nothing on this page is
available until its row says so.

## Status values

**Contract**

- `pending`: the item needs a reviewed design before work starts; the link names the current source.
- `ready`: the linked design holds a reviewed contract for the item.
- `not needed`: the Work cell states the change completely.
- `—`: delivered before this roadmap; the code and user documentation define the behaviour.

**Implementation**

- `planned`: none of the item is in the product.
- `partly implemented`: part of it is; the evidence names which part.
- `implemented`: the item is in the product.

## Maintaining the roadmap

- One row per planned item. IDs never change and are never reused; a split gets new IDs and the old one leaves the
  table. A need in another repository is written `<repository>:<ID>`.
- A row changes in the commit that changes its item: its scope, its needs, its contract or its implementation. Its
  evidence names in-repository witnesses (tests, documents, reports); a row never names the commit that changes it.
- An implemented row may leave the table once a release has shipped it and no row needs it; the repository's changelog
  or its dated history reports keep it.
- Status words live only on this page. Designs state requirements; the changelog or the dated history reports say what
  shipped.
  When this page and a design disagree on order or needs, this page wins; on behaviour, the design wins.
- The page records product facts, not how the work is run: an effort's plans, reviews and CI results stay in its own
  records.
- An item moves from Deferred to Planned on the maintainer's decision, and gets its ID then.

## Planned

Nothing is scheduled yet.

| ID | Work | Needs | Contract | Implementation | Evidence |
|---|---|---|---|---|---|

## Deferred

Not scheduled. One entry per item: what it is, and why it is not scheduled or what would make it worth doing.

- **Unexplained macOS CI cancellation**: investigate macOS unit or ThreadSanitizer jobs cancelled without useful logs,
  retaining runner diagnostics and distinguishing cancellation from a test failure; keep this as a watch item until it
  recurs with evidence that identifies a cause.
- **Daemon finalization crash**: retain the faulting native thread's stack and interpreter crash report if
  `test_daemon_exit_releases_native_resources` crashes again during interpreter shutdown; the isolated macOS Python 3.12
  failure has not been reproduced locally, so a product fix needs a usable reproducer.
- **libc++ hot-loop growth**: audit remaining per-element vector growth in hot native loops for stack-offset-dependent
  page-alias stalls, using controlled layout sweeps before changing the loops; the solver input copy and route-plan fill
  already use indexed storage, and other loops need measured exposure before a rewrite.
- **Proof checks only when needed**: move the bridge's model-derived expectations into ordinary apportionment tests and
  build/replay Lean proofs only when proof inputs change, accounting for release and distribution callers; the record,
  anchor and digest layer is already removed, and further tooling needs a concrete verification benefit.
- **Offline Core setup**: let `.superset/workspace.sh setup` prepare a wheelhouse at the CI constraints' versions and
  seed the native tests' pybind11 and googletest inputs so a worktree prepared once online can run its gates offline;
  this is useful when preparing confined workers, and needs no new downloader service.
- **GIL-build skip-guard control**: add a focused self-test that a GIL-enabled build refuses a record containing a
  `GIL-proof-only` skip, while a free-threaded build accepts the intended skips; the guard is correct today, so add the
  missing separating case when maintaining `dev/test_pytest_guard.py`.
- **One-allocation leak control**: add a negative control in `tests/py/test_native_allocation.py` that leaves exactly
  one tracked allocation live, so weakening the zero-leak bound to one fails; the present check is correct and the gap
  matters when that test is next maintained.
- **Small native fixed costs**: profile small SPF, residual SPF, path-resolution or capacity-weighted placement calls
  only when a real application profile shows their fixed costs affecting user time; old microsecond comparisons predate
  the current implementation and do not justify an optimization by themselves.
- **Free-threaded allocation fault flake**: `test_native_allocation_failure_raises_and_changes_nothing[spf.full]` in
  `tests/py/test_native_allocation.py` sometimes fails on Linux with free-threaded Python 3.14 (ten native allocations
  still live after the failed call) and passes when the same commit is run again; find what else allocates while the
  fault hook counts before reading a red run of this test as a product defect, and make the test deterministic.
