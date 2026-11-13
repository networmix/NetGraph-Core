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
- **Concurrency model controls**: `include/netgraph/core/detail/concurrency_models.hpp`, exposed through
  `bindings/python/concurrency_hooks.cpp`, runs finite-state models of the publication and phase protocols, and the
  tests that use it (`tests/cpp/contract_examples_tests.cpp` and the Python hook tests) check those models, never the
  library's own locks and atomics; `tests/cpp/concurrency_tests.cpp`, `tests/cpp/graph_concurrency_tests.cpp` and the
  ThreadSanitizer job exercise the real code. Remove the models, their bindings and their examples, unless one is shown
  to catch a defect of the real synchronization; a check that cannot fail on the product guards nothing.
- **C++ tests once per system**: every unit-test job (three systems, four Python versions and free-threaded 3.14) runs
  `make check-ci`, which rebuilds and reruns the whole C++ suite with and without the private hooks, although that suite
  does not depend on the Python version. Running it once per system, and the Python suite per version, would cut the
  runner time of a CI run; worth doing when CI time or runner queues slow the work.
- **Coverage on demand**: the coverage job (`make cov`, about 35 minutes) is the longest job of every CI run, about
  twice the next one, so it sets when a run finishes. Running it on demand or on a schedule would end a run in about 15
  minutes; keep it on every push only if a coverage drop must block a merge.
- **Ordered publication reuses its totals**: when an ordered `apply_actions` list commits, hand the carried and reserved
  totals the loop already computed to publication and swap them into the ledger after every fallible step, removing two
  passes over the touched slots (three files, +43/−14 lines; prototype on branch `prototype/ordered-publication`, proof
  and conditions in [the design page](design/routing-performance/ordered-publication.md)). On the DC-BB workload the
  whole run was 3.8 % faster at one worker and showed no difference beyond noise at its default eight, with every output
  unchanged; it waits for a landing build that meets its go/no-go rows (DC-BB faster at one worker and not slower at
  eight; `kernel.lossy.clos`, `kernel.eb.clos`, `kernel.weighted.clos`, `ledger.joint.pooled_cycle`,
  `program.stale.lossy` and 14 drift-inconclusive harness rows not slower).
- **Baseline route-program reuse (opt-in)**: a capture of a workflow step's baseline policies and a `preparation=`
  keyword, so that failure cases whose roster sources all keep their distance reuse the baseline's route programs,
  restricted and pruned to the case, instead of a new search and program build (prototype on branch
  `prototype/baseline-reuse`, its NetGraph adapter on NetGraph's `prototype/baseline-reuse-adapter`; proofs and contract
  choices in [the design page](design/routing-performance/baseline-reuse.md)). It adds public API, so it waits for a
  contract decision (adoption, capture lifetime and owner, four smaller choices); then one qualified comparison on the
  landing tree decides it: the chosen lifetime faster beyond noise on DC-BB at eight workers and not slower at one.
  Measured: 18 % faster at one worker with a per-step capture, 4.6 % at eight with a per-run capture, no difference on
  NetLab's shipped DC-BB experiment, whose cases were never counted for whether they used the reuse;
  `concurrency.delivered` was 15.0 % slower with the code present but unused, unexplained.
- **Pop-order sentence in `search.hpp`**: `include/netgraph/core/search.hpp` says both search queues "pop in
  (distance, node) order and give identical results". No public result depends on that order
  ([proof](design/routing-performance/equivalences.md#pop-order-fifo-queues-bit-lanes-and-max-flow-pair-symmetry)),
  but any patch that changes it (removing the per-bucket sort, a unit-cost FIFO queue, or a search that pops by distance
  plus a potential: A* for targeted `spf`, KSP spurs or tunnel setup) must first rewrite the sentence, for example to
  "the order among equal distances is unspecified". Nothing to do until such a patch shows its own gain; the sort
  removal alone gave none.
- **Max-flow batch store**: keep `evaluate_pairs` batch results within one NetGraph run spanning its three capacity
  analyses, since NSFNET's 39 repeated batches recur across those analyses and never inside one. Priced at 13.41 % raw /
  15.87 % corrected of NSFNET's one-worker Core time against a 13.04 % bar, with retention, copies and allocator effects
  unmeasured; worth a prototype only if the store's lifetime cost is charged and it is timed on NetGraph's command line
  at 1 and 8 workers.
- **Max-flow pair reduction (opt-in)**: solve (s, t) once and derive (t, s) on graphs whose reverse capacities, costs
  and masks equal the forward ones, with the reverse cut from one residual closure when the solve reached the maximum
  with the cut requested, `None` when the offer stopped it below or the cut was not requested, and `[]` for an empty
  roster ([proof](design/routing-performance/equivalences.md#the-reverse-pair)). Exact only for the eleven result fields
  NetGraph's export reads, so it needs an opt-in evaluator or output mode; its partial share on NSFNET at one worker is
  14.28 % raw / 16.88 % corrected of Core time, the rest unpriced.
- **Sink stop in controller solves**: stop each tier's residual Dijkstra once the sink is settled and raise unsettled
  potentials by the sink's distance; every public result of every solve from zero is unchanged
  ([proof](design/routing-performance/equivalences.md#sink-stop-in-the-flow-solver)). On `backbone_clos` at one worker a
  free search plus potential updates is 20.8-21.0 % raw / 21.5-21.7 % corrected of Core time, above the bar, but the
  changed search path is unpriced; on NSFNET it is under the bar. Deciding: a result-producing solver timed on
  `backbone_clos` at 1 and 28 workers.
- **A\* spur searches for `ksp`**: one reverse search per call as a potential, strict cost bounds, and the
  canonical-suffix shortcut only under `FEWER_HOPS` with its validity check; exact, but it must remove about 31 % (raw)
  of the KSP phase at one worker and no primary workload calls `ksp`. Worth building when a real workload makes many
  `ksp` calls.
- **Static potential for tunnel setup**: one reverse-distance potential per T1 setup plus a dead-node memory per request
  value; exact, but it must remove 34.6 % raw / 33.5 % corrected of the setup phase at one worker and no primary
  workload places T1 tunnels. Worth building when one does.
- **One native call per failure batch**: let Core run NetGraph's per-case loop for a batch of failure masks, so less
  Python work runs under the interpreter lock at eight workers. On DC-BB the lock-held caller work plus the attached
  binding work inside Core spans is 0.447 s, 3.90 % of a diagnostic run (3.59 % on the sensitivity basis), above the
  3 % bar only if a batch removes at least 45 % (64 %) of the attached work; the batch's own conversion of masks,
  demands and results is unpriced. Deciding: a Core entry point with a NetGraph adapter timed on DC-BB at 8 workers.
- **Other unknown performance leads**: fewer passes over support and programs; LOSSLESS admission with fewer sweeps;
  reduced-result certificates for failure pairs; persistent worker pools; a columnar output form; unit-cost and
  wide-cost search queues; tunnel-count sizing by cut steps; exact divisor selection with integer comparisons. Each is
  exact or unproved as stated, priced or unpriced, and has a named deciding experiment in
  [the candidates page](design/routing-performance/candidates.md); none is worth scheduling without a workload that
  spends its time there.
