# How to measure here

What the next person who times a Core change on these workloads needs to know: how the reference was built, how a
comparison is run and read, which rows move by themselves, and what each instrument measured and missed. The general
A/B/A tool and its statistics are in [`dev/perf/README.md`](../../../dev/perf/README.md); workloads and terms are in
the [README](README.md).

## The reference build

- Core `9cd7642` (code of release 0.13.0) built as a wheel from a clean export (`git archive`), offline, without build
  isolation, in one pinned build environment: scikit-build-core 1.1.2, pybind11 3.1.0, CMake 4.4.4, Ninja 1.13.2;
  CPython 3.14; macOS target 15.0. Its extension's SHA-256 is
  `208487c7aed158fb418791c22e23d09332533595053c0e4810923b31244f7c84`; rebuilt in the same environment it comes out
  byte-identical. The recipe, with `<env>` the pinned environment's venv, `<export>` a `git archive` export of the
  revision and `<wheels>` the output directory:

  ```sh
  python3.14 -m venv <env>
  <env>/bin/python -m pip install cmake==4.4.4 ninja==1.13.2 numpy==2.5.4 packaging==26.3 pathspec==1.1.1 \
      pybind11==3.1.0 scikit_build_core==1.1.2
  export CMAKE_BUILD_PARALLEL_LEVEL=4 CC="$(xcrun --find clang)" CXX="$(xcrun --find clang++)" \
      MACOSX_DEPLOYMENT_TARGET=15.0
  <env>/bin/python -m pip wheel <export> --no-deps --no-build-isolation -w <wheels>
  ```

- Every candidate is a wheel built the same way from its own clean export; builds are compared by their extension's
  hash. The repository's editable development build carries the private test hooks, so its extension differs and is
  never a reference; neither is a PyPI wheel.
- NetGraph `8a641a1` (0.27.0) from a clean export imports the wheel in its own venv; NetLab `32d5e13` (0.7.0).

## The inputs

- DC-BB: NetGraph's `docs/design/perf-0.24-vs-0.26/` scenario file, SHA-256
  `cd81a0f81d8227ffc2d688ca24df287410f991500023ecde83fbca301fa7a8c3`.
- DC-BB 192 pods: as shipped in NetGraph (`docs/design/perf-0.24-vs-0.26/inputs/dcbb192-0.26.0.yml.gz`).
- NSFNET and `backbone_clos`: as bundled in NetGraph `8a641a1`.
- The synthetic backbone comes from a fixed generator command that reproduces it byte for byte; like the other
  instruments, the generator (`gen_backbone.py`) is on the branch `research/routing-performance-instruments`. Its
  parameters: `gen_backbone.py OUT
  --draws 50`, every other argument at its default: 25 metros, `--k 2`, `--partners 5`, two Clos L16 S4 PoPs each
  (1025 nodes, 12,512 edges, 125 TE_MAXFLOW_UNLIM demands); for size, the same with `--metros 10` (410 nodes, 4992
  edges, 50 demands) and `--metros 50` (2050 nodes, 25,024 edges, 250 demands).
- Its tunnel variant replaces every demand's `preset: TE_MAXFLOW_UNLIM` with
  `model: {routing: {type: distributed_te, tunnel_count: {per_head: 16}}}`, every other key at its default (lossless, no
  cost ceiling, no minimum unreserved bandwidth, `fewer_hops`, parallel edges), with the same workflow.
- Its path queries: NetGraph's `paths` analysis (`ngraph.analysis.paths`, mode `pairwise`) on its graph for 200 ordered
  pairs drawn by `random.Random(20261108).sample` from the sorted list of the 600 ordered pairs of distinct `metroN/dc1`
  nodes, each pair queried once with `SpfOptions(Direction.FORWARD)` and once with `KspOptions(8)`, other options at
  their defaults.

## A comparison

- **Protocol.** Seven fresh-process reference/candidate/reference triples in a quiet window; the median ratio
  candidate / reference and its 95 % bootstrap interval. Rows are timed from process start to exit. Before any
  candidate, the reference is compared with itself (a reference A/A) at each row's own scale, which qualifies the window
  for that row.
- **Verdict words.** *Faster beyond noise*: the interval's upper end is below 0.97. *Slower beyond noise*: its lower
  end is above 1.03. *Not slower*: its upper end is at most 1.03. An interval that crosses 1.03 is *inconclusive at
  +3 %*: a slowdown is neither established nor excluded. A row whose two reference runs drift apart by more than 3 % is
  *inconclusive*. *Two clusters*: a row's timings fall into two separated groups, which is explained before it is
  averaged.
- **Memory.** Peak resident size of each fresh process on both sides; a row fails its memory check above
  max(1.10 × reference, reference + 8 MiB). A peak is an observation, not a promise of zero growth.
- **Outputs.** Whole command-line exports are compared on values, types, keys and list order, excluding only the field
  that names Core's binary (and, for NetLab, named execution-time fields); the comparator first shows it rejects a
  poisoned export.
- **Separate gains do not add.** A gain measured on one tree does not carry to another tree or add to another change's
  gain; a combination is judged only on its own measurement.

## Rows that move by themselves

- **Code placement.** A same-source rebuild in another link order (a *layout control*) moves some harness rows with no
  code change: code placement alone moves rows by up to about 10 %, and `kernel.lossy.clos`, `kernel.eb.clos` and
  `kernel.weighted.clos` move with code placement by 3-5 %. A gain between 0.90 and 1.00 counts only if a same-source
  layout control beside it does not move the row as much; the same holds for a slowdown before it is blamed on the code.
- **Rows to watch** on any change to the ledger or policy code: those three; `ledger.joint.pooled_cycle` and
  `program.stale.lossy` (clustered, drifting or inconclusive in the proposals' windows); `concurrency.delivered`, 15.0 %
  slower with the opt-in reuse code present but unused while its layout control was exactly 1.0000, a slowdown the
  source does not explain.
- The concurrency rows had no qualified noise floor in these windows; their reference A/A drifted.

## Harness rows a landing reruns

Both proposals name a set of Core's own harness rows (`dev/perf/`; several names carry a suffix the harness builds,
such as `.scaled` or `.one`) that a landing compares against the same tree without the change.

- **Ordered publication, 102 rows**, chosen by call chain: the rows whose timed call reaches the changed publication
  code, plus scope controls. The working-totals caller: `batch.ordered.clos`; the 85 rows on the default commit path:
  `kernel.proportional.clos`, `kernel.eb.clos`, `kernel.lossy.clos`, `kernel.proportional.grid`, `kernel.eb.grid`,
  `kernel.weighted.clos`, `capacity.costonly_eb.failures`, `capacity.costonly_prop.failures`, `maxflow.exact.totals`,
  `maxflow.exact.detailed`, `sens.exact`, `sens.exact.w2`, `sens.exact.cut`, `sens.exact.partial`, `sens.exact.branch`,
  `sens.exact.reroute`, `sens.exact.capped`, `sens.exact.cut.w2`, `sens.exact.partial.w2`, `sens.exact.branch.w2`,
  `sens.exact.reroute.w2`, `sens.exact.capped.w2`, `sens.policy_all_edges`, `sens.policy_all_edges.w2`,
  `pooled.ecmp.clos`, `pooled.lossy.clos`, `pooled.from_dag.clos`, `batch.scalar.clos`, `batch.simultaneous.clos`,
  `program.ring.lossy`, `batch.simultaneous.ring`, `global.equal.random`, `sens.equal_global`, `global.weighted.random`,
  `sens.weighted_global`, `program.increment.optimized`, `program.increment.weighted`,
  `program.increment.weighted.wide`, `program.replay.weighted`, `program.forward.global`, `program.replay.lossy`,
  `program.stale.lossy`, `program.replace.costonly`, `pooled.anycast.many`, `pooled.anycast.one`, `pooled.anycast.two`,
  `pooled.anycast.four`, `pooled.single_path.one`, `pooled.choice.collapse_available`, `pooled.choice.tree`,
  `pooled.choice.bottleneck_available`, `pooled.weight.explicit`, `pooled.weight.snapshot`,
  `pooled.weight.explicit.one`, `pooled.weight.snapshot.one`, `pooled.from_dag.weighted`, `pooled.controller.equal`,
  `global.support.equal`, `program.forward.pooled`, `program.repair.pooled`, `program.stale.pooled`,
  `ledger.aggregate_default.many`, `ledger.aggregate_default.one`, `ledger.aggregate_default.edges`,
  `ledger.aggregate_default.lossy`, `ledger.aggregate_default.lossy_one`, `ledger.aggregate_default.lossy_repeat`,
  `result.options.none`, `result.options.by_source`, `result.options.program`, `result.options.both`,
  `pooled.lossy.repeat`, `pooled.lossy.repeat.one`, `program.stale_then_place.pooled`,
  `program.stale_then_place.pooled.one`, `ledger.lifecycle.one`, `ledger.joint.pooled_cycle`, `t1.setup.fabric`,
  `t1.forward.fabric`, `t2.exact16.headroom`, `t2.exact16.cost`, `t2.upto256.headroom`, `t2.size.mesh`, `msd.loop`,
  `msd.loop.epochs`; the 16 scope controls: `pooled.from_dag_build.clos`, `pooled.from_dag_build.one`,
  `ledger.record_read.many`, `ledger.record_read.one`, `ledger.record_read.orphan`, `ledger.record_read.repair`,
  `ledger.copy.many`, `ledger.copy.one`, `concurrency.place_zero`, `concurrency.delivered`, `concurrency.program`,
  `concurrency.config`, `concurrency.two_ledgers_serial`, `concurrency.two_ledgers_parallel`,
  `concurrency.two_ledgers_serial.scaled`, `concurrency.two_ledgers_parallel.scaled`. Add `batch.ordered.one_owner`, a
  one-owner ordered list.

- **Baseline reuse, 87 rows**, chosen as the rows that run the code the opt-in prototype shares with the ordinary path,
  each smoked once on both sides with an equal correctness digest: `kernel.proportional.clos`, `kernel.eb.clos`,
  `kernel.lossy.clos`, `kernel.proportional.grid`, `kernel.eb.grid`, `kernel.weighted.clos`,
  `capacity.costonly_eb.failures`, `capacity.costonly_prop.failures`, `maxflow.exact.totals`,
  `maxflow.exact.detailed`, `sens.policy_all_edges`, `sens.policy_all_edges.w2`, `pooled.ecmp.clos`,
  `pooled.lossy.clos`, `pooled.from_dag.clos`, `pooled.from_dag_build.clos`, `batch.scalar.clos`,
  `batch.ordered.clos`, `batch.simultaneous.clos`, `program.ring.lossy`, `batch.simultaneous.ring`,
  `global.equal.random`, `global.weighted.random`, `program.increment.optimized`, `program.increment.weighted`,
  `program.increment.weighted.wide`, `program.replay.weighted`, `program.forward.global`, `program.replay.lossy`,
  `program.stale.lossy`, `program.replace.costonly`, `pooled.anycast.many`, `pooled.anycast.one`,
  `pooled.anycast.two`, `pooled.anycast.four`, `pooled.single_path.one`, `pooled.choice.collapse_available`,
  `pooled.choice.tree`, `pooled.choice.bottleneck_available`, `pooled.weight.explicit`, `pooled.weight.snapshot`,
  `pooled.weight.explicit.one`, `pooled.weight.snapshot.one`, `pooled.from_dag.weighted`, `pooled.controller.equal`,
  `global.support.equal`, `pooled.from_dag_build.one`, `program.forward.pooled`, `program.repair.pooled`,
  `program.stale.pooled`, `ledger.aggregate_default.many`, `ledger.aggregate_default.one`,
  `ledger.aggregate_default.edges`, `ledger.aggregate_default.lossy`, `ledger.aggregate_default.lossy_one`,
  `ledger.aggregate_default.lossy_repeat`, `ledger.record_read.many`, `ledger.record_read.one`,
  `ledger.record_read.orphan`, `ledger.record_read.repair`, `pooled.lossy.repeat`, `pooled.lossy.repeat.one`,
  `program.stale_then_place.pooled`, `program.stale_then_place.pooled.one`, `ledger.copy.many`, `ledger.copy.one`,
  `ledger.lifecycle.one`, `ledger.joint.pooled_cycle`, `t1.setup.fabric`, `t1.forward.fabric`, `t2.exact16.headroom`,
  `t2.exact16.cost`, `t2.upto256.headroom`, `t2.size.mesh`, `msd.loop`, `msd.loop.epochs`, `cache.cold`, `cache.warm`,
  `concurrency.place_zero`, `concurrency.delivered`, `concurrency.program`, `concurrency.config`,
  `cache.frozen.reads`, `concurrency.two_ledgers_serial`, `concurrency.two_ledgers_parallel`,
  `concurrency.two_ledgers_serial.scaled`, `concurrency.two_ledgers_parallel.scaled`.

## What each instrument measured, and its limits

- **The replay of recorded Core calls** records a workload's Core calls with full inputs (graph, mask contents,
  policies, actions, removals, getters, object lifetimes) and replays them against a chosen build. In check mode it
  compares every returned field, the ledger's state and records and any cache's statistics with the base build's; in
  timed mode it times only the Core calls. Its *Core's work* is the union of Core-call spans on the wall axis (overlaps
  counted once). Its agreement with NetGraph's own per-call timers decides what it may time: on DC-BB at one worker 7 of
  65 call kinds qualified, among them `apply_actions` (2.709 ms per call on both) and the policy constructor, and the
  qualified kinds are 97.35-97.50 % of the union; at eight workers 4 of 65 qualified. On NSFNET the replay's per-call
  clock is about 46 % of its union, so NSFNET is decided on NetGraph's command line, not on the replay; on
  `backbone_clos` no call kind qualifies. It does not read internal state in unsampled cases, and its memory peaks do
  not isolate Core's own memory.
- **Per-call clocks** (wall and CPU timers around every Core call, destructor and getter in NetGraph) cost wall time
  themselves: 1.031× at one worker and 1.075× at eight on DC-BB. Their one-thread timer correction underestimates the
  wrapper's cost by at least 0.37 s at eight workers, so a caller-CPU figure "outside Core" depends on the correction
  used (1.68 s against 1.825-1.919 s on DC-BB).
- **The lock split** divides each Core span into the part that holds the interpreter lock (argument and result
  conversion, validation, binding and runtime work, release edges), the released interval and the retake interval. The
  retake wall brackets the lock's reacquisition and is only an upper proxy for lock waiting in an instrumented run: on
  `backbone_clos` at 28 workers it sums to 55.607 s across callers from an instrumented run of 5.742 s, 36 times the
  plain run's 0.160 s, and says nothing about the plain product. Lock waiting in the plain product was not measured.
- **Sampled stacks** give stack shares, not CPU shares. Different captures have categories with different boundaries;
  their shares are shown side by side and never merged or averaged.
- **Counts** come from instrumented copies and are given raw and corrected for the instrumentation's own cost; each
  share names its own denominator, since three instruments give three different "Core times" for one workload.
