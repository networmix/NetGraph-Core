# Ordered publication reuses its totals

A local edit, ready to land under the go/no-go conditions below. The prototype is the branch
`prototype/ordered-publication` of this repository. Workloads, terms and the evidence labels are defined in the
[README](README.md).

## The change

When an ordered list of placements commits, Core already holds the working totals it computed while placing (carried
and reserved amounts per edge). Today publication rebuilds them in two further passes over the touched slots; the edit
hands them to publication and swaps them into the ledger only after every step that can fail has passed. Scalar and
simultaneous commits keep their old path (their machine code is not byte-identical: default null pointer arguments are
added, without a new default-path loop). No public API, option or retained state is added. Bindings and public
declarations are byte-identical to the base.

Size: three files, +43/−14 lines: `include/netgraph/core/detail/ledger_access.hpp` (+3/−1), `src/flow_ledger.cpp`
(+36/−10), `src/flow_policy.cpp` (+4/−3).

## What stays equal, and why

Every documented output of the calls it touches stays identical: the ledger's projections after the list
(`carried()`, `residual()`, `reserved()`, `unreserved()`) and after each action (its `residual` and
`saturated_boundary_edges` outputs); the placement results (delivered amounts, per-edge split, change flags); the
ledger's records, their ids and order, epochs, holds and retained snapshots; sharing of objects; every error with its
text and order (including which touched element fails first); and the failed-call guarantee (a raising call leaves no
partial state and can be retried).

Why, in one line: after the loop the working totals equal the ledger's totals plus, for every action that publishes,
its new entry's per-edge amounts minus its old entry's; each policy appears once in a list, so this is exactly what the
old path added edge by edge. The proof, with the widths, the error order and the atomicity argument, is in
[equivalences, "The swap of totals"](equivalences.md#the-swap-of-totals).

- **Widths**: no new sum or product; every carried total stays in [0, B_e] and every reserved total in [0, R0_e]
  (B_e and R0_e: the ledger's service and reservation baselines of edge e), both at most the edge's capacity < 2^62,
  in int64.
- **Work**: two passes over the touched slots go, one check pass stays; no loop is added, and every loop is over a
  finite list.
- **Memory**: the ordered path allocates nothing new (the per-edge change scratch is no longer built); the ledger takes
  the loop's vectors of the same length.

Evidence of equality: 231,913 fixed call histories give results and state equal to the base (counted); 16,000 random
batches through the new and old paths agree; whole command-line exports of DC-BB, NSFNET and `backbone_clos` are
byte-equal to the reference except the field that names Core's binary (`metadata.core_extension_sha256`). One planted
fault, default holds no longer subtracted, was missed by the study's oracle and caught by Core's own
`tests/py/test_distributed_te.py::test_twoh_call`: no oracle replaces the repository suite.

## Gain

Reference and candidate built the same way from the same base; NetGraph's command line timed from start to exit.

| Row | 1 worker: candidate/reference [95 % interval] | 8 workers (DC-BB's default): the same | Evidence |
| --- | --- | --- | --- |
| DC-BB whole run | 0.9624 [0.9607, 0.9662]: 3.8 % faster beyond noise, 17.41 → 16.77 s; layout control 0.9996 | 0.9884 [0.9863, 0.9951]: not slower, no gain shown, 9.86 → 9.78 s | qualified (seconds derived) |
| Core's work on DC-BB | 0.9325 [0.9308, 0.9352]: 6.7 % faster, 9.09 → 8.47 s | 0.9408 [0.9388, 0.9415]: 5.9 % faster, 1.69 → 1.59 s | qualified (seconds derived) |
| NSFNET (no ordered lists) | not slower: 0.9966 [0.9821, 1.0054] | not slower: 1.0036 [0.9983, 1.0071] | qualified |
| `backbone_clos` (599 ordered lists) | not slower: 1.0008 [0.9975, 1.0015] | not slower at 28: 0.9980 [0.9974, 1.0005] | qualified |

Core's work is where the gain comes from, not a user gain. Start-up share, so that a faster start is not read as faster
work: on DC-BB at one worker the import takes 0.04 s and an empty NetGraph run (import, load and model) 7.04 s of the
17.41 s; both controls are within noise. A cost model written before the timing expected 1-3 %; the measured 3.8 % is
the result. The gain grows with the share of time spent in ordered publication; on `backbone_clos`, where solving
dominates a 0.4 s run, it is not visible.

Peak memory flags no risk (reference / candidate: DC-BB at 1 worker 450.6 / 449.0 MB, at 8 workers 687.2 / 627.8 MB;
`backbone_clos` at 28 workers 94.4 / 94.4 MB). These are observed peaks, not a promise of zero growth.

Not measured: NetLab's `netlab run` at its 7 workers; concurrent NetLab jobs (`--run-jobs`); other interpreters.

## Risks

- `kernel.lossy.clos` (a harness row on the unchanged commit path) ran 5.07 % slower beyond noise, but the same-source
  reordered reference alone moves it by 3.1 % with two clusters, so the slowdown is not separable from code placement;
  it is not shown harmless either.
- `program.stale.lossy`: seven alternating pairs give 1.0374 [1.0065, 1.0461], inconclusive at +3 %: a slowdown of up
  to 4.6 % is neither established nor excluded.
- `kernel.eb.clos`, `kernel.weighted.clos` and `ledger.joint.pooled_cycle` split into two clusters, and their layout
  controls drifted, so the clusters are unexplained. The one-owner ordered-list row `batch.ordered.one_owner`, 0.9737
  [0.9673, 0.9966]: gain not established.
- Fourteen harness rows were inconclusive by drift and are unmeasured: `ledger.lifecycle.one`,
  `pooled.from_dag_build.clos`, `pooled.lossy.repeat.one`, `program.increment.optimized`, `program.replay.weighted`,
  `program.stale_then_place.pooled`, `sens.exact.branch`, `sens.exact.capped`, `sens.exact.reroute`,
  `concurrency.two_ledgers_parallel.scaled`, `program.increment.weighted`, `program.replay.lossy`,
  `program.ring.lossy`, `sens.exact.capped.w2`. The concurrency rows had no qualified noise floor.

## Landing

**Before landing**, the prototype's checks move into the repository: private validation and retry, and first-touch
failure order, into `tests/cpp/flow_ledger_tests.cpp`, including a two-edge hold case in which both an old-only and a
new-only edge fail; one-owner ordered allocation atomicity into `tests/cpp/apply_actions_fault_tests.cpp`. The dead
`reserved->empty() ? 0 :` arm is removed and its precondition stated, and the two private comments in
`ledger_access.hpp` are rewritten to state the rule they guard.

**Gates**, on the intended target tree: `bash .superset/workspace.sh check`, `make sanitize-test`, `make tsan-test`,
the full project CI, NetGraph's `dev/check_core_integration.sh`; the stateful, failure, retry and owner-lifetime checks;
tunnel default-hold coverage; a full-value replay of DC-BB; whole command-line export comparisons at DC-BB and NSFNET
(1 and 8 workers) and `backbone_clos` (1 and 28).

**Timing rows a landing reruns**, against the same target tree without the edit, each with its reference A/A at its own
scale and the required layout controls: DC-BB whole run at 1 and 8 workers; `backbone_clos` at 1 and 28; DC-BB's
Core's work at 1 and 8; NSFNET at 1 and 8; the 102 harness rows that run the changed code (its caller, the default
commit path and scope controls; listed in [measuring](measuring.md#harness-rows-a-landing-reruns)), plus
`batch.ordered.one_owner`; peak memory on each.

**Go/no-go.** Go when, on the landing build against its own reference with same-source layout controls, DC-BB's whole
run at one worker is faster beyond noise and not slower at eight; these five rows are not slower (upper end at most
1.03): `kernel.lossy.clos`, `kernel.eb.clos`, `kernel.weighted.clos`, `ledger.joint.pooled_cycle`,
`program.stale.lossy`; the 14 drift rows are measured and not slower; and the gates above pass. Otherwise no-go. A gain
measured on one tree does not carry to another tree or add to another change's gain.

## The prototype branch

Branch `prototype/ordered-publication`: one commit, `78a4a2f6ac6ddab54ba9fa988088df0b20c35360`, on
`9cd7642b5ea8681f7833c139969d8ed224d03a75`, three files; its tree, `afc3f0b58b29f44f8a740da4f01bb175d816b83b`, is the
tree that was measured. On a later `main` it needs rebasing, and the gates and timing rows above are then run on the
rebased tree.
