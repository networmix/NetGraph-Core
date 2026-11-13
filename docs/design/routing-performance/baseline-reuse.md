# Baseline route-program reuse (opt-in)

An opt-in call through which failure cases reuse the route programs derived once from their workflow step's baseline. It
adds public API, so it needs the maintainer's contract decision before anything lands; with that decision, one qualified
comparison decides it. The Core prototype is the branch `prototype/baseline-reuse` of this repository; the NetGraph
adapter it was measured with is NetGraph's branch `prototype/baseline-reuse-adapter`. Workloads, terms and evidence
labels are defined in the [README](README.md).

## The change

In a failure analysis NetGraph places the same demands under many failure masks, each removing elements from the
step's baseline. Today every case runs a new reverse search, builds a new support and new route programs. With the
opt-in call NetGraph captures the baseline's policies once (per workflow step, or once per run); for each failure case
an exact test decides whether it is class 2 for each demand, and if so the case's program is the baseline's program
restricted and pruned to the surviving source-to-destination walks. Placement is always recomputed on the case's own
ledger: copying the baseline's placed amounts is wrong (an equal-cost edge lost with every distance unchanged halves a
lossless split, 20 → 10). Any other case takes today's fresh path.

The call as prototyped (the prototype's documented example; not in the product):

```python
g = graph(2, [(0, 1, 10, 1), (0, 1, 10, 1)])  # two parallel edges; graph() and config() are the docs' small helpers
base = ngc.FlowPolicy(ngc.FlowLedger(g), 0, 1, config(ngc.ShortestPaths()))
before = base.replace(100)
capture = ngc.ShortestPathsReuse([base])
mask = g.eligibility(excluded_edges=[0])
hint = capture.prepare(mask)[0]
case = ngc.FlowPolicy(ngc.FlowLedger(g), 0, 1, base.config, mask, preparation=hint)
print(before.delivered, case.replace(100).delivered, hint is not None)   # 20 10 True
```

## Where it is exact, and what the prototype accepts

Proved for every shortest-path option at the support level and, for the program, for `multipath=True`,
`parallel_edges=True`, EQUAL or static WEIGHTED; refuted for `parallel_edges=False` (two parallel edges, the baseline
keeps one, the case fails it: the restricted program delivers 0, a fresh one 10) and for `multipath=False`. An empty
case roster must report the change flags a fresh evaluation reports. The proofs and counterexamples are in
[equivalences, "Failure classes and baseline-derived
programs"](equivalences.md#failure-classes-and-baseline-derived-programs).
Nobody oracle-ran simultaneous lists or the optimised controller end to end.

The prototype and its adapter accept LOSSLESS ShortestPaths with EQUAL allocation, equal sources, FEWER_HOPS,
multipath and parallel edges, without weights, and only when the case's mask contains the capture's exclusions. Static
WEIGHTED is covered by the theorem but is neither implemented nor timed. It, LOSSY and every other configuration get no
hint and run the ordinary path, and a supplied cache keeps its ordinary lookup. The result depends only on the graph,
the capture's mask, the case's mask, the endpoints and the configuration, so it is the same for a per-step and a per-run
capture. On DC-BB every one of the 4,884 failure demand-cases is class 2 and its restricted program equals a fresh one;
on DC-BB 192 pods (4,697 nodes, 148,544 directed edges; two opposite-direction demands with baseline rosters of 192
and 1 sources; 18 workflow steps whose draws sum to 200 per demand in each step and give 2 to 200 distinct failure
placements per step) every one of its own 4,884 failure demand-cases is class 2 and its restricted support equals the
fresh one byte for byte.

## What the test costs, and what a capture holds

Per physical case the test compares the mask's node and edge flags with the capture's; per demand it copies one counter
per node and initialises its scratch, then propagates the removed support members through the affected region; the
case program's materialisation is charged separately. The dense comparison and initialisation come before the
propagation, so the test's cost grows with the graph and the support, not with the failure: on DC-BB 192 pods the
per-node counter copies (22,940,148) and restriction visits (362,744,448) are about 370 times the propagation's 104,370
queued nodes and 933,091 incoming-member visits (counted; the ratio derived).

The test terminates: each node is queued at most once and each support member is removed at most once. Its counters are
int32 and bounded by the number of edges; it adds no amount or cost sum. For D demands with roster sizes R_i on a graph
of N nodes and E edges, preparing one case is O(N + E + Σ_i (N + E + R_i log(R_i + 1))). A capture holds O(D (N + E))
besides the graph and the programs it shares; a per-run capture keyed by K distinct baselines holds the sum of K
captures for the life of the process, with no eviction. The bytes a capture keeps alive were not measured.

Size: product code and type stubs +412/−12 lines, `docs/semantics.md` +52/−1, tests +795/−3; 20 files, five new. The
prototype passed an oracle of 2,720,512 rows, 2,000 seeded interaction sequences, 12 planted bugs caught, the
repository's gates with both sanitizers, and NetGraph's downstream suite with the opt-in off.

## What NetGraph changes

The measured adapter adds 155 and removes 5 lines in `ngraph/analysis/placement.py`. Besides passing `preparation=`,
it captures the baseline policies before they are removed, derives each case's hints, threads the shared preparation
through its state and case functions, and solves each step's baseline before the failure cases unless an earlier
step's capture is reused. That is a scheduling change: NetGraph's test
`test_converged_only_builds_no_baseline_and_schedules_cases_alone` fails by design in both lifetimes. It keeps per-run
captures for the life of the process in a module-level dictionary keyed by graph and baseline masks, is switched by
the environment variable `NETGRAPH_REUSE_LIFETIME` (the prototype's default is `step`, that is on), and applies only to
ordered, converged lists of lossless equal-split shortest-path demands. The adapter branch also changes one test double
in `tests/analysis/test_failure_states.py` to accept keyword arguments.

A NetGraph landing chooses the user-facing switch, the owner and lifetime of the capture, and accepts or redesigns the
scheduling change. A consumer that keeps one process across scenarios (NetLab's pipeline reuses its process pool) needs
a release or an owned context for per-run captures and a retained-memory check over several scenarios. The measured
gains include the adapter as built.

## Gain

NetGraph's command line, start to exit, with that adapter; each row is its own lifetime and must not be paired with
another.

| Row | 1 worker: candidate/reference [95 % interval] | 8 workers (DC-BB's default): the same | Evidence |
| --- | --- | --- | --- |
| DC-BB, capture per step | 0.8171 [0.8162, 0.8189]: 18 % faster, 17.81 → 14.57 s | 0.9593 [0.9551, 0.9615]: inconclusive (reference drift 4.6 %) | qualified |
| DC-BB, capture per run | 0.8121 [0.8115, 0.8247]: inconclusive (reference drift 7.3 %) | 0.9544 [0.9529, 0.9656]: 4.6 % faster, 10.01 → 9.57 s; layout control 0.9959 | qualified |
| DC-BB 192 pods (its own input sets 1 worker; the 8-worker rows use a copy changed only in that value), per run / per step | 0.7829 [0.7820, 0.7833] / 0.7829 [0.7809, 0.7848]: 22 % faster; two separate measurements, equal to four decimals (0.782866 and 0.782892) | 0.9415 [0.9376, 0.9464] / 0.9454 [0.9390, 0.9482]: 5.9 % / 5.5 % faster; layout controls 0.9961 [0.9938, 1.0046] and, at 1 worker, 1.0040 [1.0004, 1.0044], within noise | qualified |
| NetLab's shipped DC-BB experiment (whole runner, 28 workers) | — | 1.0011 [0.9983, 1.0082]: no difference | qualified |
| Tiny synthetic graphs, by cases per capture (2 nodes with two parallel edges of capacity 16 and cost 1, one demand; 6 nodes and 6 edges, 3 roots, 2 ordered demands at offer 9) | 2 cases: 1.1611 [1.1486, 1.1697] and 1.1933 [1.1781, 1.2049], 16 % and 19 % slower; 20 and 200 cases (6-node series): 1.0914 [1.0852, 1.0995] and 1.0744 [1.0658, 1.0785], 9 % and 7 % slower; the 2-node series at 20 and 200 inconclusive by drift (1.1271, 1.1270) | — | qualified |

The DC-BB 96- and 192-pod timing inputs had native outputs, details, cost distribution and failure trace disabled;
other output configurations were not timed.

A diagnostic over capture, derivation, the full failure-case `apply_actions` and release is about 20 % faster in both
lifetimes; it does not isolate program construction. The start-up and scheduling controls are within noise. On
DC-BB 192 pods the gain stays below the 1.3× bar a secondary workload needs. In NetLab's experiment (16 to 895
placements per baseline, placement a sixth of the command) nobody counted whether its cases used the reuse.

Peak memory flags no risk (DC-BB per run at 8 workers 674.9 / 650.5 MB; DC-BB 192 pods per run at 8 workers 835.7 /
720.8 MB; NetLab 1,157.8 / 1,133.2 MB). The bytes the capture keeps alive are not separated from these peaks.

## Risks

- With the opt-in path off, the shared-code harness row `concurrency.delivered` is 15.0 % slower beyond noise (about
  5.9 → 6.8 ns per getter call) and its layout control is exactly 1.0000, so the cause is not code placement; the
  source does not show it.
- `kernel.eb.clos` is 4.9 % slower but not separable from placement; `kernel.weighted.clos` and
  `ledger.joint.pooled_cycle` are unmeasured (drift).
- The break-even (how many cases a capture must serve to pay back) is unknown on every row; a controlled sweep of
  cases per baseline would decide it.

## The maintainer's contract choices

If nobody decides, nothing lands. The recommended answer of each is the default this page assumes.

| Choice | Example | Options and cost | Recommended |
| --- | --- | --- | --- |
| 1. Adopt the opt-in at all | the call above | **yes**: a new public type (`ShortestPathsReuse`), a keyword, about 460 product and documentation lines, the NetGraph adapter, and the deciding comparison; **no**: nothing changes, the measured gains are forgone | yes, if the deciding comparison passes for the chosen lifetime |
| 2. Lifetime and ownership | NetGraph keeps one capture per step, or one per run for the whole process | **per step**: simple ownership, 18 % at one worker measured; **per run**: 4.6 % at eight workers measured, but captures live as long as the process and a persistent consumer needs a release point | choose one, then run the deciding comparison on it |
| 3. Equality rule for mixed lists | a list holding a tunnel-count policy built on a support and an ECMP policy, both given the same mask, is refused because their effective masks differ: `ValueError: baseline policies must share one Eligibility` | **document** the effective-mask rule (no code); or **keep the caller's mask** beside it (a code change and new tests) | document it |
| 4. Failure guarantees of the new calls | an allocation failure inside `capture.prepare(mask)` leaves no partial group and no changed ledger; a caller out of memory gets `MemoryError` and can retry | **accept** the guarantee as written (add the new calls to the documented checked paths); or **design another** binding or ownership path first | accept |
| 5. The pybind11 limit | an allocation failure at the last step of `ShortestPathsReuse([base])` escapes as an outer `MemoryError`, as it already does for `ShortestPaths()`: the new type inherits `docs/semantics.md`'s "Limit 1: pybind11's own object machinery" | **disclose** it in the contract (no code); or **redesign** the factory with its own proof | disclose |
| 6. Old-call error text | an invalid old `FlowPolicy(...)` call's automatic `TypeError` text gains `, *, preparation: _netgraph_core.ShortestPathsPreparation = None` | **accept** with the keyword; or a **separate factory** to keep the old text (more API for a diagnostic) | accept |

## The deciding comparison

On the exact sources proposed for landing, DC-BB's whole run with the chosen capture lifetime must be faster beyond
noise at 8 workers and not slower at 1, in one qualified comparison. The study's windows were complete and met the
protocol, so they were not repeated on the same tree; this comparison decides.

**Gates and timing rows at landing.** The repository gates as for [ordered publication](ordered-publication.md#landing);
NetGraph's downstream suite with the adapter in each lifetime; the comparison above; the 87 shared-code harness rows
(listed in
[measuring](measuring.md#harness-rows-a-landing-reruns)) with `concurrency.delivered` not slower; DC-BB 192 pods at 1
and 8 workers; peak memory, and retained memory over several
scenarios for a persistent consumer. Prototype work a performance round could remove (a full template program, an
unused empty-roster program and a mask copy per case) changes the code, so its gains would need new measurements.

## The prototype branches

| Branch | Repository and base | Commit and tree | Contents |
| --- | --- | --- | --- |
| `prototype/baseline-reuse` | NetGraph-Core, on `9cd7642b5ea8681f7833c139969d8ed224d03a75` | `589e8414fb487adc78f71fd703b0666520d138c3`, tree `8f9da55980f039c637ddde41ac3816d36e7f5afd` | 20 files: the new type and keyword, `docs/semantics.md` "Preparing deletion cases", stubs, C++ and Python tests |
| `prototype/baseline-reuse-adapter` | NetGraph, on `8a641a1caa0b8f4b7d8c4df0aecf8625310cdf47` | `03fd213e5545b48523046f440c76337111a20890`, tree `dd68bc346aa0a2eb0014b35c1cd525ac79bc4b48` | the measured adapter: `ngraph/analysis/placement.py` (+155/−5; its SHA-256 is `15e0e5ec85ff467390c3771e6bde5935135dbc856000157dbf014164792a3d94`, the measured file) and one test double (+4/−4) |

Each is one commit on its base and holds the tree that was measured. On later sources they need rebasing, and the
deciding comparison then runs on the rebased trees.
