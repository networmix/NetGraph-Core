# Routing performance: where Core's time goes and what is worth changing

This directory is for maintainers. It records what a study of 2026-11-08 to 2026-11-10 established about making
NetGraph-Core's real work faster while every result, tie, error and state stays exactly as `docs/semantics.md` states
it: where the time of real workloads goes, which changes keep every result and why, what was tried and how it ended,
and two changes that are ready to land under stated conditions. Whether either change is in the product is recorded in
[`docs/roadmap.md`](../../roadmap.md), nowhere else.

## The answer

**At the worker counts users run by default, Core is a small part of the run.** On the DC-BB workload (below) at its
default 8 workers the whole run takes 10.40 s, of which scenario loading, NetGraph's own Python work, takes 6.94 s
(67 %) and Core's work 1.70 s. On NSFNET's all-pairs maximum flow the export alone takes 1.01 s of 1.87 s at one
worker (54 %). On `backbone_clos` the whole run is 0.34 s. Any Core change is bounded by these shares before its own
merits.

**A faster shortest-path algorithm cannot help on these inputs.** On DC-BB the shortest-path search is 10.1 % of
Core's time, so even a free search gives at most 1.11× on Core. Most of Core's time is in building the support and the
route programs (39.9 %), placement (25.8 %) and the ledger with its results (13.8 %). Both independent implementations
of the 2025-2026 algorithms that beat the sorting bound are slower than Dijkstra ([candidates](candidates.md)).

**At eight workers the analysis is limited on the Python side.** On DC-BB the caller's CPU outside Core calls is
1.68 s by one timer correction and 1.83-1.92 s by a later eight-thread control, against a 2.72 s analysis; inside
Core's calls the caller holds the interpreter lock for 0.19 s of CPU and releases it for 10.47 s. Whether the lock is
the critical path is not shown, so no finite eight-worker gain is inferred from a model anywhere in these pages.

**The one structural reserve in Core is reuse across failure cases.** Every one of DC-BB's 4,884 failure
demand-cases can take its route program from the failure-free baseline of its workflow step, and the work that reuse
could remove is f = 0.3485-0.3538 of the one-worker Core time before the reuse's own cost. The opt-in reuse below
takes it.

**Two changes are proposed**, each a separate landing; they were never measured together, so their gains are not
added:

- [Ordered publication reuses its totals](ordered-publication.md): a three-file local edit. DC-BB's whole run is
  3.8 % faster at one worker and shows no difference beyond noise at eight; every output is unchanged.
- [Baseline route-program reuse](baseline-reuse.md): an opt-in call. DC-BB's whole run is 18 % faster at one worker
  with a capture per workflow step and 4.6 % faster at eight workers with a capture per run (timed with native
  outputs, details, cost distribution and failure trace off); no single lifetime has been shown faster at eight
  workers and not slower at one, and it needs the maintainer's contract decision.

Nothing else passed: every other lead could not remove enough time on its workload even if its part cost nothing, was
built and gave no gain, or stays unknown with a named deciding experiment ([candidates](candidates.md)).

## The workloads

| Name | What it is | Default workers |
| --- | --- | --- |
| DC-BB (96 pods) | a data-centre fabric with a backbone (96 pods; 3,833 nodes, 91,712 directed edges), placement of two lossless equal-cost demands under 2,460 distinct cases; NetGraph's `docs/design/perf-0.24-vs-0.26/` input | 8 |
| NSFNET | NetGraph's `scenarios/nsfnet.yaml` (28 nodes): all-pairs maximum flow, 87,696 pair evaluations in 116 batches | 8 (Core's own threads) |
| `backbone_clos` | NetGraph's bundled `backbone_clos` scenario: six controller-placed demands under 1,000 failure draws, 592 distinct cases | 28 (`auto`, the CPU count) |
| DC-BB 192 pods | the same fabric at 192 pods (4,697 nodes, 148,544 directed edges), NetGraph's `dcbb192-0.26.0.yml.gz`; its shipped input sets 1 worker | 1 |
| Synthetic backbone | a generated backbone: 1025 nodes, 12,512 edges, 125 TE_MAXFLOW_UNLIM demands under 50 failure draws; also run with every demand on tunnel placement and with NetGraph's `paths` analysis (shortest and K shortest paths for 200 node pairs); generator parameters in [measuring](measuring.md#the-inputs) | `auto` (28 on the measuring machine); every quoted cost share was counted at 1 worker |
| Harness rows | Core's own timing rows (`dev/perf/`), named `<family>.<variant>.<graph>`: `kernel.lossy.clos` is the LOSSY placement kernel on a fixed Clos fabric | one thread |

DC-BB in detail: 3,600 failure draws in 18 workflow steps deduplicate to 2,460 distinct cases, 18 baselines and 2,442
failure cases; with two demands each that is 4,884 failure demand-cases (and 36 baseline ones). Its `apply_actions`
calls number 2,475: one ordered list per case plus 15 lists of the maximum-supported-demand search. Each demand joins
96 sources with one destination (and back). Costs are {1: 35,360, 5: 8,448, 10: 2,048} over its links.

NetLab sets its own worker counts: `netlab run` on TopoGen configurations (controller-placed traffic, 1,000 failure
draws) uses 7 workers; NetLab's DC-BB experiment runner uses NetGraph's default, the CPU count (28 on the measuring
machine); NetLab's autoresearch generator writes 8, DC-BB's shape. One NetLab experiment was measured (its shipped
DC-BB experiment: two sites, 1,280 nodes, 4,864 links), and no real geographic backbone ran: the synthetic backbone
stood in. Every split on these pages is a fact about these inputs, not about NetLab's geographic runs.

## Conditions of every number

Core `9cd7642` (code of release 0.13.0), NetGraph `8a641a1` (0.27.0), NetLab `32d5e13` (0.7.0); CPython 3.14.10 with
its interpreter lock; one Mac Studio (macOS, arm64, 28 CPUs, which `workers: auto` resolves to); November 2026. NetLab
supports CPython 3.11-3.13, and NetGraph and Core also support free-threaded 3.14t; on those interpreters, and where
`auto` resolves to another CPU count, the whole-run ratios are untested.

Each number carries one of four kinds of evidence:

- **qualified**: an A/B/A comparison under the timing protocol ([measuring](measuring.md)); only these decide a gain;
- **exploratory**: one run on one machine; it locates time and decides nothing;
- **counted**: an instrumented copy or a trace counted the work (searches, scans, bytes);
- **derived**: arithmetic on the others, shown where it is used.

"x % faster" means the time ratio candidate / reference is 1 − x/100. The phase times and splits below are exploratory;
the window medians of the two proposals are other runs of the same reference build and differ from them.

## Where the time goes

Whole-run phases, seconds, exploratory (one run each, no profiler):

| Workload / workers | Import | Load and model | Prepare | Run | Export | Total |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| DC-BB / 1 | 0.09 | 7.05 | 0.42 | 10.62 | 0.17 | 18.37 |
| DC-BB / 8 (default) | 0.09 | 6.94 | 0.41 | 2.72 | 0.23 | 10.40 |
| NSFNET / 1 | 0.09 | 0.05 | 0.00 | 0.69 | 1.01 | 1.87 |
| NSFNET / 8 (default) | 0.09 | 0.04 | 0.00 | 0.49 | 1.01 | 1.67 |
| `backbone_clos` / 1 | 0.07 | 0.07 | 0.01 | 0.23 | 0.03 | 0.42 |
| `backbone_clos` / 8 | 0.06 | 0.06 | 0.01 | 0.17 | 0.03 | 0.34 |
| `backbone_clos` / auto = 28 (default) | 0.06 | 0.07 | 0.01 | 0.17 | 0.03 | 0.34 |

**Core's work** is the union of the time intervals spent inside Core calls, measured by replaying the workload's
recorded Core calls; it is a diagnostic beside the whole run, never a user gain on its own. On DC-BB it is 9.35 s at
one worker and 1.70 s at eight (qualified as a noise row). Its partition comes from sampled stacks, so the shares are
stack shares, not CPU shares, and the two columns used different builds and threads, so they show no worker effect:

| DC-BB, part of Core's samples | 1 worker | 8 workers |
| --- | ---: | ---: |
| support and route-program construction | 39.9 % | 39.80 % |
| placement: admission probes, allocation, forwarding | 25.8 % | 26.22 % |
| ledger and result publication | 13.8 % | 13.00 % |
| shortest-path search (the reverse support search) | 10.1 % | 9.37 % |
| unresolved native and bindings | 5.7 % | 6.45 % |
| setup, masks, cleanup | 4.8 % | 5.16 % |

Counted on DC-BB: 4,884 searches, 446,462,698 adjacency-slot scans; 222,533,300 support edges and 222,125,828 plan
members built; 76,904 largest-fit probes and 81,788 allocations; identical at eight workers. Exact repeats are 1.6 %
of the traced native span.

**NSFNET**, per pair evaluation (exploratory): 4.44 µs, of which the flow solver 1.90 µs (residual Dijkstra 0.67,
blocking flow 0.41, potentials 0.05, closure and cut 0.14, solver setup and output 0.65) and the rest 2.54 µs (result,
record and program construction 0.64, roster and support search 0.43, binding conversion at least 0.43, ledger and
policy setup 0.36, unresolved 0.67). 39 repeated batches carry 29,484 of the 87,696 evaluations, each repeat in a
different capacity analysis.

**`backbone_clos`**, one worker: the flow solver is 68.4 % of Core-context samples (residual Dijkstra 23.6 %, blocking
flow 18.5 %, closure and cut 4.3 %, setup, output and other 21.2 %). As recorded, those parts sum to 67.6 % (derived),
not 68.4 %, and the source names no further part.

## Outside Core

These costs are NetGraph's and NetLab's own; they bound what any Core change gives the user (all exploratory).

- **Scenario loading**: 6.94 s of DC-BB's 10.40 s default run (67 %); 7.05 s of 18.37 s at one worker. Regular
  expressions, node selection and YAML parsing in Python.
- **Export of NSFNET**: 1.01 s of 1.87 s at one worker (54 %). It includes Core's getters and binding work as well as
  NetGraph's own serialisation. The getters were measured on two bases that are not one cost: one run at one worker
  recorded 1,227,789 getter events and 0.227 s excluding enum reads; a later run at eight workers recorded 1,315,488
  events, 0.289 s including 87,699 enum reads and 0.227 s excluding them. NetGraph's own object building: 393,244
  integer and 87,696 float constructions, 175,392 lists, 51 MB in all.
- **Python under the interpreter lock at 8 workers**: on DC-BB the caller's CPU outside Core calls, a proxy for
  lock-held work, is 1.68 s by one timer correction and 1.83-1.92 s by a later eight-thread control, against a 2.72 s
  analysis. Inside Core's spans the caller holds the lock for 0.19 s of CPU and releases it for 10.47 s. Adding latency
  to Core calls on the case threads reaches the wall in full. Free-threaded Python is a NetGraph-side option to
  measure, not a result.
- **NetLab**: on its shipped DC-BB experiment the placement steps are a sixth of the command (0.565 of 3.39 s).

## Not measured

- A real geographic backbone (see above).
- Lock waiting in the plain product; the lock split was measured only in instrumented runs.
- The eight-worker critical path, and so any finite eight-worker gain from a model.
- Garbage collection per Core call; Core's own memory apart from whole-process peaks; the bytes a reuse capture keeps
  alive.
- Other interpreters and CPU counts (see the conditions above).

## Terms

- **Support**: for a demand, every edge on some shortest path to its destination (a *tight* edge). **Roster**: the
  sources that take part in a demand. **Route program**: the split plan built over the part of the support the
  roster's sources reach, which placement then fills.
- **Ordered list**: a list of `apply_actions` placements committed one demand after another on one ledger.
  **Ledger**: the per-graph record of what each edge carries and holds. **Epoch**: the state one policy accumulates
  since its construction or its last remove.
- **Baseline**: the no-failure case of a workflow step. **Failure classes**, for a failure case whose mask only
  deletes elements of the baseline's mask, with `D_S` the part of the baseline support on walks from the roster
  sources: **class 1**, the failures miss `D_S`; **class 2**, every surviving roster source still reaches an allowed
  destination inside `D_S` minus the failures (then the case's support is `D_S` minus the failures pruned both ways,
  and its roster is the baseline roster minus the failed sources, with no search); **class 3**, every other case.
- **LOSSLESS** placement admits only an amount that fits on every route (found by largest-fit probes); **LOSSY**
  forwards a given offer in one pass. **ECMP**: equal split over equal-cost paths. **EQUAL** and **WEIGHTED**: route
  allocation equally or by weights. **FEWER_HOPS**: ties between equal-cost paths broken by hop count. **T1**
  (`DistributedTE`): tunnel placement by constrained shortest-path searches (CSPF). **T2**: sizing a demand by a tunnel
  count. **KSP**: the K shortest simple paths.

## Pages

| Page | For the reader who asks |
| --- | --- |
| [ordered-publication.md](ordered-publication.md) | what the ordered-publication edit changes, what it gains and what a landing must show |
| [baseline-reuse.md](baseline-reuse.md) | the opt-in reuse: its contract choices, gains, risks and the deciding comparison |
| [equivalences.md](equivalences.md) | which changes keep every result, under which assumptions, with proofs and counterexamples |
| [candidates.md](candidates.md) | "should we try X?": every lead's outcome and the literature |
| [measuring.md](measuring.md) | how to build, time and compare here, and which rows mislead |

The study's measuring instruments (the call recorder and replay, the per-call clocks, the lock split) are kept as they
ran on the branch `research/routing-performance-instruments`; the study's raw records, builds and captures were not
kept.
