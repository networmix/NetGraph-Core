# CPU optimisation implementation plan

Status: implemented in the working tree (P0–P2 and the P3 fallback review
below); not yet committed or released.
Consolidated on 27 September 2026 against NetGraph-Core `fbab119`, revised the
same day after the [queue-substitution experiment](./), then updated
with the implementation and its old/new/old measurements. This is the current
CPU plan, superseding the implementation-order guidance in the two Metal
research reports. Those reports and their raw measurements remain the historical
evidence; the earlier BFS/Dial prototypes are now understood as special cases of
one queue.

## Implementation status (27 September 2026)

What landed, as one change on top of `fbab119`:

- `StrictMultiDiGraph::from_arrays` records `max_cost()` and `cost_gcd()` in
  its existing cost scan; they follow copies and moves as members.
- `shortest_paths_core` and `shortest_paths_to_core` share one relaxation loop
  templated on the queue. `HeapQueue` is the reference binary heap on
  `(cost, -bottleneck, node)`; `BucketQueue` is the ordered bucket queue
  (`W = max_cost / gcd + 1` cyclic buckets, 64-bit occupancy bitmap, active
  bucket kept as a heap on `(-bottleneck, node)`). Dispatch is
  `W <= 65536`; `SpfQueue::Heap` or `NGRAPH_CORE_SPF_QUEUE=heap` forces the
  reference route, `Bucket`/`bucket` prefers the bucket queue when eligible.
- Scratch state is a thread-local workspace per direction with touched-node
  tracking and a hybrid reset (sparse when fewer than a quarter of the nodes
  were touched, wholesale otherwise; the same predicate chooses how outputs are
  materialised). An RAII lease restores the invariant on every exit path,
  including the reverse search's fanout exceptions. Arrays shrink when a much
  smaller graph is searched; entry arrays and buckets keep their capacity.
- Tests: `ShortestPaths.QueueDifferential_Forward` / `_Reverse` (forced heap
  against forced bucket over ten random graph families, the option matrix,
  five destinations, twice on the reused workspace; more than 25,000
  comparisons), `_ComparisonHasTeeth`, `BucketQueueKeepsCapacityOrderedPredecessors`,
  `WorkspaceIsResetAfterFanoutThrow`, `WorkspaceReuseAcrossGraphSizesAndEmptyGraphs`,
  `StrictMultiDiGraph.CostSummaryMaxAndGcd`, and a Python test that runs the
  same script under both `NGRAPH_CORE_SPF_QUEUE` values and unset.
- Documentation: README environment table and CHANGELOG (Unreleased). No
  Python API change; `_docs.py` is unchanged.

Gates passed on the final code (logs in [cpu_queue/](./)):

- Old-versus-new bit-gate ([bitgate/](bitgate/)): a C++ program that
  exercises forward and reverse SPF with fanout, `calc_max_flow` in every
  placement mode, `batch_max_flow`, `sensitivity_analysis`, `FlowPolicy`
  place/rebalance/remove sequences and KSP through the public API, built once
  against the `fbab119` sources and once against the working tree. 6,226
  output lines, identical SHA-256, also with the new build forced onto the heap.
- `bash .superset/workspace.sh check`: lint, types, 177 C++ tests, 466 Python
  tests. `make sanitize-test`: 177 C++ tests under ASan/UBSan.
- NetGraph `dev/check_core_integration.sh` against this worktree: a wheel built
  in a fresh environment, 1,278 NetGraph tests passed.
- The research harness, rebuilt against the new production sources, still
  matches its copies of the old loop across 1,351,168 comparisons on both
  routes.

One finding from the first old/new/old pass that the earlier experiment could
not show: the new code's *heap* route (the fallback for `W > 65536`) was 8–15%
slower than the old heap. Two causes were fixed: the scratch arrays were reached
through the thread-local struct, so their base pointers were reloaded around
every store (now hoisted to raw pointers for the search), and the heap
comparator was a static function that decays to a function pointer inside
`std::push_heap`/`pop_heap` (now a functor, as the old lambda was). The
measured tables below are from the final code.

A second finding concerns measurement, not Core: NetGraph's failure Monte Carlo
results differ between processes of the *same* build unless `PYTHONHASHSEED`
is fixed (set-iteration order), so the NetGraph driver pins it; with the seed
fixed, old and new digests agree ([hashseed-determinism.txt](results/netgraph/hashseed-determinism.txt)).

## Independent review (2 October 2026)

An independent reviewer (the `gpt-6-astra` model through the `codex` CLI, brief
and full report in [cpu_queue/review-astra/](review-2026-10/)) was
asked whether all easily achievable, reasonable optimisations had been
considered and verified. Its verdict: the algorithm choice is sound and the
weighted-graph and small-query gains reproduce independently, but several
small changes were missed, one exception path was wrong, and the retained
memory statement was false. Each finding was re-verified here before adoption
(harness parity on both routes, the old-versus-new bit-gate, the repository
gates and a fresh old/new/old sequence; tables below are from the final code).

Adopted:

- **Exception safety (reviewer finding 1, reproduced).** The source or target
  node's distance was written before its `touched.push_back`; a throwing
  allocation there left an unrecorded modification that the next same-sized
  search would see. Nodes are now recorded before their entries are written,
  and `WorkspaceLease` documents that rule.
- **Skip settled neighbours before scanning their parallel edges (finding
  3).** A settled node cannot be improved with non-negative costs and no longer
  accepts equal-cost or capacity-tie updates, so its group is skipped before
  any edge or bottleneck is examined. Dense unit fabrics, which the earlier
  passes recorded as neutral, now gain about 1.3x.
- **One-pass predecessor output (finding 4).** Implemented as a live per-node
  count maintained by the search (exact sizing, no over-reservation) rather
  than the reviewer's reserve-and-append variant; the second list walk is gone.
- **Inline single-item frontier (finding 5).** A frontier of one node is held
  inline and promoted into the buckets on the next push; chain graphs improve
  from about 2.1x slower to about 1.6x slower than the old heap.
- **Retention bound (finding 6, counterexample reproduced).** Bucket item
  storage accumulated the high-water marks of every bucket ever used. The queue
  now drops its item storage when it exceeds four times the last search's peak
  frontier plus 4,096 items, the growable buffers are released whenever the
  node arrays shrink, and the code comment and this plan describe the actual
  policy instead of "at most the largest graph".
- **Bitmap bookkeeping (finding 8).** `begin`/`end` touched every occupancy
  word (O(W/64) per search, measurable for tiny searches at large `W`); the
  words set during a search are now listed and only those are cleared.
- **Fixed-routing callers (finding 2).** `calc_max_flow` and
  `FlowState::place_max_flow` with `require_capacity=false` recomputed an
  identical cost-only search on every tier; the result is now computed once per
  invocation. Outputs are bit-identical (bit-gate).
- **Two redundant copies (findings 10 and 13).** `get_path_bundle` moves its
  DAG into the returned pair, and the `spf`/`spf_to` bindings move the DAG into
  the Python tuple instead of copying it.
- **Measurement corrections.** Reverse-search rows are now part of the Core
  sequence (P2's gate was previously marked done without them); the unit-cost
  rows are described as "no demonstrated regression" rather than "within
  drift"; the heap-route delta is recorded as unexplained (the reviewer's
  ablations of touched tracking, item layout and front-copy pop all showed no
  gain); and the earlier claim that a distance-only API would remove the O(N)
  output floor was wrong, since a distance vector still has N entries.

Recorded as follow-ups, not done:

- A narrower internal result for callers that only inspect `dist[dst]` or the
  DAG (max-flow, FlowPolicy, KSP spur searches), which is the only lever left
  for tiny searches on large graphs (finding 9).
- KSP's per-spur O(N+E) mask allocations and copies (finding 11), and
  `place_on_dag`'s quadratic parent grouping on wide ECMP DAGs (finding 14),
  both outside the search.
- `W_max` on a second platform, a two-level occupancy summary for very wide
  ranges, and the reviewer's boundary rows (at `W` = 65,535–65,536 full
  searches still gain about 1.4x while first-neighbour queries pay 9–17%
  against the new heap route), which argue against raising the limit from the
  full-search sweep alone (finding 8).
- A zero-copy int64 distance result in the bindings (finding 13), and the
  stale `has_min_flow` memo-key discriminator in `FlowPolicy` (finding 10).
- Refuted by the reviewer: a gcd early-out in graph construction (finding 12).

### Measured against the previous production build

Both sequences ran one fresh process per phase in the order old (A1), new,
new forced onto the heap route, old (A2), on the final code (after the
independent review). The old build is the `fbab119` sources; the harness
binaries and the NetGraph extension modules are hashed in each
`environment.txt`. Load averages were 3.8–5.9 during both sequences because
other sessions were active; the desktop was not isolated, so differences under
about 5% are within run-to-run noise here. Earlier sequences on the
pre-review code are kept in the same directories' history
(the reviewer's own runs).

Core level ([results-core](results/core/), harness `ref` = the
production function of each build, sequential, nine samples, medians in ms):

| Workload: case / nodes / edges / queries / destination | W | Old A1 / A2 | New | New, forced heap | Old / new |
|---|---:|---:|---:|---:|---:|
| random / 1,024 / 8,192 / 1 / none | 32 | 0.107 / 0.102 | 0.067 | 0.089 | 1.52–1.60 |
| random / 4,096 / 32,768 / 256 / none | 32 | 211.6 / 207.6 | 132.6 | 221.3 | 1.57–1.60 |
| random / 16,384 / 131,072 / 64 / none | 32 | 254.4 / 250.6 | 150.9 | 266.1 | 1.66–1.69 |
| random / 65,536 / 524,288 / 1 / none | 32 | 21.21 / 19.03 | 11.51 | 23.14 | 1.65–1.84 |
| random, degree 32 / 4,096 / 131,072 / 64 / none | 32 | 112.0 / 110.6 | 77.1 | 115.4 | 1.43–1.45 |
| random, masked / 16,384 / 131,072 / 64 / none | 32 | 230.0 / 229.2 | 140.9 | 236.7 | 1.63 |
| random / 262,144 / 2,097,152 / 1 / none | 32 | 136.6 / 133.8 | 96.2 | 130.3 | 1.39–1.42 |
| random / 1,024 / 8,192 / 64 / farthest | 32 | 10.81 / 10.65 | 7.09 | 11.21 | 1.50–1.52 |
| random / 16,384 / 131,072 / 64 / farthest | 32 | 251.6 / 250.5 | 149.8 | 264.0 | 1.67–1.68 |
| random / 65,536 / 524,288 / 16 / farthest | 32 | 324.3 / 330.2 | 201.2 | 336.0 | 1.61–1.64 |
| uniform cost 31 / 16,384 / 131,072 / 64 / none | 2 | 181.4 / 177.1 | 146.3 | 178.6 | 1.21–1.24 |
| unit grid / 10,000 / 39,600 / 1 / none | 2 | 0.614 / 0.618 | 0.459 | 0.635 | 1.34–1.35 |
| unit fabric / 800 / 320,000 / 64 / none | 2 | 143.9 / 144.3 | 110.2 | 111.1 | 1.31 |
| pseudo src/sink fabric / 802 / 320,200 / 64 / fixed pair | 2 | 143.4 / 143.9 | 108.0 | 110.3 | 1.33 |
| chain / 2,048 / 2,047 / 1 / none | 32 | 0.021 / 0.021 | 0.027 | 0.029 | 0.77–0.78 |
| random / 1,024 / 8,192 / 64 / first neighbour | 32 | 0.116 / 0.115 | 0.043 | 0.038 | 2.70–2.72 |
| random / 16,384 / 131,072 / 64 / first neighbour | 32 | 1.160 / 1.158 | 0.402 | 0.392 | 2.88 |
| random / 65,536 / 524,288 / 64 / first neighbour | 32 | 4.351 / 4.254 | 1.506 | 1.517 | 2.83–2.89 |
| random / 262,144 / 2,097,152 / 16 / first neighbour | 32 | 4.509 / 4.053 | 1.491 | 1.471 | 2.72–3.02 |
| reverse: random / 4,096 / 32,768 / 64 | 32 | 55.85 / 55.34 | 34.23 | 58.82 | 1.62–1.63 |
| reverse: random / 16,384 / 131,072 / 16 | 32 | 66.55 / 66.34 | 38.67 | 69.33 | 1.72 |
| reverse: unit fabric / 800 / 320,000 / 16 | 2 | 45.43 / 45.72 | 37.15 | 36.90 | 1.22–1.23 |
| reverse: chain / 2,048 / 2,047 / 1 | 32 | 0.002 / 0.003 | 0.001 | 0.001 | 2.57–2.61 |

Cost-range sweep (random 16,384 / 131,072, 16 full-source queries, costs in
`1..K`), same sequence:

| K | W | Old A1 / A2 | New | New, forced heap | Old / new |
|---:|---:|---:|---:|---:|---:|
| 1 | 2 | 44.4 / 45.3 | 36.5 | 44.9 | 1.22–1.24 |
| 31 | 32 | 64.5 / 62.6 | 37.2 | 65.9 | 1.68–1.73 |
| 1,023 | 1,024 | 59.2 / 58.9 | 35.3 | 64.5 | 1.67–1.68 |
| 32,767 | 32,768 | 51.3 / 51.2 | 42.7 | 53.3 | 1.20 |
| 1,048,575 | heap fallback | 49.7 / 49.5 | 52.8 | 51.6 | 0.94 |
| 33,554,431 | heap fallback | 51.0 / 49.6 | 51.2 | 52.8 | 0.97–1.00 |

NetGraph level ([results-netgraph](results/netgraph/),
[bench_netgraph.py](bench_netgraph.py) run by NetGraph's interpreter
against each build with `PYTHONHASHSEED=0`, five samples, medians in ms, every
workload's result digest identical across all four phases):

| Workload | Old A1 / A2 | New | New, forced heap | Old / new |
|---|---:|---:|---:|---:|
| full-source `spf`, 2-tier Clos 200×200 (400 sources) | 214.2 / 212.7 | 163.4 | 164.6 | 1.30–1.31 |
| full-source `spf`, weighted grid 100×100 (50 sources) | 220.8 / 220.3 | 181.7 | 245.4 | 1.21–1.22 |
| bound-context `max_flow`, Clos 200×200, combine, ×50 | 52.6 / 51.9 | 47.1 | 46.8 | 1.10–1.12 |
| `max_flow` with 100 link-exclusion sets, Clos 100×100 | 44.3 / 43.0 | 40.4 | 39.3 | 1.06–1.10 |
| bound-context `max_flow`, weighted grid, combine, ×20 | 202.9 / 202.1 | 174.6 | 219.1 | 1.16 |
| `max_flow`, weighted grid, pairwise corners, ×5 | 204.0 / 200.6 | 208.2 | 212.9 | 0.96–0.98 |
| demand placement, backbone_clos matrix, ×50 | 69.2 / 67.4 | 69.0 | 68.7 | 0.98–1.00 |
| placement Monte Carlo, backbone_clos, 300 iterations | 74.5 / 75.3 | 67.8 | 71.4 | 1.10–1.11 |
| max-flow Monte Carlo, backbone_clos, 300 iterations | 29.3 / 29.1 | 27.6 | 29.0 | 1.06 |
| KSP k=8, 3-tier Clos, 56 leaf pairs | 54.7 / 53.4 | 46.3 | 46.2 | 1.15–1.18 |
| `max_flow`, 3-tier Clos, 56 pod pairs, combine | 167.5 / 166.1 | 168.2 | 171.8 | 0.99–1.00 |
| placement, 3-tier Clos, 64 WCMP + 8 ECMP demands | 366.4 / 363.3 | 376.3 | 368.9 | 0.97 |
| placement with 20 link-failure sets, 3-tier Clos | 7,646 / 7,548 | 7,640 | 7,571 | 0.99–1.00 |
| scenario square_mesh.yaml, 100 iterations, serial | 10.3 / 10.3 | 9.9 | 9.9 | 1.04 |
| scenario nsfnet.yaml, 200 iterations, serial | 958.2 / 939.1 | 929.6 | 912.2 | 1.01–1.03 |
| scenario backbone_clos.yml, 100 iterations, serial | 82.3 / 81.1 | 78.3 | 77.6 | 1.04–1.05 |

Reading of the two sequences:

- **Weighted graphs** (`W` between 32 and 32,768) gain 1.4–1.8x at the Core
  level, forward and reverse, and 1.16–1.22x end to end on the weighted grid
  workloads in NetGraph, where graph preparation, Python and placement share
  the time.
- **Unit-cost fabrics** (`W = 2`), the class NetGraph's Clos scenarios and
  pseudo source/sink max-flow graphs belong to, now gain 1.2–1.35x at the Core
  level (settled-neighbour skip and one-pass output, adopted from the review)
  and 1.06–1.3x through NetGraph's SPF, max-flow, KSP and Monte Carlo paths.
  The wide-ECMP placement rows (0.97–1.00) are not a demonstrated regression:
  they lie inside the run-to-run band and are dominated by placement, whose
  quadratic parent grouping is a recorded follow-up.
- **Small destination-limited searches** gain 2.7–3.0x from the workspace
  regardless of the queue (the forced-heap column matches the bucket column).
- **Fallback class** (`W > 65,536`): the two fallback rows are 0.94 and
  0.97–1.00 of the old heap; the forced-heap column on eligible graphs is
  0.82–0.97, which is the same unexplained route penalty. The reviewer's
  ablations (touched tracking, item layout, front-copy pop) showed no gain,
  so its cause remains open; it is recorded, not rationalised.
- **Chain graphs** are now 0.77x (about +6 µs on a 21 µs query) on the bucket
  route after the inline single-item frontier, from 0.47x before.
- **Correctness**: every NetGraph workload digest was identical across old,
  new and the forced heap; the Core parity, differential tests and the
  old-versus-new bit-gate remain the semantic evidence, since these digests
  only summarise distances and parent counts.

## What the experiment established (27 September 2026)

Harness: [cpu_queue/](./) is plain C++ (no Objective-C++ or Metal) and
compiles the production `shortest_paths.cpp` beside a copy of its forward loop
with the queue and workspace injected ([spf_variants.hpp](spf_variants.hpp),
[harness.cpp](harness.cpp)). Results are in
[cpu_queue/results](results/first/) and [cpu_queue/results-short](results/short/);
the earlier research results are untouched.

**Parity.** 1,351,168 comparisons against the production `shortest_paths`,
bit-exact on distances, `parent_offsets`, `parents` and `via_edges`, over: twelve
masked random seeds, six graphs with zero-cost cycles, int64-offset costs (which
exercise the heap fallback), a unit grid, uniform cost 31, costs to 100,000, a
fabric with NetGraph-style zero-cost pseudo source/sink attachments, a chain, a
dense fabric and a 12-node multigraph with self-loops, parallel edges and
residuals at multiples of `kMinCap`; times `{multipath, single-path}` ×
`{multi-edge, single-edge}` × `{Deterministic, PreferHigherResidual}` ×
`{require_capacity}` × `{residual}` × `{masks}` × destinations `{none, first
neighbour, farthest, self, arbitrary}`, each run twice on a reused workspace.
Zero mismatches for every queue/workspace variant. The checker has teeth:
dropping the bottleneck key from the bucket comparator produces 167,636
mismatches. ASan/UBSan runs of the same matrix are clean
([plain](results/first/correctness.txt),
[sanitized](results/first/correctness-sanitized.txt),
[teeth](results/first/correctness-teeth.txt)).

**Timing.** Sequential, single thread, one process per variant, production
reference first and last (A1/A2), nine samples after two checked warmups,
medians in milliseconds for the whole query set, every timed result verified
against the reference outside the timed interval. `bucket` is the ordered bucket
queue; `fresh` allocates scratch per query as production does; `sparse` reuses a
workspace and resets only touched nodes; `hybrid` picks full or sparse reset
from the touched count after the search. Load average was 2.9–4.1 on 14 cores
throughout; the desktop was not isolated.

| Workload: case / nodes / edges / queries / destination | W | Production A1 / A2 | bucket, fresh | bucket, sparse | bucket, hybrid |
|---|---:|---:|---:|---:|---:|
| random / 1,024 / 8,192 / 1 / none | 32 | 0.101 / 0.104 | 0.078 | 0.069 | — |
| random / 4,096 / 32,768 / 256 / none | 32 | 198.7 / 205.8 | 117.7 | 122.1 | — |
| random / 16,384 / 131,072 / 64 / none | 32 | 238.8 / 245.6 | 133.4 | 140.3 | 142.5 † |
| random / 65,536 / 524,288 / 1 / none | 32 | 18.38 / 19.63 | 10.96 | 10.68 | 10.80 † |
| random, degree 32 / 4,096 / 131,072 / 64 / none | 32 | 106.7 / 109.0 | 64.1 | 65.8 | — |
| random, masked / 16,384 / 131,072 / 64 / none | 32 | 219.5 / 225.3 | 130.6 | 138.4 | — |
| random / 262,144 / 2,097,152 / 1 / none | 32 | 119.7 / 125.2 | 83.6 | 87.0 | 84.5 † |
| random / 1,024 / 8,192 / 64 / farthest | 32 | 10.13 / 10.43 | 6.26 | 6.55 | — |
| random / 16,384 / 131,072 / 64 / farthest | 32 | 240.0 / 246.3 | 134.8 | 140.8 | — |
| random / 65,536 / 524,288 / 16 / farthest | 32 | 289.4 / 299.9 | 162.7 | 176.1 | — |
| uniform cost 31 / 16,384 / 131,072 / 64 / none | 2 | 170.3 / 175.4 | 142.4 | 140.7 | — |
| unit grid / 10,000 / 39,600 / 1 / none | 2 | 0.614 / 0.622 | 0.466 | 0.497 | — |
| unit fabric / 800 / 320,000 / 64 / none | 2 | 142.1 / 143.6 | 126.1 | 127.9 | 142.1 † |
| pseudo src/sink fabric / 802 / 320,200 / 64 / fixed pair | 2 | 140.6 / 142.2 | 125.5 | 127.2 | — |
| chain / 2,048 / 2,047 / 1 / none | 32 | 0.021 / 0.021 | 0.045 | 0.045 | 0.041 † |
| random / 1,024 / 8,192 / 64 / first neighbour | 32 | 0.105 / 0.120 | 0.111 | 0.043 | — |
| random / 16,384 / 131,072 / 64 / first neighbour | 32 | 1.044 / 1.171 | 1.190 | 0.393 | 0.390 † |
| random / 65,536 / 524,288 / 64 / first neighbour | 32 | 4.034 / 4.207 | 4.332 | 1.489 | 1.484 † |
| random / 262,144 / 2,097,152 / 16 / first neighbour | 32 | 3.872 / 3.997 | 4.090 | 1.449 | 1.442 † |

† from the [second, shorter sequence](results/short/), whose own
production A1/A2 were 240.6 / 245.4, 18.51 / 18.77, 137.8 / 141.6, 0.021 / 0.023,
1.047 / 1.104, 4.206 / 4.305, 3.913 / 19.193 and 119.5 / 124.2 for the rows marked.

Cost-range sweep, random 16,384 nodes / 131,072 edges, 16 full-source queries,
costs uniform in `1..K` ([sweep files](results/first/)):

| K | W | Production | bucket, sparse | Ratio |
|---:|---:|---:|---:|---:|
| 1 | 2 | 43.7 | 35.2 | 1.24 |
| 31 | 32 | 62.1 | 35.8 | 1.74 |
| 1,023 | 1,024 | 58.2 | 34.1 | 1.71 |
| 32,767 | 32,768 | 50.0 | 35.0 | 1.43 |
| 1,048,575 | 1,048,576 | 48.6 | 46.8 | 1.04 |
| 33,554,431 | fallback | 48.9 | 56.5 (harness heap) | — |

Interpretation to carry into implementation:

- **Weighted graphs, `W` ≤ 32,768: 1.4–1.8x** over the production heap under
  the full contract, sequentially, whether the query is full-source, farthest
  destination or masked. Production A1/A2 drift was 2–7% on these rows (up to
  14% on sub-millisecond rows), so these clear the decision gate on this
  machine. The 1,024-node single query improves 1.3–1.5x.
- **Uniform-cost graphs: 1.1–1.3x, and inconclusive on the dense fabric.** With
  `W = 2` the queue is a small share of the work; edge scanning dominates the
  800-node fabric with 320,000 edges (0.97–1.13x across the two sequences). The
  unit grid improves 1.32x, uniform cost 31 improves 1.2x.
- **Destination-limited queries that touch few nodes are bounded by O(N)
  per-query overhead, not by the queue.** Production spends about 1.6 µs at
  1,024 nodes, 16 µs at 16,384, 63 µs at 65,536 and 242 µs at 262,144 per
  first-neighbour query; the queue variant alone does not change this. Reusing
  the workspace with a full reset changes nothing (consistent with the earlier
  thread-local measurement recorded in the perf notes), while resetting only
  touched nodes removes about 60% (2.4–3.0x). The remaining floor is the
  N-sized `dist` and `parent_offsets` outputs that the API returns.
- **Hybrid reset is the right default:** it matches sparse on small queries and
  full reset on full-graph queries, where sparse tracking costs 3–5%.
- **Chain graphs regress 1.7–2.1x (+15–24 µs on a 21 µs query)** in this
  pre-implementation experiment (later reduced to 1.3x by the inline
  single-item frontier adopted from the review). Every pop
  advances the active bucket, so the per-pop bookkeeping (occupancy bitmap,
  heap rebuild of a one-element bucket) exceeds a one-element binary heap.
  Replacing the modulo by a power-of-two mask did not remove it. This is the
  smallest absolute regression in the matrix and no production caller runs on
  a path graph, but it is repeatable and must be recorded, not hidden.
- **Sorting each bucket once instead of heaping it is slower** on weighted
  graphs (up to 8%) and equal on uniform graphs, except the unit grid (0.393 vs
  0.466). It is not worth a second queue.
- **The harness's own heap copy runs 7–14% slower than production** on the
  random workloads. The cause is not isolated (candidates: the touched-node
  branch per relaxation, runtime mode checks, scratch arrays reached through a
  struct instead of locals). Production work must therefore be measured against
  production, not against the harness heap, and P1 must check that its
  workspace integration does not carry this penalty. The bucket gain is far
  larger than the penalty.
- **Decomposition of the earlier prototype gain.** The research Dial package
  reached 3.6x over native on the 16k/64 batch with ten workers and a relaxed
  contract (no capacity ordering, canonical incoming-CSR DAG, cached eligibility
  bitmap). The contract-preserving queue alone gives 1.8x sequentially. Roughly
  half of the prototype's advantage came from dropping the contract and is not
  available to production; the other half is.
- **One unattributed transient.** In two of twelve fresh-allocation phases the
  262,144-node first-neighbour set took 19 ms instead of 4 ms (once in the
  harness heap copy, once in the production reference A2); none of fourteen
  reused-workspace phases did. A page-fault probe in a fresh process showed
  neither faults nor slowness. It is recorded as environmental variance, not as
  a property of either implementation.
- Timings are sequential and cover synthetic graphs; ten-worker batch timings
  from the research are not comparable row for row. No end-to-end flow
  improvement is measured yet.

## What the first consolidation overlooked

1. **Every production caller is destination-limited.** `calc_max_flow`,
   `FlowState::place_max_flow`, `FlowPolicy::get_path_bundle` and the KSP spur
   searches all pass `dst`; only NetGraph's cost/DAG-cache paths call `spf`
   without one. The earlier dispatch table sent destination-limited searches
   to the existing implementation, which would have excluded the callers that
   dominate placement time. Early-exit parity is now verified, so
   destination-limited queries are in scope from the first enablement.
2. **Zero-cost edges are the normal case downstream.** NetGraph attaches pseudo
   source/sink nodes with cost-0 augmentation edges, so "any zero-cost graph
   stays on the heap" would have excluded every max-flow analysis. Zero-cost
   pushes land in the active bucket's heap and are verified against the native
   settled-node acyclicity rule; the dispatch condition is only `W`.
3. **The dominant fixed cost of small queries is O(N) scratch initialisation**,
   which no queue change touches; touched-node reset addresses it and was not
   in the previous plan (which specified a complete reset).
4. **Uniform non-unit costs** are covered by the gcd stride rather than by a
   separate BFS; they were previously flagged as untested.
5. **Reverse SPF** (`shortest_paths_to_core`) has the same queue and takes the
   same substitution mechanically; it is not yet measured. KSP's private
   `dijkstra_single` is a cost-only search on a different key and is out of
   scope for the first change.
6. **Ordered parity is not an open research question**; it follows from queue
   monotonicity. The plan no longer needs "investigate per-layer ordering" or
   "a LIFO vector is not sufficient" as open items.

## Contract and fallback rules

Keep `Cost=int64_t`, `Cap=double`, the unreachable sentinel and graph cost-sum
bound. Preserve validation, empty graphs, invalid/blocked endpoints, masks,
reverse DAG orientation, parallel-edge IDs, multipath selection, single-path
capacity ties, destination early exit and `fanout_edges` behaviour. The queue
key stays `(distance, -bottleneck_capacity, node_id)` with the existing
`kEpsilon` stale-entry and re-push semantics; boolean eligibility cannot replace
the capacity values used for ordering. SPF excludes `rem < kMinCap` when
capacity gating applies; a supplied residual forces that gate.

Dispatch is a graph property, decided before the search starts:

| Condition | Route |
|---|---|
| `W = max_cost / gcd + 1 ≤ W_max` (any mix of zero and positive costs, any selection mode, any destination, forward or reverse once P2 lands) | Ordered bucket queue |
| `W > W_max` (wide or int64-offset costs) | Existing heap |
| Forced reference (tests, benchmarks) | Existing heap |

`W_max` starts at 65,536 (2^16): the sweep shows 1.4x at `W = 32,768` and
break-even near `W = 2^20`. Memory per workspace is `W` empty vector headers
plus a `W`-bit occupancy bitmap, about 1.6 MB at the limit. Widening the range
(two-level occupancy summary, or a radix heap whose last level is ordered the
same way) is later work with its own sweep. The graph-level `max_cost` is an
upper bound for every masked query, so no per-query rescan is needed; a
masked-out wide edge merely wastes empty buckets.

## Implementation sequence

### P0 — Parity oracle and forced route (done)

The harness is the parity oracle. Before production work, add the equivalent
differential test to [tests/cpp](../../../tests/cpp/shortest_paths_tests.cpp):
the forced-heap route against the dispatched route over random graphs, the
option matrix and the destination set above, bit-exact on all four output
arrays, plus a guard test proving the comparison rejects a wrong ordering.
Forcing must be reachable from tests and benchmarks without a public API
(an internal function parameter or a build/env switch used only there).

Keep the measured originals and hashes intact; write new results to new
directories. Dynamic profiling of NetGraph callers (full versus
destination-limited counts, touched fraction, weight range) remains useful for
the downstream gate but no longer blocks the first enablement, because both
call classes are covered.

**Exit:** differential test with teeth in the C++ suite; forced route available.

### P1 — Queue and workspace in the forward search (done)

Work in [shortest_paths.cpp](../../../src/shortest_paths.cpp) and
[strict_multidigraph.cpp](../../../src/strict_multidigraph.cpp):

- Record `max_cost` and the cost gcd in `from_arrays`'s existing cost scan as
  members, so they follow copies and moves; no process-global cache keyed by a
  graph address.
- Replace the priority queue in `shortest_paths_core` with the ordered bucket
  queue when `W ≤ W_max`; keep the heap code path intact for the fallback and
  the forced route. Keep the relaxation loop shared between routes (one loop,
  two queue types), as the harness does, so parity is structural.
- Give the search a reusable workspace (distances, bottlenecks, predecessor
  heads/tails, settled flags, entry arrays, queue buckets) with touched-node
  tracking and the hybrid reset. Concurrency: `batch_max_flow` and
  `sensitivity_analysis` run on `std::async` threads and Python callers release
  the GIL, so the workspace must be per thread (thread-local) or caller-owned;
  it must hold no pointer into a graph, residual or mask. Bound retained memory:
  one workspace per thread sized to the largest graph seen, and state the bound
  (or add a shrink rule) in the code comment and in `_docs.py` if it is
  observable. Returned vectors must own their data.
- Early exit leaves items in buckets; clear only occupied buckets.

**Exit:** the existing positive-cost SPF output hash and the flow output hash
from the August review are unchanged; the differential test passes; the full
suite and `make sanitize-test` pass; A/B/A against `main` on the harness
workloads reproduces the table above within drift, including the recorded
chain regression; no new implicit threads.

### P2 — Reverse search and KSP (done for the reverse search)

Apply the same substitution to `shortest_paths_to_core` (same key, in-adjacency)
and measure it with reverse workloads added to the harness. Leave KSP's
`dijkstra_single` alone unless profiling shows it matters.

**Exit:** parity for reverse results including `fanout_edges`, and a measured
result either way.

### P3 — Limits, regressions and dispatch review (fallback reviewed; second machine class pending)

Confirm `W_max` on a second machine class (Linux/x86 CI-class hardware), and
decide the chain-like regime: the inline single-item frontier cut the recorded
regression to about 6 µs on a 21 µs query; either accept that residual on
frontier-of-one searches or reduce the per-pop advance cost and re-measure. Do
not choose the queue per query retrospectively; the rule stays a graph property.

**Exit:** dispatch rule frozen with a held-out workload pass.

### P4 — Downstream validation, then batching and GPU reconsideration (validation done; batching/GPU unchanged)

Run NetGraph's profiles and representative demand/flow runs with the chosen
Core build, including residual updates, mask changes, cache hit/miss cases,
all placement modes and path-cost constraints. For max-flow behaviour, check
feasibility, conservation and optimality/cut certificates with nonuniform
costs, masks and cost-factor cases; SPF distance agreement alone does not
satisfy this gate.

Explicit CPU batching stays conditional on demonstrated independent-query
callers: deterministic result order, bounded worker count, per-worker scratch
and bounded output memory. Sequential flow augmentations with changing
residuals cannot be batched as if independent. Any public batch or
distance-only extension must preserve current defaults and update bindings and
`_docs.py` together; a distance-only result would also remove the O(N) output
floor measured above, which is the only remaining lever for tiny queries on
large graphs. It is not a prerequisite for shipping this change.

Revisit GPU only against the accepted CPU implementation and matching result
contracts; this plan does not require Metal, GPU frontier queues or automatic
CPU/GPU selection.

## Acceptance and validation

Correctness gates precede performance claims. Cover forward/reverse results,
both edge tie-break policies, multipath/single-path, multi-edge/single-edge,
parallel edges, self-loops, zero-cost cycles, isolated vertices, empty graphs,
invalid and blocked endpoints, all-masked graphs and malformed spans. Include
costs above `2^53`, near the constructor's bound, `W` at and just above
`W_max`, long chains, duplicate/stale queue entries, uniform non-unit weights,
and capacity values at/around `kMinCap` and the tie epsilon. Mutate masks and
residual contents in place between calls, including changes that leave
eligibility unchanged. Verify graph lifetime, workspace reuse across graphs of
different sizes, concurrent callers, repeated determinism, fanout validation,
target early exit and capped enumeration.

Relevant existing suites include [C++ SPF](../../../tests/cpp/shortest_paths_tests.cpp),
[masking](../../../tests/cpp/masking_tests.cpp),
[Python edge selection](../../../tests/py/test_spf_edge_select.py),
[path enumeration](../../../tests/py/test_paths_resolve.py),
[reverse SPF](../../../tests/py/test_spf_to.py),
[memoisation](../../../tests/py/test_flow_policy_spf_memo.py),
[thread safety](../../../tests/py/test_thread_safety.py), and
[lifetime safety](../../../tests/py/test_lifetime_safety.py).

For each change, run correctness-checked old/new/old in fresh processes on one
quiet machine, same build settings, with builds/tests excluded from the timed
run. Record raw samples, hashes, workload parameters, load, fallback counts and
peak memory. Keep warm reused-query timing, first-query timing and
changing-input end-to-end timing separate. Ship a class only when the new median
improves beyond both A1/A2 drift and sample variability, with no unrecorded
regression in its supported classes or reference fallback. A single M4 result
is insufficient to set `W_max` for every supported platform.

Reproduce the experiment from the repository root:

```sh
bash dev/perf/spf_queue/build.sh
build/perf/cpu-queue correctness
bash dev/perf/spf_queue/run_aba.sh
python3 dev/perf/spf_queue/summarize.py
bash dev/perf/spf_queue/build.sh sanitize
ASAN_OPTIONS=detect_leaks=0 build/perf/cpu-queue-sanitize correctness
```

For implementation changes, run:

```sh
bash .superset/workspace.sh setup  # when the worktree environment needs setup
bash .superset/workspace.sh check  # rebuild extension, lint, types, C++ and Python
make sanitize-test               # required for C++ changes
bash /Users/networmix/ws/NetGraph/dev/check_core_integration.sh \
  /absolute/path/to/the/chosen/NetGraph-Core-worktree
```

Use the actual worktree's absolute path for downstream verification. Local
checks do not replace CI. The prior research passed 170 C++ and 465 Python
tests plus sanitizers; that is baseline evidence, not completion of any
implementation gate.

**Immediate next step:** review and commit the working-tree change (a
version bump and CHANGELOG date at release time), then confirm `W_max` on a
Linux/x86 CI-class machine before treating the chain-graph and wide-cost
fallback numbers as portable. Batching and GPU work remain later, conditional
items.
