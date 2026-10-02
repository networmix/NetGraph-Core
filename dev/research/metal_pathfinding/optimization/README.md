# GPU optimization follow-up

27 September 2026, NetGraph-Core `fbab119`, Apple M4 Max, 32 GPU cores,
14 CPU cores, 36 GiB unified memory. Standalone research; no backend registration
or production algorithm changes.

**Yes, the first investigation missed important optimizations.** Its poor batch
results were partly an implementation limitation, especially serial CPU
predecessor construction. They are not evidence of a GPU performance ceiling.
The experiments also show why changing to a frontier or bucket algorithm alone
does not guarantee improvement. Work organization, result construction and the
CPU comparison algorithm matter together.

**The CPU algorithm comparison also changes the decision.** A follow-up CPU
Dial-style bucket implementation brings the 16k weighted batch from the native
CPU's 27.7 ms to 7.7 ms. GPU results across two passes are 3.9–7.5 ms. On the
65k single-query case the stronger CPU takes 4.4 ms and GPU takes 1.2–1.4 ms.
Improve CPU specialization first; pursue GPU acceleration selectively.

The consolidated [CPU optimisation plan](../../../perf/spf_queue/PLAN.md) is the implementation
source of truth. After the [ordered bucket queue experiment](../../../perf/spf_queue/) it
specifies one queue substitution inside the existing search (BFS and Dial are
its special cases), verified bit-exact against the production result contract,
with a hybrid-reset workspace, a width-based dispatch rule and fallback, and
downstream acceptance requirements. That change is now implemented in
`src/shortest_paths.cpp`; the plan records its old/new/old measurements at the
Core level and through NetGraph.

## What changed in the evidence

Fresh processes ran sequentially: CPU A1, original GPU A1, optimized GPU B,
original GPU A2, CPU A2. Seven measured samples plus two warmups per row;
every result was checked outside the timed interval. CPU batches use ten
reusable workers. Times include host distances and positive-cost predecessor
DAG construction, but exclude graph preparation and pipeline creation.

| Workload | Native CPU A1 / A2, ms | Original GPU A1 / A2, ms | Best tested optimized GPU, ms | Observation |
|---|---:|---:|---:|---|
| Weighted random, 1,024 nodes / 8,192 edges, one query | 0.086 / 0.097 | 0.915 / 0.972 | 0.599 | CPU still wins small individual queries |
| Weighted random, 4,096 / 32,768, 256 queries | 21.391 / 21.281 | 47.812 / 49.310 | 7.003 | Interleaved 64-bit pull makes batching promising |
| Weighted random, 16,384 / 131,072, 64 queries | 27.693 / 27.721 | 57.748 / 53.106 | 7.473 | Earlier batch pessimism is no longer supported |
| Weighted random, 65,536 / 524,288, one query | 19.316 / 19.240 | 5.184 / 5.854 | 1.373 | Guarded 32-bit pull plus GPU DAG extraction helps |
| Unit-cost fabric, 800 / 320,000, 64 queries | 17.042 / 16.525 | 69.896 / 70.228 | 8.185 | CPU BFS plus optimized DAG is already 8.296 / 7.885 ms |
| Unit-cost grid, 10,000 / 39,600, one query | 0.570 / 0.582 | 3.836 / 5.874 | 3.378 | CPU BFS plus DAG is 0.089 / 0.088 ms |
| Weighted chain, 2,048 / 2,047, one query | 0.017 / 0.017 | 4.834 / 4.908 | 4.847 | Almost no useful parallel work per frontier |

These are **exploratory measurements, not accepted performance claims**. The
desktop was not isolated: pre-phase CPU idle readings ranged from 73.28% to
87.88%, with nonzero background GPU activity. The original GPU grid result
drifted substantially between A1 and A2. A quiet-machine repeat is still required
by this repository's performance standard. Nothing was stopped or altered in
the user's other sessions to obtain these measurements.

The best method is selected after measurement, not by a working runtime
dispatcher. Original GPU columns select between the original pull and persistent
group methods; the original shared-memory method was not rerun. The complete
[18-case table and all 216 optimized rows](results/aba/summary.md),
[raw samples](results/aba/optimized-b.csv), and
[environment/source hashes](results/aba/environment.txt) preserve the scope.

## Stronger CPU comparisons

The first table compares against the current native backend. A second A/B/A
compares GPU with CPU heap and Dial-style cyclic-bucket searches using cached
eligibility and one-pass DAG assembly. Each CPU worker completes both search and
DAG for a query, retaining locality. This is a comparable experimental contract
to the GPU: exact distances and canonical-equivalent positive-cost DAG edge sets,
without CPU discovery-order parity. Both devices exclude eligibility preparation.
The weighted random graphs have integer costs 1–31; the CPU Dial prototype is
explicitly limited to maximum edge weight at most 1,024, and retains int64 sums.

| Workload | CPU cached heap + DAG A1 / A2, ms | CPU Dial + DAG A1 / A2, ms | Best GPU first / second B, ms |
|---|---:|---:|---:|
| Random 1,024 / 8,192, one query | 0.053 / 0.059 | 0.028 / 0.030 | 0.599 / 0.634 |
| Random 4,096 / 32,768, 256 queries | 13.030 / 13.024 | 7.056 / 6.925 | 7.003 / 4.001 |
| Random 16,384 / 131,072, 64 queries | 17.500 / 17.268 | 7.696 / 7.677 | 7.473 / 3.851 |
| Random 65,536 / 524,288, one query | 11.474 / 11.530 | 4.384 / 4.372 | 1.373 / 1.244 |
| Random 4,096 / 131,072, 64 queries | 6.672 / 6.662 | 3.854 / 3.865 | 2.051 / 2.215 |
| Unit-cost fabric 800 / 320,000, 64 queries | 5.125 / 5.391 | 5.421 / 5.537 | 8.185 / 8.372 |

The second B repeated all 216 optimized GPU rows, unchanged binary, between
CPU alternative A1/A2. [All second-pass results](results/cpu-algorithms-aba/summary.md)
and [raw GPU samples](results/cpu-algorithms-aba/b.csv) are retained. The winning
GPU method changed on the 4k/256 and 16k/64 batches: 32-bit pull plus GPU DAG
beat the first pass's winning 64-bit pull. Background CPU idle reached 59.82%
before one phase. We have not attributed this variability to a specific cause;
it prevents a precise crossover rule or a stable speedup claim.

This strengthens the practical recommendation: CPU optimizations are worthwhile
independently of GPU support. GPU remains promising for large single-source and
some batch workloads, but comparisons only against the existing CPU heap backend
overstate the advantage over a better CPU implementation. On the dense fabric,
the fused cached CPU heap + DAG is faster than either GPU pass. On the grid,
the earlier CPU BFS result remains faster than the heap and Dial alternatives.

## Optimizations we missed

1. **Result construction on the CPU.** The first prototype scanned every edge
   twice and built every query's DAG serially, while the CPU search comparator
   already had ten workers. A cached-eligibility, one-pass incoming-CSR DAG
   builder, parallel across queries, reduced the 16k/64-query persistent-group
   case from approximately 53–58 ms to 13.239 ms without changing the GPU search
   algorithm. This variant bundles three changes; it does not isolate the
   contribution of threading from one-pass construction or cached eligibility.
2. **Batch memory layout.** Making query the contiguous dimension reduced the
   same workload from 15.959 ms with the query-major 64-bit pull to 7.473 ms.
   Device time fell from 10.305 to 1.950 ms. Both versions use the same CPU DAG
   builder and identical relaxation counts. This improvement does not require
   narrowing distances or a new graph algorithm.
3. **Device-side DAG extraction.** Counting tight incoming edges, host prefix
   sums, then GPU filling helped the large single-query and higher-degree batch
   cases. It is not universally better: query-major frontier plus GPU DAG took
   15.944 ms on the 800-node fabric batch, versus 9.324 ms with CPU DAG assembly.
   Layout and output volume require their own decision, separate from SSSP.
4. **Exact 32-bit specialization.** Narrow distances only when
   `max_edge_cost * node_count < UINT32_MAX`; otherwise run the 64-bit pull.
   This conservatively proves that the sentinel cannot collide with finite
   shortest distances. No floating-point conversion is used. It helped some
   cases, but 64-bit pull beat 32-bit pull on the 16k batch, so narrower is not
   automatically faster. Eligibility still uses CPU double-precision residuals.
5. **Command and host overhead.** Multiple dispatches share a compute encoder
   with explicit buffer barriers; convergence flags are aggregated by SIMD
   group. These changes were tested as a bundle. An initial screening defect
   repeatedly invoked the Objective-C buffer accessor inside the conversion
   loop. Hoisting it removed large batch host overhead before final timings.
   The retained [screening notes](results/screen-known-limitations.txt) distinguish
   those preliminary measurements from the final matrix.

Apple's compute documentation supports batching dispatches and respecting
resource synchronization. Its Metal 4 command allocators and reusable command
infrastructure offer another way to reduce encoding overhead. Metal 4 itself
was **not implemented or timed** here.
[Metal calculations](https://developer.apple.com/documentation/metal/performing-calculations-on-a-gpu),
[Metal 4 command management](https://developer.apple.com/videos/play/wwdc2025/254/).

## Should we change algorithms?

**Weighted, general graphs: retain both pull and frontier candidates.**
Davidson et al. study Workfront Sweep, Near-Far and Bucketing as different
balances between redundant edge work and scheduling overhead. That is the right
comparison space for this problem; ordinary parallel Dijkstra is not an obvious
GPU default. Their NVIDIA results are not M4 speed predictions.
[2014 SSSP paper](https://escholarship.org/uc/item/8qr166v2).

The new active-frontier implementation supports one thread or one SIMD group
per active vertex, indirect dispatch and queue deduplication. It reduced the
16k/64-query relaxation work from 201.327 million adjacency entries to about
21–22 million. Yet its best full time was 10.144 ms, versus 7.473 ms for
interleaved pull. On the 10k grid it reduced edge scans from 7.920 million to
39,600, but still took 4.040 ms versus pull's 3.378 ms. It retained 200 rounds
and 25 host waits. Less edge work is real; the latency win is not automatic.

**Distance buckets remain worth investigating, but our simple version is not a
win.** The delta-4 experiment reduced the 16k batch to 8.451 million edge visits,
but needed 72 rounds and 16.825 ms. Delta-16 and delta-64 were also slower than
the best pull. This implementation scans all pending vertex slots to find and
select each bucket, uses a batch-wide minimum and bulk synchronization, and has
neither a light/heavy edge split nor a sophisticated asynchronous queue.
These results reject this implementation as the default, not delta-stepping.

ADDS specifically addresses bucket management, asynchronous work scheduling
and dynamic delta selection. It is a plausible next algorithmic investigation,
but a correct Metal port needs a progress-safe queue and synchronization design;
copying CUDA persistent-worker assumptions would not establish correctness.
Its published RTX 2080 Ti speedups do not establish gains on Apple hardware.
[ADDS, PPoPP 2021](https://www.cs.utexas.edu/~lin/papers/ppopp21.pdf).

**Uniform positive costs: test BFS on the CPU as well as the GPU.** The new CPU
BFS comparator was checked against native Dijkstra and uses the same cached
eligibility and optimized DAG construction as the GPU variants. It is a major
improvement on grids and removes most of the apparent GPU advantage on the
dense fabric. Direction-optimizing BFS can switch between sparse frontiers and
pulling from unvisited vertices. Its early exit after finding one parent applies
to distance discovery; all ECMP predecessors still require extraction.
[Direction-optimizing BFS](https://people.eecs.berkeley.edu/~krste/papers/beamer-sc2012.pdf).

The CPU plan now preserves predecessor order by construction: the ordered
bucket queue reproduces the heap's pop order exactly, and the unit-cost case is
its `W = 2` instance. Under that full contract the sequential gain on uniform-cost
graphs is 1.1–1.3x and inconclusive on the dense fabric, where edge scanning
dominates; the larger research BFS numbers came from the relaxed contract.

**Stable, road-like graphs with many repeated queries: consider preprocessing.**
PHAST uses contraction-hierarchy preprocessing followed by an ordered sweep.
The paper explicitly requires enough queries to amortize preprocessing and
focuses on low-highway-dimension graphs. For NetGraph, changing masks, residual
eligibility and weights mean shortcut validity or customization must be handled;
that applicability assessment is our inference, not a measured PHAST result.
It is not the first choice for arbitrary changing network simulations.
[PHAST](https://www.microsoft.com/en-us/research/wp-content/uploads/2011/01/phast_ipdps.pdf).

**Small or destination-limited searches: keep the CPU path.** These tests are
full-source searches; they do not compare against target-early-exit Dijkstra.
Bidirectional search or A* with a proven admissible lower bound can be relevant
to point-to-point workloads, but neither was implemented or timed here. Likewise,
CPU radix queues for wider integer costs remain an unmeasured stronger comparator;
the bounded-weight Dial alternative is measured above.

## Recommended implementation order

Follow the [CPU plan](../../../perf/spf_queue/PLAN.md): add the differential parity test with a
forced heap route, land the ordered bucket queue and hybrid-reset workspace in
the forward search as one change gated on the existing output hashes, extend to
the reverse search, freeze the width limit, then validate downstream. Batching
and a new distance-only API are conditional on actual callers, not prerequisites
for single-query CPU gains.

After those gates, revisit the GPU candidates against the accepted CPU baseline:
interleaved exact 64-bit pull, guarded 32-bit paths, and separately selected result
construction. Adaptive frontier/pull selection and improved bucket scheduling
remain later GPU research. The existing evidence does not establish full-backend
compatibility or justify automatic CPU/GPU dispatch.

## Validation and remaining limits

- Final A/B/A: 384 rows, 124,704 checked query results including warmups.
- Correctness: 12 seeded masked graphs and 12 boundary fixtures, all 14 variants,
  11,760 checked query results. Fixtures exercise forward/reverse traversal,
  parallel edges, self-loops, blocked sources, unreachable nodes, residuals just
  below and at `kMinCap`, distances near the 32-bit boundary, guarded fallback,
  and integer costs greater than `2^53`.
- Normal stress: 39,424 checked results. ASan/UBSan plus Metal API and GPU shader
  validation: another 11,760 correctness and 8,960 stress results. No reported
  sanitizer or validation errors. These repetitions are not independent graphs.
- The separate CPU comparison repeats masked positive-cost checks against native
  Dijkstra, and runs ASan/UBSan over both its masked fixtures and all seven
  performance workloads. See the [audit](results/audit.json) for exact counts.
- Across both final timing sequences and the listed correctness/stress runs:
  290,192 checked query results, including CPU comparator rows and warmups.
- Final repository check rebuilt the Python extension and passed lint, types,
  170 C++ tests and 465 Python tests. `make sanitize-test` passed all 170 C++
  tests. An initial formatting failure in the new summarizer was corrected and
  both subsequent workspace checks passed. These are local checks; CI was not
  run. No production tracked files changed, and no commit or push was made.
- All distances compare exactly. Positive-cost predecessor edge sets compare
  canonically to native CPU results; CPU discovery order is not reproduced.
- Graphs are synthetic and modest in size. No realistic demand traces, cold
  production latency, concurrent calls, per-query mask batches, cost-factor flows,
  zero-cost DAG parity or complete backend compatibility were established.
- GPU buffers and cached eligibility are prepared before timing. `upload_ms`
  includes the superset of buffers needed by all experimental variants, not a
  minimal production upload. Printed `setup_ms` is inherited phase-one pipeline
  setup and excludes the additional optimization library; it is not total cold
  initialization. GPU timestamps and distance timing are diagnostic subsets,
  not substitutes for full wall time. Do not subtract separate medians to claim
  an exact component cost.
- Queue atomics, bucket scans, vertex initialization and DAG extraction are not
  included in the `edge_visits` counter. The counter measures SSSP adjacency work.
- The prototype requires SIMD width 32 and 256-thread groups, checked at runtime,
  and imposes experimental size limits. It is not a portable library implementation.

Reproduce from the repository root:

```sh
bash dev/research/metal_pathfinding/optimization/build.sh
build/metal_pathfinding/optimized all correctness 3
build/metal_pathfinding/optimized all stress 20
bash dev/research/metal_pathfinding/optimization/run_aba.sh
python3 dev/research/metal_pathfinding/optimization/summarize.py
bash dev/research/metal_pathfinding/optimization/build.sh sanitize
MTL_DEBUG_LAYER=1 MTL_SHADER_VALIDATION=1 ASAN_OPTIONS=detect_leaks=0 \
  build/metal_pathfinding/optimized-sanitize all correctness 3
MTL_DEBUG_LAYER=1 MTL_SHADER_VALIDATION=1 ASAN_OPTIONS=detect_leaks=0 \
  build/metal_pathfinding/optimized-sanitize all stress 3
bash dev/research/metal_pathfinding/optimization/build_cpu.sh
build/metal_pathfinding/cpu-algorithms 3 correctness
bash dev/research/metal_pathfinding/optimization/run_cpu_aba.sh
python3 dev/research/metal_pathfinding/optimization/summarize_cpu.py
bash dev/research/metal_pathfinding/optimization/build_cpu.sh sanitize
ASAN_OPTIONS=detect_leaks=0 build/metal_pathfinding/cpu-algorithms-sanitize 3 correctness
ASAN_OPTIONS=detect_leaks=0 build/metal_pathfinding/cpu-algorithms-sanitize 2
make sanitize-test
bash .superset/workspace.sh check
```

Run timings sequentially without builds or tests in parallel. The compiler emits
one known warning because the reused phase-one `main` is renamed for inclusion;
that unused helper lacks an explicit return after its success path. It is not
called by this harness. The original experiment and shader hashes are unchanged.
