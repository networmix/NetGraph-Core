# Metal pathfinding investigation

**Follow-up:** the [optimization research](optimization/README.md) tests active
frontiers, buckets, batch layout and improved result construction. It finds
substantially better batch results, so the initial timings below must not be
treated as a GPU performance ceiling. Both rounds remain exploratory because
the desktop was not isolated for measurement.

The [consolidated CPU optimisation plan](../../perf/spf_queue/PLAN.md) now governs
implementation order and acceptance. The backend integration notes below remain
relevant to a later Metal follow-up.

Research against NetGraph-Core `fbab119` on 26–27 September 2026, using an
Apple M4 Max (14 CPU cores, 32 GPU cores, 36 GB unified memory).
This directory contains standalone experiments. It does not register a backend,
change the Python API, or modify production algorithms.

**Answer:** adding a Metal implementation is technically feasible, and the host
integration is fairly small. A useful, compatible replacement for `Backend.cpu()`
is a substantially larger task. The candidate interface is **distances on an
immutable graph, with optional batches**, particularly for large low-hop
full-source searches, and explicit CPU fallback. It is premature
to make Metal the default SPF implementation or promise faster flow simulation.

**What the current backend gets right.** `StrictMultiDiGraph` already has immutable
forward and reverse CSR, stable internal edge IDs, nonnegative int64 costs, and
a constructor bound on the sum of costs that prevents path-cost overflow.
These are useful GPU inputs. `Backend` and `Algorithms` provide a dispatch seam;
`FlowPolicy` calls that seam. The CPU implementation is native C++ Dijkstra,
with flat predecessor storage and Python GIL release. Its existing memoization
and reverse SPF also avoid computations that a GPU would otherwise repeat.
The appropriate comparison is this implementation, including parallel CPU
queries, not a Python Dijkstra loop.

**Where the abstraction needs work.** These are verified properties of the
current source, not newly claimed correctness defects:

| Area | Current behavior | Consequence for Metal |
|---|---|---|
| Graph ownership | `GraphHandle` contains only `shared_ptr<const StrictMultiDiGraph>`. CPU construction is non-owning; Python separately keeps the graph alive. | No place for a typed device representation or backend identity. Add owned device state and preserve topology lifetime. |
| Query interface | `spf` / `spf_to` are synchronous, single-query calls returning host vectors and a `PredDAG`. There is no batched SPF or distances-only API. | Cannot amortize dispatch across independent queries through the current interface. Every call produces a host DAG even if the caller only needs distances. |
| Cost semantics | Internal distances are int64; Python offers int64 explicitly and defaults to float64 conversion. | A float32 GPU implementation would change path selection. Keep exact integer arithmetic. |
| Predecessors | Positive-cost ECMP keeps all equal alternatives; zero-cost edges use settle order to remain acyclic. Parallel-edge selection and capacity ties are configurable. | A distance kernel plus tight-edge extraction is insufficient to reproduce the full contract. |
| Result order | Predecessor arrays follow CPU discovery order. Path enumeration consumes that order, including when capped by `max_paths`. | The prototype's deterministic edge-ID order has equivalent positive-cost edge sets, but is not byte-identical or necessarily enumeration-identical. |
| Residuals and masks | Borrowed per-call spans; capacities and residuals are double; a residual span forces capacity gating. Python currently copies residual/mask input. | Compute exact eligibility on the CPU or reproduce these semantics explicitly. Static topology upload does not remove changing-mask work. |
| Max-flow | `calc_max_flow` calls the free `shortest_paths` function directly between residual updates. | Overriding `Backend::spf` does not accelerate `Backend::max_flow`, batch max-flow, or sensitivity when those delegate to the existing CPU functions. |
| KSP | Yen-style enumeration has its own Dijkstra helper and directly calls free SPF for spurs. | Overriding backend SPF does not make KSP run on the GPU. It is also sequentially dependent across accepted paths. |
| FlowPolicy | Calls `ctx_.algorithms->spf`, generally with `dst` and changing flow state; caches some SPF results. | It can use a new backend, but single-call latency, early exit, cache hits, and residual invalidation matter. Batching across dependent placements would change semantics. |
| Downstream selection | NetGraph's context builder explicitly constructs `Backend.cpu()`. | Shipping a factory alone would not enable acceleration downstream. |

Source pointers: [backend.hpp](../../../include/netgraph/core/backend.hpp),
[algorithms.hpp](../../../include/netgraph/core/algorithms.hpp),
[shortest_paths.hpp](../../../include/netgraph/core/shortest_paths.hpp),
[shortest_paths.cpp](../../../src/shortest_paths.cpp),
[flow_policy.cpp](../../../src/flow_policy.cpp), [max_flow.cpp](../../../src/max_flow.cpp),
[k_shortest_paths.cpp](../../../src/k_shortest_paths.cpp),
[bindings](../../../bindings/python/module.cpp).
Downstream inspection was read-only: `ngraph/analysis/context.py`,
`ngraph/analysis/placement.py`, and `dev/perf/profiles.py` in the sibling NetGraph
checkout. No downstream behavior or integration speedup was tested.

**What the external research establishes.** Apple supplies a direct C++ Metal
interface, so there is no need to introduce Swift or a tensor framework just to
call compute kernels. An isolated Objective-C++ translation unit also works;
that is what this experiment builds using Command Line Tools and the system
Foundation/Metal frameworks. [Apple Metal-cpp](https://developer.apple.com/metal/cpp/)

The current Metal language specification supports 64-bit `long`/`ulong`, but
does not support `double`. Its 64-bit atomics are limited: void-returning
`atomic_min_explicit` / `atomic_max_explicit` exist; a general 64-bit atomic
fetch-min/load API is not available. This is a reason to choose a suitable
algorithm, not a reason to abandon exact costs. Apple lists 64-bit atomics from
Apple9 hardware, whereas 64-bit integer math has broader availability. Capability
checks must cover the target hardware family.
[MSL specification, sections 2.1 and 6.16.4.6](https://developer.apple.com/metal/Metal-Shading-Language-Specification.pdf),
[Apple feature tables](https://developer.apple.com/metal/Metal-Feature-Set-Tables.pdf).
The [local compiler probes](results/capabilities.txt) confirm all of these
distinctions on this M4 Max, including successful void-returning 64-bit minimum.
Compilation alone is not a performance test of that atomic operation.

MPSGraph is a graph of tensor operations, not an off-the-shelf implementation of
NetGraph's routing semantics. MLX offers custom Metal kernels, which is useful
for prototyping, but still leaves the graph algorithm and correctness contract
to us. I found no drop-in implementation of this API in these libraries.
[Apple MPSGraph](https://developer.apple.com/documentation/metalperformanceshadersgraph?changes=__5_7),
[MLX custom Metal kernels](https://ml-explore.github.io/mlx/build/html/dev/custom_metal_kernels.html).

Gunrock provides useful GPU graph-processing designs based on frontiers and
load balancing, but its advertised backends are CUDA and HIP, not Metal.
cuGraph likewise depends on NVIDIA's CUDA stack. Neither is a library switch
for an Apple GPU. Their existence does suggest a more ambitious implementation
direction: frontier compaction and work-efficient SSSP rather than repeatedly
scanning every edge. That direction was researched, not implemented or timed.
[Gunrock](https://github.com/gunrock/gunrock),
[cuGraph installation](https://docs.nvidia.com/cugraph/26.10/installation/getting_cugraph/).

**The experiments.** All three kernels compute exact int64 distances with
double-buffered Bellman–Ford-style pull relaxation over incoming CSR. Reverse
routing uses outgoing CSR. They do not race on distance writes and do not need
64-bit atomics. Per-edge eligibility is prepared on the host from finite double
residuals, node masks, and edge masks. The cost of this preparation and buffer
allocation is recorded as `upload_ms`; each warm query reuses that snapshot.
All prototypes currently model capacity-gated, multipath, multi-edge, full-source
queries. General selection options and changing eligibility are not benchmarked.
Every query in a batch shares the same mask/residual snapshot; this is not a
benchmark of independent failure scenarios with different masks per query.

| Method | Execution | Intended experiment |
|---|---|---|
| `metal_pull` | One thread per query/vertex, separate dependent dispatches; up to eight rounds per command buffer before checking convergence. | Spread one large query across the GPU; measure repeated launch/synchronization costs. |
| `metal_group` | One threadgroup per query, convergence inside the kernel, global-memory ping-pong distance arrays. | Amortize dispatch across many independent sources without a grid-wide barrier. |
| `metal_shared` | One threadgroup per query, with distance arrays in threadgroup memory. | Reduce memory/barrier overhead for small graphs; enabled only when runtime scratch-memory limits permit. |

Each returns host distance vectors. The measured full result then scans edges
on the CPU and constructs a positive-cost ECMP DAG. **That full timing includes
DAG construction**, but its order differs from the current CPU implementation.
`distance_wall_ms` ends before DAG construction; `gpu_ms` contains command-buffer
GPU timestamps, not end-to-end API latency. Neither should be presented as the
speed of a compatible `Algorithms.spf` replacement.

The pull algorithms perform `O(Q * H * (V + E))` work, where `H` is
the number of relaxation rounds (up to `V`), versus heap Dijkstra's usual
`O((V + E) log V)` per source. This makes hop depth and degree as important as
vertex count. These kernels are feasibility baselines, not an optimized
frontier/delta-stepping implementation or a bound on what Metal can achieve.
Two int64 distance buffers require `16 * V * Q` bytes: 64 MiB at 65,536 vertices
and 64 sources, before host outputs and DAGs. A dense int64 matrix at the same
vertex count would itself require 32 GiB. An all-pairs approach therefore needs
an explicit memory/work analysis and batching on this 36 GB machine.

The CPU controls use the actual current `shortest_paths` / `shortest_paths_to`
sources, built with `-O3 -flto -funroll-loops`. Controls include serial full SPF,
an independent distances-only Dijkstra, and both full SPF and distances-only
Dijkstra through ten reusable worker threads. Pool startup is outside timing;
task submission and result allocation
are inside. Profiling is disabled. The experiment does not exercise Python
binding overhead or an integrated flow-placement workload.

The full suite covers 27 graph/batch configurations: 128–65,536 vertices,
up to 524,288 directed edges, batches of 1–256 queries, nonuniform random costs,
parallel edges/self-loops, leaf–spine fabrics, 32x32 / 100x100 / 200x200 grids,
directed chains, and masked reverse queries with costs above 32 bits. The Clos
and grid shapes were selected after inspecting NetGraph's own benchmark catalogue;
these remain generated fixtures, not a replay of a user's production topology.

**Correctness evidence.** Every timed sample, including warmups, compares every
distance and the canonical `(child, parent, edge-id)` predecessor multiset to the
current C++ result. The independent CPU Dijkstra is also checked against that
result. Timing excludes verification. Twelve seeded small masked multigraphs
exercise all three GPU methods in both normal and sanitizer/Metal validation
runs. Hand-derived witnesses additionally check:

- Diamond ECMP, equal-cost parallel edges, positive self-loop exclusion, and an isolated node.
- Reverse distances and forward-oriented predecessor edges.
- Masked source, masked transit node, and combined residual/edge masks.
- A residual just below `kMinCap` versus exactly at the threshold.
- Int64 costs above 2^53, including a one-unit difference between candidate paths.

The witnesses also deliberately reproduce unsupported simplifications:

- Float32 equates 16,777,216 and 16,777,217 and rounds a residual immediately below
  `kMinCap` up to the threshold. It also erases a `1` versus `1 + 1e-8` capacity tie.
- On a zero-cost `0 -> 1 -> 2 -> 1` graph, distances agree but naive tight-edge
  reconstruction contains a cycle; the CPU DAG correctly retains only two edges.
- Full-suite CSVs record how many query DAGs retain exactly the CPU's array order.
  Canonical equality is not silently reported as byte-for-byte parity.

An initial long run exposed a convergence-flag race in the research kernels:
after the end-of-round barrier, a fast lane could reset the flag for the next
round before every lane had read it. An additional barrier after the reads
fixes that ordering. The rejected run is retained under
[results/rejected-race](results/rejected-race/), and is excluded from conclusions.
After correction, each of two 64-query stress cases passed 100 timed repetitions
plus two warmups for each of three kernels (39,168 complete query comparisons).
An additional 20 timed repetitions plus two warmups per case/kernel passed with
ASan/UBSan and Metal API/shader validation (8,448 query comparisons).
These checks found and then exercised the fix; they are not a proof of every
possible shader interleaving. See [stress.csv](results/stress.csv),
[stress-sanitized.csv](results/stress-sanitized.csv),
[validation log](results/stress-validation.txt), and
[witnesses](results/witnesses-sanitized.txt).

**Performance acceptance remains separate from correctness.** The script uses
fresh processes in CPU / GPU / CPU order on this Mac, nine timed samples per row
after two warmups, and records each raw sample, pipeline startup, upload cost,
CPU load, and GPU utilization. Builds and test suites do not overlap those runs.
The desktop has active background work. These observations must be labelled
exploratory; the repository's quiet-machine requirement is not yet satisfied.
No production speedup or universal CPU/GPU crossover is established by them.
See [A/B/A data](results/aba/) and the generated
[timing table](results/aba/summary.md), including CPU-control drift.

The final run's second CPU-load snapshots were 90.7%, 90.75%, and 73.37% idle
before A1, B, and A2. Some tiny-query CPU controls drifted substantially
(up to 67% between phase medians), so an automatic dispatch threshold cannot
be justified from this run. Selected raw observations below are milliseconds
for distances **plus DAG construction**; CPU batches use ten workers:

| Generated graph | Queries | CPU A1 | Metal B | CPU A2 |
|---|---:|---:|---:|---:|
| Random, 1,024 vertices / 8,192 edges | 1 | 0.145 | 0.707 | 0.115 |
| Random, 65,536 vertices / 524,288 edges | 1 | 20.762 | 4.922 | 22.345 |
| Random, 16,384 vertices / 131,072 edges | 64 | 28.290 | 54.004 | 28.429 |
| Clos, 800 vertices / 320,000 edges | 64 | 17.462 | 72.452 | 17.634 |
| Grid, 40,000 vertices / 159,200 edges | 1 | 2.913 | 8.408 | 2.937 |

These observations identify follow-up candidates rather than establish a
validated speedup. Large low-hop full-source searches deserve a quiet repeat;
small graphs and deep grids/chains do not provide a general reason to offload
with these algorithms. Dropping the DAG changes the experiment materially:
at 16,384 random vertices and 64 queries, CPU distances-only controls were
14.898 / 15.059 ms, with 9.176 ms for the GPU distance stage. Conversely,
the 800-vertex Clos batch was 2.672 / 2.661 ms on the CPU versus 5.681 ms for
the GPU distance stage. Batching alone is not a sufficient selection rule.
Pipeline startup in this repeat was 24.135 ms; the shader cache had already
been exercised, so this is not a cold-install startup measurement.

**Implementation sequence, revised after the CPU experiments.** Follow the
[CPU plan](../../perf/spf_queue/PLAN.md) before selecting a GPU backend. Its
[queue-substitution experiment](../../perf/spf_queue/) showed that an ordered
bucket queue inside the existing search reproduces the native result contract
bit-exactly and runs 1.3–1.9x faster on weighted graphs; BFS and
bounded-integer buckets are its special cases, with a width limit and heap
fallback. It is implemented in the working tree and measured old/new/old at the
Core level and through NetGraph. A batch or distances-only API is conditional
on demonstrated callers.
Reuse reverse SPF and caches before increasing the query count, and measure
representative NetGraph workloads end to end.

Then, if a quiet benchmark justifies it, add an opt-in `Backend.metal()` or an
explicit batch accelerator. Store owned CSR buffers and a backend identity in
the graph handle; keep the graph immutable and the masks versioned or copied per
request. A macOS-only Objective-C++ file and optional CMake Metal/Foundation
linkage are sufficient for the host layer. Keep per-request scratch state safe
under concurrent GIL-released callers. Runtime device/pipeline failures need a
documented policy, and non-Apple builds must retain the CPU-only path.
A first adapter can hold `make_cpu_backend()` and delegate its unimplemented
virtual methods there; it does not require rewriting every algorithm at once.
That is a straightforward way to add a selectable factory, but it must expose
which operations actually execute on the GPU rather than imply that delegating
max-flow or KSP became accelerated.

For the first compatible SPF implementation, retain int64 arithmetic and exact
CPU eligibility; keep zero-cost DAGs, single-path capacity tie-breaking,
fanout-edge semantics, and any unimplemented selection modes on CPU. Preserve
existing predecessor order before declaring parity for the existing result API,
as required by the consolidated CPU plan. A new distances-only batch API avoids most DAG-specific compatibility
work. It must still distinguish supported options from silently ignored ones.

The [follow-up](optimization/README.md) has since tested frontier/bucket methods,
GPU predecessor extraction and guarded 32-bit paths. Production adoption still
requires a compatible result contract and a measured win over the improved CPU.
Directly moving the existing heap-based Dijkstra loop to a GPU is unlikely to
expose useful parallel work. Dense all-pairs tensor formulations change the
memory/work complexity and are not the straightforward default for sparse
network graphs. These are engineering judgments, not measured claims about
unimplemented algorithms.

Keep max-flow, KSP, and placement on the CPU initially. Accelerating their SPF
internals requires additional dispatch plumbing and end-to-end validation of
flow amounts, feasibility, conservation, cuts, costs, masks, and cost factors;
the present distance experiment provides no such validation. Bindings and
`python/netgraph_core/_docs.py` must change together if public APIs are added.
Metal runtime tests need a real Apple GPU CI runner. Local success does not
replace CI.

**Reproduction.** Run from the repository root on a Mac with a Metal device:

```bash
bash dev/research/metal_pathfinding/build.sh
build/metal_pathfinding/capabilities
build/metal_pathfinding/experiment gpu dev/research/metal_pathfinding/sssp.metal witnesses 1
build/metal_pathfinding/experiment gpu dev/research/metal_pathfinding/sssp.metal correctness 2
build/metal_pathfinding/experiment gpu dev/research/metal_pathfinding/sssp.metal stress 100

# Close out other compute work before making performance claims.
# Use a fresh output directory to preserve the recorded run.
bash dev/research/metal_pathfinding/run_aba.sh /tmp/netgraph-metal-aba
venv/bin/python dev/research/metal_pathfinding/summarize.py /tmp/netgraph-metal-aba

bash dev/research/metal_pathfinding/build.sh sanitize
MTL_DEBUG_LAYER=1 MTL_SHADER_VALIDATION=1 ASAN_OPTIONS=detect_leaks=0 \
  build/metal_pathfinding/experiment-sanitize gpu \
  dev/research/metal_pathfinding/sssp.metal stress 20
```

Repository validation also ran `bash .superset/workspace.sh check` (extension
rebuilt; lint/types clean; 170 C++ tests and 465 Python tests passed) and
`make sanitize-test` (170 C++ tests passed). No production C++ source, binding,
or downstream code was modified. No CI run, commit, or push was made.
The [workspace check log](results/workspace-check.txt),
[core sanitizer log](results/sanitize-test.txt), and
[prototype CPU sanitizer results](results/cpu-sanitized.csv) are retained.
