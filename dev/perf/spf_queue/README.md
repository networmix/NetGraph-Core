# Ordered bucket queue experiment (CPU only)

Plain C++ harness (no Objective-C++ or Metal) behind the
[CPU optimisation plan](PLAN.md). It compiles the production
`src/shortest_paths.cpp` next to a copy of its forward search loop in which only
the priority queue and the scratch workspace are pluggable, then checks every
variant bit-exactly against the production function and times them.

- `spf_variants.hpp`: `HeapQueue` (production key), `BucketHeapQueueT` (monotone
  cyclic buckets, active bucket kept as a heap on `(-bottleneck, node)`; the
  `<true>` instantiation masks instead of taking a modulo), `BucketSortQueue`
  (sorted active bucket plus a side heap), `Workspace` (fresh, full reset,
  sparse reset, hybrid reset) and `spf_variant`, the shared loop.
- `harness.cpp`: workload generators (the research generator verbatim, plus
  uniform, wide, zero-cost-cycle and pseudo source/sink cases), the parity
  matrix, timing and the cost-range sweep. Variant names are
  `{ref, heap, bucket, bucket2, bsort}-{fresh, full, sparse, hybrid}`.
- `run_aba.sh`: one fresh process per variant, production reference first and
  last, plus the sweep. `summarize.py` prints Markdown tables.
- `results/first/`: the first sequence (nine samples), sweep files, load readings,
  source hashes and the parity outputs (plain, sanitized and the teeth check
  built with `-DSPFX_BAD_ORDER`, which must fail). `results/short/`: the second
  sequence for the hybrid and power-of-two variants.

```sh
bash dev/perf/spf_queue/build.sh
build/perf/cpu-queue correctness
build/perf/cpu-queue time bucket-hybrid 9
build/perf/cpu-queue sweep bucket-sparse 9 2097152
bash dev/perf/spf_queue/run_aba.sh
python3 dev/perf/spf_queue/summarize.py
```

Timings are sequential and single-threaded; they are not comparable row for row
with the ten-worker batch timings in the research reports.

## Production verification (after the change landed in `src/shortest_paths.cpp`)

- `bitgate/dump.cpp`: exercises forward/reverse SPF, max-flow, batch max-flow,
  sensitivity, `FlowPolicy` sequences and KSP through the public C++ API; built
  once against the previous sources (`git archive` of the base commit) and once
  against the working tree. `bitgate/sha256.txt` records the identical digests
  of `old.txt` and `new.txt` (and of the new build forced onto the heap route).
- `run_core_aba.sh` + `results/core/` (+ `summarize_core.py`): old/new/old on
  the harness workloads, where `ref` is the production function of the sources
  each binary was built from (`cpu-queue-old` is built from a `git archive` of
  the base commit); `b-new-heap` forces the heap route in the new binary as an
  in-binary control.
- `bench_netgraph.py` + `run_netgraph_aba.sh` + `results/netgraph/`: old/new/old
  on NetGraph workloads run through NetGraph's own API and interpreter (full
  SPF, bound-context max-flow with pseudo nodes, failure exclusions, demand
  placement with WCMP and ECMP policies, failure Monte Carlo, KSP, and three of
  the shipped scenarios end to end with `parallelism: 1`). Each workload
  returns a checksum that must agree across builds; `summarize_netgraph.py`
  tabulates medians and checks the checksums. The driver pins
  `PYTHONHASHSEED=0`: NetGraph's Monte Carlo digests otherwise differ between
  two processes of the same build (`results/netgraph/hashseed-determinism.txt`).
- `review-2026-10/`: brief, full report, patches, experiment sources and
  raw results of the independent review (2 October 2026); what was adopted
  from it is listed in the plan. Its full-source experimental copies, session
  log and 10 MB bit-gate dumps were not kept (they are reproducible from the
  patches against the base commit). `bash build.sh old` builds the
  old-source harness binary from `git archive` of the base commit.
- Timing modes now include reverse-search rows (`dst` = `reverse`), which are
  production-only (the harness variants implement the forward loop).
