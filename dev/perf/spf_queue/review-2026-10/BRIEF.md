# Critical review brief: SPF optimisation in NetGraph-Core

You are an independent, sceptical performance reviewer. Repository root is the
current directory (NetGraph-Core, C++20 core + pybind11). The working tree has
an UNCOMMITTED change on top of commit `fbab119` that replaces the Dijkstra
frontier in `src/shortest_paths.cpp` with a monotone bucket queue (active bucket
kept as a heap on `(-bottleneck, node)` so pop order is identical to the old
binary heap), plus a thread-local scratch workspace with touched-node tracking
and a hybrid reset. Read, in this order:

1. `dev/perf/spf_queue/PLAN.md` (the plan: decision,
   experiment, implementation status, old/new/old measurements, trade-offs).
2. `dev/perf/spf_queue/review-2026-10/change.diff`
   (the exact production + test diff) and the live files `src/shortest_paths.cpp`,
   `src/strict_multidigraph.cpp`, `include/netgraph/core/shortest_paths.hpp`.
3. `dev/perf/spf_queue/README.md` and the harness
   `spf_variants.hpp` / `harness.cpp` (plain C++; `bash dev/perf/spf_queue/build.sh`
   builds `build/perf/cpu-queue`; `build/perf/cpu-queue-old`
   is the same harness built from the previous sources; `time ref 7`,
   `short ref 7`, `sweep ref 7 2097152`, `correctness` are the modes).
4. Callers: `src/max_flow.cpp`, `src/flow_state.cpp` (place_max_flow),
   `src/flow_policy.cpp` (get_path_bundle, memo), `src/k_shortest_paths.cpp`
   (its own `dijkstra_single` and spur searches), `src/cpu_backend.cpp`.
5. Context memos: `dev/research/metal_pathfinding/optimization/README.md` (GPU research)
   and the measured dead ends listed in CPU_PLAN.md ("do not retry").

## Your task

Think hard and answer: **were all easily achievable and reasonable
optimisations considered and verified?** Specifically:

- Cheap wins inside the new code that were missed: the bucket queue (bucket
  advance, bitmap, per-pop overhead, item layout, the chain-graph regression,
  width limit), the heap route (the recorded 5–9% fallback regression and its
  real cause), the workspace (reset policy, output materialisation, the O(N)
  output floor), predecessor-list construction and the PredDAG conversion,
  early-exit handling, the reverse search, the cost/gcd metadata.
- Cheap wins around the search that were not touched: how max-flow,
  place_max_flow and FlowPolicy call SPF (repeated full searches, residual
  recomputation, memoisation, early termination), KSP's private
  `dijkstra_single`, graph construction (CSR layout, edge ordering by cost),
  the Python binding path (copies, GIL), batch/thread usage.
- Whether the measurements support the claims: A/B/A design, noise, the
  `W_max = 65536` choice, the sweep, the NetGraph-level workloads
  (`cpu_queue/bench_netgraph.py`) and whether they represent realistic
  NetGraph use; anything the author may have rationalised away.
- Correctness or robustness risks that an optimisation reviewer should flag:
  thread_local lifetime and memory bound, exception paths, reentrancy,
  int64 overflow in bucket arithmetic, graphs with zero edges, Python threads.

## Rules

- Be concrete and sceptical. Cite file paths and line numbers. Do not restate
  the plan; find what it misses or gets wrong.
- Verify, don't speculate: build the harness and run experiments where a claim
  can be checked in minutes (old binary vs new binary, forced heap via
  `NGRAPH_CORE_SPF_QUEUE=heap`, fresh processes, several samples, record load
  with `uptime`). Treat differences under ~5% as noise unless reproduced.
- Do NOT modify tracked files. Put every experiment, patch (as `.patch` files)
  or script under `dev/perf/spf_queue/review-2026-10/`.
  You may build into `build/`. Do not run `git commit`, `git stash`, `git
  checkout`, or anything that changes the working tree outside that directory.
  The Python venv is `venv/`; NetGraph lives in `/Users/networmix/ws/NetGraph`
  (read-only for you). Do not run the full test suites; they already passed.
- Other sessions share this machine; keep runs short and sequential.

## Report (write it to `dev/perf/spf_queue/review-2026-10/REPORT.md`, in English)

1. A ranked list of findings. For each: title; what exactly is missed or wrong
   (with evidence); expected gain and the workload class it applies to; cost
   and risk; verdict as one of VERIFIED-GAIN (you measured it), PLAUSIBLE-
   UNTESTED (why you could not verify), REFUTED (you measured no gain),
   ALREADY-COVERED (the plan handles it; say where).
2. For every VERIFIED-GAIN: the exact patch (file + diff) and the numbers
   (old / new / your variant, samples, load).
3. Measurement critique: what in the recorded A/B/A would you not trust, and
   what would you rerun.
4. A short verdict: is the shipped approach the right streamlined choice, and
   what (if anything) should be done before committing it.
