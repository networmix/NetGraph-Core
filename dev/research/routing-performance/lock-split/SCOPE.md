# SPLIT — scope for coordinator reading

Fill six GIL cost rows and P1/8 phase shares, with equal exports. **Stage 1 only:**
stop here as required by `briefs/SPLIT.md`; no build or workload has run.
E = `/Users/networmix/ws/work/2026-11-core-routing-research.noindex`,
T = E/tasks/split, S = E/scratch/split (TASK_SCRATCH was unset). Write only T/S.

| Release points (bindings/python) | Reach; source inventory and saved-call evidence in `stage1-evidence.json` |
| --- | --- |
| native_call.hpp:23/37: SaveThread/RestoreThread in Detached; native_call:50 and detached_argument_work:65 | All workloads. Count every detached argument/computation interval, including several in one public call. |
| flow_policy.cpp:146–149: DetachedDelete | P1/P3 Python-owned policy destruction (P3 saved shortened smoke reaches constructors). Internal C++ evaluator policies do not use this holder. |
| concurrency_hooks.cpp:141/144/148: gil_scoped_release, including call_guard | Private probe wait/release/counts; excluded by wheel's PRIVATE_HOOKS=OFF. Complete production/private caller inventory retained. |
| Saved reached callers | P1/1 and /8: graph/masks, ledger/policy construction, remove, apply_actions, **evaluate_pairs during prepare**. P2/8: graph/masks/evaluate_pairs. Shortened P3 smoke: P1's groups plus FlowPolicy.program. to_units calls do not prove its vector detach branch. Fresh six-row smokes confirm reach; other callers in NOTES. |

| Patch, accounting and checks | Proposed method; details and derived cases in `NOTES.md` / `LEDGER.md` |
| --- | --- |
| Scratch build | Archive Core 9cd7642; modify native_call.hpp/module.cpp only. Pinned offline buildenv, xcrun CC/CXX, target 15.0, four jobs; copy reference venv, install --no-deps --no-index. Preserve finalization behavior. |
| Clock reads | After SaveThread: wall r0, thread CPU c0. Before RestoreThread: c1, r1; after it: r2, c2. Inline thread-local totals: released W=r1−r0/C=c1−c0; retake W=r2−r1/C=c2−c1; count, destructor count, invalid-clock flag. No allocation/I/O in destructor. |
| Export counters | Non-detaching private snapshot getter, outside T4 outer span clocks; store deltas per call, thread/phase entry and exit. TLS shared across translation units; lifetime index avoids recycled ids. Retain T4 getters, constructors and scoped decrefs. |
| Held residual and closure | H_W=span_W−released_W−retake_W; H_C=span_C−released_C−retake_C. Independently integrate outside intervals: outside + three raw parts = thread lifetime, per phase; reconcile in/out-of-span releases with thread totals. Reject missing endpoints, overlaps, invalid counters/material negatives. |
| Meaning and overhead | Held wall includes descheduling; held CPU includes attached bindings/runtime, not only conversion. Retake wall includes attach CPU, so is an upper proxy for GIL wait. Null/full-wrapper and clock-sequence controls (1/8/28 callers) give raw/corrected values, empirical spread and overhead bucket; never subtract summed overhead from elapsed wall. |
| Core workers | Keep T4 phase process-minus-Python CPU brackets separately. Never sum overlapping per-call process deltas. P2 evaluate_pairs worker CPU is not caller CPU. |
| Exports and tool controls | Full matching plain export per row; only instrumented extension hash may differ, checked against actual binaries. Capture permits no exception. Controls cover multiple/zero detaches, destructor, exception, thread exit, native worker, poisoned export/counter and blocked-stack attribution. |
| Size / alternative | About 60–100 C++ and 150–250 incremental Python lines reusing T4. Sampling alone cannot resolve microsecond call edges. No product change or gain claim; substantial growth returns for rescope. |

| P1/8 capture | Method |
| --- | --- |
| Command | Reference wheel: `sample CHILD_PID 15 1 -file S/p1-8-sample-01/sample.txt`; ready/GO handshake, complete workflows for ≥17 s. All thread roots; export each retained result after capture and compare to plain. |
| Categories | ORIENT precedence: validated static_support search/support branches; support/program construction; placement; ledger/results; setup/masks/cleanup; unresolved. Explicit blocked leaves win. Revalidate offsets under reference wheel before use. |
| Denominator | Publish all-root exclusive counts, blocked waits inside/outside Core, unresolved and non-Core counts; conserve roots and top-of-stack counts. New phase shares exclude recognized blocked leaves. Beside carried P1/1 (9,134 Core samples, different wheel/denominator); no CPU-share or scaling inference. |

| Proposed coordinator job (heavy slot, one job at a time) | Minutes / timeout |
| --- | --- |
| `measure.py --workload dcbb96 --workers 1` and `--workers 8` | 1–3 and 1–2 / 600 s each |
| `measure.py --workload nsfnet --workers 1` and `--workers 8` | 1–3 each / 600 s each |
| `measure.py --workload backbone --workers 1` and `--workers 28` | 1–2 each / 600 s each |
| `capture.py --workload dcbb96 --workers 8 --seconds 15` | 1–2 / 180 s |

Full command prefix/output paths are proposed in NOTES; stage 2 implements and short-run
checks every command, then stops with the exact handoff. Build 2–6 min and controls 2–5 min,
separate bounded slot jobs. Stage 3 writes SPLIT.md: run/all-phase totals, five largest call
kinds by run-span wall, closure, overhead and large/small assessment. All rows exploratory:
one instrumented process/capture each, no A/B/A, significance or speedup claim.
