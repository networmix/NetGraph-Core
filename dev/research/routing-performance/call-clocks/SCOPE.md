# T4-CLOCKS: scope for coordinator reading

Stage-1 proposal, accepted with the overrides in [RULINGS.md](RULINGS.md).
The records below preserve the original proposal; RULINGS governs implementation.
Stage-2 status and deviations are in NOTES.md and HANDOFF.md.
The purpose is to distinguish caller CPU, other Python-thread CPU, native-worker
CPU and waiting at several workers, then measure response to added latency.
All proposed measurements are exploratory unless every requirement of RESEARCH §3
and PERFORMANCE §3 is met. No finite speedup or free-Core endpoint follows from a
positive-delay experiment.

All paths below are relative to
`/Users/networmix/ws/work/2026-11-core-routing-research.noindex/`.
Writes are confined to `tasks/clocks/` and `scratch/clocks/`.

## Reference and reuse

| Input or mechanism | Plan and boundary |
| --- | --- |
| Frozen design | RESEARCH SHA-256 `13b584437fecbc0e66c988187a2b9288dbb5f08f76f6070423be53d5ec7c7155`; STATE `a5057a51b0f94c1210f7e4aceceeb6827e1f899f1db5cc4332bea041e8e75d73`; both freshly checked. PLAN r1 T4 and brief T4-CLOCKS govern. |
| Runtime | `scratch/consumer/ref/venv/bin/python`; NetGraph source `worktrees/netgraph-orient` at `8a641a1`; Core extension must hash to `208487c7aed158fb418791c22e23d09332533595053c0e4810923b31244f7c84`. Assert paths, revision, ABI, GIL and hash in each child, including captures. |
| Workloads | Unchanged `scratch/orient-add/workloads/{dcbb96,nsfnet,backbone}.yml`, checked against their manifest. P1/P2 workers 8; P3 8 and explicit 28. MSD/preparation's explicitly serial calls stay serial. Count calls by phase and actual operation, not scenario workers alone. |
| Frozen T1 | Reuse/adapt `scratch/clocks/t1-recorder-r1/events.py:public_api` and its install rules, with `values.py:FIELDS` as the field inventory. It includes enum properties which ORIENT-ADD skips. Do not use T1's full-input serialization, weakref notification as a destructor timer, or its evolving live tree. Verify the supplied SHA256SUMS before copying. |
| Lifetime timing | Adapt ORIENT-ADD's `span_timers.py` and `lifetime_tracker.py` plus lexical decorations from `consumer-lifetime.patch`, in an owned consumer copy. A weakref callback observes collection; it does not surround native destruction. Holding/release scopes do, with stated lifetime perturbation and bookkeeping bias. Extend scopes to the P2 controller path and export/cleanup after inspecting their lifetimes. |
| Existing answer | STATE §3 and `scratch/orient-add/{analysis-03,p3-analysis-01}.json` give exploratory caller clocks for P1/P3 and one-worker P2 getters/profiles, not P2 native-worker CPU or a latency response. No new release-behavior probe is claimed at this scope stop. |
| Profiles | Carry the one-worker summaries from `scratch/orient-add/native-{nsfnet,backbone}-analysis-03/`. Their extension hash is `6db37900d83800b0b5c227104889f47259c5761dd640c1d0d8684fac61009ff1`, different from the pinned reference. Their PC ranges cannot be used for new captures without new same-binary witnesses. |

## Clock coverage and accounting

| Quantity | Where read; interpretation |
| --- | --- |
| Caller wall / CPU | Integer `perf_counter_ns` / `thread_time_ns` immediately around each original public function, method, constructor body and property read, including enum getters; `finally` closes raising calls. Also surround real last-owner release scopes. `__new__` shell allocation is outside constructor-body spans. |
| Process CPU | `process_time_ns` at the same call/release endpoints and phase endpoints. It includes every native worker, including workers which terminate before the call returns. Per-call process intervals overlap across callers and must never be summed. |
| Every Python thread | Register main and every `threading.Thread` lifetime, using lifetime identities rather than recyclable thread ids; intercept bootstrap, not only executor case functions. Each records its own start/end `thread_time_ns`. At every call endpoint read cumulative CPU of all registered live Python threads through Darwin `pthread_mach_thread_np` + `thread_info(THREAD_BASIC_INFO)`; archive exited threads' final counter. Use `ctypes.PyDLL` so clock inspection itself does not release the GIL. BASIC_INFO's user/system seconds and microseconds give a known conversion; validate granularity/lag against self clocks. SDK header availability is established; runtime access and precision are not yet established. |
| Thread births/exits | Initial census before collection; register new threads before their target, retain completed counters, take final counters before unregistering. Phase cuts use cumulative snapshots, not whole lifetimes extending beyond a phase. Record startup/teardown tails and mixed-clock discrepancy. Unexpected/unregistered Python threads or failed Mach reads invalidate residual attribution; do not substitute zero. |
| Core-worker CPU | For a common interval, `R = Δprocess_CPU − Σ ΔPython_thread_CPU`. First report this as **non-Python-thread CPU residual**, with clock/skew uncertainty. Attribute it to Core workers only with the thread census, native-stack evidence and null controls excluding other native work. P2's `src/evaluators.cpp:150-185` launches/join workers; it is not eight Python batch callers. P1 is a useful residual control, not an assumed zero oracle. |
| Per-call residual | Store Python vector and process endpoints for every span, with process/wall brackets around each vector read. Each residual is activity of the whole process during that interval, not exclusive CPU owned by that caller. Primary totals use phase differences once; any Core-only residual uses a partition of the global span union, never overlapping call sums. |
| Wall union | Merge intervals per thread for its outside budget, and globally for elapsed Core-boundary coverage. Report global union, per-thread sums and analysis/pipeline elapsed separately. Internal native work may overlap another caller's outside Python work. |
| Outside / not-running | Per thread and phase: `outside_wall = lifetime_wall − span_union_wall`; `outside_CPU = lifetime_CPU − exclusive_span_CPU`. Report `outside_wall − outside_CPU` and total `wall − CPU`, plus inside-span counterparts. These include queues/futures, native joins, GIL reacquisition and OS scheduling. GIL waiting alone stays unknown. |
| Proxy bias | Outside CPU is active caller work outside wrapped Core boundaries, not a demonstrated serial floor: binding work holding the GIL is inside spans; Python/runtime, retention and timer work are outside; some outside native work can release the GIL. Destructor scopes include container/decref costs and delay release. Live results/owners survive until their real export/cleanup phase. |

Snapshots are costly on getter-heavy P2. Store integer endpoint counters and event
rows in memory, write/compress only after the measured phase. No trace payload,
result canonicalization or evidence hashing occurs inside a timed phase. Report
memory and instrumented/plain wall separately; diagnostic wall is not the user baseline.

## Timer correction and controls (stage 2, not yet run)

| Control | Rule / failure outcome |
| --- | --- |
| Null boundary | Five matched direct/wrapped blocks with the complete clock/vector/retention stack; 30,000 calls each and a short multi-thread version at the workload's registry size. Keep distributions, cost inside spans and full wrapper cost separately for wall/CPU. Read clock resolution. |
| Self vs foreign clock | Main plus short-lived CPU-bound and sleeping Python threads; compare Mach deltas to each thread's self clock, including birth/exit/reused identities. An active noncaller must appear in the Python sum; a sleeper must not contribute wall as CPU. Bound skew/lag empirically, not merely from advertised resolution. |
| CPU residual | Null/sleep and a joined native-worker control from the scratch binding: CPU in an unregistered native worker must remain in residual after that worker exits. Busy noncaller Python CPU must be subtracted. Failure leaves Core-worker CPU unknown, not an invented estimate. |
| Correction | For comparable call groups, subtract median null-span CPU from span CPU; subtract `N × (full_wrapper_CPU − null_span_CPU)` from outside CPU. Keep raw and corrected values, repeat spread and snapshot uncertainty. Retention/decref control is separate. No subtraction of summed multi-thread overhead from elapsed wall or global union. No clipping negative residuals/corrected values to zero. |
| Closure | Per thread, independent complement integration of sorted event intervals must equal lifetime wall; CPU spans plus gaps must equal lifetime CPU. Nested events use exclusive accounting. Compare summed self CPU plus residual to process CPU with the measured bound, including phase-edge tails. A mismatch refuses the report row. |
| Near misses | Scratch controls poison an export field/type/order; omit an event/thread-final CPU; double-count overlapping process spans; misread a blocking stack as solver work; reuse PC ranges under the wrong extension hash. Each must fail its relevant checker. Do not add broad product tests for an experiment-only helper. |

## Sensitivity proposal and decisions needed

Use one scratch wheel from a clean export of `worktrees/core-ref`, adding only a
small private `_clocks_delay(ns)` binding. Its body detaches the GIL and sleeps
for the requested nanoseconds; δ=0 still executes the same detach/attach path.
One thin wrapper invokes it immediately before each saved original public call,
once at the outermost intercepted boundary. Constructors, methods and all getters
receive it. No sleep is inserted in each `native_call`: a public call can use
several detach scopes and many getters have none. No sleep is added per pair
inside `evaluate_pairs`; one batch remains one public call.

This inserts another detach/attach transition before the original call. The
δ=0 scratch control measures its disturbance relative to plain reference; the
δ response is conditional on this shim. Measure actual helper elapsed wall on
separate controls: requested δ is not actual scheduler latency. It does not
consume δ of CPU. A progress control must show a Python thread running while
another thread sleeps in the helper. Do not claim negative-delay behavior.

**C1 — perturbation boundary needs coordinator acceptance.** Proposed scope is
all public calls, in load/prepare/run/export, with separate analysis and pipeline
responses. Native destructors remain covered by clocks but do not get an extra
delay: a scoped release is a group of decrefs, not one destructor call. Final
cleanup is reported separately. If E5 intends a δ for each actual native
destructor too, the helper is insufficient; approve a larger binding-deleter
patch first. Do not label a group-release delay “per destructor”.

| Run design | Predeclared proposal |
| --- | --- |
| Comparisons | P1/P2 at 8; A is the same scratch wheel with δ=0, B δ=0.1 or 0.2 ms. Plain reference vs scratch δ=0 is a separate perturbation/export control. Timing lane has no per-call clock/vector/retention recording; fixed simple call counts verify exposure. Clock lane is separate. |
| Repeats/order | Seven A/B/A triples per nonzero δ per workload, fresh processes. Alternate which δ's triple comes first across repeat indices (0.1 first at indices 0,2,4,6). Each B has its own A1/A2; never pool a distant zero run as reference. All exports retained and checked outside the timed path. |
| Block lengths | At least 2 s per side. Fresh CLI processes back to back inside a block; per-side iteration counts calibrated on the zero/control rows, fixed before nonzero measurement. Allow B=1 when a single delayed process is already long. Keep analysis wall, whole process wall, call counts and RSS per process. No shortened workload or truncated getters. |
| Statistic | `a_i = mean(A1_i,A2_i)` using per-process block medians; `r_i = B_i/a_i`; report median ratio and 95% percentile bootstrap interval over whole triples, 10,000 resamples, seed 20261003 (Core `aba.py`). Also report secant slope `b_i=(B_i−a_i)/δ_ms`, median and whole-triple bootstrap interval in seconds per ms added per call. Separate slopes at 0.1 and 0.2; do not force a line through both. |
| Detectability | Before B runs, seven qualified zero/zero/zero triples at each workload's own scale. Let `e=max(abs(CI_low−1),abs(CI_high−1))`; report empirical response floor `e × median(A)` seconds and `e × median(A)/δ_ms` slope units, beside the 3% beyond-noise rule. This is a calibration resolution, not an 80%-power guarantee. A numeric minimum detectable response is currently unknown. Fix counts/drift rule from calibration; unresolved A/A, >3% reference drift or clusters invalidate a qualified response. No retry-until-significant. |
| Interpretation | Added latency response only. A slope interval including zero is inconclusive at that exposure/resolution. Neither a zero nor a positive slope establishes the free-Core endpoint or any finite gain from an optimization. Report boundary/GIL shim bias and conditional endpoint assumptions beside it. |

**C2 — requested exposure conflicts with a short qualified heavy window.** Frozen
T1's scope reports 1,930,103 P2 public calls across all phases. At serial Python
call sites δ=0.2 ms requests about 386 s of sleep per full process; seven B runs
at each δ request about 4,053 s (67.6 min) before scheduler oversleep, plain A
runs, calibration or export checking. P2 run alone also has hundreds of thousands
of getters, so excluding export does not automatically solve this.

The brief caps each heavy job at approximately ten minutes and asks resumable
parts, whereas PERFORMANCE §3 requires the qualified series inside one heavy
hold. Do not silently replace all-call exposure by only `evaluate_pairs`, omit
getters, or call a multi-window series RESEARCH §3-qualified. Recommendation:
retain full exposure, run one bounded delayed child per heavy part, and report
the seven-triple slopes/intervals as **exploratory across parts**. Coordinator
must accept that qualification or authorize a longer reserved hold. If choosing
the latter, plan at least 90–120 min for P2 alone and revise the job/time budget;
measure oversleep in a pilot before committing the window. No jobs queued yet.

## E6: new captures and carried one-worker evidence

| Capture / analysis | Method and boundaries |
| --- | --- |
| New rows | P2/8 and P3/8,28, one 15 s `sample` capture each at 1 ms, 17 s or more of complete repeated workflows after readiness. Use reference wheel, no clock/delay wrappers; export after capture. Child timeout 120 s, whole job 150 s. `--no-sample` short control precedes handoff; only coordinator attaches sampler. Explicit worker override must handle 8 and 28, not ORIENT-ADD's current special case for 1. |
| Every thread | Adapt ORIENT-ADD parser to all roots; remove its one-main-thread assertion. Preserve thread roots and exclusive self weights (`node inclusive − immediate children`) and conserve totals per root. Independently cross-check the sample's top-of-stack totals. C++ evaluator worker ancestry differs from Python caller ancestry: recognize `run_tasks`/worker paths, not only `evaluate_pairs` on the main stack. |
| Split | Solver residual searches; potentials; blocking flow; closure/cut; solver setup/output/other; roster/support searches; ledger/policy fixed costs; result/record/program construction; binding conversion/publication; unresolved. P3's corresponding fixed setup is per demand/case, not per all-pairs evaluation. A sampled share is not removable time. |
| Waits | Separate explicit pthread/futex/condition/join/GIL waits and readiness wait, even beneath a Core ancestor. Unknown scheduler/runtime stacks stay unresolved. Report raw all-thread stack counts and identified active-native-context shares; `sample` is not CPU-weighted thread accounting and blocked stacks must not inflate the solver denominator. No scaling of an all-thread share by caller wall into microseconds/evaluation. |
| Binary attribution | Generate disassembly of the actual pinned reference and new source/call-site witnesses for inlined phases and stripped worker paths. Reject hash mismatch. Recognized named frames may be used; unresolved PCs stay unresolved until a separating witness exists. Do not port old numeric offsets by name. |
| One-worker comparison | Carry ORIENT-ADD's exclusive counts/shares with original denominator, extension hash and exploratory label. Show new shares beside them with different build/worker/denominator caveat. Existing main-thread one-worker samples did not measure a multi-worker CPU share. |
| P2 getters | Independent reference P2/8 run: full getter clocks by field during export; plain and instrumented exports use identical sinks. Include enum fields from frozen T1. Report count, raw/corrected CPU and wall, null overhead, diagnostic export wall and plain wall; run workers=8 does not imply eight exporter threads. |

## Equality, build and handoff

| Item | Required check / cost |
| --- | --- |
| Full exports | Fresh plain reference controls at matching workers, same scenario and sink. Strict recursive scalar/type/key/list-order comparison for every clock, sample, scratch δ and null-control workload export. No `metadata.workers` exception when counts match. Only scratch build identity `metadata.core_extension_sha256` may differ: preserve/check its presence/type and verify its actual expected hash separately. Any additional producer-identity field requires an explicit named rule first. Comparator must reject poisoned value, type, missing/extra key and list order. Compare after timing, not in the timed phase. |
| Scratch build | Offline clean source export; pinned `scratch/wheels/buildenv`, `--no-build-isolation`, xcrun CC/CXX, deployment target 15.0, four jobs. Owned venv copy under `scratch/clocks/`, install `--no-deps --no-index`. Record compiler/SDK/flags, patch and extension hashes; assert imports inside copied venv. Preserve symlinks and check no install target resolves into the reference. No product edits, cache migration, transition rules or roadmap changes. |
| Expected size | About 450–650 new/adapted Python lines across collector/driver/analyzer/controls plus a small delay/native-worker-control binding patch; unchanged parser/attribution helpers carried separately. Four stage-1 records and two short reports. Smallest alternative (old main-thread timers + one `Detached` sleep) misses native CPU, enum/getter coverage and per-public-call exposure. Foreign-thread snapshots and all-root attribution address these concrete omissions. Stop and rescope if this grows several times the estimate. |
| Stage 1 | Read-only inspection and records; zero heavy minutes. Stage 2 awaits coordinator word. |
| Stage 2 | Controls/smokes roughly 5–10 heavy minutes; clean build/venv roughly 5–10 more, each job ≤600 s. Short checks bounded and untimed; lint/format/AST with pinned consumer tools. Native helper checks proportionate to scratch-only delay/worker behavior, no full product release gates. Then stop with exact measured commands. |
| Stage 3 clocks/captures | Coordinator unconfined: clocks/getters and controls about 5–10 min; three captures about 2–5 min including checking. Analyze already-recorded data locally. P2 clock overhead may require separate bounded phases; never truncate calls to make a budget. |
| Stage 3 sensitivity | Full P2 exposure ≥67.6 min requested sleep for nonzero B runs alone; provisionally 90–120+ min including oversleep/calibration. P1 cost depends on phase counts and overlap, provisional 20–40 min. C2 must resolve hold/qualification and budget before scheduling. At most one waiting/running job from this task. |
| Coordinator needed | Accept clock/coverage plan and C1's δ boundary; resolve C2's window/qualification; later run measured commands and samples outside confinement. Runtime denial or unreliable foreign-thread clocks returns here as an explicit blocker, never a permission workaround. |
| Next stop | After script/build/smoke validation, provide exact bounded commands and new output paths; no measured series starts before that handoff. After reports coordinator freezes and arranges the specified fresh Opus reader; this worker starts no subagents. |
