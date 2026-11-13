# T1-REPLAY: scope for coordinator reading

Status: stage 1 authorized by [RULINGS.md](RULINGS.md). Untimed inventories are complete.
The coordinator's eight rulings override the original proposals below; implementation now proceeds to its planted-mismatch stop.

The tool serves the researchers comparing Core builds: a wrong cut, edge allocation or ledger
record must fail before a speed result can influence a decision. It records the pinned consumer's
Core boundary once, replays the actual calls and lifetimes, and uses a separate exhaustive-output
pass for differential checking. No product file or shared environment changes.

Paths below are relative to `/Users/networmix/ws/work/2026-11-core-routing-research.noindex/`
unless written absolute. Task files live in `tasks/replay/`; generated evidence in `scratch/replay/`.

## Inputs and scope boundary

| Input | Boundary and evidence |
| --- | --- |
| Reference | Wheel `scratch/wheels/ref/netgraph_core-0.13.0-cp314-cp314-macosx_15_0_arm64.whl`; SHA-256 verified as `3a79bfcd7de29612ac5126fd22fd7080a88ef5466202c827041b3e4ca8bb8550`. Use `scratch/consumer/ref/venv/bin/python`, never the editable build for timing. |
| Consumer | `worktrees/netgraph-orient`, HEAD verified `8a641a1caa0b8f4b7d8c4df0aecf8625310cdf47`. Reference interpreter imports it successfully by explicit source path. |
| Contract | `worktrees/core-ref`, HEAD verified `9cd7642b5ea8681f7833c139969d8ed224d03a75`; `docs/semantics.md`, public stub and bindings. Private hooks are excluded. |
| Design | `design/RESEARCH.md` SHA-256 verified `13b584437fecbc0e66c988187a2b9288dbb5f08f76f6070423be53d5ec7c7155`; accepted `plan/PLAN.md` T1. |
| P1 | `scratch/orient-add/workloads/dcbb96.yml`, SHA-256 `cd81a0f81d8227ffc2d688ca24df287410f991500023ecde83fbca301fa7a8c3`; 2460 cases, MSD plus placement, 1/8 Python workers. |
| P2 | `scratch/orient-add/workloads/nsfnet.yml`, SHA-256 `64f19e3a1577f9fcaae83f9118722129de49e6fe72715509328e5d72c41a6f04`; 116 run batches / 87,696 pairs, 1/8 native workers. Preparation batches are additional calls. Export getters belong to the replay. |
| P3 | `scratch/orient-add/workloads/backbone.yml`, SHA-256 `5fbadce0c78f6d9d2fc0117ed2578c50be0dafc7d9a4dd7b130d8f1304c598a7`; 592 distinct cases from 1000 draws, MSD plus placement, 1/28 Python workers. |
| Existing release answer | `scratch/replay/tmp/scope-controls-01.json`: reference wheel returns ECMP delivery 5, per-edge `[3,2,3,2]`; controller delivery 12, cost 24, cut `[0,1]`, curve `[(2,12)]` on the unequal-capacity diamond in LEDGER. GIL is enabled. This probes the chosen release, not just source. |

The workload call inventory is below. Counts distinguish import/load/prepare/run/export/cleanup.
A call-count trace contains no instrumentation-generated Core reads. Python wrapper collection is
counted separately from native destruction: evaluator-internal ledgers, policies and results'
native shared owners are not Python calls.

## Inventory

Fresh reference-wheel traces at the pinned worker counts are complete. [INVENTORY.md](INVENTORY.md)
contains every observed call kind by phase and every observed wrapper lifetime count. Raw files:
`scratch/replay/tmp/scope-p{1,2,3}-02.json`, matching bounded logs and exports; validation:
`scope-validation-01.json` and `scope-validation-02.json` in the same directory.

| Workload | Total calls / kinds | Main stateful calls and scope differences |
| --- | ---: | --- |
| P1 | 215,859 / 52 | 4922 policy constructors, 2461 ledgers, 2475 action batches, 9872 removals; one empty preparation evaluator batch; import defaults, 91,716 prepare conversions and all export getters included. |
| P2 | 1,930,103 / 42 | 119 evaluator batches: 3 preparation pairs plus 87,696 run pairs in 116 batches. Export result/options getters and enum names are included; run also reads routing 175,514 times and each of CentralizedTE source-allocation/support/tunnel-count 87,699 times. |
| P3 | 114,618 / 50 | 3558 policy constructors, 593 ledgers, 599 action batches, 7152 removals; one empty preparation batch; 3552 retained program reads and cost-curve reads included. |

Both count passes have identical non-enum counts; the second adds public enum property reads.
Their full exports are exact by JSON type/key/list order/scalar value. Against the saved consumer
exports under PyPI, only `$.metadata.core_extension_sha256` differs. This is a recording probe
control, not full-output replay equivalence. Export metadata confirms run workers 8/8/28 and
serial MSD. Native reference SHA-256 is
`208487c7aed158fb418791c22e23d09332533595053c0e4810923b31244f7c84`;
PyPI native SHA-256 is `6db37900d83800b0b5c227104889f47259c5761dd640c1d0d8684fac61009ff1`.

### Inputs and lifetimes carried by those calls

| Call family | Full input payload and retained objects |
| --- | --- |
| `to_units`, `as_float`, `Graph` | Physical scalars/arrays and exact resolution; node count; original input-order source/destination/capacity/cost arrays. Graph views keep their owner alive. Preserve the EdgeId permutation. |
| Models/config/options | Constructor arguments, including omitted versus explicit `None`, enum type/name/value, every boolean, weights, exact fraction ceilings, supports, outputs, evaluator offer/design. Models and query options may outlive the workflow and are shared between policies. |
| `Graph.eligibility`, `Eligibility.intersect` | Actual excluded node/edge lists, graph reference, both operands and returned alias relationship. Equal exclusion payloads do not imply one mask identity. Baseline mask and graph span all cases. |
| `FlowLedger` | Graph reference and both baselines, including `None`. One private ledger per case; MSD retains one roster across probes. Preserve remove/replace ordering and monotonically allocated local policy/record ids. |
| `FlowPolicy`, `remove`, policy getters | Ledger, normalized or original endpoints as actually passed, config, eligibility, cache; later getters/removals refer to the same policy. Destruction may leave orphan records and is not `remove()`. |
| `Replace`, `Forward`, `Place`, `apply_actions` | Ordered action list; actual policy refs, offers, designs, black holes, mode, outputs. Actions keep policies alive. Input containers and result-list elements have distinct reference lifetimes. Only observed action kinds are timed. |
| `evaluate_pairs` | Graph, ordered complete endpoint pairs, common options, per-pair masks and workers. Preparation's workers=1 stays serial even on a scenario-worker row. Results survive until their actual last consumer read/export. |
| Result/program/record getters | Object references and method arguments, not reconstructed totals. Programs are replayed through their producer calls: the public program lacks a serialization of split tables (`semantics.md`, Evaluators / Limits). Never rebuild a program from its visible tunnels alone. |
| Array/owner release | Record last observed Python references, array-owner dependencies and call/job position with weakrefs, without a recorder-owned strong registry. Replay drops its own references at those events. Native internal destruction remains inside the producing call or final owner release. |

No public cache constructor appears in the fresh P1–P3 traces, including consumer import, preparation, run and export.
The generic check reader still covers a supplied cache on a small control; hidden evaluator-private
caches have no public statistics to read. P1–P3 do not constitute coverage of all Core models,
searches, errors or cache concurrency; later prototype-specific oracles remain required.

## Concurrency and timed scope

| Row | Scheduling and ordering |
| --- | --- |
| P1 / 1, 8 | Record logical failure-case jobs at `schedule_cases`, with step barriers and case input order; MSD is a serial job with its complete persistent roster. Replay the same jobs in a `ThreadPoolExecutor`, not one future per individual Core call. |
| P2 / 1, 8 | Python invokes each batch in sequence; explicit `evaluate_pairs(workers=1 or 8)` launches Core's own workers. Never put these batches in eight Python workers. Preserve pair order; prepare batches retain their recorded explicit worker count. |
| P3 / 1, 28 | Same Python case scheduler as P1. Resolve `auto` to 28 on this machine and record both configured and effective counts; reject a silent host-dependent change. MSD remains serial. |

For every row, a per-call span begins immediately before calling the saved original bound function
and ends after its Python return/conversion; includes GIL detach/reattach, native locks, native
work and binding result creation. Graph/model/result immutable reads keep the GIL; native work
and ledger-gate waits release it; policy destruction also detaches (`bindings/python/native_call.hpp`,
`flow_policy.cpp:144`). This is the same public boundary as the consumer timers, not pure detached
C++ time. Constructor interception observes `__init__`; Python shell allocation in `__new__`
is outside that span on both sides and must be reported as excluded boundary cost, not silently
counted as measured construction. Final-owner release spans include their decref/container cost. Concurrent wall-span sums are not elapsed replay wall time.

Keep three quantities: total replay elapsed wall (dispatch and scheduling included), sum of Core
span wall/CPU by kind, and union of span intervals on the wall axis (overlaps counted once).
The union is the proposed Core-only wall denominator for P1/P3; serial P2 reduces to its span sum.
Report dispatch gaps and the null-dispatch control beside it. E7 qualifies per-kind agreement;
large gaps or changed scheduling make a row unsuitable for a Core-only decision. No subtraction of
summed wrapper overhead from multi-worker elapsed wall. This denominator is accepted by ruling 1; E0 qualifies union and elapsed separately.

Module-global immutable defaults that survive cleanup are marked retained-at-end; do not invent a destructor event for them. Match their lifetime in the check lanes and report process-finalization work outside replay coverage.

The global trace records call entry and exit order plus thread/job-local sequence and dependencies.
The replay preserves dependencies, job boundaries and serial operations, allowing independent jobs
to interleave naturally. It does not lock all calls into the recorder's incidental global ordering.
Changed by ruling 14: every returned field is compared in every job of both lanes. Full state
snapshots are taken after each call on main, all P2/P3 jobs and the fixed P1 sampled jobs, before that job's next
call; there is no primary shared mutable ledger/cache across cases. The fixed sample is job
index modulo 64 in {0,1}, widened by the first case of any otherwise uncovered workflow group:
79 P1 cases (including g0001-000000); P3 has full state on all 592 cases, and P2 is all on main. Exact P1 job lists and
workflow/direction checks are in scratch/replay/state-sample-r14-01.json. A defect affecting only
internal state in an unsampled P1 case is outside replay coverage unless a returned field changes.

## Recording and build selection

| Choice | Proposed representation / refusal |
| --- | --- |
| Files | Versioned manifest, JSONL events and typed NumPy payload files. Every event has phase, logical job, thread-local index, entry/exit ordinal, operation, args/kwargs, returned object refs or exception, and release events. Returned containers are unpacked into element slots so a registry-held list cannot delay their native releases. Footer has event/job counts and completion status. A partial, empty, unknown-kind or unresolved-reference recording fails. |
| Values | Tagged `None`, bool, arbitrary Python int, binary64 hex/bits including infinity, enum, exact fraction, tuple, list, dict, NumPy scalar/array. Preserve array dtype, shape, strides/order and mutability when relevant to the recorded input; unsupported encountered types fail with call index and path. Arrays carry full bytes, masks carry full exclusions. Hashes can deduplicate payload bytes; references never replace missing payloads. |
| Objects | Separate logical occurrence ids and payload ids: equal masks remain separate, `intersect` aliases remain aliases, arrays retain native owners. Native-produced values stay as dependency refs to the producing call; store enough observable fields to diagnose them, never pickle a pybind object. |
| Check streams | Run reference and candidate in separate fresh interpreters. Store the first full state of each owner and exact per-call sparse deltas thereafter (ruling 4); compare by logical job/call/field, not completion order. Release diagnostic snapshots immediately; preserve only consumer-retained snapshots. Optional-output lane has matching lifetimes on both builds. |
| Exact state encoding | FlowRecord edge-flow lists compact to `$int_pairs` int64 only after exact list/tuple/int, row-length and domain checks (edge < 2**31, amount < 2**62). Order is unchanged; native ndarray outputs retain their original metadata. Record deltas use complete canonical FlowIds, explicit order and additions/removals/changes; duplicate ids fail. Independent state reconstruction and boundary controls verify this representation. Iterator validation changed under ruling 13 without reducing fields or snapshot frequency. |
| Build provenance | Explicit interpreter/source path; assert `netgraph_core.__file__`, `_netgraph_core.__file__`, version, native SHA-256, Python ABI/GIL, consumer revision and input digests in each process. A candidate is a coordinator-provided wheel venv. Never import two native builds into one process. |
| PyPI | Record with the frozen locally built reference, then replay the same representation under the existing PyPI consumer venv as a check-only cross-build control. No assumption that two wheels labelled 0.13.0 have identical binaries. |
| Size estimate, not measurement | P1 inputs/events 20–100 MiB, P2 150–500 MiB (about 1.23 million export reads), P3 10–50 MiB. Full-output check observations can be much larger: four ledger projections alone are `4*8*91712 = 2,934,784` bytes per P1 snapshot before compression. Use streaming and payload deduplication; report actual bytes in stage 1. At 4 GiB per recording or 20 GiB transient check evidence, stop and rescope rather than silently truncate. |

## Comparator: exact coverage

Run two untimed lanes from the same recorded inputs: (a) original output choices and lifetimes,
(b) all seven `OutputOptions` flags true for producing calls. The second lane compares against
reference with those same flags; it must not feed changed optional-field presence back into the
consumer's recorded control flow. Small controls cover all-off and each individual flag, plus
independent `program`/`by_source` selection. An all-on pass alone does not prove all combinations.

| Value | Every field read and compared |
| --- | --- |
| `PlacementResult` | `admitted`, `by_source`, `by_tunnel`, `call_delivered`, `call_edge_drops`, `call_source_drops`, `call_undelivered`, `delivered`, `delivery_delta`, `edge_drops`, `edge_flows`, `edge_reservations`, `marginal_cost_curve`, `min_cut`, `offered`, `plan`, `program`, `requested`, `reservations_changed`, `residual`, `routes_changed`, `saturated_boundary_edges`, `source_drops`, `termination`, `total_cost`, `unadmitted`, `weights_changed`. |
| `SourceFlow`; `TunnelFlow` | Source: `source`, `offered`, `admitted`, `delivered`, `source_drops`, `total_edge_drops`. Tunnel: `entry_id`, `head`, `status`, `assigned`, `delivered`, `total_edge_drops`; row/container order retained. |
| Program; tunnel | `sources`, `destinations`, `source_weights`, `tunnels`; each tunnel's `entry_id`, `head`, `local_index`, `destination`, `edges`, `path_cost`, `status`, `requested_bandwidth`, `reserved_bandwidth`, `planned_bandwidth`, `weight`. Check documented sharing/lifetime relationships, not addresses. |
| `TunnelPlan` | `design_by_head`, `design_demand`, `headroom_penalty`, `heads_with_capacity`, `insufficient_resolution`, `max_bandwidth_per_tunnel`, `max_fractional_bandwidth_per_tunnel`, `max_lossless_offer`, `planned_bandwidth_per_tunnel`, `total_planned_bandwidth`, `total_reserved_bandwidth`, `total_tunnel_cost`, `tunnels_per_head`. Primary workloads normally return `None`; tunnel controls exercise non-null values. |
| Ledger, policy after each call | Each live ledger's graph/uid relation, capacity, `carried`, `residual`, `reserved`, `unreserved`, full ordered `records()` with every `FlowRecord` field; policy `delivered`, resolved `config`, `program`, `records()`. Changed by ruling 9: individual `record(flow_id)` is checked on recorded consumer calls and in each check run's tiny-stream and two-source preflight, rather than repeated in every snapshot. Candidates changing record derivation or getters require a full lookup cross-check on P1's MSD roster before their rows count. Changed by ruling 14: read main's and each fixed sampled job's live mutable owners again after every call, including getters; never assume a getter preserved state to avoid checking it. No timed snapshots. Release events are checked as state transitions too. |
| `FlowRecord` and route | `flow_id`, `policy_id`, `sources`, `destinations`, `admission`, `entry_id`, `route`, `min_route_cost`, `offered`, `delivered`, `edge_flows`, `edge_drops`, `edge_reservations`; DAG `graph_uid`, `parent_offsets`, `parents`, `via_edges`. No checking only record counts or delivery totals. |
| Graph, mask, configuration | Every public getter/field in the pinned stub for an encountered type. Graph input permutation and CSR, resolution, graph uid relationships; mask allowed-node/edge vectors and uid relationships; constructor options and resolved policy configuration. Explicit type readers fail on an unhandled field/type. |
| Cache | `graph_uid`, `stats().hits`, `misses`, `entries`, `bytes`, `live_values` after each relevant call/release on main and sampled jobs; lifecycle operations in controls. Quiesce a shared-cache control before asserting a multi-field snapshot, because the contract promises individually atomic counters only. |
| Search/sensitivity if encountered | `SpfResult.costs`, `dag`, `target_cost`; ordered KSP `(cost, Path)` list; Path `graph_uid`, `hops`, `cost`, `source`, `destination`, equivalent DAG; sensitivity `baseline`, ordered entries `edge`, `delivery_reduction`. Primary inventory decides whether these need timed support. |
| Failure | Exception type (module/class) and full text, then complete exposed mutable post-state; same-call retry control. Unexpected exceptions or an unreadable field fail. No private state hooks. |

| Field class for coordinator plants | Equality / rejection rule |
| --- | --- |
| Scalar integers, booleans, enums | Exact type and value; bool is not int, enum class/name/value all match; retain integers beyond 64 bits. |
| Float / fraction | Exact binary64 bits (signed zero included), tagged infinities; exact rational numerator/denominator. No tolerance. Unexpected NaN fails rather than comparing equal accidentally. |
| Sequence and structure | Exact key sets and typed values; tuple/list distinguished, lengths and order exact, missing/extra element or field fails. `None`, empty and absent are different. No sorting cuts/curves/DAGs to hide ordering differences. |
| Arrays | Dtype, rank/shape, every element, read-only/view-owner contract; report first differing index. Values are never reduced to a checksum verdict. |
| Identity / sharing | Canonicalize only process-unique graph/ledger/mask uid numbers with a bijection to logical owners; preserve 0, equality, inequality and non-reuse. FlowId's ledger component follows that map; record counter and policy-local id remain exact. Never normalize arbitrary integers or every FlowId to ordinal position. |
| State / error | Same structural comparator reaches ledger projections, records and cache counters, including after failure and removal; error class/full text exact. |

Preflight uses a real tiny reference result twice: change only its `min_cut` while delivery/cost
stay fixed, then change only one ledger carried unit (and separately a record in the controls).
The actual comparison entrypoint must reject both with nonzero status, logical call index, field
path, reference value and candidate value. Coordinator plants in each class above follow stage 1.
The intended success path emits nonempty per-job coverage and a final checked-call summary.

## Timed mode and experiments (stage 2)

Input decoding, payload IO, external literal/array preparation, dispatch plans, saved unbound
native functions/descriptors, worker startup and null calibration occur before measured blocks.
Containers whose members are Core-produced objects are completed between spans; those objects
cannot exist before their producer calls. Pure external containers are already built.
Native constructor calls and their dependent objects are part of the measured stream: a policy
produced by a timed constructor necessarily becomes a later action's input. Such dependencies
cannot themselves be preconstructed without removing the constructor or changing its ledger.
Prepare fixed payloads and argument containers in advance, fill dependency slots between spans,
and time the original constructor/call only. No reset/evaluate shortcut replaces the real sequence. Never prebind methods to native result instances: that would extend their lifetimes.
This is the proposed interpretation of “build every input before the timed loop”; the coordinator accepted this interpretation in ruling 2.

One iteration is the full recorded Core sequence from consumer import through export and cleanup;
import/load constructor calls are tagged and included where actually reached, without re-importing Python modules inside an iteration. Workload setup and
Core graph/model construction are reported separately from run/export, with one total over the
real reuse lifetime. Each repetition gets fresh mutable owners and the same graph/model sharing
boundaries as the recording. Result getters and final-owner releases execute in their original
order; results are not retained wholesale past their recorded lifetime. Auditing, canonicalization
and comparison run in different untimed processes, never as a timed iteration's tail.

| E0 row | Initial replay block target | Initial fresh CLI processes per consumer block |
| --- | --- | --- |
| P1 / 1 | 2 s minimum; one full sequence will likely exceed it | 1 |
| P1 / 8 | 2 s minimum (rulings 17/23: 1.1 count margin) | 1 |
| P2 / 1 | 2 s minimum | 2 |
| P2 / 8 | 2 s minimum (rulings 17/23: 1.1 count margin) | 3 |
| P3 / 1 | 2 s minimum | 6 |
| P3 / 28 | 2 s minimum (rulings 17/23: 1.1 count margin) | 12 |

These counts are initial estimates from exploratory STATE pipeline times, not qualified lengths.
In the coordinator's window calibrate each row on reference, fix an integer repeat/process count
before its seven triples, and reject any block below 2 s of its declared measurement. Increase to
the next longer block on noisy rows, within the two-hour total; retain failed shorter blocks.
Planned CLI: `PYTHONPATH=<pinned-source> <wheel-venv>/bin/python -B -m ngraph --quiet run <input> --results <unique-output>`. For one-worker variants, prepare and name a scratch YAML differing only in the relevant workflow worker fields; keep both original and derived input hashes.

Consumer blocks sum complete fresh CLI-process start-to-exit times; input generation is outside,
all load/prepare/run/export/materialization is inside, each process gets a fresh output path.
They use the reference wheel, not STATE's PyPI environment. Add import and empty-run subprocess
controls and compare complete exports with only explicitly identified producer-identity fields
excluded by name, never all metadata. Poison one export before timing.

| Experiment | Exact output and decision boundary |
| --- | --- |
| E0, 12 rows | Six workload/worker combinations each for replay and consumer; seven fresh-process A/A/A triples per row. Save raw block/iteration times, counts, summed spans, union and elapsed wall, process CPU, peak RSS, build/import provenance, clock resolution, before/after machine/swap/power state, null-dispatch overhead. |
| Statistics | Use `aba.py`'s `B / mean(A1,A2)` per triple; median ratio and 95% whole-triple bootstrap interval (10,000 resamples, fixed seed 20261108), drift `median(abs(A1/A2-1))`; report all blocks and clusters. Drift >3%, short blocks, sleep/power/swap problems keep a row inconclusive. Memory: largest per-side fresh-process peak, fail above `max(1.10*reference, reference+8 MiB)`. |
| Smallest detectable gain | Report an explicitly operational A/A noise bound: for A/A upper ratio bound `U`, a uniform speed gain must exceed `U/0.97` to cross the beyond-noise bar; for a native bar `B`, it must reach `B*U`. These are shifted observed intervals, not a prospective statistical-power claim. Show spread/CI and noise-only bound as well; use a longer block if uncertainty defeats the applicable 3% or 1.15 native decision bar. |
| E7 | Same reference wheel, input sequence, phase, actual call/getter/release counts, outputs, lifetime and worker count on both replay and consumer. Pair the six workload/worker rows; give count, raw total, per-call mean and uncertainty for every kind. Fresh matched native-timer runs use the recorder's wrappers with payload recording disabled; prior ORIENT-ADD times are diagnostic only. |
| E7 qualification | Proposed agreement band 5% per kind, fixed before timing; bootstrap matched block ratios. Require interval within `[0.95,1.05]`, identical counts/inputs and overhead small enough to resolve that band. Keep getters grouped by exact field; report unresolved tiny kinds when timer/null cost is comparable. No subtraction of overlapping wall spans. An unqualified kind cannot support a replay decision; use consumer comparison. |
| Destructor E7 | Weakref collection alone cannot time native destruction. Ruling 3 excludes a release adapter: time replay reference drops, report releases as an unqualified group against ORIENT-ADD diagnostics. Destruction gains require a consumer comparison. |

All full recordings/checks use `HEAVY_OWNER=replay` and bounded commands; worker-run times remain
exploratory. Only the coordinator runs qualified E0/E7 unconfined, under one
`HEAVY_OWNER=replay-timing` hold for at most two hours. Scripts print the full plan and expected
budget before starting; each child has a timeout; new output names preserve evidence.

## Size, dependencies and stops

| Stage / mechanism | Cost, smaller alternative and what it misses |
| --- | --- |
| Recorder + tagged values + object/dependency events | About 250–350 lines. Reusing hash-only Tracer is smaller but cannot replay masks, payloads or lifetimes; wrapping only `apply_actions`/evaluators misses preparation, removals, constructors and export getters. |
| Check reader + driver + controls | About 250–350 lines plus 80–140 control lines. A delivery-only comparator misses the concrete equal-delivery/different-cut and orphan-record defects the research must reject. Explicit small readers are preferred to a general object serializer. |
| Timed mode + E0/E7 orchestration | About 180–250 lines, reusing the same event stream and `aba.py` statistics where possible. Separate per-workload benchmark scripts would duplicate lifetime/order logic. |
| Total estimate | 680–950 tool lines plus controls, several hours in stages. Inventory probe is disposable scratch evidence. Stop for rescoping if implementation needs roughly twice this size or unbounded retention to work. |
| Existing environments | Reference venv and consumer PyPI venv are sufficient now. `ruff` is present in `worktrees/core-ref/venv/bin/ruff`; it is absent from the frozen consumer-reference venv. Use Core's ruff with caches off, without modifying either venv. |
| Scratch consumer | None: ruling 3 prohibits a release adapter or scratch consumer copy. No new environment or network dependency. |
| Stage 1 | Record all P1–P3, reference/reference lane (a) at 1/scenario workers, lane (b) at scenario workers, PyPI/reference at a/1 only (ruling 14), every returned field everywhere plus sampled full-state coverage, and size reports, comparator preflight. Initial limits: 180 s per full recording, 600 s per full check; use bounded job groups if needed. Stop for coordinator plants. |
| Stage 2 | Implement timed mode and E0/E7 scripts; only short untimed dry runs here. Stop for coordinator's qualified window/raw outputs, then analyze E0/E7. Fresh readers are coordinator-managed after freeze. |

Stage 1 is accepted under ruling 15: the coordinator rejected 36 planted mismatches and passed
four controls. Stage 2 code is isolated under timing/; the accepted root sources and matrix
streams remain unchanged. TIMING.md gives the boundary, metadata-only CLI launcher, pilot,
single-hold two-hour window and analysis commands. Only short unqualified controls run here;
primary qualification and E0.md/E7.md await the coordinator's window and raw data.

Rulings 19–23 supersede earlier broad timing/state statements: only P1 may make a replay union
decision, only if changed code reaches E7-approved kinds and preserves call sequence, visible
lifetimes, cross-thread behaviour and Python object counts. P2/P3 timing is diagnostic and
consumer comparison decides; CLI rows determine memory risk. P2/P3 check mode retains full
state, independently of these timing restrictions. Pilot counts include 10% margin; only pilot
exports retain gzip archives, while compared window exports retain size/SHA receipts.
