# Independent review of the SPF optimisation

Reviewed 2 October 2026 against `fbab119e794c843838f4f3f96d75cb7f32f06845` and the supplied uncommitted diff. Line references below refer to the live, unmodified production files, not the experimental copies.

**No: the main algorithm choice is sensible, but several small, worthwhile changes were neither considered nor verified.** I measured gains in neighbour processing, predecessor output construction, singleton-frontier handling, and fixed-routing callers. More importantly, I reproduced a new exception-safety defect and disproved the stated retained-memory bound. The weighted-graph speedup itself survives independent measurement.

All experimental sources, patches, scripts and results are in this directory; executables/objects are under `build/perf/`. All 97 tracked files were checked byte-for-byte against the review-start snapshot; none changed ([verification](results/tracked-verification.txt)). No full repository suite, commit, checkout, stash, or downstream source modification was performed.

## 1. Ranked findings

### 1. A failed touched-list allocation corrupts the next SPF call

**Verdict: VERIFIED-GAIN — correctness restoration, not a latency claim. Priority: fix before committing.**

`src/shortest_paths.cpp:726–728` writes the source distance and bottleneck **before** the potentially allocating `touched.push_back`. Reverse SPF does the same at `:932–934`. If that allocation throws, the lease runs, but sparse cleanup at `:433–436`/`:475–479` has no record of the modified node. The next same-sized search can suppress a valid relaxation against the stale zero distance. Sparse output materialisation then reports the reachable node as unreachable.

[Reproducer](robustness.cpp), [public old-build reproducer](fault-public.cpp), and [exact fix](exception-fix.patch), affecting those two blocks in `src/shortest_paths.cpp`:

| Eight-node, one-edge graph; inject one four-byte allocation failure, then search from/to the other endpoint | Old | New | Fix |
|---|---:|---:|---:|
| Forward distance to the reachable node | 1 | `INT64_MAX` | 1 |
| Reverse distance from the reachable node | 1 | `INT64_MAX` | 1 |

Raw outputs: [old](results/fault-old.txt), [new](results/fault.txt), [fix](results/fault-fixed.txt). The patch records the source/target first, then performs the nonthrowing writes. It does not change successful-call work or outputs. Cost/risk: two statement reorderings; very low. The existing fanout-throw test at `tests/cpp/shortest_paths_tests.cpp:641–668` does not exercise this failure point. The comment claiming restoration on *every* exit at `src/shortest_paths.cpp:497` is presently false.

### 2. Do not recompute an identical SPF inside fixed-routing max-flow calls

**Verdict: VERIFIED-GAIN.**

`src/max_flow.cpp:123–151` and `src/flow_state.cpp:537–552` repeatedly call SPF. With `require_capacity=false`, they pass no residual; graph, masks, endpoints and selection remain identical. The next placement still needs current residuals, but the next SPF does not. This is a particularly simple caller optimisation, independent of the queue. It is unnecessary when `shortest_path=true`, which already stops after one successful tier.

[Exact patch](fixed-routing.patch) caches the result **within the invocation**, in both files, only re-searching when capacity-aware routing is enabled. It does not cache changing residual-aware routes. On eight weighted 4,096-node/32,768-edge queries, seven samples per fresh process:

| Workload | Old A1 / A2, ms | New A1 / A2, ms | Reuse variant, ms | Reduction from new |
|---|---:|---:|---:|---:|
| `place_max_flow`, fixed routing | 5.757 / 5.646 | 3.857 / 3.784 | 2.134 | 44–45% |
| Same, node/edge exclusions | 4.255 / 4.150 | 3.090 / 2.951 | 1.693 | 43–45% |
| `calc_max_flow`, exclusions | 6.032 / 6.033 | 4.942 / 4.846 | 3.531 | 27–29% |
| `calc_max_flow`, no exclusions | 7.105 / 6.903 | 5.139 / 4.972 | 4.486 | Unstable samples; do not use for a precise claim |

The unmasked `calc_max_flow` variant samples span 3.851–5.648 ms; I am explicitly not treating its median as equally strong evidence. The other gains also appeared in the screening pass. The accepted rerun fixes all query endpoints in the masks before warmup; the earlier screening files are retained and labelled `screen-flow-*`.

Evidence: [driver](flow-driver.cpp), [build/run script](build-flow.sh), [old A1](results/flow-old1.csv), [new A1](results/flow-new1.csv), [variant](results/flow-fixed.csv), [new A2](results/flow-new2.csv), [old A2](results/flow-old2.csv), [load](results/load-flow.txt): 1-minute load 4.46–4.76. Output hashes agree; capacity and conservation checks run outside timing. See section 2 for stronger flow validation.

Expected applicability: fixed IP/IGP routing through these helpers, with repeated placement enabled. No promised improvement for residual-aware default max-flow or the recommended single-pass IP-ECMP combination. Cost/risk: small local lifetime extension; review the old result's extra lifetime during replacement on the residual-aware route. Keep this as a separate caller patch if the SPF change is to remain narrowly scoped.

### 3. Skip settled neighbours before collecting their parallel edges

**Verdict: VERIFIED-GAIN.**

At `src/shortest_paths.cpp:583–640` and `:832–881`, the code collects eligible edges and computes their bottleneck even when the neighbour is settled. With nonnegative costs, a settled neighbour cannot receive a shorter path; equal-cost predecessor changes and single-path capacity improvements are already prohibited for settled neighbours at `:644–655`/`:884–895`. Skipping its contiguous neighbour run early preserves the existing outputs.

[Exact patch](settled.patch): two condition changes in `src/shortest_paths.cpp`. This is not a different tie-break, cached eligibility policy, or new queue.

| Dense unit fabric, 800 nodes / 320,000 edges, 16 full-source queries | Old A1 / A2 | New A1 / A2 | Skip-settled |
|---|---:|---:|---:|
| Pass A, ms | 36.784 / 36.827 | 34.775 / 35.326 | 30.552 |
| Pass B, ms | 36.005 / 36.967 | 35.225 / 36.691 | 29.811 |
| Reverse search, 16 targets, ms | 48.257 / 46.690 | 46.923 | 37.042 |

Forward improvement is 12–19% against the new code; reverse improvement is 21%. Weighted random graphs do **not** establish a gain: pass B is 36.487 versus 36.437/36.404 ms. Do not advertise this as a universal improvement. The extra branch could be unfavourable in other graph families.

Evidence: `results/{old-a1,new-a1,settled,new-a2,old-a2}.csv`, the corresponding `*-b*` files, and `results/reverse-*.csv`; [first load log](results/load.txt), [second load log](results/load-second.txt). Loads during these forward phases were 3.53–4.42; reverse phases 4.67. Full forward old-loop comparison: 1,351,168 comparisons, zero failures. Reverse hashes agree with the separately built old implementation; targeted sanitizers and the cross-caller bitgate also pass.

Cost/risk: very small implementation, low semantic risk after nonnegative-weight proof and parity checks. This directly challenges the suggestion at `CPU_PLAN.md:150–155` that edge-scanning dominance makes unit fabrics merely neutral: some of that edge processing is avoidable.

### 4. Forward PredDAG conversion need not walk every surviving list twice

**Verdict: VERIFIED-GAIN.**

`src/shortest_paths.cpp:679–700` walks each predecessor list to count entries, prefixes the counts, sizes/initialises two output arrays, and walks the same lists again to fill them. On dense equal-cost DAGs this is measurable work, not an irreducible result-size floor.

[Exact patch](append.patch) keeps the sparse path unchanged. For the dense case it reserves the number of recorded entries as an upper bound, appends surviving entries once in node/list order, and writes each final offset immediately. It preserves parent ordering and all four output arrays.

| Same dense forward workload as finding 3 | Old A1 / A2, ms | New A1 / A2, ms | One-pass conversion, ms |
|---|---:|---:|---:|
| Pass B | 36.005 / 36.967 | 35.225 / 36.691 | 29.755 |
| Pass C | 37.396 / 37.468 | 35.800 / 35.380 | 31.159 |

This is a reproduced 12–19% reduction. Random-graph changes remain below the review's significance threshold. Combined with settled-neighbour skipping, the dense row is 26.327 ms versus 35.800/35.380 ms, about 26% lower; that combined timing is one phase, not a separately replicated universal claim.

Evidence: [patch](append.patch), [combined patch](combined.patch), `results/{append,append-c,combined}.csv`, [pass C load](results/load-third.txt). Pass B load was 3.53; pass C was busier, 7.93–8.46. The individual improvement survives both passes. [Forward comparison](results/correctness-append.txt): 1,351,168 comparisons, zero failures.

Expected applicability: dense multipath outputs with many live predecessor entries. Cost/risk: small conversion change, but reserving *all* recorded entries can retain excess space in the returned vectors when many predecessors were superseded. This must be weighed against the exact-sized current output; it is not a free memory improvement. Reverse conversion is a scatter by successor and requires ascending-tail order (`:981–1015`), so this exact patch cannot simply be copied there.

### 5. The chain regression is partly avoidable without another search algorithm

**Verdict: VERIFIED-GAIN.**

For a frontier with one item, `BucketQueue` still divides/modulos its key, sets/clears occupancy, advances the cursor and touches a vector/heap (`src/shortest_paths.cpp:329–362`). The plan tried modulo masking and accepted the regression; that does not test avoiding the bucket machinery altogether for an isolated item.

[Exact patch](inline.patch) stores a lone newly pushed item inline, promoting it into the existing buckets if a second item arrives. After ordinary bucket operation it only re-enters this mode when the frontier empties. Keys, queue selection and pop order are unchanged. This is not retrospective dispatch to a second Dijkstra implementation.

| Chain: 2,048 nodes, 256 queries with varying sources | Old A1 / A2, ms | New A1 / A2, ms | Inline item, ms |
|---|---:|---:|---:|
| Pass B | 2.929 / 3.114 | 5.232 / 5.462 | 3.541 |
| Pass C | 3.039 / 3.115 | 5.607 / 5.569 | 3.618 |

Reduction from the new bucket route: 32–35%. It still trails the old implementation by about 14–21% in these passes. These are 256 varying-length chain queries, **not** the plan's single full-chain query; do not quote the whole batch as single-query latency. Other focus rows show no reproduced material gain or regression.

Evidence: `results/{inline,inline-c}.csv`, the B/C controls above, and [correctness](results/correctness-inline.txt): 1,351,168 comparisons, zero failures. Same load ranges as finding 4. A smaller experiment retaining just a sole item's bucket location ([singleton.patch](singleton.patch)) gave only about 9% in the first chain pass; bypassing the bucket storage is the useful version.

Cost/risk: modest queue-state complexity; promotion after zero-cost pushes and cyclic cursor transitions deserve explicit tests. Applicability includes narrow-frontier portions of ordinary graphs, not just graphs that are globally chains. There is no workload evidence behind the plan's categorical statement that no production caller runs on a path graph (`CPU_PLAN.md:287–292`).

### 6. Retained memory is not bounded by one largest-graph workspace

**Verdict: PLAUSIBLE-UNTESTED for the reclamation optimisation; the counterexample itself is verified. Priority: correct the bound before committing.**

Only the fixed node arrays go through `ensure_sized` (`src/shortest_paths.cpp:402–405,424–430,465–472`). Entry arrays, touched/selection buffers, heap storage, bucket headers and each individual bucket's item capacity retain their high-water marks. Both directional workspaces coexist per thread (`:494–495`).

The stronger problem is the **sum of independent bucket high-water marks**. [Memory probe](robustness.cpp), [output](results/memory.txt): with 63 graphs, each having only 1,027 nodes and 1,026 edges and a maximum live frontier of 1,024 items, retained bucket-item capacity rises from **16,400 bytes to 1,032,208 bytes**. A subsequent two-node graph shrinks `dist.capacity()` to two but retains all 1,032,208 bucket bytes.

Thus `:260–264` is misleading and `CPU_PLAN.md:368–369` accounts only for empty headers/bitmap. A conservative bound is proportional to `W_max × historical maximum queued items`, plus node/entry/heap buffers, **per direction per live thread**, rather than just the largest simultaneous frontier. The multiplier can be much more important than the 1.6 MB header number. Entries and heap capacity also survive a graph-size shrink.

Expected benefit: memory reclamation after heterogeneous graphs or worker-pool traffic; no measured latency gain. Cost/risk: add a total retained-byte budget or reclaim bucket/item storage on substantial graph-size shrink; avoid destroying hot capacities every query. I did not benchmark such a policy because choosing a memory budget requires a workload trade-off; pretending that `W_max` already supplies it would be incorrect. Thread exit does free the TLS vectors; long-lived Python/executor threads are the relevant retention case.

### 7. The fallback regression is real in some rows, but its “workspace cost” explanation is not established

**Verdict: PLAUSIBLE-UNTESTED for the remaining causal attribution. Specific cheap remedies below are REFUTED as material gains.**

`CPU_PLAN.md:156–162` converts an unresolved difference into “the cost of the reusable workspace.” The new route does extra touched-list work, copies dense distances (`src/shortest_paths.cpp:675`), and resets arrays after copying (`:438–443`); these are candidates, not a measured attribution. The three previously attempted mitigations were not repeated.

I isolated three other possibilities, using the same production loop:

| Pass A/B, 16 full queries | New forced heap | Experiment | Outcome |
|---|---:|---:|---|
| Random 16k, remove touched tracking and force dense handling | 65.364 / 66.960 ms | 64.784 ms | Under 5%; no established full-search gain |
| Wide ~2^20 costs, same ablation | 52.299 / 53.113 ms | 51.040 ms | Under 5%; insufficient to attribute the gap |
| Near-target queries, same ablation | 0.399 / 0.408 ms | 1.233 ms | Destroys the intended workspace benefit |
| Wide ~2^20, explicit heap item/comparator | 52.299 / 53.113 ms | 55.679 ms | No gain |
| Wide ~2^20, copy heap front before pop, like old priority_queue | 52.248 ms | 52.007 ms | Noise |

Patches: [notouch](notouch.patch), [heapitem](heapitem.patch), [frontpop](frontpop.patch). The 24-byte heap item versus 16-byte bucket item is a genuine layout difference; merely replacing the tuple with an explicitly ordered struct did not help. Do not narrow the double bottleneck or int64 distance to pack a key.

The quarter-touched threshold at `src/shortest_paths.cpp:399` is not established as the crossover: the reported workload matrix mostly samples very sparse or nearly full searches. Reset and output materialisation need not have the same optimal threshold; reverse sparse conversion additionally sorts touched nodes (`:986–987`). Intermediate touched fractions, locality and graph-size churn remain unmeasured. This is a plausible tuning opportunity, not a verified gain, and the no-touch ablation is not a threshold sweep.

Fresh old/new controls vary: pass B wide ~2^20 is 48.972/50.562 ms old versus 53.175/52.534 ms new; pass A is 50.820/52.036 versus 53.541/53.120. The boundary run at the same broad cost range is 25.107/25.151 versus 27.498 ms for eight queries. These support a workload-dependent penalty, sometimes clearly over 5%, not a fixed universal 5–9% tax. Cost/risk of further attribution: targeted search/finish/reset profiling on a quieter machine. Expected recoverable percentage remains unknown. Neither “noise” nor “inevitable” is an adequate root cause.

### 8. `W_max=65536` is conservative and useful, not a verified optimal or portable threshold

**Verdict: PLAUSIBLE-UNTESTED for the remaining bitmap/cleanup refinements.** The bounded-dispatch design and deferred wider-queue alternatives are already covered by the plan.

The recorded sweep jumps from `W=32768` to approximately `2^20`; it never measures the chosen boundary (`cpu_queue/harness.cpp:339–341`). Passing `2097152` as harness `maxW` does **not** change the production `ref` route's compile-time limit: `Runner::run` immediately calls `run_ref` (`:147`), while `src/shortest_paths.cpp:270,368–393` decides actual dispatch. The CSV's harness `fallback` field therefore is not a production fallback counter.

I added the missing near-boundary rows, seven samples, 16,384 nodes/131,072 edges:

| Actual W | Old A1 / A2, ms | New auto, ms | New forced heap, ms |
|---:|---:|---:|---:|
| 32,768; eight full queries | 25.959 / 26.272 | 17.628 | 27.216 |
| 65,535; eight full queries | 25.433 / 25.761 | 18.147 | 26.350 |
| 65,536; eight full queries | 25.641 / 26.609 | 18.499 | 27.175 |
| 65,537; eight full queries | 25.659 / 25.454 | 26.491 | 26.865 |
| 65,535; 64 near queries | 1.131 / 1.123 | 0.455 | 0.390 |
| 65,536; 64 near queries | 1.171 / 1.149 | 0.435 | 0.399 |

[Raw boundary files](results/boundary-new1.csv), corresponding old/heap controls, and [load](results/load-second.txt), 4.42–4.67. The boundary supports full-search acceleration on this M4; near queries pay a 9–17% bucket overhead relative to the new heap while still beating old substantially. This is evidence against treating W alone as a universal profitability predictor, not a request to add a complex dispatcher immediately.

The occupancy bitmap is already implemented. However, `settle_active` can scan up to `ceil(W/64)` words (`:349`); `end()` also scans **every word**, including zero ones (`:317`), and `begin()` clears every word (`:312`). “Clear only occupied buckets” does not mean O(number of occupied buckets) total cleanup. A second-level nonempty-word bitmap or touched-word list addresses this overhead; its benefit at small W is unproven. Both workspaces also call both queues' `end()`, so a heap query after a wide bucket query still scans that old bitmap.

The plan already defers a hierarchical bitmap/radix alternative and x86 verification (`CPU_PLAN.md:367–373,433–438`); these are not omitted discoveries. Do not widen the limit solely from the old full-search sweep, especially before bounding retained item capacity and measuring cold construction, sparse frontiers, masks and graph changes. No second machine was available in this review.

### 9. An O(N) result floor exists, but the proposed distance-only remedy is misstated; early exit also has a compatibility trap

**Verdict: PLAUSIBLE-UNTESTED for a target-specific internal result path.**

`src/shortest_paths.cpp:671–687` necessarily fills N distances and N+1 offsets under the current API. Sparse reset correctly removes scratch initialisation, not these outputs. But a **distance-only vector still has N elements**; it cannot remove the O(N) floor as claimed at `CPU_PLAN.md:454–457`. A scalar target-distance result, compact target subgraph/path, caller-owned output, or a lazy/sparse result contract is needed. The internal callers commonly inspect only `dist[dst]` (`src/max_flow.cpp:142`, `src/flow_policy.cpp:184`), and `place_max_flow` does not use distances at all (`src/flow_state.cpp:543–549`). They are concrete candidates for a narrower internal result interface, without changing Python defaults. Removing the distances alone still leaves the dense DAG offsets.

Also, the current early-exit code is not simply “process exactly through the target's distance.” Stale-entry `continue` branches at `:566–569` bypass the stop check at `:661`. A stale equal-distance entry behind the target can cause a more expensive node to be expanded. [Six-node reproducer](early-exit.cpp): target cost 1, yet a node reachable only by expanding a distance-2 node gets distance 3. [Old](results/early-old.txt) and [new](results/early-new.txt) both do this.

Moving the stop test to the loop head would save work, but changes the returned tentative off-target distances/DAG and fails strict old-result parity. Likewise, breaking immediately when the destination is popped or skipping all relaxations beyond its distance is not a drop-in optimisation. Use an explicit narrower internal contract if taking this route. Expected gain: large only when output/unused exploration dominates; not measured. Cost/risk: moderate API/contract work. Bidirectional search/A* and canonical incoming-CSR reconstruction are not equally cheap substitutes under the current ordered-output contract.

### 10. FlowPolicy already avoids many searches; it still has avoidable copies and post-search rejection

**Verdict: PLAUSIBLE-UNTESTED for the remaining copy and bounded-search opportunities.** Existing memoisation and shared seeding already cover much of this area.

Memoisation is not missing: `src/flow_policy.cpp:148–205` matches exact state stamps or residual bytes, and `:434–444` computes the initial bundle once. The header documents a prior 94% repeat rate for the measured EqualBalanced cycle and 0% for Proportional (`include/netgraph/core/flow_policy.hpp:215–228`). This is historical profiling evidence, not a distribution I remeasured here. Residual-aware memoisation must retain exact capacities: unchanged eligibility does not preserve bottleneck ordering or zero-cost DAG choices.

There is nevertheless a definite redundant copy: after obtaining an owned `dag`, `get_path_bundle` returns `std::make_pair(dag, dst_cost)` at `src/flow_policy.cpp:238`. Moving that local would avoid another full DAG copy on both hit and miss. A hit still needs its own owned result because the cache retains its DAG. Also, `has_min_flow` remains a memo-key discriminator at `:159–160`, although residual awareness now depends solely on `require_capacity_` at `:94–101`, and EB skips the threshold-derived edge mask. Revisit that stale key/comment with a call-count reproducer before changing it.

Cost ceilings are applied after SPF at `:207–232`; when a bound is already known, an internal bounded-target search could reject earlier. Respect `best_path_cost_` updates and equal-cost/zero-edge behaviour. Expected gain: one DAG copy per bundle return, fewer equivalent-key misses, and less work on rejected high-cost routes; workload speedup unmeasured. Cost/risk: move is trivial; memo-key/ceiling changes require focused semantic tests. Broad new memoisation is not justified by the existing measured Proportional miss rate. I did not build an isolated memo-hit/bundle-return benchmark for this private helper; the source proves the copy, not a measurable end-to-end gain.

### 11. KSP's main missed cheap cost is not necessarily its private Dijkstra queue

**Verdict: PLAUSIBLE-UNTESTED for mask/result work.** Deferring private-queue replacement is already justified in the plan.

`dijkstra_single` is used for the initial path (`src/k_shortest_paths.cpp:204`); its comparator only compares cost (`:61–63`), so substituting SPF's capacity/node key can change first-path selection. Most subsequent spur work calls the newly optimised public SPF (`:283`). Replacing the private queue without profiling may have little impact on k>1.

Every spur builds byte masks of N/E elements (`:247–267`), then allocates bool arrays and copies them again (`:269–279`). Reuse bool scratch directly, reset only the prefix/edge exclusions that changed, and preserve the base masks. KSP already caps enumeration and rejects candidate families against the cost ceiling (`:289–305`), but the ceiling is checked *after* the complete SPF/DAG result. A bounded internal search could help rejected spurs. Expected benefit: fewer O(N+E) allocations/copies per spur, especially on large graphs and small reachable regions; no percentage established. Cost/risk: moderate scratch-lifetime and ordering work. The one k=8 Clos benchmark does not profile these components or cover long WAN paths and changing k.

### 12. Graph layout and metadata are not free, but the obvious gcd micro-optimisation is not a demonstrated win

**Verdict: REFUTED for a material gcd-guard gain in these measurements.** A separate adjacency-local attribute layout remains an untested possibility.

[One-line gcd patch](gcd.patch) stops calling `std::gcd` once the running gcd reaches 1 (`src/strict_multidigraph.cpp:66`). I timed the **whole constructor**, not just that operation: 16,384 nodes, 131,072 edges, nine samples, old/new/guard/new/old. Random-input construction was 9.060/9.137 ms old, 9.164/9.289 new, 9.122 guarded. Already cost-ordered input was 1.344/1.374 old, 1.430/1.448 new, 1.364 guarded; the repeat was 1.448/1.456 new versus 1.411 guarded. The smaller repeat does not reproduce a >5% gain. Do not sell this as a significant speedup. Output graph hashes agree. Raw files: `results/build-*.csv`, `results/build-repeat-*.csv`; load logs are the B/C logs.

The graph already has forward and reverse CSR. Global `(cost,src,dst)` edge ordering (`src/strict_multidigraph.cpp:80–108`) means a CSR row gathers cost/capacity through edge IDs, and the comment that *all* parallel edges are consecutive at `src/shortest_paths.cpp:580` is too strong across cost tiers. Cost-first rows are useful, but not the same as CSR-local cost/capacity arrays. Reordering edge IDs or neighbour visitation risks documented edge-ID tie-breaking and predecessor order. An auxiliary CSR-aligned attribute layout could preserve IDs/order but costs memory and still needs edge-ID gathers for residuals/masks. No layout experiment establishes that trade-off here.

The max/gcd metadata itself is correct and cheap to access, travels with graph copies, handles zero costs, and avoids repeated eligibility scans. A masked-out wide edge can force the whole query onto the heap; computing a tighter per-query profile is not free and needs amortisation over a stable mask. It is not automatically a cheap win.

### 13. Python copies and threading deserve separate measurements; “release the GIL” is already done

**Verdict: PLAUSIBLE-UNTESTED for owned-output transfer and worker reuse.** GIL release and bounded existing batches are already implemented.

`bindings/python/module.cpp:185–206,228–255` copies residuals/masks before releasing the GIL. The copied inputs give the native call a snapshot; blindly borrowing writable NumPy memory would remove that protection while another Python thread runs. The GIL is already released around SPF, reverse SPF, max-flow, KSP and placement. `CpuBackend` forwards spans directly (`src/cpu_backend.cpp:24–37,40–51`); removing its small duplicated validation is not a credible main optimisation.

Outputs do have extra traffic: the binding allocates/copies or converts the distance vector and passes `res.second` as an lvalue into `py::make_tuple` (`bindings/python/module.cpp:208–217,257–266`). The installed pybind11 caster maps an lvalue under automatic policy to a copy (`venv/lib/python3.14/site-packages/pybind11/include/pybind11/detail/type_caster_base.h:1650–1659`). Move the owned DAG; consider a capsule owning the int64 vector for a zero-copy int64 result. Float64's sentinel conversion still needs a pass. No separately rebuilt pybind extension variant was timed in this review, so the existing-module NetGraph rerun cannot verify these transfer changes. Expected benefit: output-heavy and tiny-search Python calls; not measured, and public ownership/mutability/default-dtype tests are needed.

`batch_max_flow` already dynamically schedules bounded `std::async` workers (`src/max_flow.cpp:354–379`), and sensitivity has its own bounded workers (`:427–454`). TLS scratch is reusable within those workers but is destroyed when each batch's workers exit. Persistent workers may amortise first-use allocation; nested Python pools plus native batches can oversubscribe the machine. The NetGraph benchmark explicitly uses serial Monte Carlo/scenarios, so it cannot establish throughput or persistent-pool memory behaviour. Sequential augmentations sharing residuals cannot be batched as independent searches.

### 14. There is another potentially larger dense-placement hotspot outside SPF

**Verdict: PLAUSIBLE-UNTESTED.**

`src/flow_state.cpp:174–199` finds unique DAG parents by linear search and then rescans the entire parent list once per unique parent. For a node with P distinct predecessors, this is quadratic parent processing before the placement solve. Preserve first-appearance and edge accumulation order with a reusable indexed grouping/counting structure; do not sort indiscriminately because floating-point accumulation and group order are observable.

Expected applicability: wide ECMP/Clos DAGs where queue changes hardly move end-to-end placement time. Expected asymptotic improvement is O(P²) to O(P), but I did not measure a latency gain or implement it: maintaining arbitrary caller-supplied PredDAG order and all three placement semantics needs a separate patch and feasibility/conservation/cut checks. Cost/risk: larger than the four small verified performance changes above. This is a concrete profiling target, not evidence that the existing solver can simply be replaced.

### 15. Arithmetic/empty-graph/Python-thread risks are mostly covered; reentrancy is a separate assumption

**Verdict: ALREADY-COVERED for ordinary arithmetic, edgeless graphs and independent threads.** Same-thread reentrancy remains a separate, unverified use case.

The constructor enforces nonnegative costs and total cost below `2^62` (`src/strict_multidigraph.cpp:40–66`). Therefore `max_cost/gcd+1` cannot overflow int64; cursor advancement corresponds to an actual queued nonnegative distance, and its stride step is bounded by an edge-cost window. I found no ordinary signed-overflow counterexample in the bucket arithmetic. All-zero costs give W=1; no-edge graphs need no edge access despite empty residual/mask sizes comparing equal to E=0.

[Targeted sanitizer probe](safety.cpp) covered exact W=65,535/65,536/65,537, zero costs, costs near the constructor's aggregate bound, N=0 and N>0/E=0, both directions, option combinations, reused workspaces and concurrent C++ threads. It passed under ASan/UBSan. This does not substitute for Python-thread throughput measurements. Independent immutable-graph calls receive independent TLS state; Python input copies and GIL release remain unchanged.

There is no in-use guard or nested lease stack for `tls_forward`/`tls_reverse` (`src/shortest_paths.cpp:494–501`). A same-direction recursive call on one thread would reset the outer call's workspace. I found no supported callback/reentry path in these synchronous loops; ordinary Python threads are not such a path. Document the non-reentrant assumption or add an in-use fallback if callbacks/custom allocation hooks become supported. Do not present a hypothetical recursive call as a reproduced production failure.

## 2. Experiments, exact patches and validation coverage

The reproducible performance patches are linked in findings 2–5; the exception fix is linked in finding 1. Each `.patch` contains the full unified diff against the unmodified live source. Experimental `.cpp` copies let them be built without applying anything to the working tree.

- [build-review.sh](build-review.sh) creates old sources/headers using read-only `git show fbab119:...` and builds the focused harness variants. Same Apple clang 21, C++20, `-O3 -funroll-loops -fno-math-errno -fno-trapping-math` as the supplied harness. The supplied `build.sh` was also run successfully.
- [driver.cpp](driver.cpp) reuses the supplied generators/query setup and adds focused, reverse and boundary rows. Random16/wide rows contain 16 queries; Clos contains 16; near contains 64; chain contains 256 varying sources. Every phase is a fresh process; two warmups precede seven measured samples. Raw CSVs include all samples.
- [run-review.py](run-review.py), [run-second.py](run-second.py), [run-third.py](run-third.py) record phase order/load and keep old and new controls on both sides of the focus experiments. [environment](results/environment.txt) records compiler/host/source/binary hashes; [final manifest](results/final-sha256.txt) covers later experimental builds and sources. Pass C is deliberately reported as noisier, not quietly substituted for pass B.
- Baseline, skip-settled, inline-item and one-pass-output versions each passed the supplied forward comparison matrix: **1,351,168 comparisons each**, against the preserved old-loop variants, not just distances. These are repeated combinations on a finite set of graphs, not millions of independent graphs. Files: `results/correctness-{baseline,settled,inline,append}.txt`.
- [validate-review.sh](validate-review.sh) builds settled+one-pass output together with the fixed-routing caller patch under **ASan/UBSan**, with sanitizer failures fatal. Its targeted fixture checks nonuniform-cost flow feasibility/conservation, all single-edge exclusions, exhaustive small-graph cut optimality and reported cut capacity, plus FlowPolicy cost factors 1 and 2. [Results](results/validate-flow-variant.txt), the build/run log (not retained).
- The same sanitized combined/caller build passes the existing cross-caller dump: **6,226 lines byte-identical** to the old binary, including forward/reverse/fanout, placement modes, FlowPolicy sequences, batch/sensitivity and KSP. [Hashes](results/bitgate-sha256.txt). Native batch/sensitivity thread budgets were set to 1 for this check to keep load bounded. The inline queue was separately checked by the full forward matrix; it was not included in that combined sanitized build.
- [safety.cpp](safety.cpp) adds the width/int64/edgeless/concurrency probes described above. [Output](results/safety.txt).
- [run-netgraph-review.py](run-netgraph-review.py) reruns four NetGraph workloads with five samples, `PYTHONHASHSEED=0`, bytecode writes disabled, explicitly selected old/new modules. [Module paths/hashes/load](results/netgraph-environment.txt) verify the same extension binaries used by the recorded experiment. No downstream rebuild was needed or performed; these runs measure the existing uncommitted change, **not** my C++ variants.

NetGraph rerun, medians in ms, identical recorded workload digests in all phases:

| Workload | Old A1 / A2 | New | New forced heap |
|---|---:|---:|---:|
| Full SPF, Clos 200×200 | 212.167 / 212.584 | 215.551 | 215.680 |
| Full SPF, weighted 100×100 grid | 219.249 / 219.160 | 169.151 | 235.446 |
| Bound-context max-flow, weighted grid combine | 201.715 / 204.482 | 160.344 | 216.545 |
| Backbone Clos placement | 68.333 / 68.736 | 67.171 | 67.445 |

The weighted benefit is independently reproduced. Sub-5% Clos/placement changes are noise under this review's rule. The broader 16-row scenario suite was read but not rerun.

## 3. What the measurements do and do not establish

Already-measured dead ends were not retried: sorting active buckets once, power-of-two masking as a cure for chains, re-hoisting TLS array pointers, re-fixing the function-pointer comparator, and moving growable buffers into locals. None is presented here as a new recommendation.

1. **The strongest original claims are sound.** A matched old/new/forced-heap/old sequence, result checks and fresh processes are much better than comparing the slower experimental heap copy with its sibling bucket implementation. The weighted Core gain and approximately 1.25–1.30× weighted NetGraph gains survived this review. The small-query workspace gain also survives. GPU ten-worker tables and relaxed canonical-DAG prototypes are not comparable replacements; the plan correctly separates them.

2. **A/B/A is not machine isolation.** Original Core load ranged 6.5–12.1, and original NetGraph load 3.9–6.5. Medians from consecutive samples within one process do not remove thermal/frequency, scheduler or allocator drift. My own load logs show 3.5–4.7 during the first two passes and roughly 8 during the third. These are useful large-effect confirmations on a shared machine, not the plan's promised quiet-machine acceptance evidence. Do not derive confidence intervals from nine correlated samples as if they were nine isolated process runs.

3. **The “no unit workload regressed beyond A1/A2 drift” statement is literally too strong.** The recorded backbone placement row is 68.5/67.9 old versus 71.4 new; three-tier placement is 370.7/366.7 versus 381.3 (`CPU_PLAN.md:126–139,150–155`). Both exceed the old bracketing drift while being around the declared noise band. The defensible conclusion is “not a demonstrated regression,” not “within A1/A2 drift.” The first pre-comparator-fix pass is useful context but is not an independent replication of the final binary.

4. **The sweep samples the wrong dimensions for a universal cutoff.** One positive uniform-random cost distribution, one N/degree and full-source queries cannot establish a universal W limit. Add exact boundary widths, stride>1, all-zero costs, sparse/high-diameter graphs, outlier edges, changing masks, first-query allocation and memory high-water behaviour. This review supplies boundary and some safety rows, not that entire performance matrix. Linux/x86 remains open. A masked wide edge can cause a fallback, rather than “merely wasting empty buckets.”

5. **All supplied Core timing rows are warmed.** Expected `run_ref` calls occur before timing (`harness.cpp:263`), followed by two warmups (`:266`). This warms the new TLS workspace before its first measured call. All graph construction is outside timing. The timing loop retains the entire batch's output vectors until checking (`:267–271`), unlike callers that immediately consume/discard each result. This matters for allocation/cache pressure and the large-output conversion variant. Add cold first-call and immediate-consumption measurements rather than replacing the warmed numbers.

6. **The NetGraph workloads are useful software-path coverage, not a representative demand distribution.** They include real APIs, shipped scenarios, pseudo nodes, failures, WCMP/ECMP and KSP. However, the positive result is dominated by a synthetic planar 100×100 grid with costs sampled from 1–20 (`bench_netgraph.py:85–110`). It is not evidence that a typical metric-weighted WAN simulation speeds up 25%; the included NSFNET scenario was already near neutral. Bound contexts are created outside timing (`:174–175,187–188,200–205`), while pairwise grid builds a context inside (`:231–237`); their different gains do not isolate queue effectiveness. Timed Python functions also hash/normalise results and perform NumPy reductions. Serial scenarios/Monte Carlo do not represent parallel worker-pool reuse or oversubscription.

7. **Digest agreement has limits.** SPF checksums contain only a distance sum and parent count (`bench_netgraph.py:178–196`), so they cannot detect reordered/wrong parents, swapped distances or wrong via-edges. `norm` drops fields named `id` recursively (`:38–51`). Only the last inner repetition's checksum is compared (`:376–381`). Keep the full bitgate/differential evidence separate from these weaker end-to-end digests. Hash-seed pinning is a legitimate correction, but does not turn this into a semantic max-flow certificate. The targeted certificate added here provides a small independent check, not complete application acceptance.

8. **Reverse P2 was not finished according to its own gate.** `CPU_PLAN.md:426–431` requires a measured reverse result; the plan elsewhere admits none. This review's reverse measurements show 14.188/13.807 → 8.539 ms on random weighted 4k graphs, and 0.402/0.409 → 0.741 ms on chains. Add these to the evidence rather than labelling an unmeasured stage done.

9. **Rerun priorities:** first the exception reproducer after the production fix; then isolated new/variant/new and old/new/old on dense fabrics, narrow frontiers and wide-cost fallback, with shorter interleaved phase blocks and source hashes. Include cold/warm and immediate-consumption results, realistic caller counts/touched fractions/cost ranges, output size and total per-thread retained bytes. Profile search versus output versus reset to settle the fallback attribution. Finally run the normal project/integration/sanitizer gates after applying any production patches, and CI plus a second CPU architecture. Local measurements and previously passed suites do not replace those checks on the final code.

## 4. Commit recommendation

**Keep the ordered bucket queue, exact pop-order contract and sparse workspace strategy. They are a good streamlined baseline. Do not commit this exact tree unchanged.**

Fix the reproduced allocation-failure contamination first. Correct the memory claim and either bound aggregate retained bucket capacity or explicitly accept and document the actual retention policy. Consider the settled-neighbour skip and one-pass forward output together: their dense-fabric gain is large enough to merit the modest additional review. The inline singleton path and fixed-routing reuse are measured, contained follow-ups; they can be separate commits rather than expanding the queue patch without limit.

Retain the conservative width limit while documenting that its second-platform and cold/memory acceptance remains open. Record the remaining heap penalty as unexplained rather than an unavoidable workspace tax. There is no justification to block this CPU direction on GPU work, a general batch API, a different KSP algorithm, or a wholesale graph-layout redesign. There is also no evidence for saying that all reasonable cheap wins were exhausted.
