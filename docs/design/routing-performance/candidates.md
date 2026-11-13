# Every lead and how it ended

A maintainer asking "should we try X?" finds X here, with the number that decided it and, if it is still open, the
experiment that would decide it. Workloads and terms are defined in the [README](README.md); exactness statements are
in [equivalences](equivalences.md). The literature search is dated 2026-11-10.

## How a lead was judged

- **Bars.** A new mechanism had to be able to remove at least `1 − 1/B` of its row's Core time before it was built,
  where B is the Core-time gain it had to reach: on a primary workload B = 1.15, so 13.04 % (3/23 = 1 − 1/1.15); on a
  secondary row (a harness row, a NetGraph analysis no primary workload runs, or a synthetic workload) B = 1.3, so
  23.08 % (3/13 = 1 − 1/1.3). Above about 1000 product lines the native bar was 1.5× on any row (a screen removing one
  third of its native time); an opt-in candidate that needs a contract decision uses 13 % or 23 % at any size and is
  decided by its consumer comparison, without a native or size bar at landing. A lead moving lock-held Python work into
  Core had to move at least 3 % of the run's wall time. A small local edit was instead built and kept only if its gain
  survived its layout control with no shared-code row slower beyond noise.
- **Raw and corrected.** Shares from instrumented counts are given on two bases: raw, and corrected for the
  instrumentation's own cost. Both must clear the bar. Unless stated otherwise, these cost-screen shares are
  exploratory one-worker measurements, each using its own instrument's raw and corrected Core-time denominator, not
  the replay times or the sampled stack shares of the [README](README.md); the multi-worker shares were not
  measured.
- **"Large"** means a replacement that produces the call's full results, beyond a small edit; the study did not fund
  one once its counts were in.

Outcomes: **applicable** (in place, or proposed); **lost** (built, no gain); **closed** (could not remove enough even if
free, with the number); **does not fit** (a contract clause it cannot keep); **unknown** (with its deciding
experiment).

## Shortest-path search

| Lead | Outcome | Deciding number or experiment |
| --- | --- | --- |
| Core's Dial bucket queue | applicable, in place | any replacement of DC-BB's search removes at most 0.101 of Core time |
| Newer SSSP algorithms (2025-2026, beyond the sorting bound) | lost | both independent implementations are slower than Dijkstra; DC-BB's search 0.101 < 1/3 |
| Removing the per-bucket node sort | lost (built) | command line DC-BB 0.9976 / 0.9998, NSFNET 0.9969 / 0.9970, `backbone_clos` 0.9963 / 0.9966 at 1 / default workers; Core's work 0.9926 / 0.9996 |
| FIFO queue for unit costs | closed on DC-BB (costs {1, 5, 10}); unknown on unit-cost harness `spf` rows | seven `spf.*` rows are 36.0-95.0 % search, free; deciding: a result-producing queue (Large) |
| Wide-cost queues (beyond 65,536 buckets) | unknown | DC-BB and NSFNET never reach that region; `spf.full.wide_fallback16k`, `spf.to.wide_fallback16k`, `spf.to.width65537` do (search 89-95 %); deciding: a result-producing queue (Large) |
| Thorup's 2004 directed integer queue | unknown in the wide-cost region; undirected methods and bottom-up BFS do not fit | directed edges; every equal-cost predecessor |
| Bit lanes in `evaluate_pairs` | closed on NSFNET | a free search is at most 9.73 % of its Core time |
| Failure bit lanes | unknown | on NetGraph's `edge_sensitivity` pair the lane envelope is 16.23-21.95 % raw; native transfer unknown |
| Parallel search inside one query | lost | below about 4,000 nodes and against query-level parallelism; batched GPU unknown, no caller |
| Inline DAG recording | applicable, in place | any exact distances give Core's DAG by one O(m) pass |
| Multi-root SIMD | unknown, no caller | — |
| Many-root bit lanes (one sweep for many roots) | unknown | no caller reads lane masks; worth a look only when one appears |

## Goal-directed search, hierarchies and indexes

| Lead | Outcome | Deciding number or experiment |
| --- | --- | --- |
| A* for targeted `spf` | exact ([equivalences](equivalences.md#early-stopping-searches)); closed as an SPF-only replacement on NetGraph's `paths` analysis | the whole SPF phase there is 2.6196 % raw / 4.1677 % corrected < 23.08 %; other targeted series unknown (deciding: a caller with repeated target sets) |
| ALT landmarks | the two-landmark variant closed on `paths`; others unknown | same ceiling |
| Bidirectional search | exact with the strict rule; closed on `paths` as an SPF-only replacement | same ceiling |
| SPF-only batches sharing endpoints | closed on `paths` (≤ 4.1677 %); KSP sharing across pairs unknown | envelope 62.74-62.98 % raw / 65.99-66.09 % corrected; deciding: result-producing KSP sharing (Large) |
| Contraction hierarchies, manual hierarchies, regions | do not fit | masks per call; every equal-cost route; exactness |
| CCH, PHAST, separators | closed on Clos and as SPF-only on `paths`; as a bound source for KSP and T1 only the replacement of the existing `static_costs` preparation was priced | 13.11-13.28 % raw < 23.08 %; backbone indexes unknown |
| CCH potentials | `spf` part closed; KSP and T1 parts unknown | KSP phase 74.5-74.6 %, T1 phase 66.6-66.7 % raw ceilings; deciding: a priced bound source |
| Distance tables, hub labels, APSP | do not fit as an answer engine; as an SPF distance source closed on `paths`; elsewhere unknown | a backbone workload |
| Distance-only engines (dynamic SSSP, replacement paths, oracles) | do not fit as a drop-in; unknown as distance sources | exact new-state distances plus the O(m) pass would keep every tie; cost unchecked |
| Symmetry (orbits) | promising only on symmetric fabrics for invariant outputs; unknown, no caller | NSFNET is not symmetric and DC-BB's failures break symmetry |

## Failures and reuse

| Lead | Outcome | Deciding number or experiment |
| --- | --- | --- |
| [Baseline route-program reuse](baseline-reuse.md) (class 2) | promising under a condition; needs a contract decision | DC-BB at 1 worker, per-step capture, 0.8171; at 8 workers, per-run capture, 0.9544; deciding: the chosen lifetime faster at 8 and not slower at 1 on the landing tree |
| Reuse when failures miss the support (class 1) | closed on DC-BB and DC-BB 192 pods; applicable as Core's single-edge skip in `evaluate_edge_failures` | 0 of 4,884 failure demand-cases on each |
| Local repair of changed distances (class 3) | closed on DC-BB and DC-BB 192 pods (no class-3 case); unknown elsewhere | a caller with class-3 cases |
| Exact-repeat memoisation | closed on DC-BB | exact repeats ≤ 2.5 % |
| Max-flow batch store | unknown on NSFNET | 39 repeated batches recur across NetGraph's three capacity analyses, none inside one; priced on NSFNET at 1 worker at 13.410787 % raw / 15.870057 % corrected against 13.043478 % (raw margin 2.78 ms in the limiting process; the five processes' Core time was 756.801607-787.836773 ms raw, 633.854833-659.865019 ms corrected), retention and allocator effects unmeasured; deciding: a store kept within one NetGraph run spanning all three analyses, timed on the command line (Large) |
| Batch bounds across failure versions | closed on DC-BB (every failure is class 2); unknown on batches with overlapping versions | — |
| Lazy LOSSY re-forwarding | closed on `program.stale_then_place.pooled` (5.46-6.75 % raw); unknown on `kernel.lossy.clos` (20.56-22.18 % raw) and for LOSSLESS | native transfer of the envelope |
| Symbolic verifiers | do not fit | per-case integer results |

## K shortest paths and tunnel setup

| Lead | Outcome | Deciding number or experiment |
| --- | --- | --- |
| Yen with Lawler's partition | applicable, in place | — |
| CSPF as routers and PCEs run it (prune below the request, then a shortest-path search per request) | applicable: keep | resource-constrained shortest paths with a second additive constraint are outside the contract |
| A* spur searches for `ksp` | exact with strict bounds and the `FEWER_HOPS` suffix test; unknown | must remove 30.94-30.96 % raw / 29.40-29.53 % corrected of the KSP phase on the synthetic backbone at 1 worker's `paths` queries; deciding: a C++ spur engine (Large) |
| KSP spur variants; sBELA* (Linares López, Herman, IJCAI 2026) | unknown; loopy KSP does not fit | Core's full order key, EdgeId ties and parallel identity are not proved for sBELA*; deciding: that check, then its cost on the saved queries |
| Repeated-k and shared-endpoint path queries | unknown | no growing-k series exists |
| Static potential and dead-node memory for tunnel setup | exact per request value; unknown | must remove 34.58-34.67 % raw / 33.50-33.58 % corrected of the T1 setup phase on the synthetic backbone at 1 worker with tunnels; deciding: a CSPF engine (Large) |
| Decremental tunnel distances | unknown | after the tunnel-setup engine |
| Tunnel-setup reuse across failure cases, narrow form (raw setup-body reuse) | exact under the four conditions of [equivalences](equivalences.md#reuse-of-the-failure-free-te-solve-and-tunnel-setup); closed | on the synthetic backbone with tunnels at 1 worker: 22.37 % raw / 23.05 % corrected < 23.08 %; a wider form that also skips input packing and reconstruction is unknown |

## Exact flows

| Lead | Outcome | Deciding number or experiment |
| --- | --- | --- |
| Successive shortest paths with potentials and blocking flows | applicable: keep the solver | the only family that yields the marginal cost curve without extra solves in 64-bit integers |
| Sink stop in the residual searches | exact ([equivalences](equivalences.md#sink-stop-in-the-flow-solver)); closed on NSFNET; unknown on `backbone_clos` | at 1 worker: NSFNET, a free Dijkstra plus potential updates is 61.886599 ms, 7.2297 % raw / 8.3567 % corrected; `backbone_clos`, ceiling 20.77-20.96 % raw / 21.49-21.67 % corrected, above the bar, net unpriced; deciding: a result-producing solver (Large) |
| Reachability-based controller roster | lost (built) | `backbone_clos` at 1 worker 0.9833 [0.9806, 0.9846], a small real effect under the 3 % bar; `backbone_clos` at 28 workers and NSFNET at 8 are within noise; NSFNET at 1 is unqualified because its reference A/A drifted, so neither a gain nor neutrality is established there |
| Practical max-flow engines (push-relabel, pseudoflow, Boykov-Kolmogorov, EIBFS) in place of Core's Dinic | lost in NSFNET's isolated max-flow phase; unknown on harness T2 rows | no published comparison on network graphs |
| Approximate flow drop-ins (interior point, almost-linear, LP) | do not fit and lose at Core's sizes | exact integer results |
| Exact combinatorial min-cost flow (Haeupler, Jiang, Saranurak 2026) | unknown | no implementation; curve and Core's allocation not established |
| Tunnel-count sizing by cut steps or other flow kernels | does not fit as the default; unknown on harness T2 rows | 29.41-68.79 % envelopes; deciding: a result-producing solver (Large) |
| Cost-free deletion loss | exact; closed on the NSFNET pair (first repair 2.94-3.47 % raw); `sens.exact*` rows unknown | a result-producing change (Large) |
| Reduced-result certificates for failure pairs (opt-in) | unknown, not pursued | On NSFNET at 1 worker, a partial kernel-only ledger: 4.11 / 4.44 ms short of the bar raw, 10.95 / 10.61 ms above it corrected; deciding: a remeasurement with capture, tests and reconstruction priced |
| Solver-only reuse of a certified prefix of demands | closed on `backbone_clos` at 1 worker (6.6384 % / 10.7464 %) and the synthetic backbone at 1 worker (14.9779 % / 18.0515 %); unknown at their default workers | — |
| [Max-flow pair reduction](equivalences.md#the-reverse-pair) (opt-in) | exact for the eleven exported fields only, with the cut derived only for a solve that reached the maximum with the cut requested (`None` otherwise, `[]` for an empty roster); unknown on NSFNET | NSFNET at 1 worker: partial 14.283527 % raw / 16.884134 % corrected, the rest unpriced; deciding: an opt-in evaluator timed on the command line (Large) |
| Cut trees, quotients, undirected and one-path indexes | do not fit | curves, cuts, per-edge results, every equal-cost route, exact integers |

## Placement, admission and the ledger

| Lead | Outcome | Deciding number or experiment |
| --- | --- | --- |
| [Ordered publication reuses its totals](ordered-publication.md) | applicable, proposed | DC-BB at 1 worker 0.9624 [0.9607, 0.9662], layout 0.9996; at 8 workers 0.9884 (within noise) |
| Evaluation-only failure lists | closed | 0.088835 of Core time < 13.04 % |
| Admission with fewer sweeps (window, lanes, delta allocation) | unknown, not pursued | admission is 25.38 % of DC-BB's one-worker Core time (a counted phase of an instrumented copy over its own 9.74 s Core time, not the sampled 25.8 % of the README) and 51.40 % of it must go; deciding: a costed admission scheme |
| Clearing only written slots | lost (built) | DC-BB at 1 worker 1.0324 [1.0277, 1.0353], borderline slower; at 8 workers 1.0117 |
| T_max without search | its harness rows lose or have nothing to remove; `kernel.weighted.clos` below the bar nominally (42.58-42.98 % raw / 19.33-22.51 % corrected), native unknown | exact only on forests and in-trees |
| Fewer passes over support and programs | unknown, not pursued | the phase is 38.68 % of DC-BB's one-worker Core time (3.767 s of the same instrumented copy's 9.74 s, not the sampled 39.9 % of the README); a small mechanism must remove 33.72 % of it, a large one 86.18 %; no mechanism with a cost was identified |
| Exact divisor (highest-averages) selection for wide weighted splits | above 16 members `kernel.weighted.clos` closed (10.96-11.57 % raw); `pooled.weight.explicit` and `.snapshot` unknown | the published linear-time code compares quotients in floating point and can return wrong seats; deciding: an exact port with integer comparisons, Core's lower-index tie rule and a near-tie plant (Large) |

## Engineering and the application boundary

| Lead | Outcome | Deciding number or experiment |
| --- | --- | --- |
| One native call per failure batch | unknown, not pursued | the identified lock-held caller work (0.259 s, 2.25 % of an 11.46 s diagnostic run at 8 workers) plus the attached binding work inside Core spans (0.189 s) is 0.447 s, 3.90 % (3.59 % on the sensitivity basis), above the 3 % bar: a batch reaches it if it removes at least 45 % (64 %) of the attached work; two instrumented runs, the batch's own conversion unpriced; deciding: a Core batch entry point with a NetGraph adapter timed on DC-BB at 8 workers (Large) |
| Persistent worker pools | unknown, not pursued | NSFNET starts 928 threads in 116 batches; `backbone_clos` already uses one pool; deciding: dispatch, reset and idle-memory cost against today's |
| Columnar output form | unknown, not pursued | garbage-collection attribution moves from 286.9 to 11.1 ms under observer changes: no causal figure; needs a new output form, a contract decision, first |
| Vertex reordering and other layouts | lost for DC-BB's search; other layouts unknown | Core's graphs fit in one L2 cache |
| Whole-library profile-guided optimisation | lost | a large loss on untrained rows |
| Free-threaded Python | outside Core | conditional on the lock being the critical path, not shown |
| Oracle libraries (NetworkX, Boost) | applicable as test oracles; unknown as speed references | — |
| The replay of recorded Core calls | applicable, done | decides only DC-BB rows whose call kinds it reproduces ([measuring](measuring.md#what-each-instrument-measured-and-its-limits)) |

## The literature that matters

| Source | Fits Core's inputs? |
| --- | --- |
| Dial's buckets; Meyer-Sanders (Δ-stepping with Δ = 1 is Dial) | yes, in place: linear in nodes, edges and the largest distance for Core's integer costs |
| Duan, Mao, Mao, Shu, Yin (STOC 2025, JACM 2026); Duan, Mao, Shu, Yin (ICALP 2026) | no: the sorting term they beat is one the buckets already avoid; one tree per search; the independent implementations lose about 2× (Rust, Valko et al.) to 2.8-27.5× (C++, Castro et al.) to a binary-heap Dijkstra |
| Thorup 1999; Thorup 2004 (JCSS 69(3)) | 1999: undirected, no. 2004: directed integer SSSP, unknown only beyond 65,536 buckets |
| Haeupler et al., universal optimality (FOCS 2024) | a property of Dijkstra with a working-set heap; integer buckets lie outside its comparison model |
| Goldberg-Harrelson (ALT); Strasser-Zeitz | a reused reverse distance is a feasible potential under masks and gates; exact for targeted `spf` with the ≤ d* stop |
| Contraction hierarchies, CCH, CRP, PHAST, hub labels | no on fabrics: every balanced separator of a fat tree has Θ(n) switches, so a third of the nodes search a clique of that size; labels unknown on Clos |
| Ramalingam-Reps (1996) | yes: keeps every tight in-edge; its first phase is the exact class-2 test |
| Tree-based dynamic SSSP, router iSPF, replacement paths, sensitivity oracles | not drop-ins: one parent or distances only; mostly unique shortest paths |
| Bilò, Chechik, Choudhary, Cohen, Schirneck (ICALP 2026) replacement-path coverings | unknown as an auxiliary mechanism for directed multi-failure cases, not an ECMP engine |
| Yen; Lawler's partition | yes, in place |
| Feng; Martins-Pascoal (canonical suffix); PeeK (K-upper bound) | exact only under `FEWER_HOPS` with the suffix validity check; bounds keep ties only when strict |
| D'Emidio, Di Stefano (ALENEX 2026); UPD (2610.20062v2) | unknown: they order by cost alone; Core's order and parallel-edge identity unproved |
| US 8,854,956 B2 | prior art for the tunnel-setup potential, without ties or measurements |
| Kovács (thesis 2019) | supports Core's flow family: fastest on networks with few terminals, and the only one giving the curve in int64 |
| Network simplex, cost scaling | no: another optimal allocation, no curve, costs beyond int64 |
| Schätzle, Meierhans, Kyng, Probst Gutenberg (2610.18328v1, practical almost-linear min-cost flow) | no: floating point, heuristic, two non-optimal outputs, slower than LEMON |
| Haeupler, Jiang, Saranurak (2610.32556v1) | theory only: no implementation; curve and allocation not established |
| Fortz-Thorup lazy ECMP re-forwarding | only for placement at a fixed offer; DC-BB is LOSSLESS |
| Reitzig-Wild SandwichSelect (Mathematical Programming 2024) | unknown for wide WCMP; the published code compares in floating point and can return wrong seats |
| Velo (NSDI 2025), Hoyan2 (SIGCOMM 2025), NetDice, Jingubang | confirm per-case failure enumeration with a support-based skip as current practice |
| pybind11 3.0.2 and 3.1.0; PEP 779 | lower binding overhead; free threading is a NetGraph-side option |
