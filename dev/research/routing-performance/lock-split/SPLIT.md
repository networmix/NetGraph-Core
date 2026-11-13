# SPLIT — held, released and retake costs; P1 at eight workers

In these instruments, the attached caller budget is small for P1, but material for P2/P3.
The six lock splits and P1/8's all-thread phase shares are measured; all exports match.

**Status: measured-exploratory.** Six instrumented processes and one 15 s capture, run by the
coordinator unconfined through the heavy slot on 2026-11-10; no A/B/A, significance, speedup
or removable-work claim. P1=dcbb96 (1/8 workers), P2=nsfnet (1/8), P3=backbone (1/28).
Inputs and NetGraph `8a641a1` are pinned. Core is `9cd7642`, CPython 3.14.10 with GIL.
Reference extension SHA-256 `208487c7aed158fb418791c22e23d09332533595053c0e4810923b31244f7c84`;
scratch instrumented extension `4df35c735285eabd0797daf2e2bd1081f85c860087f33d2ef01f0f3f30021626`.
Only the task's scratch bindings changed.

H = attached residual, R = released interval, A = retake interval; W = wall, C = caller thread
CPU. `H_W = span_W − R_W − A_W`, and independently `H_C = span_C − R_C − A_C`.
H includes argument/result conversion, validation, binding/runtime and release-edge work;
H_W also includes descheduling while attached. A_W brackets RestoreThread, including attach
work: it is an upper proxy for lock waiting, not a pure GIL-wait measurement. R_W can include
native joins/waits. Core's own workers are not in R_C. Small A_C/A_W values also differ because
of read ordering; do not compare them as perfectly coincident physical intervals.
Production releases use `Detached` in native_call.hpp (native_call/argument callbacks and
FlowPolicy's Python-holder deleter). Private probe gil_scoped_release/call_guard sites are
excluded from the wheel. Reads are W/C after SaveThread, C/W before RestoreThread, W/C
afterwards. Own-thread cumulative TLS snapshots surround T4's outer call clocks; invalid
flags, counts and thread entry/exit snapshots are retained. Full release-site inventory:
`tasks/split/stage1-evidence.json` and `SCOPE.md`.

**Analysis (`ngraph.run`), seconds.** Triples are H / R / A after approximate reader correction.
Tables round to 1 us; closure uses unrounded records.
Raw parts/endpoints remain in analysis.json, report-data-01.json and the arrays. Plain/clock
wall illustrates disturbance, not a performance comparison. Top-five call details follow in the appendix.

| Row | Run wall plain / clock | Raw caller-span W / C | Corrected W: H / R / A | Corrected C: H / R / A | Non-Python CPU bracket, raw |
| --- | ---: | ---: | ---: | ---: | ---: |
| P1/1 | 10.526 / 10.857 | 9.055806 / 9.046343 | 0.112769 / 8.917142 / 0.000945 | 0.111787 / 8.915920 / 0.001468 | [-0.000003; 0.000003] |
| P1/8 | 2.616 / 3.956 | 19.189208 / 10.776375 | 0.194289 / 10.471664 / 8.498427 | 0.188804 / 10.465430 / 0.104928 | [0.006606; 0.006612] |
| P2/1 | 0.667 / 5.331 | 0.738502 / 0.647819 | 0.191888 / 0.363001 / 0.000031 | 0.190233 / 0.362698 / 0.000049 | [-0.000002; 0.000002] |
| P2/8 | 0.499 / 5.251 | 0.568232 / 0.337224 | 0.230215 / 0.157522 / 0.000053 | 0.226192 / 0.016967 / 0.000068 | [1.149650; 1.149654] |
| P3/1 | 0.208 / 0.672 | 0.201147 / 0.191808 | 0.020336 / 0.156430 / 0.000337 | 0.019589 / 0.156328 / 0.000383 | [-0.000004; 0.000002] |
| P3/28 | 0.160 / 5.742 | 55.968703 / 0.345560 | 0.140633 / 0.196652 / 55.606795 | 0.056060 / 0.196035 / 0.077416 | [0.001260; 0.001264] |

P1 has 54,506 analysis spans, 29,572 detach intervals and 4,922 Python-holder destructor
intervals at either worker count; P2: 614,480 / 348 / 0; P3: 60,723 / 20,196 / 3,558.
P2's 116 analysis evaluate_pairs calls launch/join native workers at workers=8; its residual
is predominantly Core-worker plus pthread/runtime CPU (T4 census/native-stack evidence), not
solver-only CPU. P1/8's 6.61 ms and P3/28's 1.26 ms residuals include Python startup/teardown
tails and remain unattributed. Brackets are endpoint-read brackets, not CIs. Registered analysis
**lifetime** threads are 1/140, 1/1, 1/29 for P1/P2/P3, not simultaneous worker counts.

**All six phases**, import/load/prepare/run/export/cleanup, seconds, same corrections.
These totals exclude initial Core/NumPy imports, calibration and inter-phase driver work;
they are not whole CLI timings. C remains caller CPU; native workers are separate above.
Cleanup includes model teardown/gc.collect inside its span.

| Row | Calls / detaches / destructors | Corrected W: H / R / A | Corrected C: H / R / A |
| --- | ---: | ---: | ---: |
| P1/1 | 220,787 / 29,577 / 4,922 | 0.183814 / 8.918099 / 0.000946 | 0.180905 / 8.916877 / 0.001468 |
| P1/8 | 220,787 / 29,577 / 4,922 | 0.260574 / 10.472600 / 8.498428 | 0.253107 / 10.466366 / 0.104929 |
| P2/1 | 1,930,113 / 354 / 0 | 0.500806 / 0.363038 / 0.000031 | 0.497609 / 0.362736 / 0.000049 |
| P2/8 | 1,930,113 / 354 / 0 | 0.543203 / 0.157558 / 0.000053 | 0.530286 / 0.017003 / 0.000068 |
| P3/1 | 115,810 / 20,201 / 3,558 | 0.037236 / 0.156458 / 0.000337 | 0.036238 / 0.156356 / 0.000383 |
| P3/28 | 115,810 / 20,201 / 3,558 | 0.156910 / 0.196683 / 55.606795 | 0.072134 / 0.196066 / 0.077417 |

P2/8's corrected held CPU includes 0.226192 s in run and 0.272871 s in export;
the 0.031044 s cleanup span includes whole-model teardown/GC, not just conversion. Operators such as
hash/equality/repr/copy outside T4's intercepted API remain in outside time; no out-of-span
Detached interval was found in any measured phase.

**P1 sampled phases.** Exclusive weights are inclusive minus immediate children, assigned once
by Core ancestry; explicit blocked leaves win. All 541 observed roots, each with Core context,
are retained. Five workflows were retained/exported after capture; every export equals plain
P1/8. Retention changes lifetimes/memory. P1/1 uses wheel `6db37900…`, one main-thread root and
all 9,134 Core samples. P1/8 uses the new reference and 34,503 Core samples **excluding** 1,281
blocked leaves. Different binaries/coverage/denominators prevent a worker-effect inference;
these are stack shares, not CPU shares.

| Exclusive category | P1/1 carried count / % | P1/8 count / % |
| --- | ---: | ---: |
| Support / program construction | 3,644 / 39.89% | 13,732 / 39.80% |
| Placement | 2,353 / 25.76% | 9,045 / 26.22% |
| Ledger / result materialisation | 1,259 / 13.78% | 4,484 / 13.00% |
| Search | 920 / 10.07% | 3,234 / 9.37% |
| Unresolved native / bindings | 523 / 5.73% | 2,226 / 6.45% |
| Setup / masks / cleanup | 435 / 4.76% | 1,782 / 5.16% |

All self samples: 46,389; Core including waits: 35,784; outside: 10,605. Recognized blocked
waits: 6,815 (14.69% of all), 1,281 inside Core and 5,534 outside. Nonblocked outside: 4,810
Python/runtime + 261 readiness reads. Unresolved Core: 2,226; unknown waits can remain there.
Ten named top-of-stack counts cross-check the weights. Reference static_support return offsets
+140/+180 are revalidated; stripped target names remain a source-order/call-context inference.
No old absolute PC ranges are reused. Use **p1-8-sample-analysis-02**: -01 failed because the
coordinator wrapper contained Markdown fences; the capture was not repeated.

**Closure and overhead.** Independent wall/CPU interval complements give, per phase/thread,
`outside + H + R + A = lifetime`, with **0 ns integer closure error**. TLS totals reconcile;
invalid clocks, outside detaches and negative raw held intervals are absent. Corrected parts,
corrected outside and overhead also close (<10 ns in unrounded JSON). This proves
accounting closure, not clock accuracy or exact GIL transition instants.

| Row | Inner-null W / C, ns/call | Full-wrap CPU, us/call | Concurrent outside correction, us/call | Run span reader correction W / C, ms |
| --- | ---: | ---: | ---: | ---: |
| P1/1 | 294.9 / 152.2 | 6.645 | not used | 24.950 / 17.168 |
| P1/8 | 293.6 / 153.9 | 6.668 | 26.39–30.13 | 24.828 / 17.213 |
| P2/1 | 298.6 / 154.2 | 6.718 | not used | 183.582 / 94.838 |
| P2/8 | 293.5 / 152.8 | 6.709 | not used | 180.442 / 93.997 |
| P3/1 | 296.4 / 155.8 | 6.682 | not used | 24.044 / 15.509 |
| P3/28 | 305.4 / 164.2 | 6.691 | 80.16–97.56 | 24.623 / 16.048 |

Native clock-only sequences cost median 298–301 ns/detach. Per detach their released brackets
contribute about 186–188 ns wall / 86–87 ns CPU and retake brackets about 12 ns wall /
111–113 ns CPU; the remainder is charged to held-reader cost. Subtract these from R/A/H and
one inner-null span from H per public call. Charge full-wrapper minus inner-null outside,
never again inside spans or as a sum subtracted from global elapsed wall. Null costs include
call/loop effects and do not undo scheduling, retention or registry-lock distortion. The APIs
advertise 1 ns resolution; that is not effective accuracy. Signed corrected values are kept.
P2/8 has **one** Python caller: its saved eight-Python-thread null is inapplicable and ignored.
P1/8 and P3/28 concurrent controls are sensitivities, not guaranteed error bounds or serial
floors; P3/28's outside-CPU sensitivity even crosses zero. Raw values remain authoritative
accounting observations. In particular P3/28's raw retake sum is 55.607040 s across overlapping
callers, while instrumented run wall is 5.742 s versus plain control 0.160 s. This is substantial
instrumentation disturbance, not a measurement of the plain product's GIL bottleneck.

**Large or small at scenario workers?** The measured attached binding/runtime budget is small
for P1/8 and cannot be called negligible for P2/8 or P3/28. Corrected H_C is 1.75%, 93.00% and
17.01% of their respective **caller-span CPU**. P2's 93% omits its native workers. As an explicitly
conditional CPU proxy, divide H_C by H_C+R_C, add the P2 native-worker residual to that denominator,
and exclude retake CPU: raw→corrected fractions are about 1.9→1.8%, 21.5→16.2%, 25.6→22.2%.
Thus P2/P3 have a material attached budget in these instruments; P2's visible getters account
for much of its held CPU, while P1's detached computation dominates. Pure conversion is not separated from
validation/runtime or release-edge cost, and correction is approximate; these fractions are
neither plain-wall shares nor attainable gains. No conclusion about a finite free-Core endpoint
or a product optimization follows from them.

Evidence root S = `/Users/networmix/ws/work/2026-11-core-routing-research.noindex/scratch/split/`:
S/measured-{p1-1,p1-8,p2-1,p2-8,p3-1,p3-28}-01/{plain,clocks}/ and analysis.json contain all raw/
corrected phase/call/thread data and calibration blocks; captures are S/p1-8-sample-01/ and
S/p1-8-sample-analysis-02/. S/report-data-01.json/log independently checks exports, raw sums,
closure, sample weights/top counts and unchanged instruments. Task scripts/patch/commands are
in `/Users/networmix/ws/work/2026-11-core-routing-research.noindex/tasks/split/`; reproduce the
read-only check with reference Python `-B tasks/split/summarize.py <new-output.json>`. No workload
was rerun. This task claims no product suite, sanitizer, CI or publication acceptance.

**Call-detail appendix: five largest analysis kinds per row**, ranked by **raw span wall**, seconds.
Corrected triples use the same H/R/A order. Long retake waits can determine the ranking.
CoreObjects.release is a scoped container/decref span, not one destructor; P3's 1,186 such
spans contain 3,558 DetachedDelete intervals. Rest is an explicit subtotal for closure.

| Row / call kind | Calls | Raw span wall | Corrected W: H / R / A | Corrected C: H / R / A |
| --- | ---: | ---: | ---: | ---: |
| P1/1 `apply_actions` | 2,475 | 6.597210 | 0.007577 / 6.587865 / 0.000295 | 0.007342 / 6.587208 / 0.000565 |
| P1/1 `FlowPolicy.__init__` | 4,922 | 1.902971 | 0.021528 / 1.878343 / 0.000171 | 0.021363 / 1.878006 / 0.000312 |
| P1/1 `FlowLedger.__init__` | 2,461 | 0.276518 | 0.045989 / 0.229002 / 0.000062 | 0.045983 / 0.228890 / 0.000104 |
| P1/1 `Eligibility.intersect` | 2,460 | 0.127337 | 0.001491 / 0.124310 / 0.000072 | 0.001410 / 0.124267 / 0.000105 |
| P1/1 `FlowPolicy.remove` | 9,872 | 0.097824 | 0.002148 / 0.089596 / 0.000206 | 0.001989 / 0.089556 / 0.000233 |
| P1/1 rest | 32,316 | 0.053946 | 0.034035 / 0.008025 / 0.000139 | 0.033699 / 0.007994 / 0.000148 |
| P1/8 `apply_actions` | 2,475 | 8.773374 | 0.014449 / 7.727795 / 1.029665 | 0.013872 / 7.724140 / 0.013947 |
| P1/8 `FlowPolicy.__init__` | 4,922 | 3.966751 | 0.038402 / 2.174787 / 1.750647 | 0.037222 / 2.173324 / 0.027956 |
| P1/8 `FlowPolicy.remove` | 9,872 | 2.703630 | 0.016903 / 0.131775 / 2.549107 | 0.015592 / 0.131396 / 0.025561 |
| P1/8 `FlowLedger.__init__` | 2,461 | 1.431698 | 0.063775 / 0.275554 / 1.090911 | 0.063177 / 0.275294 / 0.013539 |
| P1/8 `Eligibility.intersect` | 2,460 | 1.370958 | 0.006658 / 0.142307 / 1.220537 | 0.006110 / 0.142075 / 0.013287 |
| P1/8 rest | 32,316 | 0.942798 | 0.054103 / 0.019445 / 0.857559 | 0.052830 / 0.019201 / 0.010638 |
| P2/1 `evaluate_pairs` | 116 | 0.399526 | 0.036582 / 0.362855 / 0.000019 | 0.036573 / 0.362555 / 0.000037 |
| P2/1 `FlowPolicyConfig.routing` | 175,514 | 0.135545 | 0.083138 / 0.000000 / 0.000000 | 0.082465 / 0.000000 / 0.000000 |
| P2/1 `CentralizedTE.source_allocation` | 87,699 | 0.043488 | 0.017302 / 0.000000 / 0.000000 | 0.017092 / 0.000000 / 0.000000 |
| P2/1 `PlacementResult.program` | 87,696 | 0.037044 | 0.010859 / 0.000000 / 0.000000 | 0.010636 / 0.000000 / 0.000000 |
| P2/1 `PlacementResult.marginal_cost_curve` | 87,696 | 0.037042 | 0.010857 / 0.000000 / 0.000000 | 0.010641 / 0.000000 / 0.000000 |
| P2/1 rest | 175,759 | 0.085857 | 0.033151 / 0.000145 / 0.000011 | 0.032828 / 0.000143 / 0.000012 |
| P2/8 `evaluate_pairs` | 116 | 0.196040 | 0.038591 / 0.157341 / 0.000039 | 0.038574 / 0.016789 / 0.000052 |
| P2/8 `FlowPolicyConfig.routing` | 175,514 | 0.151981 | 0.100471 / 0.000000 / 0.000000 | 0.099202 / 0.000000 / 0.000000 |
| P2/8 `PlacementResult.marginal_cost_curve` | 87,696 | 0.048594 | 0.022857 / 0.000000 / 0.000000 | 0.022261 / 0.000000 / 0.000000 |
| P2/8 `CentralizedTE.source_allocation` | 87,699 | 0.044232 | 0.018494 / 0.000000 / 0.000000 | 0.017942 / 0.000000 / 0.000000 |
| P2/8 `PlacementResult.program` | 87,696 | 0.039171 | 0.013434 / 0.000000 / 0.000000 | 0.012848 / 0.000000 / 0.000000 |
| P2/8 rest | 175,759 | 0.088214 | 0.036368 / 0.000181 / 0.000014 | 0.035365 / 0.000178 / 0.000016 |
| P3/1 `apply_actions` | 599 | 0.140051 | 0.001871 / 0.137796 / 0.000027 | 0.001841 / 0.137768 / 0.000034 |
| P3/1 `FlowPolicy.__init__` | 3,558 | 0.020914 | 0.003796 / 0.014940 / 0.000058 | 0.003694 / 0.014907 / 0.000080 |
| P3/1 `FlowPolicy.remove` | 7,152 | 0.006764 | 0.000994 / 0.001397 / 0.000112 | 0.000928 / 0.001380 / 0.000120 |
| P3/1 `CoreObjects.release` | 1,186 | 0.005674 | 0.003768 / 0.000446 / 0.000043 | 0.003738 / 0.000437 / 0.000046 |
| P3/1 `FlowPolicyConfig.routing` | 7,111 | 0.004087 | 0.001980 / 0.000000 / 0.000000 | 0.001866 / 0.000000 / 0.000000 |
| P3/1 rest | 41,117 | 0.023656 | 0.007927 / 0.001851 / 0.000096 | 0.007522 / 0.001837 / 0.000102 |
| P3/28 `FlowPolicy.__init__` | 3,558 | 32.955777 | 0.011771 / 0.022290 / 32.919559 | 0.011104 / 0.022166 / 0.047759 |
| P3/28 `apply_actions` | 599 | 11.281891 | 0.003842 / 0.165876 / 11.111811 | 0.003747 / 0.165835 / 0.012681 |
| P3/28 `CoreObjects.release` | 1,186 | 10.567155 | 0.013322 / 0.001383 / 10.551018 | 0.012711 / 0.001256 / 0.014452 |
| P3/28 `FlowLedger.__init__` | 593 | 0.572457 | 0.001449 / 0.001695 / 0.568953 | 0.001403 / 0.001675 / 0.001024 |
| P3/28 `FlowPolicy.remove` | 7,152 | 0.397211 | 0.082092 / 0.003681 / 0.307102 | 0.007861 / 0.003498 / 0.000887 |
| P3/28 rest | 47,635 | 0.194212 | 0.028156 / 0.001728 / 0.148353 | 0.019234 / 0.001605 / 0.000613 |

