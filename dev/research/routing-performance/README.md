# Routing-performance study: the measuring instruments, as used

These are the scripts and reports of the instruments behind `docs/design/routing-performance/` (study of 2026-11-08 to
2026-11-10, Core `9cd7642`, NetGraph `8a641a1`), kept on this branch so the next measurement can start from them. They
are not maintained and not part of the product: they are copied as they ran. Absolute paths in them point into the
study's work directory, which no longer exists; its raw records, builds and captures were not kept. Each instrument
needs its paths and environment set up again before it runs; `docs/design/routing-performance/measuring.md` gives the
reference build recipe and the protocol.

| Directory | What it is |
| --- | --- |
| `lock-split/` | the lock split: `split.patch` puts wall and thread-CPU clocks around the binding's lock release and retake (`bindings/python/native_call.hpp`); `measure.py`, `run_workload.py`, `collector.py`, `analyze.py`, `calibrate.py` run a workload on that build and split each Core span into held, released and retake parts; `capture.py` and `analyze_sample.py` take and classify a sampled capture; `summarize.py` rebuilds the report's tables; `SPLIT.md` is the report, `SCOPE.md` and `HANDOFF.md` the method and the exact commands |
| `call-clocks/` | per-call wall and CPU clocks around every Core call, destructor and getter in NetGraph (`collector.py`, `run_workload.py`), their calibration and overhead controls, the added-latency sensitivity build (`delay.patch`, `sensitivity.py`) and the reports `E5.md` and `E6.md` (in Russian) |
| `replay/` | the recorder and replay of a workload's Core calls with check and timed modes, and its reports |
| `profiling/` | whole-run phase timing, native sampling and the classifier of sampled stacks into Core's phases (`analyze_native_sample.py`), lifetime and getter counters |
| `inputs/` | `gen_backbone.py`, the generator of the synthetic backbone, and `manifest.json`, the hashes of the three workload files |

The prototypes these instruments measured are the branches `prototype/ordered-publication` and
`prototype/baseline-reuse`.
