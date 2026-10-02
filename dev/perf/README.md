# Performance tooling

Developer-only tools for measuring NetGraph-Core. Nothing here ships in the
wheel or the sdist, and nothing here is run by CI; the repository gates
(`make check-ci`, `make sanitize-test`) remain the acceptance tests.

- [`spf_queue/`](spf_queue/): the shortest-path frontier queue and workspace
  work (plan, parity harness, old-versus-new bit-gate, NetGraph-level
  benchmark, recorded A/B/A results and the independent review). Start with
  [`spf_queue/PLAN.md`](spf_queue/PLAN.md).
- [`../benchmark_profiling_overhead.py`](../benchmark_profiling_overhead.py):
  SPF and placement micro-benchmark used for the PGO training run
  (`make install-pgo`) and for measuring profiling overhead.

Conventions for performance claims (see `AGENTS.md`): verify correctness
first (bit-identical outputs where the contract is exact, certificates for
flows), then measure old/new/old in fresh processes on one quiet machine,
record load, hashes and raw samples, and treat differences inside the
run-to-run band as inconclusive. The research record behind the queue work,
including the Apple GPU investigation that preceded it, is in
[`../research/`](../research/).
