# Reference-wheel native attribution

The new captures use extension SHA-256
`208487c7aed158fb418791c22e23d09332533595053c0e4810923b31244f7c84`.
ORIENT-ADD's old PC ranges belong to another binary and are not applied.
`analyze_sample.py` imports only its `raw_pc` syntax reader and tree parser.

Disassembly was generated read-only with:

```bash
/usr/bin/otool -tV /Users/networmix/ws/work/2026-11-core-routing-research.noindex/scratch/consumer/ref/venv/lib/python3.14/site-packages/_netgraph_core.cpython-314-darwin.so > /Users/networmix/ws/work/2026-11-core-routing-research.noindex/scratch/clocks/reference-arm64.txt
```

Source is `worktrees/core-ref/src/min_cost_flow.cpp`, clean HEAD `9cd7642`.
The symbol `detail::solve_min_cost_flow` starts at `0xc251c`; its call at
`0xc2d68` enters `0xc4e24`. Arguments set closure and potentials true, matching
`MinCostFlowSolver::run` at source lines 275–307. The latter is mostly inlined.

| Half-open reference PC range | Phase | Separating instruction/source witnesses |
| --- | --- | --- |
| `[0xc4ef0,0xc51cc)` | Residual Dijkstra initialization and queue setup | Loads network arc offsets/destinations/costs; `0xc4f08` increments dijkstra_runs; `0xc4f28` fills distances with max int64; source lines 446–471. |
| `[0xc5cc8,0xc6140)` | Residual Dijkstra queue processing and relaxations | Bucket minimum scan/migration at `0xc5d18`–`0xc5df0`; return to reached-target read `0xc51c0`; settled-vector growth at `0xc60c0`–`0xc613c`; source lines 472–498. Excludes marginal-curve growth beginning at `0xc6140`. |
| `[0xc51cc,0xc52e4)` | Potentials / settle | `0xc51f4` adds dist[v] to pi[v], `0xc5200` writes ever_settled[v]; subsequent shadow checkpoint update follows source lines 504–518. |
| `[0xc5338,0xc5c44)` | Blocking flow | Reach generation at `0xc5368`; admissibility comparison `0xc5494`; Dinic phase increment `0xc5588`; DFS stack read `0xc5870`; forward/reverse residual updates `0xc5bdc`/`0xc5bf0`; source lines 533–641. `0xc5c44` begins post-drain progress/curve handling and is excluded. |

The checker asserts twelve actual instructions, including the enclosing call,
before applying these ranges. Controls distinguish a Dijkstra PC, a blocking-flow
PC, an explicit wait beneath a native ancestor, and malformed tree weights.
These controls support attribution logic; only coordinator captures supply sample
counts. Source and complete disassembly remain available for independent review.

These are partial ranges, not a promise that every solver instruction is named.
Tier-cost calculation, curve materialization, closure/cut, solver setup/output,
unrecognized stripped evaluator-worker/binding frames and cold branches remain
unresolved unless a named stack identifies them. Named roster searches remain
separate from solver searches. No unresolved bucket is counted as zero work.

Every root is retained and its self weights conserved. Explicit blocking leaves
take priority over a Core ancestor. All-thread sample counts are not CPU time;
shares must name whether their denominator includes blocked, unresolved and
outside-Core stacks. Multi-worker shares are not scaled by caller wall into
microseconds per pair.
