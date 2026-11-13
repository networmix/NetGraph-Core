Current disposition: superseded by ruling 18; window-02 stands. See E0.md/E7.md for final verdicts.

# Window-02: qualification stops at the machine-state condition

The standard analysis completed successfully, but its statistical flags do not establish the
quiet-machine condition throughout the window. A non-system renderer exceeded the allowed
half-core threshold before a measured block. The code checks background CPU only at the pilot's
opening; later snapshots are recorded but not included in the analysis qualification flags.
No causal timing bias is inferred from `ps`'s lifetime-averaged CPU percentage. The observed
condition nevertheless needs the coordinator's ruling before these data support decisions.

This is a qualification defect/condition failure, not evidence of numerically wrong Core results.
As requested, stop here instead of publishing complete qualified E0/E7 tables or changing code.

| Reproducible witness | Value |
| --- | --- |
| Raw file | /Users/networmix/ws/work/2026-11-core-routing-research.noindex/scratch/replay/timing-window-02/run/window.json |
| Selector | blocks: row=1, experiment=cli, triple=0, slot=A1; before.other_cpu |
| Time | 2026-11-08T12:49:54+00:00 |
| Process | PID 3579, Superset Helper (Renderer), 71.9% CPU |
| Requirement | rules/PERFORMANCE.md §3: no other process above half a core; the exception is at most two system background daemons at about one core each |
| Analysis result | p1 / workers=8 / cli: qualified_noise=true |
| Opening-only CPU gate | /Users/networmix/ws/work/2026-11-core-routing-research.noindex/tasks/replay/timing/experiment.py:332 |
| Later block validity gate | /Users/networmix/ws/work/2026-11-core-routing-research.noindex/tasks/replay/timing/experiment.py:570 (power, sleep, deadline) |
| Final noise flag | /Users/networmix/ws/work/2026-11-core-routing-research.noindex/tasks/replay/timing/experiment.py:774 (smoke, short block, swap, direction, statistical noise; no background-CPU gate) |
| All non-allowlisted witnesses | /Users/networmix/ws/work/2026-11-core-routing-research.noindex/scratch/replay/timing-analysis-02/quiet-window-witnesses.json |

The starting snapshot was quiet (load 1.044, no listed other CPU consumers), all recorded power
states remained AC, the sleep marker did not change and swap-in/out counters remained zero.
That does not remove later renderer/Claude/browser activity from the saved snapshots. CPU flags
also occur in other rows; mediaanalysisd reaches 207.6%. Do not silently remove affected blocks,
relax the machine condition or infer how much their background load changed the measurements.

The raw window is structurally complete: six rows, 126 triples, 378 A1/B/A2 slots. The tool equals
snapshot-r2. Native per-kind analysis provisionally reports 13 passes out of 360 row/phase/kind
entries, before the machine-condition review. Those flags are not final E7 permission.
P2/8 matched also correctly fails the 2 s minimum in four slots (1.987488484, 1.991562908,
1.997795121 and 1.996617540 s); all its kinds remain unqualified. Every E0 A/A memory verdict
passes; all six matched consumer/replay memory comparisons fail, reflecting different process
surfaces rather than a candidate regression. Complete numerical rows remain in analysis.json.

Raw and analysis evidence are retained. No code was changed and no window was rerun. The
coordinator must decide the affected-row/window disposition and the CPU-condition handling.

## Superseded disposition: ruling 18

The coordinator rechecked these witnesses and accepted window-02. PERFORMANCE's quiet condition
gates the opening; later state is recorded, and a failure is rerun once only where the load
sensitivity explains it. Thus the earlier stop was appropriate pending clarification, but its
window-wide qualification hold is no longer current. E0.md/E7.md now contain the accepted
original analysis verdicts with witnessed triples/processes and the prescribed sensitivity.
No CPU gate was added to timing/. The exact rule and the size decision are in RULINGS.md 18/16.
