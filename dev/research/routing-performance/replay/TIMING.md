# Timed replay and the E0/E7 window

The entry point is `timing/worker.py replay RECORDING OUT --workers N --repeats R`.
Use a new output directory for every process. The ordinary output is unqualified timing data;
only the coordinator's complete reference window can qualify noise or fidelity.

The stage-1 sources and `matrix-r14-01` remain unchanged. The timing implementation lives in
`timing/`, outside the stage-1 source-hash set, so candidate checks still use the accepted streams.

## Measurement boundary

The plan reads every event and array payload before its block. It prepares external scalar,
enum, fraction and constant container values, native dispatch functions, array templates and
workers before measurement. Fresh external arrays are allocated for each iteration before the
block; their shape, strides and mutability use the accepted decoder. Containers depending on
Core-produced values are completed between spans. Their construction is not native time.

Every native constructor and dependency-producing call runs at its original point. Constructor
`__new__` is outside both replay and consumer clocks; `__init__` is inside. No bound native
instance methods are cached. The replay preserves recorded releases and retains only the recorded
module defaults through the end of the block; their finalization stays outside coverage.
Explicit drops are timed in `CoreObjects.release` and remain unqualified for E7.

The same payload-disabled recorder subclass supplies interception and clocks on both paths. It
uses the accepted descriptor installation, scheduler and weak-reference registry, without input
encoding or result snapshots. Replay does not patch Core classes; it wraps the saved unbound
functions. The matched consumer patches the same API boundaries inside NetGraph. P1/P3 use
Python case workers; P2 has one Python caller and the requested native evaluator workers.

Each block reports the union of wall intervals, first-to-last elapsed wall, their difference,
phase unions, per-kind wall and caller-thread CPU sums/counts, process CPU, GC counters and peak
RSS. Caller-thread CPU excludes P2's native worker CPU; fresh-process CPU includes all threads.
The thread-pool startup is outside replay blocks. GC remains enabled with its default thresholds.
Null controls retain the clock and dispatch costs; nothing is subtracted from measured time.

`--dry` additionally compares compiled arguments with the accepted decoder, checks unchanged
inputs, records full check-mode observations and checks release boundaries. Those audits and
payload writes occur only in explicitly unqualified dry runs, never in a qualified block.

## Prepare the six rows

The generated variants change only workflow worker fields; the generator reparses them,
restores the original worker fields and compares the remaining document exactly. Original and
derived SHA-256 values and the changed paths are recorded. Empty-workflow controls are separate
inputs, with their own hashes.

```bash
EFFORT=/Users/networmix/ws/work/2026-11-core-routing-research.noindex
PY="$EFFORT/scratch/consumer/ref/venv/bin/python"
T="$EFFORT/tasks/replay/timing"

HEAVY_OWNER=replay bash "$EFFORT/rules/tools/heavy.sh" \
  "$PY" -B "$T/worker.py" prepare \
  "$EFFORT/scratch/orient-add/workloads/manifest.json" \
  "$EFFORT/scratch/replay/timing-inputs-NEW"
```

An explicitly identified candidate interpreter may run the same timed plan with
`--expected-prefix`, `--expected-sha256` and `--expected-version`, as the check-mode build guard
requires. Its untimed checks must already pass against the retained stage-1 references. E0/E7
qualification itself uses only the frozen reference wheel.

## Coordinator window

Run unconfined in one quiet reservation. The script never acquires additional slots inside the
reservation. Pilot children and measured children each have a timeout at most 600 s; the pilot
has an approximately ten-minute cap. The two-hour deadline starts before the pilot, not after it.

```bash
HEAVY_OWNER=replay-timing bash "$EFFORT/rules/tools/heavy.sh" \
  bash "$T/window.sh" \
  "$EFFORT/scratch/replay/timing-inputs-NEW/report.json" \
  "$EFFORT/scratch/replay/timing-window-NEW"
```

The pilot measures one reference sequence on each path and fixes repeat counts before the
triples. Targets are 2 s at every worker count (ruling 17); 2 s is also the hard minimum.
Counts use 2 s × 1.1 before ceiling to leave margin (ruling 23).
It prints the full plan and a projection including 20% reserve. It refuses the measured campaign
if the remaining two-hour budget or projected 20 GiB evidence budget does not fit. No bound or
statistical threshold is raised automatically. Each replay block is one fresh process with full
sequences repeated inside it; each consumer block sums complete fresh CLI processes.

E0 runs seven A/A/A triples for replay and seven for the consumer on each of the six rows.
E7 runs seven consumer/replay/consumer triples on the same rows, with identical whole-sequence
counts and inputs. The source input and recording hashes must agree, and derived documents are
checked independently. Actual API counts, Python workers and native evaluator workers are saved.

The fresh CLI launcher adds only before/after GC metadata around the real `ngraph` module entry.
Its cost is included in complete process wall time, never subtracted. The pilot also runs the
unmodified `python -m ngraph` command, compares its entire export and records the timing ratio
as a diagnostic. Import-only and empty-workflow controls report startup costs separately.
Exports are compared by typed values and key/list order outside the child. Pilot exports keep
lossless gzip copies. Other in-window exports keep only a receipt with path, size, SHA-256 and
the comparison reference, then their plaintext is removed after its full comparison. No matrix
stream is deleted. The GC launcher/control and per-block evidence-size scan remain by choice.

Machine state is saved before and after each measured block: load, CPU consumers, mains power,
swap-in/out counters and the kernel sleep marker. Missing state prevents a qualified run; a power
switch, sleep or deadline breach voids the window. Swap changes keep affected rows inconclusive.
Raw child reports, commands, state readings and partial results survive a failure. Pilot exports
remain archived; already-compared measured exports retain size/SHA receipts rather than copies.

## Analysis and limits

```bash
"$PY" -B "$T/experiment.py" analyse \
  "$EFFORT/scratch/replay/timing-window-NEW/run/window.json" \
  "$EFFORT/scratch/replay/timing-analysis-NEW"
```

The analysis uses `aba.py`'s whole-triple bootstrap (10,000 resamples, seed 20261108), reference
drift, same-side clusters and balanced fresh-process RSS allowance. It reports E0's union and
elapsed noise bounds separately; opposite beyond-noise directions are inconclusive. The smallest detectable gain is the
noise-only bound U (ruling 21); bar×U and U/0.97 are the true gains a candidate needs to pass at that block
length.

E7 names every phase/kind, including missing kinds as unqualified. A kind qualifies only with
matching call/release counts, stable references, its entire CI within [0.95,1.05], and in-span null
cost below 5% of its per-call mean. Releases never qualify. An unqualified kind requires consumer
comparison; no whole-union fidelity claim hides it. Memory verdicts are reported independently.

The short smoke mode uses a small real scenario and one triple. It always remains unqualified.
Full primary block counts, machine access outside the sandbox and the two-hour projection are
not established by those controls. The coordinator runs the window; E0.md/E7.md follow its data.

## Current decision boundary (rulings 19–23)

The Native decision number is the full row union ratio with elapsed beside it, never a selected
kind's sum or ratio. P1 may decide only where every phase/kind changed code reaches is approved
on that row, with approved-kind union share and row clock share reported. Existing per-kind
sums give exact serial shares and only bounds for overlapping P1/8 intervals; E7.md states the
bounds and the non-identifiability. P2/P3 replay rows are diagnostic; their candidates are decided
on the consumer comparison because clock costs materially dilute their ratios.

Call-sequence, Python-visible lifetime, cross-thread timing/cache-order or Python object-count
changes are consumer-only. Check mode still checks the candidate build with its new path off.
Matching GC settings does not match GC activity. This includes opt-in REUSE-P1 forms.

Replay memory has a plan/harness floor: P1/1 771.63 MiB, P1/8 961.03 MiB, P2/1 5980.88 MiB,
P2/8 5981.81 MiB, P3/1 440.16 MiB and P3/28 446.03 MiB in window-02. Read memory risk from the
CLI rows, not a replay-plan allowance. E0.md gives each floor and allowed increment.

Recorded load counts symmetrically for passes and failures. A changed filtered verdict selects
one later rerun; >3 loaded triples makes the whole group not judgeable. P3/1 matched, including
apply_actions, loses approval until the next window that times it; no window opens for it alone.
P1 check mode uses the fixed 79-case state sample. P2 and P3 keep full state on every job.
