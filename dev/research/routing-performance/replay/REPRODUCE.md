# Recording and checking Core workloads

Use the frozen reference interpreter to record. Run checks in separate wheel interpreters; the
native import path, version and binary SHA-256 are verified before work starts. No environment
is installed or modified. All output paths must be new. The campaign launcher uses the common
heavy slot and writes each command, timeout outcome and final summary to its own log.

From `/Users/networmix/ws/work/2026-11-core-routing-research.noindex/tasks/replay`:

```bash
EFFORT=/Users/networmix/ws/work/2026-11-core-routing-research.noindex
PY="$EFFORT/scratch/consumer/ref/venv/bin/python"
export PYTHONDONTWRITEBYTECODE=1
export XDG_CACHE_HOME="$EFFORT/scratch/replay/cache"
export MPLCONFIGDIR="$EFFORT/scratch/replay/mplconfig"

"$PY" -B campaign.py record backbone "$EFFORT/scratch/replay/record-p3-NEW"
"$PY" -B campaign.py record nsfnet "$EFFORT/scratch/replay/record-p2-NEW"
"$PY" -B campaign.py record dcbb96 "$EFFORT/scratch/replay/record-p1-NEW"

"$PY" -B campaign.py check "$EFFORT/scratch/replay/record-p3-NEW" \
  "$EFFORT/scratch/replay/check-p3-a28-NEW" --workers 28 --lane a --pypi
"$PY" -B campaign.py check "$EFFORT/scratch/replay/record-p3-NEW" \
  "$EFFORT/scratch/replay/check-p3-a1-NEW" --workers 1 --lane a
"$PY" -B campaign.py check "$EFFORT/scratch/replay/record-p3-NEW" \
  "$EFFORT/scratch/replay/check-p3-b28-NEW" --workers 28 --lane b
```

For P1 use workers 1/8 and P2 1/8; lane b runs at the scenario count only. A campaign launches two
fresh reference processes and compares their observations; `--pypi` adds a fresh PyPI process and
its comparison. Each recording has a 180 s limit, each check/comparison a 600 s limit. Use
`--shards N` to partition observations by case jobs when a full check would exceed its limit.
Every shard executes the same entire native dependency sequence, and only snapshots its selected
case jobs; main-thread observations belong to shard 0. Summing the shards' `observed` counts must
recover the recording's coverage. A timeout is a failed/incomplete run, not a clean check.

A single check for an explicitly identified candidate wheel:

```bash
HEAVY_OWNER=replay bash "$EFFORT/rules/tools/heavy.sh" \
  "$PY" -B "$EFFORT/scratch/orient-add/scripts/bounded.py" 600 \
  "$EFFORT/scratch/replay/candidate-NEW.log" \
  /absolute/candidate/venv/bin/python -B replay.py \
  --expected-prefix /absolute/candidate/venv \
  --expected-sha256 FULL_NATIVE_SHA256 --expected-version WHEEL_VERSION \
  check "$EFFORT/scratch/replay/record-p3-NEW" \
  "$EFFORT/scratch/replay/candidate-NEW" --workers 28 --lane a
```

Only P1 candidate checks use `--state-sample`. P2/P3 retain full state on every job. P1 uses `--shard 0/2` and `--shard 1/2` in separate
fresh output directories, against references with the same settings. P2/P3 use `--shard 0/1`.
Tiny controls intentionally retain full state on all their jobs.

The CLI's `compare REFERENCE_OBSERVATIONS CANDIDATE_OBSERVATIONS` command checks all recorded
fields exactly, loading array payloads by value even when their content-address keys agree.
It refuses an empty or incomplete stream. Build identity may differ between sides, but recording,
mode, worker count, shard, tool source hashes, event counts and coverage must agree.

## Files and field paths

A recording directory contains `manifest.json`, full array payloads in `payloads/*.npz`, and one
`JOB.jsonl.gz` stream per logical job. `main` contains import/load/prepare, serial work, case-group
barriers, export and cleanup. A call event records `op`, typed `args`/`kwargs`, output references,
constructor target and exception. Release events drop the recorded references. Containers retain
their Python type; `range` retains start/stop/step; arrays retain dtype, shape, strides and
mutability. A graph/program/mask/policy input references its actual producer call. An array input
also has its full recorded content. Masks are created from their full exclusion arguments.

Observation streams use `call` for the recording's event index, `op`, `result`, `error`, `state`
and `inputs_unchanged`. `state[OWNER].full` is the first complete public state. Later entries are
exact recursive deltas: `fields`, `items`, `replace`, `removed`, and array `indices` plus new
`values`. A `null` delta means no change. Values retain tuple/list distinctions, enum class/name/
value, arbitrary-width integers, binary64 bit strings and rational numerator/denominator. Graph,
ledger and mask uids are bijectively mapped to logical owners; local record/policy ids are exact.
Known ledger uids in FlowId input arguments are translated back into the selected process's uid.

An array value is a payload reference with dtype/shape. To plant a wrong array value, write a new
payload and update that one observation's reference; do not damage a shared payload used by other
observations. `preflight.py` demonstrates this on `FlowLedger.carried`, and separately changes a
real nonempty `PlacementResult.min_cut`. It invokes the actual compare CLI and requires nonzero
exit plus job, call, field path and both values. The reference stream stays untouched.

```bash
"$PY" -B controls.py "$EFFORT/scratch/replay/control-record-NEW"
"$PY" -B replay.py check "$EFFORT/scratch/replay/control-record-NEW" \
  "$EFFORT/scratch/replay/control-ref-NEW" --workers 2 --lane a
"$PY" -B preflight.py "$EFFORT/scratch/replay/control-ref-NEW" \
  "$EFFORT/scratch/replay/preflight-NEW"
```

The tiny controls also exercise all-off/individual output flags, all-on tunnel results, orphan
records, failure and retry, cache lifecycle, search results, wide exact counters/cost and rational
T2 bandwidth. These are short untimed probes. Whole-workload recording/check commands use the
heavy slot. Product suites, sanitizers and acceptance timing are not part of this stage.

## Release boundary

The recorder uses weak references and never retains observed native wrappers. The replay keeps
one slot per recorded wrapper and drops it at the recorded event. Native pybind holders may
release a child before the parent's own weakref callback fires. For adjacent release events the
replay drops every slot, then asserts that all those wrappers are dead before the next Core call,
parallel group or job end. It does not insert a consumer release adapter. Module-retained defaults
remain marked retained at the end; interpreter-finalization work is outside replay coverage.
Releases will remain an unqualified E7 group, as the accepted rulings require.

## Sequential execution and lookup preflight

Changed by ruling 9: bulk record snapshots retain every public `FlowRecord` field, while the
observer-added `record_lookup` subtree is omitted. Actual consumer `record()` calls are replayed
and read fully. Before every check, the tiny real-Core stream and complete single-/two-source
lookup controls run; a failure prevents a clean observation manifest. A candidate changing record
derivation or either getter additionally needs the P1 MSD-roster lookup cross-check in its brief.

At most one replay job may be waiting or running in the heavy slot. Do not launch campaigns in
parallel. The all-state matrix was stopped for cost; the authorized replacement uses every
returned field everywhere and a fixed full-state sample. `RUN-COST.md` records its bounded
return-cost measurement, two-child P1 partition and nine-row budget including reserve.

```bash
"$PY" -B campaign.py matrix "$EFFORT/scratch/replay" \
  "$EFFORT/scratch/replay/matrix-r14-NEW" --p1-shards 2
```

The matrix passes `--prune` and `--state-sample` only for P1, uses PyPI only on a/1, and makes fresh reference pairs
for every row. Each child executes the entire native sequence. Returned-value observations are
partitioned by job index modulo the child count, with main on child 0; the union must recover
all recorded call-kind counts and jobs without overlap. Full state is recorded after every call
only on main and fixed sample jobs (index modulo 64 in {0,1}, widened by the first case in any
otherwise uncovered case group). Every output manifest lists the global state sample. P1 has
79 sampled cases spanning all 18 case workflow steps and both directions. P3 retains full state
on main and all 592 cases; P2 has no case jobs and its main state is complete. Only P1 has the
internal-only unsampled-state limitation. Later candidates use these same coverage choices.

A change only to internal state of an unsampled case is not seen unless a later returned field
changes. The consumer comparison checks fields NetGraph reads. This limitation must accompany
any check verdict. A failed check or 600 s timeout stops the driver; no unchanged retry is made.

`$int_pairs` stores FlowRecord's public list of two-int tuples as an exact int64 matrix. Strict
container/element types, row lengths, order and domains are checked before conversion with
iterator validation and `numpy.fromiter`. Public ndarray outputs keep their original array
representation. Record deltas use complete canonical FlowIds, explicit order and additions,
removals and changes. `statecheck.py` independently reconstructs full states; hashes never
establish equality.

The driver checks aggregate evidence size between children. All required side comparisons must
pass before pruning. It lists literal paths/bytes/SHA-256 first, rechecks each file, deletes only
listed streams/payloads, verifies absence, and keeps logs, comparison outputs, manifests and
deletion receipts. Child 0 and child 1 streams are always kept for coordinator plants. With the
authorized two-child P1 matrix all new child streams therefore remain retained.

## Timed mode and E0/E7

The accepted stage-1 root sources are preserved in the snapshots. Timing has a separate entry point under
`timing/`; its exact boundaries, prepared inputs, coordinator window command and analysis are
in [TIMING.md](TIMING.md). The current six-row input manifest is
`/Users/networmix/ws/work/2026-11-core-routing-research.noindex/scratch/replay/timing-inputs-01/report.json`.
The coordinator's completed window-02 and its E0/E7 qualification are documented in E0.md
and E7.md. Later code controls are unqualified and do not replace that window's numbers.

Current-source provenance for checks: changing campaign.py changes the deliberately broad
root tool hash set, even though replay.py/values.py are unchanged. P1 and P2 retained references
therefore use their matching snapshot-r3/task/replay.py entry point when making candidate
observations. P3's replacement full references use the current root sources, frozen next as
snapshot-r4. Match every observation manifest's tool_sha256; never waive that guard or compare
new-source observations against an old-source reference. The saved same-source E0/E7 window
continues to describe the original timing code; future counts/archive edits need their own pilot.

The replacement P3 full-state references are
`/Users/networmix/ws/work/2026-11-core-routing-research.noindex/scratch/replay/check-p3-full-r22-01/{a1,a28,b28}/s00-ref1`.
The matching ref2 streams and a1 PyPI streams are retained beside them. Each observes main plus
all 592 cases (114,618 calls, 50 kinds). P3's older sampled matrix-r14 observations are historical
evidence and are not the full-state references for later candidates.
