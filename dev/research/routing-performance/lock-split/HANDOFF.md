# SPLIT — stage 2 coordinator commands

Run these commands **unconfined**, sequentially. Stage 2 has built the scratch wheel
and checked six shortened plain/instrumented rows, their exports and closure, plus
the capture launcher without sampling. None of the six full measured rows or the
actual P1/8 capture has run. All output paths below must be new.

The rows are exploratory, one instrumented process each. `measure.py` also runs a
matching plain reference and null calibration, then analyzes the retained arrays in
a separate process. No A/B/A or qualified gain claim is made. Each child has a
240-second timeout inside the outer 600-second slot job. Do not nest heavy holds.
The parallel null is a Python-caller stress control. P2/8 has one Python caller;
its saved eight-Python-thread control is inapplicable and is ignored in SPLIT.md.

```bash
set -euo pipefail
E=/Users/networmix/ws/work/2026-11-core-routing-research.noindex
T="$E/tasks/split"
S="$E/scratch/split"
PY="$E/scratch/consumer/ref/venv/bin/python"
BOUND="$E/scratch/orient-add/scripts/bounded.py"
export TASK_SCRATCH="$S" PYTHONDONTWRITEBYTECODE=1
export TMPDIR="$S/tmp" XDG_CACHE_HOME="$S/cache" MPLCONFIGDIR="$S/mplconfig"

test "$(git -C "$E/worktrees/core-ref" rev-parse HEAD)" = 9cd7642b5ea8681f7833c139969d8ed224d03a75
test "$(git -C "$E/worktrees/netgraph-orient" rev-parse HEAD)" = 8a641a1caa0b8f4b7d8c4df0aecf8625310cdf47
git -C "$E/worktrees/core-ref" diff --quiet HEAD
git -C "$E/worktrees/netgraph-orient" diff --quiet HEAD

HEAVY_OWNER=split bash "$E/rules/tools/heavy.sh" \
  "$PY" -B "$BOUND" 600 "$S/measured-p1-1-01.log" \
  "$PY" -B "$T/measure.py" --workload dcbb96 --workers 1 --out "$S/measured-p1-1-01"

HEAVY_OWNER=split bash "$E/rules/tools/heavy.sh" \
  "$PY" -B "$BOUND" 600 "$S/measured-p1-8-01.log" \
  "$PY" -B "$T/measure.py" --workload dcbb96 --workers 8 --out "$S/measured-p1-8-01"

HEAVY_OWNER=split bash "$E/rules/tools/heavy.sh" \
  "$PY" -B "$BOUND" 600 "$S/measured-p2-1-01.log" \
  "$PY" -B "$T/measure.py" --workload nsfnet --workers 1 --out "$S/measured-p2-1-01"

HEAVY_OWNER=split bash "$E/rules/tools/heavy.sh" \
  "$PY" -B "$BOUND" 600 "$S/measured-p2-8-01.log" \
  "$PY" -B "$T/measure.py" --workload nsfnet --workers 8 --out "$S/measured-p2-8-01"

HEAVY_OWNER=split bash "$E/rules/tools/heavy.sh" \
  "$PY" -B "$BOUND" 600 "$S/measured-p3-1-01.log" \
  "$PY" -B "$T/measure.py" --workload backbone --workers 1 --out "$S/measured-p3-1-01"

HEAVY_OWNER=split bash "$E/rules/tools/heavy.sh" \
  "$PY" -B "$BOUND" 600 "$S/measured-p3-28-01.log" \
  "$PY" -B "$T/measure.py" --workload backbone --workers 28 --out "$S/measured-p3-28-01"

HEAVY_OWNER=split bash "$E/rules/tools/heavy.sh" \
  "$PY" -B "$BOUND" 180 "$S/p1-8-sample-01.log" \
  "$PY" -B "$T/capture.py" --workload dcbb96 --workers 8 --seconds 15 \
  --expected "$S/measured-p1-8-01/plain" --out "$S/p1-8-sample-01"
```

Budget: allow about 1–3 minutes per row and 1–2 minutes for the capture. These
are reservations, not measured user timings. A failed/timeout job ends with its
logs retained; do not overwrite it or silently continue to the next row.

Each measured row writes `plain/`, `clocks/` (including `events.i64`,
`python-vectors.i64`, `detached.i64`, complete export and summary),
`parallel-null.json`, `analysis.json` and `index.json`. `analysis.json` includes
per-phase/per-thread and per-call splits, the top five kinds by run-span wall,
native-worker residual brackets and raw/corrected closure. Corrections remain
diagnostic, not bounds on instrumented scheduling. Signed corrected values stay
signed. Cleanup's held residual includes whole-model Python teardown and GC.

The capture repeats complete full workflows for at least 17 seconds, retains
their results, and exports/checks **every** result after sampling. Its `index.json`
must have `sample_exit=0`, `worker_exit=0`, `full_export_equal=true` and a positive
`exports_compared`. Retention changes lifetimes/memory; sampling biases remain.

After the capture, the sample analysis can run without another measured workload:

```bash
HEAVY_OWNER=split bash "$E/rules/tools/heavy.sh" \
  "$PY" -B "$BOUND" 120 "$S/p1-8-sample-analysis-01.log" \
  "$PY" -B "$T/analyze_sample.py" "$S/p1-8-sample-01" \
  "$S/measured-p1-8-01/plain" "$S/p1-8-sample-analysis-01"
```

It checks all roots, exclusive-weight conservation, named top-of-stack counts,
export equality and live reference-wheel static_support instruction witnesses.
The six category percentages use Core-context samples excluding recognized
blocked leaves. Unknown work remains unresolved. The carried P1/1 column uses a
different binary and denominator; neither column is a CPU share or worker speedup.

Return the output directories/logs to this session for `SPLIT.md`. No further
decision is currently needed; this is the stage-2 measured-run handoff.
