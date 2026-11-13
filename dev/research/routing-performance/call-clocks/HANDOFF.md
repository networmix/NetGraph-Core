# T4-CLOCKS: coordinator commands after stage 2

Run unconfined, sequentially, in new output paths. Do not rerun an existing path.
Stage 2 does not run the measured clock series, timing window or samplers.
The source/size and final check status are summarized in EXIT.md.

```bash
EFFORT=/Users/networmix/ws/work/2026-11-core-routing-research.noindex
TASK="$EFFORT/tasks/clocks"
S="$EFFORT/scratch/clocks"
PY="$EFFORT/scratch/consumer/ref/venv/bin/python"
REF=208487c7aed158fb418791c22e23d09332533595053c0e4810923b31244f7c84
BOUND="$EFFORT/scratch/orient-add/scripts/bounded.py"
test "$(git -C "$EFFORT/worktrees/core-ref" rev-parse HEAD)" = 9cd7642b5ea8681f7833c139969d8ed224d03a75
test "$(git -C "$EFFORT/worktrees/netgraph-orient" rev-parse HEAD)" = 8a641a1caa0b8f4b7d8c4df0aecf8625310cdf47
git -C "$EFFORT/worktrees/core-ref" diff --quiet HEAD
git -C "$EFFORT/worktrees/netgraph-orient" diff --quiet HEAD

# One registered-eight-thread null control at full calibration length.
HEAVY_OWNER=clocks bash "$EFFORT/rules/tools/heavy.sh" \
  "$PY" -B "$BOUND" 180 "$S/measured-null-8-01.log" \
  "$PY" -B "$TASK/calibrate_parallel.py" "$S/measured-null-8-01.json"

# All phases, all public boundaries, scoped native releases; no added latency.
HEAVY_OWNER=clocks bash "$EFFORT/rules/tools/heavy.sh" \
  "$PY" -B "$BOUND" 600 "$S/measured-p1-clocks-01.log" \
  "$PY" -B "$TASK/run_workload.py" --workload dcbb96 --workers 8 \
  --lane clocks --native-sha "$REF" --expected "$S/pilot-01/plain" \
  --out "$S/measured-p1-clocks-01"

HEAVY_OWNER=clocks bash "$EFFORT/rules/tools/heavy.sh" \
  "$PY" -B "$BOUND" 600 "$S/measured-p2-clocks-01.log" \
  "$PY" -B "$TASK/run_workload.py" --workload nsfnet --workers 8 \
  --lane clocks --native-sha "$REF" --expected "$S/smokes-02/nsfnet-8-full-plain" \
  --out "$S/measured-p2-clocks-01"

# Analyze without timing; P2 export getter rows are in phases.export.operations.
"$PY" -B "$TASK/analyze_clocks.py" "$S/measured-p1-clocks-01" \
  "$S/measured-p1-analysis-01.json" --parallel-calibration "$S/measured-null-8-01.json"
"$PY" -B "$TASK/analyze_clocks.py" "$S/measured-p2-clocks-01" \
  "$S/measured-p2-analysis-01.json" --parallel-calibration "$S/measured-null-8-01.json"
```

The null control's overhead is diagnostic. The analyzer retains raw values,
single-thread null correction, and sensitivity to the eight-thread correction.
Those are not guaranteed error bounds. Snapshot/lifecycle and collector-record
locks are instrumentation bias; caller CPU outside spans is not a proven serial
floor. A residual close to clock/startup tails stays unresolved.

Reserve one quiet, on-mains window for the sensitivity series and hold other
tasks. The pilot's conservative projection is under 45 minutes. It included
redundant YAML override preparation removed from the final full P1/8 path; its
process wall is a reservation estimate, not the final reference timing. The coordinator
must inspect machine snapshots and clusters, and confirm no sleep/competing work,
before qualifying any statistics. Full process timing includes source/setup,
imports, load/prepare, run, export and cleanup; JSON export comparisons occur in
the parent after the timed process exits. Analysis and parallel schedule wall
are also reported separately. δ is applied only on P1 case-execution threads.

```bash
HEAVY_OWNER=clocks-timing bash "$EFFORT/rules/tools/heavy.sh" \
  "$PY" -B "$BOUND" 2800 "$S/measured-sensitivity-01.log" \
  "$PY" -B "$TASK/sensitivity.py" --pilot "$S/pilot-01" \
  --window-minutes 45 --out "$S/measured-sensitivity-01"
```

This runner performs seven A/A/A triples followed by seven A/B/A triples at each
nonzero δ; 63 fresh processes. `delta=0` invokes the same detached helper.
The timing window's own deadline bounds each child; the outer timeout includes
shutdown/statistics margin. No ordinary worker job uses this longer hold.

```bash
HEAVY_OWNER=clocks bash "$EFFORT/rules/tools/heavy.sh" \
  "$PY" -B "$BOUND" 150 "$S/measured-p2-native-8-01.log" \
  "$PY" -B "$TASK/capture.py" --workload nsfnet --workers 8 \
  --expected "$S/smokes-02/nsfnet-8-full-plain" --out "$S/measured-p2-native-8-01"

HEAVY_OWNER=clocks bash "$EFFORT/rules/tools/heavy.sh" \
  "$PY" -B "$BOUND" 150 "$S/measured-p3-native-8-01.log" \
  "$PY" -B "$TASK/capture.py" --workload backbone --workers 8 \
  --expected "$S/smokes-02/backbone-8-full-plain" --out "$S/measured-p3-native-8-01"

HEAVY_OWNER=clocks bash "$EFFORT/rules/tools/heavy.sh" \
  "$PY" -B "$BOUND" 150 "$S/measured-p3-native-28-01.log" \
  "$PY" -B "$TASK/capture.py" --workload backbone --workers 28 \
  --expected "$S/smokes-02/backbone-28-full-plain" --out "$S/measured-p3-native-28-01"

"$PY" -B "$TASK/analyze_sample.py" "$S/measured-p2-native-8-01" \
  "$S/smokes-02/nsfnet-8-full-plain" "$S/measured-p2-native-analysis-01"
"$PY" -B "$TASK/analyze_sample.py" "$S/measured-p3-native-8-01" \
  "$S/smokes-02/backbone-8-full-plain" "$S/measured-p3-native-8-analysis-01"
"$PY" -B "$TASK/analyze_sample.py" "$S/measured-p3-native-28-01" \
  "$S/smokes-02/backbone-28-full-plain" "$S/measured-p3-native-28-analysis-01"
```

After these outputs return, this session prepares E5.md and E6.md. One-worker
profiles are carried with their original binary hashes and denominators.
New all-thread captures retain explicit waits and unresolved shares.
