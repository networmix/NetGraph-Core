"""P3 backbone: accepted lifetime collector plus same-process pipeline clocks."""

import argparse
import gc
import concurrent.futures.thread as thread_pool
import json
import os
import sys
import time
from contextlib import contextmanager
from pathlib import Path

from collections import Counter
import lifetime_tracker
from lifetime_tracker import LifetimeTimers

BASE = Path(__file__).resolve().parents[1]
REPO = BASE / "consumer-lifetime"
sys.path.insert(0, str(REPO))
phases = {}
import_start = time.perf_counter()
import netgraph_core as core  # noqa: E402 -- select the pinned consumer before importing
import ngraph  # noqa: E402 -- select the pinned consumer before importing
import yaml  # noqa: E402 -- select the pinned consumer before importing
from ngraph.scenario import Scenario  # noqa: E402 -- select the pinned consumer before importing

phases["import"] = time.perf_counter() - import_start

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("workload", choices=("dcbb96", "backbone"))
parser.add_argument("workers", choices=("1", "2", "8", "default"))
parser.add_argument("out", type=Path)
args = parser.parse_args()
args.out.mkdir()
assert Path(ngraph.__file__).is_relative_to(REPO)
text = (BASE / f"workloads/{args.workload}.yml").read_text()
if args.workers != "default":
    document = yaml.safe_load(text)
    for step in document["workflow"]:
        if step["type"] == "TrafficMatrixPlacement":
            step["workers"] = int(args.workers)
    text = yaml.safe_dump(document, sort_keys=False)
tick = time.perf_counter()
scenario = Scenario.from_yaml(text)
phases["load"] = time.perf_counter() - tick
tick = time.perf_counter()
prepared = ngraph.prepare(scenario)
phases["prepare"] = time.perf_counter() - tick
instrumentation_setup_start = time.perf_counter()
collector = LifetimeTimers(core)
lifetime_tracker.COLLECTOR = collector
collector.install(core)
calibration = [collector.calibrate() for _ in range(3)]
instrumentation_setup_s = time.perf_counter() - instrumentation_setup_start
original_worker = thread_pool._worker


def measured_worker(*a, **kw):
    record = collector.begin_thread()
    try:
        return original_worker(*a, **kw)
    finally:
        collector.end_thread(record)


thread_pool._worker = measured_worker
steps = []


@contextmanager
def hook(step):
    tick = time.perf_counter()
    yield
    steps.append({"name": step.name, "wall_s": time.perf_counter() - tick})


load_before = os.getloadavg()
collector.active = True
main = collector.begin_thread()
process_cpu_start = time.process_time()
start = time.perf_counter()
with collector.scope("ngraph.run"):
    result = ngraph.run(prepared, step_hook=hook)
wall = time.perf_counter() - start
phases["run"] = wall
process_cpu = time.process_time() - process_cpu_start
collector.end_thread(main)
collector.active = False
load_after = os.getloadavg()
threads = collector.summarize()
cal = {k: sorted(c[k] for c in calibration)[1] for k in calibration[0]}
ncalls = sum(t["calls"] for t in threads)
raw_cpu = sum(t["cpu_s"] for t in threads)
raw_native = sum(t["native_span_cpu_s"] for t in threads)
raw_native_wall = sum(t["native_span_wall_s"] for t in threads)
overhead = ncalls * cal["cost_per_wrap_cpu_ns"] / 1e9
native_overhead = ncalls * cal["recorded_null_cpu_ns"] / 1e9
lock = raw_cpu - raw_native - (overhead - native_overhead)
native = raw_native - native_overhead
# Lower bounds on the new wall: serial outside CPU, worker capacity including
# that CPU, and coordinator-only work. A work/span bound, not Amdahl equality.
free_floor = lock
twice_floor = lock + native / 2 if args.workers == "1" else None

summary = {
    "exploratory": True,
    "workload": args.workload,
    "workers": args.workers,
    "analysis_wall_s": wall,
    "instrumentation_setup_excluded_s": instrumentation_setup_s,
    "total_thread_cpu_s": raw_cpu,
    "process_cpu_s": process_cpu,
    "uncovered_process_cpu_s": process_cpu - raw_cpu,
    "native_span_wall_s": raw_native_wall,
    "native_span_cpu_s": raw_native,
    "outside_span_cpu_s": raw_cpu - raw_native,
    "corrected_outside_cpu_s": lock,
    "corrected_native_cpu_s": native,
    "wrapper_calls": ncalls,
    "estimated_wrap_overhead_cpu_s": overhead,
    "estimated_wrap_overhead_wall_serial_s": ncalls
    * cal["cost_per_wrap_wall_ns"]
    / 1e9,
    "free_core_wall_lower_bound_s": free_floor,
    "free_core_speedup_ceiling": wall / free_floor,
    "twice_core_wall_lower_bound_s": twice_floor,
    "twice_core_speedup_ceiling": wall / twice_floor
    if twice_floor is not None
    else None,
    "calibration": calibration,
    "threads": threads,
    "steps": steps,
    "loadavg_before": load_before,
    "loadavg_after": load_after,
    "installed_boundaries": collector.installed,
    "ngraph_path": ngraph.__file__,
    "bias": "CPU outside spans includes Python/runtime plus retention/release bookkeeping; scoped wrapper destruction is timed; attached binding conversion remains in native spans. Output/prepared objects alive after run are excluded from analysis destruction. Bounds condition on unchanged outside work; no latency prediction",
    "summary": "analysis thread clocks and Core spans recorded; last result exported",
}
getter_events = [
    e
    for thread in collector.threads
    for e in thread["events"]
    if e[0] in collector.properties
]
release_events = [
    e
    for thread in collector.threads
    for e in thread["events"]
    if e[0] == "CoreObjects.release"
]
getter_cpu = sum(e[4] - e[3] for e in getter_events) / 1e9
release_cpu = sum(e[4] - e[3] for e in release_events) / 1e9
summary["getter_calls"] = len(getter_events)
summary["getter_cpu_s"] = getter_cpu
summary["getter_wall_s"] = sum(e[2] - e[1] for e in getter_events) / 1e9
summary["destruction_scope_calls"] = len(release_events)
summary["destruction_scope_cpu_s"] = release_cpu
summary["destruction_scope_wall_s"] = sum(e[2] - e[1] for e in release_events) / 1e9
summary["outside_without_getters_and_destruction_raw_s"] = (
    raw_cpu - raw_native + getter_cpu + release_cpu
)
summary["destroyed_wrappers"] = dict(
    sum((Counter(row["destroyed_wrappers"]) for row in collector.releases), Counter())
)
summary["release_rows"] = collector.releases

tick = time.perf_counter()
with (args.out / "results.json").open("x") as stream:
    ngraph.write_json(result, stream)
phases["export"] = time.perf_counter() - tick
# Thread records contain only times/names, never native result references.
# Clear the user objects at the ordinary post-export cleanup boundary.
tick = time.perf_counter()
del result, prepared, scenario
gc.collect()
phases["cleanup"] = time.perf_counter() - tick
summary["phases_s"] = phases
summary["pipeline_wall_s"] = sum(phases.values())
summary["pipeline_other_phases_s"] = summary["pipeline_wall_s"] - wall
summary["pipeline_clock_scope"] = (
    "sum of import/load/prepare/run/export/cleanup; excludes input read/YAML override, harness startup, instrumentation setup/calibration and evidence analysis/writes; export wrappers inactive but installed"
)
with (args.out / "summary.json").open("x") as stream:
    json.dump(summary, stream, indent=2)
with (args.out / "events.jsonl").open("x") as stream:
    for i, t in enumerate(collector.threads):
        for name, w0, w1, c0, c1 in t["events"]:
            stream.write(
                json.dumps(
                    {
                        "thread": i,
                        "operation": name,
                        "wall_start_ns": w0,
                        "wall_stop_ns": w1,
                        "cpu_start_ns": c0,
                        "cpu_stop_ns": c1,
                    }
                )
                + "\n"
            )
    stream.write(
        json.dumps({"summary": "all Core call spans recorded", "calls": ncalls}) + "\n"
    )
print(
    json.dumps(
        {
            k: v
            for k, v in summary.items()
            if k
            not in ("threads", "installed_boundaries", "calibration", "release_rows")
        }
    ),
    flush=True,
)
