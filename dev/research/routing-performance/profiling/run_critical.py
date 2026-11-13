"""Whole analysis phase CPU partition at one/two/eight workers; fresh process."""

import argparse
import concurrent.futures.thread as thread_pool
import json
import os
import sys
import time
from contextlib import contextmanager
from pathlib import Path

from span_timers import SpanTimers

BASE = Path(__file__).resolve().parents[1]
REPO = BASE.parents[1] / "worktrees/netgraph-orient"
sys.path.insert(0, str(REPO))
import netgraph_core as core  # noqa: E402 -- select the pinned consumer before importing
import ngraph  # noqa: E402 -- select the pinned consumer before importing
import yaml  # noqa: E402 -- select the pinned consumer before importing
from ngraph.scenario import Scenario  # noqa: E402 -- select the pinned consumer before importing

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("workers", type=int, choices=(1, 2, 8))
parser.add_argument("out", type=Path)
args = parser.parse_args()
args.out.mkdir()
assert Path(ngraph.__file__).is_relative_to(REPO)
document = yaml.safe_load((BASE / "workloads/dcbb96.yml").read_text())
for step in document["workflow"]:
    if step["type"] == "TrafficMatrixPlacement":
        step["workers"] = args.workers
scenario = Scenario.from_yaml(yaml.safe_dump(document, sort_keys=False))
prepared = ngraph.prepare(scenario)
collector = SpanTimers()
collector.install(core)
calibration = [collector.calibrate() for _ in range(3)]
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
start = time.perf_counter()
result = ngraph.run(prepared, step_hook=hook)
wall = time.perf_counter() - start
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
twice_floor = lock + native / 2 if args.workers == 1 else None
summary = {
    "exploratory": True,
    "workers": args.workers,
    "analysis_wall_s": wall,
    "total_thread_cpu_s": raw_cpu,
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
    "bias": "CPU outside spans approximates GIL-held work and includes implicit object destruction; spans include binding and GIL reacquisition; free/twice bounds condition on unchanged outside work and parallel capacity, no latency prediction",
    "summary": "analysis thread clocks and Core spans recorded; last result exported",
}
with (args.out / "results.json").open("x") as stream:
    ngraph.write_json(result, stream)
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
            if k not in ("threads", "installed_boundaries", "calibration")
        }
    ),
    flush=True,
)
