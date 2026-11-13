"""Holding/release bookkeeping control: every native object stays alive."""

import json
import time
from pathlib import Path

import netgraph_core as core

from lifetime_tracker import LifetimeTimers

objects = [core.OutputOptions() for _ in range(1000)]
collector = LifetimeTimers(core)
collector.active = True
record = collector.begin_thread()
start = time.thread_time()
for _ in range(100):
    with collector.scope("no-destruction control"):
        for obj in objects:
            collector.retain(obj)
elapsed = time.thread_time() - start
collector.end_thread(record)
collector.active = False
assert all(not row["destroyed_wrappers"] for row in collector.releases)
clear_cpu = sum(row["cpu_ns"] for row in collector.releases) / 1e9
held = sum(row["held_wrappers"] for row in collector.releases)
report = {
    "held_references": held,
    "total_thread_cpu_s": elapsed,
    "dict_clear_cpu_s": clear_cpu,
    "dict_clear_ns_per_reference": clear_cpu * 1e9 / held,
    "holding_bookkeeping_ns_per_reference": (elapsed - clear_cpu) * 1e9 / held,
    "summary": "all native wrappers survived; estimated pure holding/weakref/counter and dict-decref overhead, diagnostic scale only",
}
p = Path(__file__).resolve().parents[1] / "lifetime-overhead-01.json"
with p.open("x") as stream:
    json.dump(report, stream, indent=2)
print(json.dumps(report), flush=True)
