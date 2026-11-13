"""Cheap separating clock and descriptor controls for the research harness."""

import json
import time
from pathlib import Path

import netgraph_core as core

from evidence import equal
from original_analyze_native_sample import parse
from span_timers import SpanTimers

out = Path(__file__).resolve().parents[1] / "timer-controls-01.json"
collector = SpanTimers()
collector.install(core)
collector.active = True
record = collector.begin_thread()
w0, c0 = time.perf_counter(), time.thread_time()
time.sleep(0.02)
sleep_wall, sleep_cpu = time.perf_counter() - w0, time.thread_time() - c0
w0, c0 = time.perf_counter(), time.thread_time()
while time.perf_counter() - w0 < 0.02:
    pass
busy_wall, busy_cpu = time.perf_counter() - w0, time.thread_time() - c0
options = core.OutputOptions()
assert options.by_source is False
assert core.to_units(1.0, 0.5) == 2
collector.end_thread(record)
collector.active = False
rows = collector.summarize()
assert sleep_cpu < sleep_wall / 4
assert busy_cpu > busy_wall / 2
assert rows[0]["calls"] == 3
assert set(rows[0]["call_counts"]) == {
    "OutputOptions.__init__",
    "OutputOptions.by_source",
    "to_units",
}
# A hand-derived tree has self weights 2/1/2; malformed inclusive children fail.
tree = "Call graph:\n5 root\n  3 branch\n    2 leaf\nTotal number in stack\n"
nodes = parse(tree)
assert [n.count - sum(nodes[c].count for c in n.children) for n in nodes] == [2, 1, 2]
try:
    parse(tree.replace("  3 branch", "  6 branch"))
except AssertionError:
    pass
else:
    raise AssertionError("invalid tree accepted")
equal({"values": [1, 2]}, {"values": [1, 2]})
report = {
    "sleep_wall_s": sleep_wall,
    "sleep_thread_cpu_s": sleep_cpu,
    "busy_wall_s": busy_wall,
    "busy_thread_cpu_s": busy_cpu,
    "descriptor_calls": rows,
    "summary": "thread CPU separates sleep from busy work; property/constructor/function counted; sample tree self weights and malformed control passed",
}
with out.open("x") as stream:
    json.dump(report, stream, indent=2)
print(json.dumps(report), flush=True)
