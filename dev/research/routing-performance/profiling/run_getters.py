"""Export Core properties: grouped wall timers and null-call overhead control."""

import argparse
import io
import json
import statistics
import sys
import time
from collections import defaultdict
from pathlib import Path

BASE = Path(__file__).resolve().parents[1]
REPO = BASE.parents[1] / "worktrees/netgraph-orient"
sys.path.insert(0, str(REPO))
import netgraph_core as core  # noqa: E402 -- select the pinned consumer before importing
import ngraph  # noqa: E402 -- select the pinned consumer before importing
import yaml  # noqa: E402 -- select the pinned consumer before importing
from ngraph.scenario import Scenario  # noqa: E402 -- select the pinned consumer before importing

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("out", type=Path)
args = parser.parse_args()
args.out.mkdir()
document = yaml.safe_load((BASE / "workloads/nsfnet.yml").read_text())
for step in document["workflow"]:
    if step["type"] == "MaxFlow":
        step["workers"] = 1
prepared = ngraph.prepare(Scenario.from_yaml(yaml.safe_dump(document, sort_keys=False)))
result = ngraph.run(prepared)
# Plain and timed exports use the same prepared results and fresh identical sinks.
tick = time.perf_counter()
plain = io.StringIO()
ngraph.write_json(result, plain)
plain_wall = time.perf_counter() - tick
rows = defaultdict(lambda: [0, 0])


def timed(name, original):
    def wrapper(*a, **kw):
        start = time.perf_counter_ns()
        value = original(*a, **kw)
        elapsed = time.perf_counter_ns() - start
        row = rows[name]
        row[0] += 1
        row[1] += elapsed
        return value

    return wrapper


calibration = []
for repeat in range(5):

    def null():
        return None

    wrapper = timed("calibration", null)
    count = 100000
    tick = time.perf_counter_ns()
    for _ in range(count):
        null()
    direct = time.perf_counter_ns() - tick
    rows.clear()
    tick = time.perf_counter_ns()
    for _ in range(count):
        wrapper()
    full = time.perf_counter_ns() - tick
    calibration.append(
        {
            "full_wrap_ns": (full - direct) / count,
            "inside_timer_ns": rows["calibration"][1] / count,
        }
    )
rows.clear()
installed = []
for name, kind in list(vars(core).items()):
    if type(kind).__name__ != "pybind11_type" or hasattr(kind, "__members__"):
        continue
    for field, member in list(vars(kind).items()):
        if isinstance(member, property):
            label = name + "." + field
            setattr(
                kind,
                field,
                property(
                    timed(label, member.fget), member.fset, member.fdel, member.__doc__
                ),
            )
            installed.append(label)
tick = time.perf_counter()
with (args.out / "results.json").open("x") as stream:
    ngraph.write_json(result, stream)
export_wall = time.perf_counter() - tick
raw = {k: list(v) for k, v in rows.items()}
assert json.loads(plain.getvalue()) == json.loads(
    (args.out / "results.json").read_text()
)
inside = statistics.median(c["inside_timer_ns"] for c in calibration)
full = statistics.median(c["full_wrap_ns"] for c in calibration)
report_rows = [
    {
        "getter": k,
        "calls": v[0],
        "raw_s": v[1] / 1e9,
        "corrected_s": max(0, (v[1] - v[0] * inside) / 1e9),
    }
    for k, v in sorted(raw.items())
]
summary = {
    "exploratory": True,
    "plain_export_stringio_s": plain_wall,
    "timed_export_file_s": export_wall,
    "getter_calls": sum(v[0] for v in raw.values()),
    "raw_getter_s": sum(v[1] for v in raw.values()) / 1e9,
    "corrected_getter_s": sum(r["corrected_s"] for r in report_rows),
    "full_wrapper_overhead_s": sum(v[0] for v in raw.values()) * full / 1e9,
    "calibration": calibration,
    "rows": report_rows,
    "installed": installed,
    "summary": "plain and getter-timed full exports equal; null timer overhead subtracted",
}
(args.out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
print(
    json.dumps({k: v for k, v in summary.items() if k not in ("rows", "installed")}),
    flush=True,
)
