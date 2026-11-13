"""Strict export checks and compact machine-readable additional survey tables."""

import hashlib
import json
import platform
import pstats
import subprocess
import sys
from collections import Counter, defaultdict
from pathlib import Path

from evidence import equal, normalize

BASE = Path(__file__).resolve().parents[1]
EFFORT = BASE.parents[1]
OLD = EFFORT / "scratch/orient-consumers"
checks = []


def read(path):
    return json.loads(path.read_text())


def compare(expected_path, actual_path):
    equal(normalize(read(expected_path)), normalize(read(actual_path)))
    checks.append(
        {
            "expected": str(expected_path),
            "actual": str(actual_path),
            "outcome": "all values/types/keys/order equal except metadata.workers",
        }
    )


for p in sorted(BASE.glob("backbone-*/results.json")):
    compare(BASE / "backbone-wall-1-01/results.json", p)
for w in (1, 2, 8):
    compare(
        OLD / "dcbb96-1-02/results.json", BASE / f"dcbb-critical-{w}-01/results.json"
    )
compare(OLD / "nsfnet-1-02/results.json", BASE / "nsfnet-getters-01/results.json")
for workload in ("backbone", "nsfnet"):
    p = BASE / f"native-dry-{workload}-01/1/results.json"
    if p.exists():
        compare(
            BASE / "backbone-wall-1-01/results.json"
            if workload == "backbone"
            else OLD / "nsfnet-1-02/results.json",
            p,
        )
# Comparator controls cover the guarantees used in all comparisons.
original = {"a": [1, 2], "metadata": {"workers": 1}}
mutations = [
    {"a": [1, 3], "metadata": {"workers": 1}},
    {"a": [1.0, 2], "metadata": {"workers": 1}},
    {"a": [1, 2]},
    {"a": [2, 1], "metadata": {"workers": 1}},
    {},
    {"metadata": {"workers": 1}, "a": [1, 2]},
]
for mutation in mutations:
    try:
        equal(normalize(original), normalize(mutation))
    except AssertionError:
        checks.append({"control": mutation, "outcome": "rejected"})
    else:
        raise AssertionError("poison accepted")
equal(normalize(original), normalize({"a": [1, 2], "metadata": {"workers": 8}}))
checks.append({"control": "metadata.workers override", "outcome": "accepted"})
profiles = {}
for phase in ("import", "load", "prepare", "run", "export", "cleanup"):
    stats = pstats.Stats(str(BASE / f"backbone-profile-1-01/{phase}.pstats"))
    profiles[phase] = {
        "total_self_s": sum(v[2] for v in stats.stats.values()),
        "known_python_result_self_s": sum(
            v[2]
            for k, v in stats.stats.items()
            if k[0].endswith(
                ("/analysis/costs.py", "/analysis/totals.py", "/results.py")
            )
        ),
        "core_builtins": [
            {"function": k[2], "calls": v[1], "self_s": v[2]}
            for k, v in stats.stats.items()
            if "_netgraph_core" in k[2]
        ],
        "top_self": [
            {
                "file": k[0],
                "line": k[1],
                "function": k[2],
                "calls": v[1],
                "self_s": v[2],
            }
            for k, v in sorted(
                stats.stats.items(), key=lambda x: x[1][2], reverse=True
            )[:15]
        ],
    }
critical = {}
for w in (1, 2, 8):
    x = read(BASE / f"dcbb-critical-{w}-01/summary.json")
    events = [
        json.loads(line)
        for line in (BASE / f"dcbb-critical-{w}-01/events.jsonl")
        .read_text()
        .splitlines()
    ]
    observed = Counter()
    previous = {}
    for event in events[:-1]:
        assert event["wall_stop_ns"] >= event["wall_start_ns"]
        assert event["cpu_stop_ns"] >= event["cpu_start_ns"]
        key = event["thread"]
        assert event["wall_start_ns"] >= previous.get(key, 0), (
            "nested or overlapping spans on same thread"
        )
        previous[key] = event["wall_stop_ns"]
        observed[key] += 1
    assert sum(observed.values()) == x["wrapper_calls"] == events[-1]["calls"]
    assert x["wrapper_calls"] == 49584
    checks.append(
        {
            "critical_workers": w,
            "outcome": "nonnegative, nonoverlapping thread spans; all 49584 events accounted",
        }
    )
    slots = defaultdict(lambda: defaultdict(float))
    for t in x["threads"][1:]:
        slot = t["name"].rsplit("_", 1)[1]
        for field in ("cpu_s", "native_span_cpu_s", "native_span_wall_s", "calls"):
            slots[slot][field] += t[field]
        slots[slot]["lifetimes"] += 1
        slots[slot]["case_count"] += t["call_counts"].get("FlowLedger.__init__", 0)
        slots[slot]["outside_span_cpu_s"] += t["cpu_s"] - t["native_span_cpu_s"]
    x["worker_slots_across_18_pools"] = dict(sorted(slots.items()))
    prior = read(OLD / f"dcbb96-wall-{'1' if w == 1 else 'default'}-01/summary.json")
    if w in (1, 8):
        nonanalysis = sum(prior["phases_s"].values()) - prior["phases_s"]["run"]
        x["pipeline_nonanalysis_prior_s"] = nonanalysis
        x["derived_pipeline_current_analysis_plus_prior_other_s"] = (
            nonanalysis + x["analysis_wall_s"]
        )
        x["derived_pipeline_free_core_speedup_ceiling"] = (
            nonanalysis + x["analysis_wall_s"]
        ) / (nonanalysis + x["corrected_outside_cpu_s"])
        if w == 1:
            x["derived_pipeline_twice_core_speedup_ceiling"] = (
                nonanalysis + x["analysis_wall_s"]
            ) / (nonanalysis + x["twice_core_wall_lower_bound_s"])
    if w > 1:
        for key in list(x):
            if "twice_core" in key:
                del x[key]
        x["finite_core_speedup_status"] = (
            "not identified from aggregate multi-worker sums"
        )
    critical[str(w)] = x
traces = [
    json.loads(line)
    for line in (BASE / "backbone-trace-1-01/calls.jsonl").read_text().splitlines()
    if '"operation"' in line
]
trace_total = sum(c["native_s"] for c in traces)
seen = set()
repeat_s = 0
for c in traces:
    k = (c["operation"], c["full"])
    if k in seen:
        repeat_s += c["native_s"]
    seen.add(k)
# Read actual default-worker metadata from full export, independent of schema defaults.
default = read(BASE / "backbone-wall-default-01/results.json")
default_workers = [
    s["result"].get("metadata", {}).get("workers") for s in default["steps"]
]
identities = {}
for repo in ("netgraph-orient", "core-orient"):
    p = EFFORT / "worktrees" / repo
    identities[repo] = {
        "commit": subprocess.check_output(
            ["git", "-C", str(p), "rev-parse", "HEAD"], text=True
        ).strip(),
        "tracked_status": subprocess.check_output(
            ["git", "-C", str(p), "status", "--porcelain", "--untracked-files=no"],
            text=True,
        ),
    }
import _netgraph_core as extension  # noqa: E402 -- identity after analysis

identities["extension"] = {
    "path": extension.__file__,
    "sha256": hashlib.sha256(Path(extension.__file__).read_bytes()).hexdigest(),
}
identities["python"] = sys.version
identities["machine"] = platform.platform()
lifetimes = {}
for label, workload, workers in [
    ("dcbb-lifetime-1-01", "dcbb96", 1),
    ("dcbb-lifetime-2-01", "dcbb96", 2),
    ("dcbb-lifetime-8-01", "dcbb96", 8),
    ("backbone-lifetime-1-01", "backbone", 1),
]:
    folder = BASE / label
    if not (folder / "summary.json").exists():
        continue
    x = read(folder / "summary.json")
    compare(
        OLD / "dcbb96-1-02/results.json"
        if workload == "dcbb96"
        else BASE / "backbone-wall-1-01/results.json",
        folder / "results.json",
    )
    assert x["destroyed_wrappers"]["FlowPolicy"] == (
        4922 if workload == "dcbb96" else 3558
    )
    assert x["destroyed_wrappers"]["FlowLedger"] == (
        2461 if workload == "dcbb96" else 593
    )
    assert abs(x["uncovered_process_cpu_s"]) < 0.01
    checks.append(
        {
            "lifetime": label,
            "outcome": "all policies/ledgers destroyed inside measured release; full export matches; process CPU and summed threads agree within 10 ms",
        }
    )
    events = [
        json.loads(line) for line in (folder / "events.jsonl").read_text().splitlines()
    ]
    previous = {}
    operations = defaultdict(lambda: {"calls": 0, "wall_s": 0.0, "cpu_s": 0.0})
    for event in events[:-1]:
        thread = event["thread"]
        assert event["wall_start_ns"] >= previous.get(thread, 0)
        previous[thread] = event["wall_stop_ns"]
        assert event["cpu_stop_ns"] >= event["cpu_start_ns"]
        op = operations[event["operation"]]
        op["calls"] += 1
        op["wall_s"] += (event["wall_stop_ns"] - event["wall_start_ns"]) / 1e9
        op["cpu_s"] += (event["cpu_stop_ns"] - event["cpu_start_ns"]) / 1e9
    assert (
        sum(op["calls"] for op in operations.values())
        == x["wrapper_calls"]
        == events[-1]["calls"]
    )
    outside_wall = sum(
        thread["wall_s"] - thread["native_span_wall_s"] for thread in x["threads"]
    )
    cal = {
        k: sorted(row[k] for row in x["calibration"])[1] for k in x["calibration"][0]
    }
    outside_wall_corrected = (
        outside_wall
        - x["wrapper_calls"]
        * (cal["cost_per_wrap_wall_ns"] - cal["recorded_null_wall_ns"])
        / 1e9
    )
    x["outside_span_wall_s"] = outside_wall
    x["corrected_outside_wall_s"] = outside_wall_corrected
    x["outside_wall_minus_cpu_raw_s"] = outside_wall - x["outside_span_cpu_s"]
    x["outside_wall_minus_cpu_corrected_s"] = (
        outside_wall_corrected - x["corrected_outside_cpu_s"]
    )
    for thread in x["threads"]:
        thread["outside_span_wall_s"] = thread["wall_s"] - thread["native_span_wall_s"]
        thread["outside_span_cpu_s"] = thread["cpu_s"] - thread["native_span_cpu_s"]
        thread["outside_wait_or_deschedule_raw_s"] = (
            thread["outside_span_wall_s"] - thread["outside_span_cpu_s"]
        )
    if workers > 1:
        for key in list(x):
            if "twice_core" in key:
                del x[key]
        x["finite_core_speedup_status"] = (
            "not identified from aggregate multi-worker sums; no interpolation from free-Core endpoint"
        )

    x["operations"] = dict(operations)
    x["held_wrapper_references"] = sum(
        row["held_wrappers"] for row in x["release_rows"]
    )
    calibration = read(BASE / "lifetime-overhead-01.json")
    x["estimated_holding_bookkeeping_s"] = (
        x["held_wrapper_references"]
        * calibration["holding_bookkeeping_ns_per_reference"]
        / 1e9
    )
    x["estimated_dict_clear_decref_s"] = (
        x["held_wrapper_references"] * calibration["dict_clear_ns_per_reference"] / 1e9
    )
    slots = defaultdict(lambda: defaultdict(float))
    for thread in x["threads"][1:]:
        slot = thread["name"].rsplit("_", 1)[1]
        for field in ("cpu_s", "native_span_cpu_s", "native_span_wall_s", "calls"):
            slots[slot][field] += thread[field]
        slots[slot]["cases"] += thread["call_counts"].get("FlowLedger.__init__", 0)
        slots[slot]["outside_span_cpu_s"] += thread["outside_span_cpu_s"]
        slots[slot]["outside_span_wall_s"] += thread["outside_span_wall_s"]
        slots[slot]["outside_wait_or_deschedule_raw_s"] += thread[
            "outside_wait_or_deschedule_raw_s"
        ]
    x["worker_slots_across_pools"] = dict(sorted(slots.items()))
    if workload == "dcbb96" and workers in (1, 8):
        prior = read(
            OLD / f"dcbb96-wall-{'1' if workers == 1 else 'default'}-01/summary.json"
        )
        rest = sum(prior["phases_s"].values()) - prior["phases_s"]["run"]
        x["pipeline_other_phases_prior_s"] = rest
        x["conditional_pipeline_free_core_ceiling"] = (rest + x["analysis_wall_s"]) / (
            rest + x["corrected_outside_cpu_s"]
        )
        if workers == 1:
            x["conditional_pipeline_twice_core_ceiling"] = (
                rest + x["analysis_wall_s"]
            ) / (rest + x["twice_core_wall_lower_bound_s"])
    lifetimes[label] = x

native_samples = {}
for workload in ("nsfnet", "backbone"):
    native = read(BASE / f"native-{workload}-analysis-03/summary.json")
    compare(
        OLD / "nsfnet-1-02/results.json"
        if workload == "nsfnet"
        else BASE / "backbone-wall-1-01/results.json",
        BASE / f"native-{workload}-outside-01/1/results.json",
    )
    assert sum(row["samples"] for row in native["categories"]) == native["samples"]
    assert native["samples"] > 10000
    native_samples[workload] = native
    checks.append(
        {
            "native_capture": workload,
            "outcome": "sample >10000; complete exclusive partition; native phase witnesses verified; named top-of-stack counts and full export match",
        }
    )

summary = {
    "exploratory": True,
    "identities": identities,
    "checks": checks,
    "backbone": {
        "inventory": read(BASE / "backbone-profile-1-01/inventory.json"),
        "wall_one": read(BASE / "backbone-wall-1-01/summary.json"),
        "wall_default": read(BASE / "backbone-wall-default-01/summary.json"),
        "actual_default_workers_by_step": default_workers,
        "profile_one": read(BASE / "backbone-profile-1-01/summary.json"),
        "profile_phases": profiles,
        "native_one": read(BASE / "backbone-native-1-01/native-timers.json"),
        "native_control": read(BASE / "backbone-native-1-01/summary.json"),
        "repeats": read(BASE / "backbone-trace-1-01/repeats.json"),
        "all_trace": {
            "calls": len(traces),
            "distinct": len(seen),
            "native_s": trace_total,
            "repeat_native_s": repeat_s,
            "repeat_share": repeat_s / trace_total,
        },
    },
    "getters": read(BASE / "nsfnet-getters-01/summary.json"),
    "critical_public_boundaries_only_superseded": critical,
    "lifetime_corrected": lifetimes,
    "native_sample_status": "two captures analysed; matched native batch denominator used",
    "native_samples": native_samples,
    "summary": "complete export/clock/poison checks passed; both native captures partitioned; destruction/getters and outside wall/CPU covered; multi-worker finite-speedup estimate omitted",
}
out = Path(sys.argv[1])
with out.open("x") as stream:
    json.dump(summary, stream, indent=2)
print(
    json.dumps(
        {
            "summary": summary["summary"],
            "checks": len(checks),
            "default_workers": default_workers,
        }
    ),
    flush=True,
)
