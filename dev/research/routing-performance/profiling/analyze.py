"""Compare every exported value and publish fresh survey tables from evidence."""

import copy
import json
import pstats
import sys
from pathlib import Path

base = Path(__file__).resolve().parents[1]
workloads = ["readme", "square", "nsfnet", "topogen", "dcbb96"]


def normalize(value, parent=None):
    # This is the complete allowed change for an explicit workers override.
    if isinstance(value, dict):
        return {
            key: (
                None
                if key == "workers" and parent == "metadata"
                else normalize(item, key)
            )
            for key, item in value.items()
        }
    if isinstance(value, list):
        return [normalize(item) for item in value]
    return value


def equal(a, b, path="root"):
    assert type(a) is type(b), (path, type(a), type(b))
    if isinstance(a, dict):
        assert list(a) == list(b), (path, "keys/order")
        for key in a:
            equal(a[key], b[key], f"{path}.{key}")
    elif isinstance(a, list):
        assert len(a) == len(b), (path, "length")
        for index, (x, y) in enumerate(zip(a, b, strict=True)):
            equal(x, y, f"{path}[{index}]")
    else:
        assert a == b, (path, a, b)


rows, checks = {}, []
for workload in workloads:
    one = base / f"{workload}-1-02"
    default = base / (
        "dcbb96-default-01" if workload == "dcbb96" else f"{workload}-default-02"
    )
    expected = normalize(json.loads((one / "results.json").read_text()))
    comparisons = [default, base / f"{workload}-trace-02"]
    comparisons += [
        p
        for p in (
            base / f"{workload}-wall-1-01",
            base / f"{workload}-wall-default-01",
            base / f"{workload}-native-1-01",
        )
        if (p / "results.json").exists()
    ]
    for folder in comparisons:
        candidate = normalize(json.loads((folder / "results.json").read_text()))
        equal(expected, candidate)
        checks.append(
            {
                "workload": workload,
                "comparison": folder.name,
                "outcome": "all exported values, types, keys and list order agree except metadata.workers",
            }
        )
    profiles = {}
    for phase in ("import", "load", "prepare", "run", "export"):
        stats = pstats.Stats(str(one / f"{phase}.pstats"))
        native, result_self = [], 0.0
        for key, value in stats.stats.items():
            if "_netgraph_core" in key[2]:
                native.append(
                    {"function": key[2], "calls": value[1], "self_s": value[2]}
                )
            if key[0].endswith(
                ("/analysis/costs.py", "/analysis/totals.py", "/results.py")
            ):
                result_self += value[2]
        profiles[phase] = {
            "core_builtin_rows": native,
            "known_python_result_self_s_lower_bound": result_self,
            "top_self_rows": [
                {
                    "file": k[0],
                    "line": k[1],
                    "function": k[2],
                    "calls": v[1],
                    "self_s": v[2],
                }
                for k, v in sorted(
                    stats.stats.items(), key=lambda row: row[1][2], reverse=True
                )[:20]
            ],
        }
    trace = base / f"{workload}-trace-02"
    calls = [
        json.loads(line) for line in (trace / "calls.jsonl").read_text().splitlines()
    ]
    calls = [call for call in calls if "operation" in call]
    # Static-support request keys derived from constructor inputs and current Core code.
    support_keys, support_repeats, seen, support_s = [], 0, set(), 0.0
    per_stage = {}
    for c in calls:
        group = per_stage.setdefault(c["stage"], {})
        group[c["operation"]] = group.get(c["operation"], 0) + 1
        if (
            c["operation"] == "FlowPolicy.__init__"
            and c["input"]["config"]["routing"]["type"] == "ShortestPaths"
        ):
            full = c["input"]
            key = json.dumps(
                [
                    full["ledger"]["graph"],
                    full["destinations"],
                    "REVERSE",
                    full["eligibility"],
                ],
                sort_keys=True,
            )
            support_keys.append(key)
            if key in seen:
                support_repeats += 1
                support_s += c["native_s"]
            seen.add(key)
    aggregate_seen = set()
    repeat_s = total_s = 0.0
    for c in calls:
        key = (c["operation"], c["full"])
        total_s += c["native_s"]
        if key in aggregate_seen:
            repeat_s += c["native_s"]
        aggregate_seen.add(key)
    rows[workload] = {
        "questions": [
            {"step": step["name"], "question": step["question"]}
            for step in expected["steps"]
        ],
        "profile_one": json.loads((one / "summary.json").read_text()),
        "profile_default": json.loads((default / "summary.json").read_text()),
        "inventory": json.loads((one / "inventory.json").read_text()),
        "profiles": profiles,
        "repeats": json.loads((trace / "repeats.json").read_text()),
        "stage_call_counts": per_stage,
        "static_support_requests_modelled_from_constructor_inputs": {
            "calls": len(support_keys),
            "distinct_keys": len(seen),
            "repeats": support_repeats,
            "constructor_time_on_repeated_keys_s": support_s,
            "warning": "constructor span includes more than static search; no search-kernel timing",
        },
        "all_traced_operations": {
            "calls": len(calls),
            "distinct_inputs": len(aggregate_seen),
            "native_span_s": total_s,
            "repeat_native_span_s": repeat_s,
            "repeat_native_span_share": repeat_s / total_s if total_s else 0,
        },
    }
    for workers in ("1", "default"):
        path = base / f"{workload}-wall-{workers}-01/summary.json"
        if path.exists():
            rows[workload][f"wall_{workers}"] = json.loads(path.read_text())
    native_path = base / f"{workload}-native-1-01/native-timers.json"
    if native_path.exists():
        rows[workload]["native_timers"] = json.loads(native_path.read_text())

# Demonstrate that the comparator rejects a changed computed value.
poison = copy.deepcopy(expected)
poison["steps"][0]["result"]["best_rows"][0]["result"]["delivered_units"] += 1
try:
    equal(expected, poison)
except AssertionError:
    checks.append({"control": "changed computed result", "outcome": "rejected"})
else:
    raise AssertionError("comparison accepted poisoned result")
with (base / (sys.argv[1] if len(sys.argv) > 1 else "analysis-01.json")).open(
    "x"
) as stream:
    stream.write(
        json.dumps(
            {
                "workloads": rows,
                "checks": checks,
                "summary": "all requested export comparisons passed; poison rejected",
            },
            indent=2,
        )
        + "\n"
    )
print(
    json.dumps(
        {"summary": "analysis completed", "workloads": len(rows), "checks": len(checks)}
    )
)
