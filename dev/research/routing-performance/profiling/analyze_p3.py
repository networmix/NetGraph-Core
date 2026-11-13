"""P3 1/8/auto lifecycle clocks, plain controls and conditional free-Core endpoints."""

import argparse
import json
import statistics
from collections import Counter
from pathlib import Path

from evidence import equal, normalize

BASE = Path(__file__).resolve().parents[1]
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("out", type=Path)
args = parser.parse_args()
expected = normalize(json.loads((BASE / "backbone-wall-1-01/results.json").read_text()))
rows = []
checks = []
for workers in (1, 8, 28):
    diag = BASE / f"backbone-p3-lifetime-{workers}-01"
    plain = BASE / f"backbone-p3-wall-{workers}-01"
    value = json.loads((diag / "summary.json").read_text())
    control = json.loads((plain / "summary.json").read_text())
    for folder in (diag, plain):
        full = json.loads((folder / "results.json").read_text())
        equal(expected, normalize(full))
        placement = full["steps"][2]["result"]
        assert placement["metadata"]["workers"] == workers
        checks.append(
            {
                "export": str(folder / "results.json"),
                "outcome": "all values/types/keys/order match accepted complete backbone export except metadata.workers; actual worker count verified",
            }
        )
    assert value["destroyed_wrappers"]["FlowPolicy"] == 3558
    assert value["destroyed_wrappers"]["FlowLedger"] == 593
    assert value["wrapper_calls"] == 60723
    assert abs(value["process_cpu_s"] - value["total_thread_cpu_s"]) < 0.01
    events = [
        json.loads(line) for line in (diag / "events.jsonl").read_text().splitlines()
    ]
    previous = {}
    calls = Counter()
    for event in events[:-1]:
        thread = event["thread"]
        assert event["wall_start_ns"] >= previous.get(thread, 0)
        assert event["wall_stop_ns"] >= event["wall_start_ns"]
        assert event["cpu_stop_ns"] >= event["cpu_start_ns"]
        previous[thread] = event["wall_stop_ns"]
        calls[thread] += 1
    assert sum(calls.values()) == value["wrapper_calls"] == events[-1]["calls"]
    for thread in value["threads"]:
        assert calls[thread["thread_index"]] == thread["calls"]
        thread["outside_wall_s"] = thread["wall_s"] - thread["native_span_wall_s"]
        thread["outside_cpu_s"] = thread["cpu_s"] - thread["native_span_cpu_s"]
        thread["outside_wait_or_deschedule_raw_s"] = (
            thread["outside_wall_s"] - thread["outside_cpu_s"]
        )
    cal = {
        key: statistics.median(c[key] for c in value["calibration"])
        for key in value["calibration"][0]
    }
    count = value["wrapper_calls"]
    outside_wall = sum(thread["outside_wall_s"] for thread in value["threads"])
    outside_cpu = sum(thread["outside_cpu_s"] for thread in value["threads"])
    corrected_wall = (
        outside_wall
        - count * (cal["cost_per_wrap_wall_ns"] - cal["recorded_null_wall_ns"]) / 1e9
    )
    corrected_cpu = (
        outside_cpu
        - count * (cal["cost_per_wrap_cpu_ns"] - cal["recorded_null_cpu_ns"]) / 1e9
    )
    assert abs(corrected_cpu - value["corrected_outside_cpu_s"]) < 1e-12
    assert corrected_cpu > 0
    assert abs(sum(value["phases_s"].values()) - value["pipeline_wall_s"]) < 1e-12
    diag_other = value["pipeline_wall_s"] - value["analysis_wall_s"]
    plain_pipe = sum(control["phases_s"].values())
    plain_run = control["phases_s"]["run"]
    plain_other = plain_pipe - plain_run
    holding = json.loads((BASE / "lifetime-overhead-01.json").read_text())
    held_refs = sum(row["held_wrappers"] for row in value["release_rows"])
    row = {
        "workers": workers,
        "exploratory": True,
        "instrumented": value,
        "plain_control": control,
        "outside_wall_raw_s": outside_wall,
        "outside_cpu_raw_s": outside_cpu,
        "outside_wall_corrected_s": corrected_wall,
        "outside_cpu_corrected_s": corrected_cpu,
        "outside_wait_or_deschedule_raw_s": outside_wall - outside_cpu,
        "outside_wait_or_deschedule_corrected_s": corrected_wall - corrected_cpu,
        "whole_wrapper_overhead_wall_s": count * cal["cost_per_wrap_wall_ns"] / 1e9,
        "whole_wrapper_overhead_cpu_s": count * cal["cost_per_wrap_cpu_ns"] / 1e9,
        "holding_bookkeeping_estimated_s": held_refs
        * holding["holding_bookkeeping_ns_per_reference"]
        / 1e9,
        "free_core_run_floor_s": corrected_cpu,
        "instrumented_pipeline_other_s": diag_other,
        "instrumented_free_core_run_ceiling": value["analysis_wall_s"] / corrected_cpu,
        "instrumented_free_core_pipeline_floor_s": diag_other + corrected_cpu,
        "instrumented_free_core_pipeline_ceiling": value["pipeline_wall_s"]
        / (diag_other + corrected_cpu),
        "plain_pipeline_wall_s": plain_pipe,
        "plain_analysis_wall_s": plain_run,
        "plain_pipeline_other_s": plain_other,
        "outside_cpu_limit_for_1_1_plain_s": plain_pipe / 1.1 - plain_other,
        "outside_cpu_headroom_for_1_1_plain_s": plain_pipe / 1.1
        - plain_other
        - corrected_cpu,
        "plain_free_core_run_ceiling": plain_run / corrected_cpu,
        "plain_free_core_pipeline_floor_s": plain_other + corrected_cpu,
        "plain_free_core_pipeline_ceiling": plain_pipe / (plain_other + corrected_cpu),
        "finite_core_speedup": "not inferred",
        "scope": "free endpoint conditional on unchanged outside work; diag pipeline includes timer overhead and inactive export wrappers; paired plain phase costs give a separate user-pipeline estimate; no subtraction of summed overhead from multi-worker elapsed wall",
    }
    rows.append(row)
    checks.append(
        {
            "workers": workers,
            "outcome": "60723 nonoverlapping call/release events conserved, all 3558 policies/593 ledgers destroyed inside spans, process/thread CPU within 10 ms, phase sum and null-wrap correction agree",
        }
    )
auto = rows[-1]
result = {
    "task": "ORIENT-ADD-P3",
    "exploratory": True,
    "rows": rows,
    "checks": checks,
    "auto_p3_threshold": 1.1,
    "auto_plain_conditional_ceiling_passes_threshold": auto[
        "plain_free_core_pipeline_ceiling"
    ]
    >= 1.1,
    "auto_instrumented_conditional_ceiling_passes_threshold": auto[
        "instrumented_free_core_pipeline_ceiling"
    ]
    >= 1.1,
    "summary": "full P3 1/8/auto lifecycle spans and measured pipeline phases checked against matched plain controls; only conditional free-Core endpoints derived",
}
with args.out.open("x") as stream:
    json.dump(result, stream, indent=2)
print(
    json.dumps(
        {
            "summary": result["summary"],
            "checks": len(checks),
            "rows": [
                {
                    k: r[k]
                    for k in (
                        "workers",
                        "plain_pipeline_wall_s",
                        "plain_analysis_wall_s",
                        "outside_cpu_corrected_s",
                        "instrumented_free_core_pipeline_ceiling",
                        "plain_free_core_pipeline_ceiling",
                    )
                }
                for r in rows
            ],
        }
    ),
    flush=True,
)
