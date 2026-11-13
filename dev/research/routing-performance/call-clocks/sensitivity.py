"""Coordinator timing window: seven A/A/A, then seven A/B/A at each P1 delta."""

import argparse
import datetime
import json
import os
from pathlib import Path
import statistics
import subprocess
import sys
import time

from pilot import SHA
from run_workload import ROOT, equal

sys.path.insert(0, str(ROOT / "worktrees/core-ref/dev/perf"))
from aba import bootstrap_interval  # noqa: E402


def machine():
    result = {
        "utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "load": os.getloadavg(),
    }
    for key, command in [
        ("power", ["/usr/bin/pmset", "-g", "batt"]),
        ("swap", ["/usr/sbin/sysctl", "vm.swapusage"]),
        ("cpu", ["/bin/ps", "-axo", "pid,pcpu,comm"]),
    ]:
        result[key] = subprocess.check_output(command, text=True)
    return result


def stats(triples, delta, metric="wall_s"):
    ratios, slopes, drift = [], [], []
    for a, b, c in triples:
        ref = (a[metric] + c[metric]) / 2
        ratios.append(b[metric] / ref)
        slopes.append((b[metric] - ref) / (delta / 1e6) if delta else b[metric] - ref)
        drift.append(abs(a[metric] / c[metric] - 1))
    return {
        "median_ratio": statistics.median(ratios),
        "ratio_ci95": bootstrap_interval(ratios, 10000, 20261003),
        "median_slope_s_per_ms"
        if delta
        else "median_AA_difference_s": statistics.median(slopes),
        "slope_ci95" if delta else "AA_difference_ci95_s": bootstrap_interval(
            slopes, 10000, 20261003
        ),
        "reference_drifts": drift,
        "max_reference_drift": max(drift),
        "median_reference_s": statistics.median(
            (a[metric] + c[metric]) / 2 for a, _, c in triples
        ),
    }


def run(args):
    pilot = json.loads((args.pilot / "summary.json").read_text())
    assert not pilot["over_45_minutes"], "pilot exceeds authorized projection"
    assert 0 < args.window_minutes <= 45
    args.out.mkdir()
    initial = machine()
    assert "AC Power" in initial["power"] and initial["load"][0] < 4, initial
    print(
        json.dumps(
            {
                "expected_window_s": pilot["projected_window_s"],
                "maximum_s": args.window_minutes * 60,
                "processes": 63,
            }
        ),
        flush=True,
    )
    deadline = time.monotonic() + args.window_minutes * 60
    expected = json.loads((args.pilot / "plain/results.json").read_text())
    rows, groups, states = [], {0: [], 100000: [], 200000: []}, [initial]
    order = [("cal", i, 0) for i in range(7)]
    for i in range(7):
        order.extend(
            ("sensitivity", i, d)
            for d in ((100000, 200000) if i % 2 == 0 else (200000, 100000))
        )
    with (args.out / "blocks.jsonl").open("x") as raw:
        for stage, repeat, delta in order:
            before = machine()
            triple = []
            for slot, exposure in (("A1", 0), ("B", delta), ("A2", 0)):
                label = f"{stage}-{repeat}-{delta}-{slot}"
                target = args.out / label
                command = [
                    str(ROOT / "scratch/clocks/venv-delay/bin/python"),
                    "-B",
                    str(Path(__file__).with_name("run_workload.py")),
                    "--workload",
                    "dcbb96",
                    "--workers",
                    "8",
                    "--lane",
                    "delay",
                    "--delta-ns",
                    str(exposure),
                    "--native-sha",
                    SHA,
                    "--out",
                    str(target),
                ]
                remaining = deadline - time.monotonic()
                assert remaining > 0, "window expired"
                with (args.out / (label + ".log")).open("x") as log:
                    tick = time.monotonic()
                    child = subprocess.run(
                        command,
                        stdout=log,
                        stderr=subprocess.STDOUT,
                        timeout=min(180, remaining),
                    )
                    wall = time.monotonic() - tick
                    log.write(f"\nSUMMARY exit={child.returncode} elapsed_s={wall}\n")
                assert child.returncode == 0
                # User-path process is over before equality/audits begin.
                actual = json.loads((target / "results.json").read_text())
                assert actual["metadata"]["core_extension_sha256"] == SHA
                equal(expected, actual, identity=True)
                summary = json.loads((target / "summary.json").read_text())
                row = {
                    "stage": stage,
                    "repeat": repeat,
                    "delta": delta,
                    "slot": slot,
                    "wall_s": wall,
                    "analysis_s": summary["phases"]["run"]["wall_s"],
                    "parallel_s": sum(
                        r["wall_s"] for r in summary["parallel_schedules"]
                    ),
                    "calls": summary["parallel_call_count"],
                    "peak_rss_bytes": summary["peak_rss_bytes"],
                    "full_export_equal": True,
                    "native_sha256": SHA,
                }
                assert wall >= 2, "block shorter than RESEARCH minimum"
                assert row["calls"] == pilot["rows"][1]["call_count"]
                raw.write(json.dumps(row) + "\n")
                raw.flush()
                triple.append(row)
                rows.append(row)
            after = machine()
            assert "AC Power" in after["power"] and time.monotonic() < deadline
            states.append({"before": before, "after": after})
            groups[delta].append(triple)
            if stage == "cal" and repeat == 6:
                calibration = stats(groups[0], 0)
                print(
                    json.dumps(
                        {
                            "calibration_before_B": calibration,
                            "empirical_response_floor_s": max(
                                abs(v - 1) for v in calibration["ratio_ci95"]
                            )
                            * calibration["median_reference_s"],
                        }
                    ),
                    flush=True,
                )
                if (
                    calibration["max_reference_drift"] > 0.03
                    or calibration["ratio_ci95"][0] < 0.97
                    or calibration["ratio_ci95"][1] > 1.03
                ):
                    report = {
                        "calibration": calibration,
                        "machine": states,
                        "summary": "INCONCLUSIVE: A/A calibration failed; no nonzero delta scheduled",
                    }
                    (args.out / "summary.json").write_text(
                        json.dumps(report, indent=2) + "\n"
                    )
                    raw.write(json.dumps({"summary": report["summary"]}) + "\n")
                    return
            print(
                json.dumps(
                    {
                        "completed": len(rows),
                        "remaining_seconds": deadline - time.monotonic(),
                    }
                ),
                flush=True,
            )
        raw.write(
            json.dumps(
                {"summary": "all 63 fresh-process blocks completed and exports equal"}
            )
            + "\n"
        )
    analysis = {str(delta): stats(triples, delta) for delta, triples in groups.items()}
    phase_analysis = {
        metric: {
            str(delta): stats(triples, delta, metric)
            for delta, triples in groups.items()
        }
        for metric in ("analysis_s", "parallel_s")
    }
    control = analysis["0"]
    response_floor = (
        max(abs(v - 1) for v in control["ratio_ci95"]) * control["median_reference_s"]
    )
    for delta in (100000, 200000):
        row = analysis[str(delta)]
        row["AA_response_floor_s"] = response_floor
        row["AA_slope_floor_s_per_ms"] = response_floor / (delta / 1e6)
        row["statistical_checks_pass"] = (
            max(row["max_reference_drift"], control["max_reference_drift"]) <= 0.03
            and control["ratio_ci95"][0] >= 0.97
            and control["ratio_ci95"][1] <= 1.03
        )
        triples = groups[delta]
        ref_rss = max(r["peak_rss_bytes"] for t in triples for r in (t[0], t[2]))
        b_rss = max(t[1]["peak_rss_bytes"] for t in triples)
        row["peak_rss_A_B_bytes"] = [ref_rss, b_rss]
        row["memory_not_over_budget"] = b_rss <= max(
            ref_rss * 1.10, ref_rss + 8 * 1024**2
        )
    report = {
        "exploratory_until_window_verified": True,
        "statistics": analysis,
        "phase_statistics": phase_analysis,
        "machine": states,
        "seed": 20261003,
        "resamples": 10000,
        "scope": "public calls only on P1 case-execution Python threads; zero traverses same GIL-free helper",
        "conditions": "coordinator must additionally confirm quiet reservation, no machine sleep and inspect per-block machine load/clusters before calling statistics qualified",
        "summary": "P1 added-latency series complete; no free-Core endpoint or optimization gain inferred",
    }
    (args.out / "summary.json").write_text(json.dumps(report, indent=2) + "\n")
    print(report["summary"], flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pilot", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--window-minutes", type=int, default=45)
    run(parser.parse_args())
