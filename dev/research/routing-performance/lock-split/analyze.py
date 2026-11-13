"""Reconcile outer spans and native TLS deltas without mixing wall and CPU clocks."""

import argparse
import json
from pathlib import Path
import statistics

import numpy as np

from reuse import load

FIELDS = (
    "released_wall",
    "released_cpu",
    "retake_wall",
    "retake_cpu",
    "detaches",
    "destructors",
)


def delta(before, after):
    assert len(before) == len(after) == 7 and before[6] == after[6] == 0
    values = [int(y) - int(x) for x, y in zip(before[:6], after[:6], strict=True)]
    assert min(values) >= 0, "native counters decreased"
    assert values[5] <= values[4], "destructor count exceeds detaches"
    assert values[4] or not any(values), "elapsed native work without a detach count"
    return values


def analyze(source, parallel_calibration=None):
    s = json.loads((source / "summary.json").read_text())
    report = load("analyze_clocks").analyze(source, parallel_calibration)
    events = np.fromfile(source / "events.i64", dtype=np.int64).reshape(-1, 17)
    snapshots = np.fromfile(source / "detached.i64", dtype=np.int64).reshape(-1, 14)
    assert len(events) == len(snapshots) > 0
    native_cal = (
        np.asarray(s["native_clock_calibration"], dtype=np.float64)
        / s["native_clock_calibration_calls"]
    )
    native_median = np.median(native_cal, axis=0)
    null_w = statistics.median(row["null_span_wall_ns"] for row in s["calibration"])
    null_c = statistics.median(row["null_span_cpu_ns"] for row in s["calibration"])

    def parts(rows, values):
        native = np.sum(values, axis=0).tolist() if len(values) else [0] * 6
        span_w, span_c = (
            int(np.sum(rows[:, b] - rows[:, a])) for a, b in ((3, 4), (5, 6))
        )
        held = [span_w - native[0] - native[2], span_c - native[1] - native[3]]
        # Approximate clock-only subtraction, not a bound on scheduling bias.
        correction = [native[4] * v for v in native_median]
        raw = [held[0], held[1], *native[:4]]
        overhead = [
            len(rows) * null_w + correction[4] - correction[0] - correction[2],
            len(rows) * null_c + correction[5] - correction[1] - correction[3],
            *correction[:4],
        ]
        return {
            "calls": len(rows),
            "detaches": native[4],
            "destructors": native[5],
            "span_wall_s": span_w / 1e9,
            "span_cpu_s": span_c / 1e9,
            "raw_s": dict(
                zip(
                    ("held_wall", "held_cpu", *FIELDS[:4]),
                    (v / 1e9 for v in raw),
                    strict=True,
                )
            ),
            "corrected_s": dict(
                zip(
                    ("held_wall", "held_cpu", *FIELDS[:4]),
                    ((v - c) / 1e9 for v, c in zip(raw, overhead, strict=True)),
                    strict=True,
                )
            ),
            "clock_overhead_s": {
                "wall": (len(rows) * null_w + correction[4]) / 1e9,
                "cpu": (len(rows) * null_c + correction[5]) / 1e9,
            },
        }

    all_deltas = []
    for before_after in snapshots:
        all_deltas.append(delta(before_after[:7], before_after[7:]))
    all_deltas = np.asarray(all_deltas, dtype=np.int64)
    for index, (name, phase) in enumerate(report["phases"].items()):
        selected = events[:, 1] == index
        rows, ns, values = events[selected], snapshots[selected], all_deltas[selected]
        phase["split"] = parts(rows, values)
        # Endpoint-read allowance from observed clock-only sequences, not a GIL-wait bound.
        allowance = max(1000, int(native_cal[:, 4:].max() * 4))
        assert np.all(
            rows[:, 4] - rows[:, 3] - values[:, 0] - values[:, 2] >= -allowance
        )
        assert np.all(
            rows[:, 6] - rows[:, 5] - values[:, 1] - values[:, 3] >= -allowance
        )
        for thread in phase["threads"]:
            tid = thread["thread"]
            original = s["threads"][tid]
            start, stop = (
                (s["phases"][name]["split_start"], s["phases"][name]["split_stop"])
                if tid == 0
                else (original["split_start"], original["split_stop"])
            )
            mine = rows[:, 2] == tid
            cursor, gaps = start, np.zeros(6, dtype=np.int64)
            for pair in ns[mine]:
                gaps += delta(cursor, pair[:7])
                cursor = pair[7:]
            gaps += delta(cursor, stop)
            total = delta(start, stop)
            assert np.array_equal(np.sum(values[mine], axis=0) + gaps, total), (
                "TLS closure"
            )
            thread["split"] = parts(rows[mine], values[mine])
            thread["outside_detached_ns_counts"] = dict(
                zip(FIELDS, (int(v) for v in gaps), strict=True)
            )
            raw = thread["split"]["raw_s"]
            thread["total_clock_overhead_s"] = {}
            for clock in ("wall", "cpu"):
                inside = (
                    raw["held_" + clock]
                    + raw["released_" + clock]
                    + raw["retake_" + clock]
                )
                lifetime = (
                    thread["lifetime_wall_s"] if clock == "wall" else thread["cpu_s"]
                )
                assert abs(inside + thread["outside_" + clock + "_s"] - lifetime) < 1e-8
                corrected = thread["split"]["corrected_s"]
                overhead = (
                    thread["split"]["clock_overhead_s"][clock]
                    + thread["outside_" + clock + "_s"]
                    - thread["corrected_outside_" + clock + "_s"]
                )
                thread["total_clock_overhead_s"][clock] = overhead
                corrected_inside = sum(
                    corrected[part + "_" + clock]
                    for part in ("held", "released", "retake")
                )
                assert (
                    abs(
                        corrected_inside
                        + thread["corrected_outside_" + clock + "_s"]
                        + overhead
                        - lifetime
                    )
                    < 1e-8
                )
        phase["call_splits"] = {
            label: parts(rows[rows[:, 0] == operation], values[rows[:, 0] == operation])
            for operation, label in enumerate(s["names"])
            if np.any(rows[:, 0] == operation)
        }
        phase["top_five_by_raw_wall"] = sorted(
            phase["call_splits"],
            key=lambda label: phase["call_splits"][label]["span_wall_s"],
            reverse=True,
        )[:5]
    report["all_phase_split"] = parts(events, all_deltas)
    report["native_clock_calibration_per_detach_ns"] = {
        "median": native_median.tolist(),
        "min": np.min(native_cal, axis=0).tolist(),
        "max": np.max(native_cal, axis=0).tolist(),
    }
    report["clock_resolution_s"] = {
        name: __import__("time").get_clock_info(name).resolution
        for name in ("perf_counter", "thread_time")
    }
    report["summary"] = (
        "PASS: outer wall/CPU complements and native TLS in-span/outside-span closure; raw and approximate clock-corrected splits retained"
    )
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    parser.add_argument("out", type=Path)
    parser.add_argument("--parallel-calibration", type=Path)
    args = parser.parse_args()
    report = analyze(args.source, args.parallel_calibration)
    with args.out.open("x") as stream:
        json.dump(report, stream, indent=2)
        stream.write("\n")
    print(report["summary"])
