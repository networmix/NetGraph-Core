"""Conserve clock intervals; never sum process CPU across overlapping calls."""

import argparse
import json
from pathlib import Path
import statistics

import numpy as np


def union(intervals):
    total, left, right = 0, None, None
    for a, b in sorted(intervals):
        assert b >= a
        if right is None or a > right:
            if right is not None:
                total += right - left
            left, right = a, b
        else:
            right = max(right, b)
    return total + (0 if right is None else right - left)


def complement(intervals, begin, end):
    cursor, gaps = begin, 0
    for a, b in sorted(intervals):
        assert begin <= a <= b <= end
        gaps += max(0, a - cursor)
        cursor = max(cursor, b)
    return gaps + end - cursor


def analyze(source, parallel_calibration=None):
    summary = json.loads((source / "summary.json").read_text())
    assert summary["byteorder"] == __import__("sys").byteorder
    phase_intervals = sorted(
        (p["start_wall_ns"], p["stop_wall_ns"]) for p in summary["phases"].values()
    )
    assert all(begin <= end for begin, end in phase_intervals), (
        "negative phase duration"
    )
    assert all(
        left[1] <= right[0] for left, right in zip(phase_intervals, phase_intervals[1:])
    ), "phase intervals overlap"
    events = np.fromfile(source / "events.i64", dtype=np.int64).reshape(-1, 17)
    vectors = np.fromfile(source / "python-vectors.i64", dtype=np.int64)
    assert len(events) == summary["events"] > 0
    assert (events[:, 4] >= events[:, 3]).all()
    assert (events[:, 6] >= events[:, 5]).all()
    assert (events[:, 14:17:2] % 2 == 0).all()
    assert max(events[:, 13] + events[:, 14]) <= len(vectors)
    assert max(events[:, 15] + events[:, 16]) <= len(vectors)
    cal = summary["calibration"]
    correction_cpu = statistics.median(
        c["full_wrap_cpu_ns"] - c["null_span_cpu_ns"] for c in cal
    )
    correction_wall = statistics.median(
        c["full_wrap_wall_ns"] - c["null_span_wall_ns"] for c in cal
    )
    corrections = [correction_cpu]
    if parallel_calibration:
        other = json.loads(parallel_calibration.read_text())
        corrections.extend(
            r["full_wrap_cpu_ns"] - r["null_span_cpu_ns"] for r in other["threads"]
        )
    phases = {}
    for index, (name, phase) in enumerate(summary["phases"].items()):
        rows = events[events[:, 1] == index]
        w0, w1 = phase["start_wall_ns"], phase["stop_wall_ns"]
        threads = []
        for thread in summary["threads"]:
            tid = thread["index"]
            mine = rows[rows[:, 2] == tid]
            if tid == 0:
                begin, end = w0, w1
                c0, c1 = phase["caller_start_cpu_ns"], phase["caller_stop_cpu_ns"]
            elif thread["start_wall"] >= w0 and thread["end_wall"] <= w1:
                begin, end = thread["start_wall"], thread["end_wall"]
                c0, c1 = thread["start_self"], thread["end_self"]
            elif thread["end_wall"] < w0 or thread["start_wall"] > w1:
                assert len(mine) == 0
                continue
            else:
                raise AssertionError(
                    ("thread straddles phase; needs explicit cuts", tid, name)
                )
            wi = [(int(r[3]), int(r[4])) for r in mine]
            ci = [(int(r[5]), int(r[6])) for r in mine]
            span_wall, span_cpu = union(wi), union(ci)
            outside_wall = complement(wi, begin, end)
            outside_cpu = complement(ci, c0, c1)
            assert span_wall + outside_wall == end - begin
            assert span_cpu + outside_cpu == c1 - c0
            # No nested public calls are allowed in this collector's outermost lane.
            assert span_wall == sum(b - a for a, b in wi)
            assert span_cpu == sum(b - a for a, b in ci)
            threads.append(
                {
                    "thread": tid,
                    "name": thread["name"],
                    "calls": len(mine),
                    "lifetime_wall_s": (end - begin) / 1e9,
                    "cpu_s": (c1 - c0) / 1e9,
                    "span_wall_s": span_wall / 1e9,
                    "span_cpu_s": span_cpu / 1e9,
                    "outside_wall_s": outside_wall / 1e9,
                    "outside_cpu_s": outside_cpu / 1e9,
                    "outside_not_running_s": (outside_wall - outside_cpu) / 1e9,
                    "corrected_outside_cpu_s": (
                        outside_cpu - len(mine) * correction_cpu
                    )
                    / 1e9,
                    "outside_cpu_calibration_sensitivity_s": [
                        (outside_cpu - len(mine) * c) / 1e9
                        for c in (max(corrections), min(corrections))
                    ],
                    "corrected_outside_wall_s": (
                        outside_wall - len(mine) * correction_wall
                    )
                    / 1e9,
                }
            )
        snap0, snap1 = summary["phase_snapshots"][name]
        python_cpu = snap1[4] - snap0[4]
        residual_min = snap1[2] - snap0[3] - python_cpu
        residual_max = snap1[3] - snap0[2] - python_cpu
        grouped = {}
        for operation, label in enumerate(summary["names"]):
            selected = rows[rows[:, 0] == operation]
            if len(selected):
                grouped[label] = {
                    "calls": len(selected),
                    "wall_s": int(np.sum(selected[:, 4] - selected[:, 3])) / 1e9,
                    "cpu_s": int(np.sum(selected[:, 6] - selected[:, 5])) / 1e9,
                }
        phases[name] = {
            "wall_s": phase["wall_s"],
            "global_span_union_s": union([(int(r[3]), int(r[4])) for r in rows]) / 1e9,
            "python_mach_cpu_s": python_cpu / 1e9,
            "non_python_cpu_residual_bracket_s": [
                residual_min / 1e9,
                residual_max / 1e9,
            ],
            "self_minus_mach_cpu_s": sum(t["cpu_s"] for t in threads)
            - python_cpu / 1e9,
            "threads": threads,
            "operations": grouped,
        }
    return {
        "exploratory": True,
        "source": str(source),
        "workers": summary["workers"],
        "native_sha256": summary["native_sha256"],
        "global_span_union_all_phases_s": sum(
            phase["global_span_union_s"] for phase in phases.values()
        ),
        "phases": phases,
        "calibration": cal,
        "summary": "PASS: all event/vector bounds, per-thread wall/CPU complement closure and non-overlap; residual computed once per phase, not summed across calls",
    }


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
