"""Check the retained coordinator outputs and extract report tables; runs no workload."""

from collections import Counter
import hashlib
import json
from pathlib import Path
import re
import statistics
import sys

import numpy as np

from analyze_sample import classify, parse
from reuse import ROOT
from run_workload import equal

TASK = Path(__file__).resolve().parent
SCRATCH = ROOT / "scratch/split"
CASES = (
    ("p1-1", "dcbb96", 1),
    ("p1-8", "dcbb96", 8),
    ("p2-1", "nsfnet", 1),
    ("p2-8", "nsfnet", 8),
    ("p3-1", "backbone", 1),
    ("p3-28", "backbone", 28),
)
PARTS = (
    "held_wall",
    "held_cpu",
    "released_wall",
    "released_cpu",
    "retake_wall",
    "retake_cpu",
)


def read(path):
    return json.loads(path.read_text())


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def check_parts(rows, native, recorded):
    span_w = int(np.sum(rows[:, 4] - rows[:, 3]))
    span_c = int(np.sum(rows[:, 6] - rows[:, 5]))
    sums = np.sum(native, axis=0).tolist()
    values = (span_w - sums[0] - sums[2], span_c - sums[1] - sums[3], *sums[:4])
    for key, value in zip(PARTS, values, strict=True):
        assert abs(recorded["raw_s"][key] - value / 1e9) < 1e-10, key
    assert (
        recorded["calls"] == len(rows)
        and recorded["detaches"] == sums[4]
        and recorded["destructors"] == sums[5]
    )
    for clock, total in (("wall", span_w), ("cpu", span_c)):
        corrected = sum(
            recorded["corrected_s"][part + "_" + clock]
            for part in ("held", "released", "retake")
        )
        assert abs(corrected + recorded["clock_overhead_s"][clock] - total / 1e9) < 1e-9


def clock_row(tag, workload, workers, identity):
    source = SCRATCH / f"measured-{tag}-01"
    index, report, summary = (
        read(source / name)
        for name in ("index.json", "analysis.json", "clocks/summary.json")
    )
    plain = read(source / "plain/summary.json")
    assert (
        index["workload"] == workload
        and index["workers"] == workers
        and index["smoke"] is False
    )
    assert index["summary"].startswith("PASS:") and summary["smoke"] is False
    assert (
        summary["native_sha256"]
        == report["native_sha256"]
        == identity["scratch_sha256"]
    )
    assert plain["native_sha256"] == identity["reference_sha256"]
    expected, actual = (
        read(source / "plain/results.json"),
        read(source / "clocks/results.json"),
    )
    assert expected["metadata"]["core_extension_sha256"] == identity["reference_sha256"]
    assert actual["metadata"]["core_extension_sha256"] == identity["scratch_sha256"]
    equal(expected, actual, identity=True)
    del expected, actual
    assert (
        summary["byteorder"] == sys.byteorder
        and summary["event_width"] == 17
        and summary["detached_width"] == 14
    )
    events = np.memmap(source / "clocks/events.i64", mode="r", dtype=np.int64).reshape(
        -1, 17
    )
    snapshots = np.memmap(
        source / "clocks/detached.i64", mode="r", dtype=np.int64
    ).reshape(-1, 14)
    assert len(events) == len(snapshots) == summary["events"] > 0
    assert not np.any(snapshots[:, (6, 13)])
    native = snapshots[:, 7:13] - snapshots[:, :6]
    assert np.all(native >= 0) and np.all(native[:, 5] <= native[:, 4])
    assert np.all((native[:, 4] != 0) | np.all(native == 0, axis=1))
    check_parts(events, native, report["all_phase_split"])
    max_raw_error = max_corrected_error = 0.0
    outside_detaches = phase_threads = negative_held = 0
    minima = {"wall": None, "cpu": None}
    for phase_id, (name, phase) in enumerate(report["phases"].items()):
        chosen = events[:, 1] == phase_id
        rows, ns, nd = events[chosen], snapshots[chosen], native[chosen]
        check_parts(rows, nd, phase["split"])
        for clock, left, right, r_col, a_col in (
            ("wall", 3, 4, 0, 2),
            ("cpu", 5, 6, 1, 3),
        ):
            held = rows[:, right] - rows[:, left] - nd[:, r_col] - nd[:, a_col]
            if len(held):
                minima[clock] = (
                    int(np.min(held))
                    if minima[clock] is None
                    else min(minima[clock], int(np.min(held)))
                )
                negative_held += int(np.sum(held < 0))
        for thread in phase["threads"]:
            phase_threads += 1
            tid = thread["thread"]
            original = summary["threads"][tid]
            if tid == 0:
                p = summary["phases"][name]
                start_w, end_w, start_c, end_c = (
                    p["start_wall_ns"],
                    p["stop_wall_ns"],
                    p["caller_start_cpu_ns"],
                    p["caller_stop_cpu_ns"],
                )
                start_n, end_n = p["split_start"], p["split_stop"]
            else:
                start_w, end_w, start_c, end_c = (
                    original[key]
                    for key in ("start_wall", "end_wall", "start_self", "end_self")
                )
                start_n, end_n = original["split_start"], original["split_stop"]
            mask = rows[:, 2] == tid
            mine, pairs, deltas = rows[mask], ns[mask], nd[mask]
            check_parts(mine, deltas, thread["split"])
            cursor_n, gaps_n = np.array(start_n[:6]), np.zeros(6, dtype=np.int64)
            assert start_n[6] == end_n[6] == 0
            for pair in pairs:
                gap = pair[:6] - cursor_n
                assert np.all(gap >= 0)
                gaps_n += gap
                cursor_n = pair[7:13]
            tail = np.array(end_n[:6]) - cursor_n
            assert np.all(tail >= 0)
            gaps_n += tail
            assert np.array_equal(
                np.sum(deltas, axis=0) + gaps_n,
                np.array(end_n[:6]) - np.array(start_n[:6]),
            )
            outside_detaches += int(gaps_n[4])
            for clock, left, right, begin, end in (
                ("wall", 3, 4, start_w, end_w),
                ("cpu", 5, 6, start_c, end_c),
            ):
                assert np.all(mine[:, left] >= begin) and np.all(mine[:, right] <= end)
                assert np.all(mine[:, right] >= mine[:, left])
                assert len(mine) < 2 or np.all(mine[1:, left] >= mine[:-1, right])
                previous = np.r_[begin, mine[:, right]]
                following = np.r_[mine[:, left], end]
                outside = int(np.sum(following - previous))
                spans = int(np.sum(mine[:, right] - mine[:, left]))
                assert outside + spans == end - begin
                assert abs(outside / 1e9 - thread["outside_" + clock + "_s"]) < 1e-9
                raw = thread["split"]["raw_s"]
                error = abs(
                    sum(
                        raw[part + "_" + clock]
                        for part in ("held", "released", "retake")
                    )
                    + outside / 1e9
                    - (end - begin) / 1e9
                )
                max_raw_error = max(max_raw_error, error)
                corrected = thread["split"]["corrected_s"]
                error = abs(
                    sum(
                        corrected[part + "_" + clock]
                        for part in ("held", "released", "retake")
                    )
                    + thread["corrected_outside_" + clock + "_s"]
                    + thread["total_clock_overhead_s"][clock]
                    - (end - begin) / 1e9
                )
                max_corrected_error = max(max_corrected_error, error)
        for operation, label in enumerate(summary["names"]):
            mask = rows[:, 0] == operation
            if np.any(mask):
                check_parts(rows[mask], nd[mask], phase["call_splits"][label])
    assert max_raw_error < 1e-8 and max_corrected_error < 1e-8
    run = report["phases"]["run"]
    top = sorted(
        run["call_splits"],
        key=lambda k: run["call_splits"][k]["span_wall_s"],
        reverse=True,
    )[:5]
    assert top == run["top_five_by_raw_wall"]
    counts = report["all_phase_split"]
    cal = summary["calibration"]
    return {
        "tag": tag.upper().replace("-", "/"),
        "source": str(source),
        "workers": workers,
        "run": run,
        "all_phase_split": counts,
        "plain_run_wall_s": plain["phases"]["run"]["wall_s"],
        "plain_phase_wall_s": {k: v["wall_s"] for k, v in plain["phases"].items()},
        "top_five": {k: run["call_splits"][k] for k in top},
        "calibration": {
            "null_span_wall_ns": statistics.median(c["null_span_wall_ns"] for c in cal),
            "null_span_cpu_ns": statistics.median(c["null_span_cpu_ns"] for c in cal),
            "full_wrapper_cpu_ns": statistics.median(
                c["full_wrap_cpu_ns"] for c in cal
            ),
            "matched_parallel_cpu_correction_us": [
                min(
                    t["full_wrap_cpu_ns"] - t["null_span_cpu_ns"]
                    for t in read(source / "parallel-null.json")["threads"]
                )
                / 1000,
                max(
                    t["full_wrap_cpu_ns"] - t["null_span_cpu_ns"]
                    for t in read(source / "parallel-null.json")["threads"]
                )
                / 1000,
            ],
            "native_reads_ns": report["native_clock_calibration_per_detach_ns"],
        },
        "closure": {
            "integer_error_ns": 0,
            "max_raw_float_error_ns": max_raw_error * 1e9,
            "max_corrected_float_error_ns": max_corrected_error * 1e9,
            "phase_thread_rows": phase_threads,
            "outside_detaches": outside_detaches,
            "negative_raw_held_event_clock_pairs": negative_held,
            "minimum_raw_held_ns": minima,
        },
        "export_equal": True,
    }


def sample_row():
    source = SCRATCH / "p1-8-sample-01"
    report = read(SCRATCH / "p1-8-sample-analysis-02/summary.json")
    index, worker = read(source / "index.json"), read(source / "worker.json")
    assert (
        index["sample_exit"] == index["worker_exit"] == 0 and index["full_export_equal"]
    )
    exports = sorted(source.glob("results-*.json"))
    assert (
        len(exports) == index["exports_compared"] == worker["complete_workflows"] == 5
    )
    expected = read(SCRATCH / "measured-p1-8-01/plain/results.json")
    for path in (*exports, source / "results.json"):
        equal(expected, read(path))
    text = (source / "sample.txt").read_text()
    assert sha(source / "sample.txt") == report["sample_sha256"]
    nodes = parse(text)
    roots = [n for n in nodes if n.parent is None]
    totals, native, functions, native_roots = Counter(), Counter(), Counter(), set()
    for i, node in enumerate(nodes):
        weight = node.count - sum(nodes[c].count for c in node.children)
        if not weight:
            continue
        chain, cursor = [], i
        while cursor is not None:
            chain.append(nodes[cursor])
            cursor = nodes[cursor].parent
        chain.reverse()
        category, is_native = classify([n.frame for n in chain])
        totals[category] += weight
        if is_native:
            native[category] += weight
            native_roots.add(chain[0].line)
        functions[node.frame.split("  (in ", 1)[0]] += weight
    assert dict(totals) == report["categories"]
    assert (
        sum(totals.values())
        == sum(n.count for n in roots)
        == sum(r["samples"] for r in report["roots"])
    )
    assert len(roots) == len(report["roots"]) and len(native_roots) > 1
    denominator = sum(native.values()) - native["blocked_wait"]
    assert denominator == report["core_samples_excluding_wait"]
    assert sum(native.values()) == report["core_samples_including_wait"]
    for category, percent in report["native_phase_share_percent"].items():
        assert abs(100 * native[category] / denominator - percent) < 1e-10
    top = text.split("Sort by top of stack, same collapsed", 1)[1].split(
        "Binary Images:", 1
    )[0]
    for line in top.splitlines():
        match = re.fullmatch(r"\s*(netgraph::.+?)  \(in .+?\)\s+(\d+)\s*", line)
        if match:
            assert (
                functions[match[1]] == int(match[2]) == report["top_checks"][match[1]]
            )
    return {
        "report": report,
        "worker": worker,
        "all_samples": sum(totals.values()),
        "native_roots": len(native_roots),
        "exports_equal": len(exports),
        "summary": "PASS: raw all-root sample conservation, categories, top counts, hash and all five exports",
    }


def main(out):
    identity = read(SCRATCH / "build-identity.json")
    manifest = read(TASK / "stage2-checks.json")
    for name, record in manifest["files"].items():
        if Path(name).suffix in (".py", ".sh", ".patch"):
            assert sha(TASK / name) == record["sha256"], (
                "instrument changed after stage 2",
                name,
            )
    for label, prefix in (
        ("reference", ROOT / "scratch/consumer/ref/venv"),
        ("scratch", SCRATCH / "venv"),
    ):
        (extension,) = (prefix / "lib/python3.14/site-packages").glob(
            "_netgraph_core*.so"
        )
        assert sha(extension) == identity[label + "_sha256"]
    rows = [clock_row(*case, identity) for case in CASES]
    sample = sample_row()
    old = read(ROOT / "scratch/orient-consumers/native-stack-analysis-01/summary.json")
    assert sha(Path(old["source_sample"])) == old["sample_sha256"]
    record = {
        "identity": identity,
        "rows": rows,
        "sample": sample,
        "carried_p1_one": old,
        "summary": "PASS: unchanged instruments/binaries, all six full exports, raw/derived per-call and per-thread splits, exact integer closure, corrected closure and five sampled exports verified; no workload or capture rerun",
    }
    with out.open("x") as stream:
        json.dump(record, stream, indent=2)
        stream.write("\n")
    print(record["summary"], flush=True)


if __name__ == "__main__":
    main(Path(sys.argv[1]))
