"""Recheck saved window outputs and emit compact report evidence, without timing."""

import argparse
import datetime
import filecmp
import json
from pathlib import Path
import statistics

from analyze_clocks import analyze
from run_workload import ROOT, equal
from sensitivity import SHA, stats
from aba import side_clusters  # noqa: E402 -- sensitivity selects the pinned perf tool


def main(out):
    scratch = ROOT / "scratch/clocks"
    folder = scratch / "measured-sensitivity-01"
    report = json.loads((folder / "summary.json").read_text())
    records = [
        json.loads(line) for line in (folder / "blocks.jsonl").read_text().splitlines()
    ]
    rows = [r for r in records if "stage" in r]
    assert len(rows) == 63 and "summary" in records[-1]
    order = [("cal", i, 0) for i in range(7)]
    for i in range(7):
        order.extend(
            ("sensitivity", i, d)
            for d in ((100000, 200000) if i % 2 == 0 else (200000, 100000))
        )
    expected_order = [
        (stage, i, d, slot) for stage, i, d in order for slot in ("A1", "B", "A2")
    ]
    assert [
        (r["stage"], r["repeat"], r["delta"], r["slot"]) for r in rows
    ] == expected_order
    first = folder / "cal-0-0-A1/results.json"
    export = json.loads(first.read_text())
    assert export["metadata"]["core_extension_sha256"] == SHA
    equal(
        json.loads((scratch / "pilot-01/plain/results.json").read_text()),
        export,
        identity=True,
    )
    helpers = {0: [], 100000: [], 200000: []}
    for row in rows:
        label = f"{row['stage']}-{row['repeat']}-{row['delta']}-{row['slot']}"
        source = folder / label
        summary = json.loads((source / "summary.json").read_text())
        actual_delta = row["delta"] if row["slot"] == "B" else 0
        assert summary["delta_ns"] == actual_delta
        assert (
            summary["workers"] == 8
            and summary["lane"] == "delay"
            and not summary["smoke"]
        )
        assert summary["native_sha256"] == SHA == row["native_sha256"]
        assert summary["parallel_call_count"] == row["calls"] == 49200
        assert filecmp.cmp(first, source / "results.json", shallow=False)
        helpers[actual_delta].append(summary["helper_elapsed_ns"] / 49200)
    recomputed = {}
    for metric in ("wall_s", "analysis_s", "parallel_s"):
        saved = (
            report["statistics"]
            if metric == "wall_s"
            else report["phase_statistics"][metric]
        )
        result = {}
        for delta in (0, 100000, 200000):
            triples = []
            for i in range(7):
                triple = [
                    r
                    for r in rows
                    if r["repeat"] == i
                    and r["delta"] == delta
                    and r["stage"] == ("cal" if delta == 0 else "sensitivity")
                ]
                assert len(triple) == 3
                triples.append(triple)
            calculated = json.loads(json.dumps(stats(triples, delta, metric)))
            assert all(saved[str(delta)][k] == v for k, v in calculated.items())
            for side in (
                [t[0][metric] for t in triples] + [t[2][metric] for t in triples],
                [t[1][metric] for t in triples],
            ):
                assert side_clusters(side, 0.03) is None
            result[str(delta)] = calculated
        recomputed[metric] = result
    states = [report["machine"][0]]
    high = []
    for triple, pair in zip(order, report["machine"][1:]):
        for boundary in ("before", "after"):
            state = pair[boundary]
            states.append(state)
            for line in state["cpu"].splitlines()[1:]:
                cols = line.strip().split(None, 2)
                if len(cols) != 3:
                    continue
                cpu = float(cols[1])
                if cpu > 50 and ("Chrome" in cols[2] or "photolibraryd" in cols[2]):
                    high.append(
                        {
                            "triple": triple,
                            "boundary": boundary,
                            "utc": state["utc"],
                            "pid": int(cols[0]),
                            "cpu_percent": cpu,
                            "command": cols[2],
                        }
                    )
    assert len(states) == 43 and all("AC Power" in s["power"] for s in states)
    assert any("Chrome" in item["command"] for item in high)
    one_worker = []
    for k in ("01", "02", "03"):
        source = scratch / f"p1-one-clocks-{k}"
        saved = json.loads((scratch / f"p1-one-analysis-{k}.json").read_text())
        assert analyze(source) == saved
        assert (
            saved["workers"] == 1
            and saved["native_sha256"]
            == "208487c7aed158fb418791c22e23d09332533595053c0e4810923b31244f7c84"
        )
        assert filecmp.cmp(
            source / "results.json",
            scratch / "p1-one-plain-01/results.json",
            shallow=False,
        )
        phase = saved["phases"]["run"]
        totals = {
            key: sum(t[key] for t in phase["threads"])
            for key in (
                "cpu_s",
                "span_cpu_s",
                "outside_cpu_s",
                "outside_wall_s",
                "outside_not_running_s",
                "corrected_outside_cpu_s",
                "corrected_outside_wall_s",
            )
        }
        one_worker.append(
            {
                "run": k,
                "N_all_phases_s": saved["global_span_union_all_phases_s"],
                "analysis_wall_s": phase["wall_s"],
                "analysis_union_s": phase["global_span_union_s"],
                "cpu_residual_bracket_s": phase["non_python_cpu_residual_bracket_s"],
                "totals": totals,
                "all_fields_reproduced": True,
                "full_export_byte_identical": True,
            }
        )
    evidence = {
        "created_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "statistics_recomputed": recomputed,
        "clusters": "none under pinned aba.py >3% gap / at least two blocks per cluster definition",
        "helper_median_ns_per_call": {
            str(k): statistics.median(v) for k, v in helpers.items()
        },
        "minimum_blocks_s": {
            key: min(r[key] for r in rows)
            for key in ("wall_s", "analysis_s", "parallel_s")
        },
        "machine": {
            "snapshots": len(states),
            "first_utc": states[0]["utc"],
            "last_utc": states[-1]["utc"],
            "all_AC": True,
            "one_minute_load_range": [
                min(s["load"][0] for s in states),
                max(s["load"][0] for s in states),
            ],
            "swap_usage_records": sorted({s["swap"].strip() for s in states}),
            "selected_background_cpu": high,
            "granularity": "before/after each triple, not every process block",
            "swap_io_counters_recorded": False,
        },
        "memory": {
            k: {
                n: report["statistics"][k][n]
                for n in ("peak_rss_A_B_bytes", "memory_not_over_budget")
            }
            for k in ("100000", "200000")
        },
        "one_worker": one_worker,
        "N_range_s": [
            min(r["N_all_phases_s"] for r in one_worker),
            max(r["N_all_phases_s"] for r in one_worker),
        ],
        "N_median_s": statistics.median(r["N_all_phases_s"] for r in one_worker),
        "strict_window_qualification": "not established: Chrome CPU spike at triple boundary, per-process machine snapshots and swap-in/out counters absent",
        "summary": "PASS: 63 export/call/build/order checks and all saved bootstrap statistics reproduced; three P1/1 analyses and exports reproduced; strict window qualification is not claimed",
    }
    with out.open("x") as stream:
        json.dump(evidence, stream, indent=2)
        stream.write("\n")
    print(evidence["summary"], flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("out", type=Path)
    main(parser.parse_args().out)
