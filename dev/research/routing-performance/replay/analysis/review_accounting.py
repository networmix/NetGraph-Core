"""Read-only window-02 accounting and symmetric ruling-20 dispositions. Window-02 only: its filtered verdicts
omit the analyser's union/elapsed direction veto and take the matched short/swap vetoes from all seven
triples (REVIEW-T1 Opus N2; no verdict of window-02 changes). A later window computes ruling 20's filtered
verdicts with the analyser's own verdict code on the kept triples, vetoes included."""

from pathlib import Path
import json
import statistics
import collections
import sys

base = Path("/Users/networmix/ws/work/2026-11-core-routing-research.noindex")
root = base / "scratch/replay"
task = base / "tasks/replay"
sys.path.insert(0, str(task / "timing"))
from experiment import statistics_for  # noqa: E402

w = json.loads((root / "timing-window-02/run/window.json").read_text())
a = json.loads((root / "timing-analysis-02/analysis.json").read_text())
wit = json.loads((root / "timing-analysis-02/quiet-window-witnesses.json").read_text())[
    "witnesses"
]
excluded = collections.defaultdict(set)
for x in wit:
    excluded[x["raw_block"]["row"], x["raw_block"]["experiment"]].add(
        x["raw_block"]["triple"]
    )
rows = []
changes = []


def refs_ok(s):
    return s["reference_drift"] <= 0.03 and not any(s["clusters"].values())


for i, p in enumerate(w["plan"]["rows"]):
    original = {
        r["experiment"]: r
        for r in a["rows"]
        if r["workload"] == p["workload"] and r["workers"] == p["workers"]
    }
    group = original["matched"]
    sets = {}
    sensitivity = {}
    for mode in ("replay", "cli", "matched"):
        bb = [b for b in w["blocks"] if b["row"] == i and b["experiment"] == mode]
        triples = [
            [
                next(b for b in bb if b["triple"] == t and b["slot"] == label)
                for label in ("A1", "B", "A2")
            ]
            for t in range(7)
        ]
        bad = excluded[i, mode]
        kept = [t for n, t in enumerate(triples) if n not in bad]
        majority = len(bad) > 3
        sensitivity[mode] = {
            "witnessed_triples": sorted(bad),
            "remaining": len(kept),
            "majority": majority,
            "cases": {},
        }
        if mode == "matched":
            for name, k in group["kinds"].items():
                final = bool(k["qualified"])
                reason = "original verdict retained"
                filtered = None
                diagnostic = final
                if bad and all(name in b["kinds"] for t in kept for b in t) and kept:
                    filtered = statistics_for(
                        [
                            [{**b, "union_s": b["kinds"][name]["wall_s"]} for b in t]
                            for t in kept
                        ],
                        "union_s",
                    )
                    minmean = min(
                        b["kinds"][name]["wall_s"] / b["kinds"][name]["calls"] * 1e9
                        for t in kept
                        for b in t
                    )
                    null = max(x["in_span_ns"] for x in p["null"]["rows"])
                    diagnostic = (
                        k.get("counts_match", False)
                        and k.get("release_counts_match", False)
                        and refs_ok(filtered)
                        and 0.95 <= filtered["ci95"][0] <= filtered["ci95"][1] <= 1.05
                        and null / minmean < 0.05
                        and not group["short_blocks"]
                        and not group["swap_changed"]
                    )
                changed = bool(bad) and diagnostic != bool(k["qualified"])
                if majority:
                    final = False
                    reason = "not judgeable in this window: more than three triples witnessed; rerun once with next window"
                elif changed:
                    final = False
                    reason = "not judgeable now: verdict changes without witnessed triples; rerun once with next window"
                    changes.append((i, mode, name))
                sensitivity[mode]["cases"][name] = {
                    "original_qualified": k["qualified"],
                    "filtered": filtered,
                    "diagnostic_pass": diagnostic,
                    "changed": changed,
                    "final_qualified": final,
                    "disposition": reason,
                }
            sets = {
                name
                for name, v in sensitivity[mode]["cases"].items()
                if v["final_qualified"]
            }
        elif bad:
            filtered = {
                key: statistics_for(kept, key) for key in ("union_s", "elapsed_s")
            }
            bar = 1.15 if mode == "replay" else 1 / 0.97
            diagnostic = (
                all(refs_ok(v) and v["ci95"][1] < bar for v in filtered.values())
                and not original[mode]["short_blocks"]
                and not original[mode]["swap_changed"]
            )
            old = original[mode]["qualified_noise"] and all(
                v["ci95"][1] < bar for v in original[mode]["metrics"].values()
            )
            changed = diagnostic != old
            sensitivity[mode]["cases"]["row"] = {
                "filtered": filtered,
                "original_pass": old,
                "diagnostic_pass": diagnostic,
                "changed": changed,
                "disposition": "not judgeable; rerun once with next window"
                if majority or changed
                else "original verdict retained",
            }
            if changed:
                changes.append((i, mode, "row"))
    replay = [
        b
        for b in w["blocks"]
        if b["row"] == i and b["experiment"] == "replay" and b["slot"] == "A1"
    ]
    null = max(x["in_span_ns"] for x in p["null"]["rows"])
    U = statistics.median(b["union_s"] / b["repeats"] for b in replay)
    spans = statistics.median(
        sum(v["calls"] for v in b["kinds"].values()) / b["repeats"] for b in replay
    )
    H = spans * null / 1e9
    h = H / U

    def own(factor):
        return (1 - h) / (1 / factor - h) if 1 / factor > h else None

    bounds = []
    for b in replay:
        approved = sum(v["wall_s"] for name, v in b["kinds"].items() if name in sets)
        unapproved = sum(
            v["wall_s"] for name, v in b["kinds"].items() if name not in sets
        )
        bounds.append(
            (max(0, 1 - unapproved / b["union_s"]), min(1, approved / b["union_s"]))
        )
    cli = original["cli"]["memory"]
    memory = original["replay"]["memory"]
    matched = [
        b
        for b in w["blocks"]
        if b["row"] == i and b["experiment"] == "matched" and b["slot"] != "B"
    ]
    plain = [b for b in w["blocks"] if b["row"] == i and b["experiment"] == "cli"]
    rows.append(
        {
            "row": i,
            "workload": p["workload"],
            "workers": p["workers"],
            "null_ns": null,
            "spans_per_sequence": spans,
            "reference_union_s": U,
            "clock_model_s": H,
            "clock_share": h,
            "own_gain_at_1_15": own(1.15),
            "own_gain_at_0_97": own(1 / 0.97),
            "own_gain_at_e0_shifted_1_15": own(
                original["replay"]["metrics"]["union_s"]["native_bar_gain"]
            ),
            "replay_peak_mib": memory["rss_a1_max_bytes"] / 2**20,
            "replay_allowance_increment_mib": (
                memory["rss_allowance_bytes"] - memory["rss_a1_max_bytes"]
            )
            / 2**20,
            "cli_peak_mib": cli["rss_a1_max_bytes"] / 2**20,
            "cli_allowance_increment_mib": (
                cli["rss_allowance_bytes"] - cli["rss_a1_max_bytes"]
            )
            / 2**20,
            "instrumented_consumer_to_cli_wall": statistics.median(
                b["fresh_wall_s"] / b["repeats"] for b in matched
            )
            / statistics.median(b["fresh_wall_s"] / b["repeats"] for b in plain),
            "approved_union_share_bounds": [
                min(v[0] for v in bounds),
                max(v[1] for v in bounds),
            ],
            "final_qualified_kinds": sorted(sets),
            "sensitivity": sensitivity,
        }
    )
out = (
    Path(sys.argv[1]).resolve()
    if len(sys.argv) > 1
    else root / "r19-r20-assessment-01.json"
)
with out.open("x") as f:
    json.dump(
        {
            "rows": rows,
            "verdict_changes_after_filter": changes,
            "summary": "ruling-19 accounting and two-sided ruling-20 diagnostic; timing data unchanged",
        },
        f,
        indent=2,
    )
print(
    json.dumps(
        {
            "changes": changes,
            "rows": [
                {
                    k: r[k]
                    for k in (
                        "workload",
                        "workers",
                        "clock_share",
                        "own_gain_at_1_15",
                        "replay_peak_mib",
                        "approved_union_share_bounds",
                        "final_qualified_kinds",
                    )
                }
                for r in rows
            ],
        },
        indent=2,
    )
)
