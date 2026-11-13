"""Historical ruling-18 diagnostic; current two-sided dispositions use review_accounting.py."""

import argparse
from collections import defaultdict
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "timing"))
from experiment import statistics_for  # noqa: E402


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("window", type=Path)
    parser.add_argument("analysis", type=Path)
    parser.add_argument("witnesses", type=Path)
    parser.add_argument("out", type=Path)
    args = parser.parse_args()
    window = json.loads(args.window.read_text())
    results = json.loads(args.analysis.read_text())["rows"]
    witnessed = defaultdict(set)
    details = defaultdict(list)
    for witness in json.loads(args.witnesses.read_text())["witnesses"]:
        selector = witness["raw_block"]
        key = (selector["row"], selector["experiment"])
        witnessed[key].add(selector["triple"])
        details[key].append(witness)
    output = []
    for index, plan in enumerate(window["plan"]["rows"]):
        for mode in ("replay", "cli", "matched"):
            result = next(
                row
                for row in results
                if row["workload"] == plan["workload"]
                and row["workers"] == plan["workers"]
                and row["experiment"] == mode
            )
            blocks = [
                b
                for b in window["blocks"]
                if b["row"] == index and b["experiment"] == mode
            ]
            triples = [
                [
                    next(b for b in blocks if b["triple"] == t and b["slot"] == label)
                    for label in ("A1", "B", "A2")
                ]
                for t in range(7)
            ]
            excluded = sorted(witnessed[(index, mode)])
            kept = [triple for t, triple in enumerate(triples) if t not in excluded]
            record = {
                "row": index,
                "workload": plan["workload"],
                "workers": plan["workers"],
                "experiment": mode,
                "excluded_triples": excluded,
                "witnesses": details[(index, mode)],
                "sensitivity": {},
            }
            if mode == "matched":
                for name, kind in result["kinds"].items():
                    ratio_ok = (
                        "ci95" in kind
                        and 0.95 <= kind["ci95"][0] <= kind["ci95"][1] <= 1.05
                    )
                    eligible = (
                        kind.get("counts_match", False)
                        and kind.get("release_counts_match", False)
                        and kind.get("qualified_noise", False)
                        and kind.get("null_fraction", 1) < 0.05
                        and not result["short_blocks"]
                        and not result["swap_changed"]
                        and not ratio_ok
                    )
                    if not eligible:
                        continue
                    filtered = (
                        statistics_for(
                            [
                                [
                                    {**b, "union_s": b["kinds"][name]["wall_s"]}
                                    for b in triple
                                ]
                                for triple in kept
                            ],
                            "union_s",
                        )
                        if kept
                        else None
                    )
                    # This is an explanatory diagnostic, not a new seven-triple qualification.
                    passes = (
                        filtered is not None
                        and filtered["reference_drift"] <= 0.03
                        and not any(filtered["clusters"].values())
                        and 0.95 <= filtered["ci95"][0] <= filtered["ci95"][1] <= 1.05
                    )
                    explained = bool(excluded) and (passes or len(excluded) > 3)
                    record["sensitivity"][name] = {
                        "original_ci95": kind["ci95"],
                        "filtered": filtered,
                        "remaining_triples": len(kept),
                        "diagnostic_pass": passes,
                        "explained_by_load": explained,
                        "disposition": "unqualified now; rerun once in a later window"
                        if explained
                        else "unqualified; ratio interval remains outside band",
                    }
            else:
                # The noise-only AA upper bound must resolve the row's declared bar.
                bar = 1.15 if mode == "replay" else 1 / 0.97
                failing = not result["qualified_noise"] or any(
                    metric["ci95"][1] >= bar for metric in result["metrics"].values()
                )
                record["detectable_gain_check"] = {
                    "noise_only_upper_below_bar": all(
                        metric["ci95"][1] < bar for metric in result["metrics"].values()
                    ),
                    "bar_speed_factor": bar,
                }
                if failing:
                    filtered = (
                        {
                            metric: statistics_for(kept, metric)
                            for metric in ("union_s", "elapsed_s")
                        }
                        if kept
                        else None
                    )
                    passes = (
                        filtered is not None
                        and all(
                            stat["reference_drift"] <= 0.03
                            and not any(stat["clusters"].values())
                            and stat["ci95"][1] < bar
                            for stat in filtered.values()
                        )
                        and not result["short_blocks"]
                    )
                    explained = (
                        bool(excluded)
                        and not result["short_blocks"]
                        and (passes or len(excluded) > 3)
                    )
                    record["sensitivity"]["row"] = {
                        "filtered": filtered,
                        "diagnostic_pass": passes,
                        "explained_by_load": explained,
                        "remaining_triples": len(kept),
                        "disposition": "unqualified now; rerun once in a later window"
                        if explained
                        else "original failure remains",
                    }
            output.append(record)
    args.out.mkdir()
    report = {
        "rows": output,
        "inputs": {
            "window": str(args.window),
            "analysis": str(args.analysis),
            "witnesses": str(args.witnesses),
        },
        "qualification": "original analysis verdicts retained; filtered results only select one future rerun",
        "summary": "completed ruling-18 load sensitivity",
    }
    with (args.out / "report.json").open("x") as stream:
        json.dump(report, stream, indent=2)
    cases = [case for row in output for case in row["sensitivity"].values()]
    print(
        json.dumps(
            {
                "witnesses": sum(len(row["witnesses"]) for row in output),
                "groups": sum(bool(row["witnesses"]) for row in output),
                "filtered_cases": len(cases),
                "rerun_once": sum(case["explained_by_load"] for case in cases),
                "summary": report["summary"],
            }
        )
    )


if __name__ == "__main__":
    main()
