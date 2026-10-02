"""Summarize fresh-process A/B/A CSVs; timings are milliseconds per whole batch."""

from __future__ import annotations

import csv
import sys
from pathlib import Path


def read(path: Path) -> dict[tuple[str, int, int, int], dict[str, dict[str, str]]]:
    cases: dict[tuple[str, int, int, int], dict[str, dict[str, str]]] = {}
    with path.open() as file:
        for row in csv.DictReader(file):
            assert row["correct"] == "pass", row
            key = (row["case"], int(row["n"]), int(row["e"]), int(row["queries"]))
            cases.setdefault(key, {})[row["method"]] = row
    return cases


def main() -> None:
    root = Path(sys.argv[1])
    a1, b, a2 = (read(root / name) for name in ("a1.csv", "b.csv", "a2.csv"))
    assert a1.keys() == b.keys() == a2.keys()
    print("Exploratory observations, not a production performance acceptance result.")
    print(
        "Inspect load-*.txt and gpu-load-*.txt; speedup claims require a quiet machine."
    )
    print("The prototype matches supported DAG edge sets, not CPU predecessor order.\n")
    print("All times are medians in milliseconds per batch, with nine timed samples.")
    print("CPU A columns use serial SPF for Q=1 and the 10-worker pool otherwise.")
    print("GPU includes host distance copies and serial CPU ECMP DAG reconstruction.")
    print("Upload and pipeline setup are excluded and recorded separately in raw data.")
    print("The table picks the fastest measured GPU method for each case.\n")
    print(
        "| Case | V / E | Q | CPU A1 | GPU B | CPU A2 | GPU method | "
        "Distance-only B | CPU drift |"
    )
    print("|---|---:|---:|---:|---:|---:|---|---:|---:|")
    for key, runs in b.items():
        name, n, e, q = key
        baseline = "cpu_full" if q == 1 else "cpu_pool10"
        left = float(a1[key][baseline]["wall_ms"])
        right = float(a2[key][baseline]["wall_ms"])
        best = min(runs.values(), key=lambda row: float(row["wall_ms"]))
        dist = float(best["distance_wall_ms"])
        drift = 100 * (right / left - 1)
        print(
            f"| {name} | {n:,} / {e:,} | {q} | {left:.3f} | "
            f"{float(best['wall_ms']):.3f} | {right:.3f} | {best['method']} | "
            f"{dist:.3f} | {drift:+.1f}% |"
        )
    print(
        "\nDistance-only timings stop before DAG reconstruction; they are not SPF API timings."
    )
    print(
        "CPU distance-only and serial full-batch controls are also available in a1/a2.csv."
    )
    print("\nDistances-only comparison (both sides omit DAG construction):\n")
    print("| Case | V / E | Q | CPU A1 | GPU B | CPU A2 | GPU method |")
    print("|---|---:|---:|---:|---:|---:|---|")
    for key, runs in b.items():
        name, n, e, q = key
        baseline = "cpu_distance" if q == 1 else "cpu_pool10_distance"
        left = float(a1[key][baseline]["wall_ms"])
        right = float(a2[key][baseline]["wall_ms"])
        best = min(runs.values(), key=lambda row: float(row["distance_wall_ms"]))
        print(
            f"| {name} | {n:,} / {e:,} | {q} | {left:.3f} | "
            f"{float(best['distance_wall_ms']):.3f} | {right:.3f} | {best['method']} |"
        )


if __name__ == "__main__":
    main()
