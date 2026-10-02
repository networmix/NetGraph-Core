#!/usr/bin/env python3
"""Tabulate cpu_queue A/B/A results as Markdown (sequential single-thread medians)."""

import csv
import glob
import os
import sys

root = (
    sys.argv[1]
    if len(sys.argv) > 1
    else os.path.join(os.path.dirname(__file__), "results", "first")
)


def load(pattern):
    rows = {}
    for path in sorted(glob.glob(os.path.join(root, pattern))):
        tag = (
            os.path.basename(path)[len("time-") : -4]
            if "time-" in path
            else os.path.basename(path)[len("sweep-") : -4]
        )
        with open(path) as f:
            for r in csv.DictReader(f):
                key = (r["case"], r["n"], r["e"], r["queries"], r["dst"])
                rows.setdefault(key, {})[tag] = r
    return rows


def fmt(r):
    if r is None:
        return "—"
    s = f"{float(r['median_ms']):.3f}"
    return s + (" (heap fallback)" if r.get("fallback") == "1" else "")


time_rows = load("time-*.csv")
order = [
    "ref-a1",
    "heap-fresh-b",
    "heap-full-b",
    "heap-sparse-b",
    "bucket-fresh-b",
    "bucket-full-b",
    "bucket-sparse-b",
    "bsort-fresh-b",
    "bsort-sparse-b",
    "ref-a2",
]
seen = {t for v in time_rows.values() for t in v}
present = [t for t in order if t in seen] + sorted(t for t in seen if t not in order)
present = [t for t in present if t != "ref-a2"] + (
    ["ref-a2"] if "ref-a2" in seen else []
)
print(
    "| Workload: case / nodes / edges / queries / dst | W | "
    + " | ".join(present)
    + " |"
)
print("|---|---:|" + "---:|" * len(present))
for key, v in time_rows.items():
    any_row = next(iter(v.values()))
    label = f"{key[0]} / {key[1]} / {key[2]} / {key[3]} / {key[4]}"
    print(
        f"| {label} | {any_row['W']} | "
        + " | ".join(fmt(v.get(t)) for t in present)
        + " |"
    )

sweep_rows = load("sweep-*.csv")
if sweep_rows:
    tags = sorted({t for v in sweep_rows.values() for t in v})
    print()
    print("| Cost sweep: case / nodes / queries | W | " + " | ".join(tags) + " |")
    print("|---|---:|" + "---:|" * len(tags))
    for key, v in sweep_rows.items():
        any_row = next(iter(v.values()))
        print(
            f"| {key[0]} / {key[1]} / {key[3]} | {any_row['W']} | "
            + " | ".join(fmt(v.get(t)) for t in tags)
            + " |"
        )
