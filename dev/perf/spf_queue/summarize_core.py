#!/usr/bin/env python3
"""Tabulate the Core-level old/new/old A/B/A (results-core) as Markdown."""

import csv
import os
import sys

root = (
    sys.argv[1]
    if len(sys.argv) > 1
    else os.path.join(os.path.dirname(__file__), "results", "core")
)
phases = ["a1-old", "b-new", "b-new-heap", "a2-old"]
for kind in ("time", "sweep"):
    rows = {}
    for p in phases:
        with open(f"{root}/{kind}-{p}.csv") as f:
            for r in csv.DictReader(f):
                key = (r["case"], r["n"], r["e"], r["queries"], r["dst"], r["W"])
                rows.setdefault(key, {})[p] = float(r["median_ms"])
    print(
        f"| {kind}: case / nodes / edges / queries / destination | W | "
        + " | ".join(phases)
        + " | old/new | old/new-heap |"
    )
    print("|---|---:|" + "---:|" * (len(phases) + 2))
    for k, v in rows.items():
        lo, hi = min(v["a1-old"], v["a2-old"]), max(v["a1-old"], v["a2-old"])
        print(
            f"| {' / '.join(k[:5])} | {k[5]} | "
            + " | ".join(f"{v[p]:.3f}" for p in phases)
            + f" | {lo / v['b-new']:.2f}–{hi / v['b-new']:.2f} | {lo / v['b-new-heap']:.2f}–{hi / v['b-new-heap']:.2f} |"
        )
    print()
