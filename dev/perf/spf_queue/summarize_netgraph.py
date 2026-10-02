#!/usr/bin/env python3
"""Tabulate the NetGraph-workload A/B/A (old/new/old, fresh processes) as Markdown."""

import csv
import glob
import os
import sys

root = (
    sys.argv[1]
    if len(sys.argv) > 1
    else os.path.join(os.path.dirname(__file__), "results", "netgraph")
)
phases = ["a1-old", "b-new", "b-new-heap", "a2-old"]
rows: dict[str, dict[str, dict]] = {}
for path in sorted(glob.glob(os.path.join(root, "*.csv"))):
    tag = os.path.basename(path)[:-4]
    with open(path) as f:
        for r in csv.DictReader(f):
            rows.setdefault(r["workload"], {})[tag] = r
present = [p for p in phases if any(p in v for v in rows.values())]
print(
    "| Workload | "
    + " | ".join(f"{p} ms" for p in present)
    + " | old/new | checksums |"
)
print("|---|" + "---:|" * len(present) + "---:|---|")
mismatch = 0
for name, v in rows.items():
    med = {p: float(v[p]["median_s"]) * 1000 for p in present if p in v}
    ratio = ""
    if "a1-old" in med and "a2-old" in med and "b-new" in med:
        ratio = f"{min(med['a1-old'], med['a2-old']) / med['b-new']:.2f}–{max(med['a1-old'], med['a2-old']) / med['b-new']:.2f}"
    sums = {v[p]["checksum"] for p in present if p in v}
    ok = "same" if len(sums) == 1 else "DIFFER"
    if ok != "same":
        mismatch += 1
    print(
        f"| {name} | "
        + " | ".join(f"{med[p]:.2f}" if p in med else "—" for p in present)
        + f" | {ratio} | {ok} |"
    )
print()
print(f"checksum mismatches: {mismatch}")
