#!/usr/bin/env python3
"""Stronger CPU algorithm comparison; no production-dispatch claim."""

import csv
import json
from pathlib import Path

root = Path(__file__).resolve().parent / "results/cpu-algorithms-aba"
data = {p: list(csv.DictReader((root / f"{p}.csv").open())) for p in ("a1", "b", "a2")}
assert [len(v) for v in data.values()] == [14, 216, 14]
assert all(r["correct"] == "pass" for rs in data.values() for r in rs)


def key(row):
    return tuple(row[k] for k in ("case", "n", "e", "queries"))


lines = [
    "# Stronger CPU algorithm A/B/A",
    "",
    "All methods include exact distances and canonical-equivalent positive-cost DAG edge sets. "
    "CPU batches use 10 workers. Eligibility is cached outside timing for both devices. "
    "Seven measured samples and two checked warmups. These remain exploratory, "
    "non-isolated desktop measurements. GPU choice is an oracle across 12 tested variants.",
    "",
    "| Case / N / E / Q | Cached heap + DAG A1 / A2 ms | Dial + DAG A1 / A2 ms | Best GPU B ms | GPU method |",
    "|---|---:|---:|---:|---|",
]
summary = []
for k in dict.fromkeys(key(r) for r in data["a1"]):
    phases = [
        {r["method"]: r for r in data[p] if key(r) == k} for p in ("a1", "b", "a2")
    ]
    a, b, c = phases
    best = min(b.values(), key=lambda r: float(r["wall_ms"]))
    heap = [float(p["cpu_heap_cached_full"]["wall_ms"]) for p in (a, c)]
    dial = [float(p["cpu_dial_full"]["wall_ms"]) for p in (a, c)]
    lines.append(
        f"| {' / '.join(k)} | {heap[0]:.3f} / {heap[1]:.3f} | {dial[0]:.3f} / {dial[1]:.3f} | "
        f"{float(best['wall_ms']):.3f} | {best['method']} |"
    )
    summary.append(
        {"case": k, "heap_a1_a2": heap, "dial_a1_a2": dial, "best_gpu": best}
    )
audit = {
    "phase_rows": {p: len(v) for p, v in data.items()},
    "checked_query_results": sum(
        (int(r["samples"]) + 2) * int(r["queries"]) for rs in data.values() for r in rs
    ),
    "summary": summary,
}
(root / "summary.md").write_text("\n".join(lines) + "\n")
(root / "summary.json").write_text(json.dumps(audit, indent=2) + "\n")
print("\n".join(lines))
