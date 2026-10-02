#!/usr/bin/env python3
"""Produce auditable tables from all five sequential benchmark phases."""

import csv
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent / "results" / "aba"
PHASES = ("cpu-a1", "base-a1", "optimized-b", "base-a2", "cpu-a2")
rows = {p: list(csv.DictReader((ROOT / f"{p}.csv").open())) for p in PHASES}
assert [len(rows[p]) for p in PHASES] == [48, 36, 216, 36, 48]
assert all(r["correct"] == "pass" for rs in rows.values() for r in rs)


def key(r):
    return tuple(r[k] for k in ("case", "n", "e", "queries"))


def table(phase, case):
    return {r["method"]: r for r in rows[phase] if key(r) == case}


def elapsed(r):
    return float(r["wall_ms"])


keys = list(dict.fromkeys(key(r) for r in rows["optimized-b"]))
lines = [
    "# Exploratory A/B/A measurements",
    "",
    "Warm full-result latency in ms, including host distances and positive-cost predecessor DAG. "
    "Seven measured samples plus two checked warmups per row. CPU batches use 10 reusable workers. "
    "Best variant is selected after measurement; this is an oracle comparison, not an implemented dispatcher. "
    "Old GPU columns cover the original pull and group variants; the shared-memory variant was not rerun. "
    "The machine was not isolated; these are not accepted performance claims.",
    "",
    "| Case | N / E / Q | CPU full A1 / A2 | Old GPU A1 / A2 | Best tested optimized GPU | ms | CPU BFS full A1 / A2 |",
    "|---|---:|---:|---:|---|---:|---:|",
]
summary = []
for k in keys:
    ca, cb, a, b, c = (
        table(p, k) for p in ("cpu-a1", "cpu-a2", "base-a1", "optimized-b", "base-a2")
    )
    # Choose one old algorithm by its mean of A1/A2, then show both measurements.
    old = min(a, key=lambda name: (elapsed(a[name]) + elapsed(c[name])) / 2)
    best = min(b, key=lambda name: elapsed(b[name]))
    best_label = best + (
        " (64-bit pull fallback)" if b[best]["fallback64"] == "1" else ""
    )
    bfs = (
        f"{elapsed(ca['cpu_bfs_full']):.3f} / {elapsed(cb['cpu_bfs_full']):.3f}"
        if "cpu_bfs_full" in ca
        else "—"
    )
    lines.append(
        f"| {k[0]} | {' / '.join(k[1:])} | {elapsed(ca['cpu_full']):.3f} / {elapsed(cb['cpu_full']):.3f} | "
        f"{elapsed(a[old]):.3f} / {elapsed(c[old]):.3f} | {best_label} | {elapsed(b[best]):.3f} | {bfs} |"
    )
    summary.append(
        {
            "case": k,
            "cpu_a1": elapsed(ca["cpu_full"]),
            "cpu_a2": elapsed(cb["cpu_full"]),
            "old_method": old,
            "old_a1": elapsed(a[old]),
            "old_a2": elapsed(c[old]),
            "best": b[best],
        }
    )

lines += [
    "",
    "## All optimized variants",
    "",
    "No variants omitted. Columns: full wall / distance wall / device work, ms. "
    "Device work includes DAG kernels where selected. Edge visits count SSSP adjacency entries, "
    "excluding bucket vertex scans, queue management and DAG extraction.",
    "",
]
for k in keys:
    lines += [
        f"### {k}",
        "",
        "| Variant | Full / distance / device ms | Edge visits | Rounds | Host waits | 64-bit fallback |",
        "|---|---:|---:|---:|---:|---|",
    ]
    for r in table("optimized-b", k).values():
        times = " / ".join(
            f"{float(r[f]):.3f}" for f in ("wall_ms", "distance_ms", "gpu_ms")
        )
        lines.append(
            f"| {r['method']} | {times} | {r['edge_visits']} | {r['rounds']} | {r['waits']} | {r['fallback64']} |"
        )

checks = sum(
    (int(r["samples"]) + 2) * int(r["queries"]) for rs in rows.values() for r in rs
)
audit = {
    "phase_rows": {p: len(rs) for p, rs in rows.items()},
    "checked_query_results": checks,
    "summary": summary,
}
(ROOT / "summary.md").write_text("\n".join(lines) + "\n")
(ROOT / "summary.json").write_text(json.dumps(audit, indent=2) + "\n")
print(json.dumps({k: v for k, v in audit.items() if k != "summary"}, indent=2))
