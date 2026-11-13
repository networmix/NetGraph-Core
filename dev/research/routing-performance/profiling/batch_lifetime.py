"""Bounded sequential remainder of the coordinator-directed lifetime correction."""

import subprocess
import sys
from pathlib import Path

base = Path(__file__).resolve().parents[1]
for workload, workers, label in [
    ("dcbb96", "2", "dcbb-lifetime-2-01"),
    ("dcbb96", "8", "dcbb-lifetime-8-01"),
    ("backbone", "1", "backbone-lifetime-1-01"),
]:
    out = base / label
    subprocess.run(
        [
            sys.executable,
            "-B",
            str(base / "scripts/bounded.py"),
            "180",
            str(out) + ".log",
            sys.executable,
            "-B",
            str(base / "scripts/run_lifetime.py"),
            workload,
            workers,
            str(out),
        ],
        check=True,
        timeout=200,
    )
print("SUMMARY lifecycle batch completed")
