"""Three complete fresh P3 processes in one heavy slot; bounded individually."""

import subprocess
import sys
from pathlib import Path

base = Path(__file__).resolve().parents[1]
for workers, label in [("1", "1"), ("8", "8"), ("default", "28")]:
    out = base / f"backbone-p3-lifetime-{label}-01"
    subprocess.run(
        [
            sys.executable,
            "-B",
            str(base / "scripts/bounded.py"),
            "120",
            str(out) + ".log",
            sys.executable,
            "-B",
            str(base / "scripts/run_p3_lifetime.py"),
            "backbone",
            workers,
            str(out),
        ],
        check=True,
        timeout=150,
    )
    plain = base / f"backbone-p3-wall-{label}-01"
    subprocess.run(
        [
            sys.executable,
            "-B",
            str(base / "scripts/bounded.py"),
            "120",
            str(plain) + ".log",
            sys.executable,
            "-B",
            str(base / "scripts/run_wall.py"),
            "backbone",
            workers,
            str(plain),
        ],
        check=True,
        timeout=150,
    )
print(
    "SUMMARY P3 1/8/auto lifetime/pipeline runs and matched plain controls completed",
    flush=True,
)
