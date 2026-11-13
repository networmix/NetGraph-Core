"""Independent fresh processes, bounded and serialized in one heavy hold."""

import os
import subprocess
import sys
from pathlib import Path

base = Path(__file__).resolve().parents[1]
for script, mode, suffix, extra in [
    ("run_wall.py", "default", "wall-default-01", []),
    ("native_timers.py", "1", "native-1-01", []),
    ("run_profile.py", "1", "profile-1-01", []),
    ("run_profile.py", "1", "trace-1-01", ["--trace"]),
]:
    out = base / f"backbone-{suffix}"
    command = [
        sys.executable,
        "-B",
        str(base / "scripts/bounded.py"),
        "1200",
        str(out) + ".log",
        sys.executable,
        "-B",
        str(base / "scripts" / script),
        "backbone",
        mode,
        str(out),
        *extra,
    ]
    subprocess.run(
        command,
        check=True,
        timeout=1250,
        env={**os.environ, "PYTHONDONTWRITEBYTECODE": "1"},
    )
print("SUMMARY backbone batch completed")
