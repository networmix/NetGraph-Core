"""Bound one evidence command; preserve stdout/stderr and a final result line."""

import json
import os
import subprocess
import sys
import time
from pathlib import Path

limit, log, *command = sys.argv[1:]
start = time.monotonic()
with Path(log).open("x") as stream:
    try:
        result = subprocess.run(
            command,
            stdout=stream,
            stderr=subprocess.STDOUT,
            timeout=float(limit),
            env={**os.environ, "PYTHONDONTWRITEBYTECODE": "1"},
        )
        code, outcome = result.returncode, "completed"
    except subprocess.TimeoutExpired:
        code, outcome = 124, "caught: timeout"
    record = {
        "summary": outcome,
        "exit": code,
        "elapsed_s": time.monotonic() - start,
        "command": command,
    }
    stream.write("\nSUMMARY " + json.dumps(record) + "\n")
print(json.dumps(record), flush=True)
sys.exit(code)
