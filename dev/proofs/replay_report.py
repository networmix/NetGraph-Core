"""Bind a successful standalone kernel replay to its current model-check evidence."""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

work = Path(sys.argv[1])
check = json.loads((work / "check.json").read_text())
assert check["status"] == "PASS" and check["accepted"]
assert check["dependency_origin"] in {"cache+local", "source", "mixed"}, (
    "State PROOFS_DEPENDENCY_ORIGIN=cache+local|source|mixed truthfully before release replay"
)
report = {
    "status": "PASS",
    "checker": f"Lean {check['lean_version']} leanchecker --fresh Checks (same Lean kernel)",
    "check_sha256": hashlib.sha256((work / "check.json").read_bytes()).hexdigest(),
    "log_sha256": hashlib.sha256((work / "replay.log").read_bytes()).hexdigest(),
    "inventory_sha256": check["inventory_sha256"],
    "source_inputs": check["source_inputs"],
    "dependency_builds": check["dependency_origin"],
}
(work / "replay.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
print(f"Standalone kernel replay: {work / 'replay.json'}")
