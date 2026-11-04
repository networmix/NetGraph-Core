"""Report a successful standalone kernel replay and its checked theorem count."""

from __future__ import annotations

import json
import sys
from pathlib import Path

work = Path(sys.argv[1])
check = json.loads((work / "check.json").read_text())
assert check["status"] == "PASS" and check["theorems"]
assert check["dependency_origin"] in {"cache+local", "source", "mixed"}, (
    "State PROOFS_DEPENDENCY_ORIGIN=cache+local|source|mixed truthfully before release replay"
)
report = {
    "status": "PASS",
    "checker": f"Lean {check['lean_version']} leanchecker --fresh Checks (same Lean kernel)",
    "theorem_count": len(check["theorems"]),
    "dependency_builds": check["dependency_origin"],
}
(work / "replay.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
print(f"Standalone kernel replay: {work / 'replay.json'}")
