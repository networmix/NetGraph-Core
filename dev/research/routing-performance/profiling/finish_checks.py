"""Final read-only product check and research-script formatting validation."""

import ast
import hashlib
import json
import subprocess
import sys
from pathlib import Path

BASE = Path(__file__).resolve().parents[1]
EFFORT = BASE.parents[1]
repo = EFFORT / "worktrees/netgraph-orient"
changed = []
for line in subprocess.check_output(
    ["git", "-C", str(repo), "ls-files", "ngraph"], text=True
).splitlines():
    original = repo / line
    candidate = BASE / "consumer-lifetime" / line
    if original.read_bytes() != candidate.read_bytes():
        changed.append(line)
assert sorted(changed) == ["ngraph/analysis/msd.py", "ngraph/analysis/placement.py"], (
    changed
)
for project in ("netgraph-orient", "core-orient"):
    status = subprocess.check_output(
        [
            "git",
            "-C",
            str(EFFORT / "worktrees" / project),
            "status",
            "--porcelain",
            "--untracked-files=no",
        ],
        text=True,
    )
    assert not status, (project, status)
for script in (BASE / "scripts").glob("*.py"):
    ast.parse(script.read_text(), filename=str(script))
ruff = repo / "venv/bin/ruff"
results = []
for args in [("check",), ("format", "--check")]:
    command = [str(ruff), *args, "--no-cache", str(BASE / "scripts")]
    result = subprocess.run(command, capture_output=True, text=True, timeout=30)
    assert result.returncode == 0, result.stdout + result.stderr
    results.append({"command": command, "output": result.stdout + result.stderr})
report = {
    "scratch_consumer_changed_files": changed,
    "ruff": results,
    "script_sha256": {
        str(p.relative_to(BASE)): hashlib.sha256(p.read_bytes()).hexdigest()
        for p in sorted((BASE / "scripts").glob("*.py"))
    },
    "summary": "product tracked trees clean; scratch consumer differs in exactly 2 lifetime-decorated files; scripts parse/lint/format checks pass",
}
with Path(sys.argv[1]).open("x") as stream:
    json.dump(report, stream, indent=2)
print(json.dumps({"summary": report["summary"]}), flush=True)
