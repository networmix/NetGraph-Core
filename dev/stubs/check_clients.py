"""Type-check the client programs of tests/typing against the installed netgraph_core.

Run by ``bash dev/stubs/run.sh check`` with the interpreter of an environment holding an
installed netgraph_core (a wheel, or an installation built from the sdist) and pinned Pyright:

* the installation is the one checked: netgraph_core is imported from this interpreter's
  site-packages, its __init__.pyi equals the committed stub byte for byte, py.typed is there and
  _docs.py is not; Pyright's search paths (its verbose log) hold that site-packages and not the
  source tree's python/ directory, and the clients run from a copy outside the source tree;
* tests/typing/valid_client.py: zero Pyright diagnostics of any severity, and it runs;
* tests/typing/invalid_client.py: the diagnostics are exactly the lines marked
  ``# expect: RULE[, RULE...]`` (one diagnostic of each rule listed), at least one per typing
  case of REQUIRED below;
* a control: a copy of the valid client with one wrong line appended gets exactly that one
  diagnostic (the valid client is analysed, its zero is not vacuous);
* every Pyright run ends normally (``pyright`` below); one that does not stops the check.

Usage: python -I dev/stubs/check_clients.py LABEL WORK_DIR
"""

from __future__ import annotations

import collections
import importlib
import json
import os
import re
import shutil
import subprocess
import sys
import sysconfig
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
CLIENTS = ROOT / "tests/typing"
STUB = ROOT / "python/netgraph_core/__init__.pyi"
EXPECT = re.compile(r"#\s*expect:\s*(\w+(?:\s*,\s*\w+)*)\s*$")
# The typing cases the invalid client must keep, one per marked line (the result classes the
# module refuses to construct, wrong dtypes, wrong argument types, missing Optional checks; a
# line containing the text, marked): a marked line deleted with its mark is caught.
REQUIRED = [
    "ng.Eligibility()",
    "ng.PredDAG()",
    "ng.FlowRecord()",
    "ng.PlacementResult()",
    "ng.SourceFlow()",
    "ng.DagCacheStats()",
    "ng.SensitivityEntry()",
    "ng.SensitivityResult()",
    "ng.SpfResult()",
    "ng.ForwardingProgram()",
    "ng.Tunnel()",
    "ng.TunnelFlow()",
    "ng.Pinned([dag])",
    "npt.NDArray[np.float64] = g.capacity_units()",
    "text: str = g.capacity_units().item(0)",
    'ng.Graph("wrong"',
    "ng.Graph(2, [0], [1], [1.5]",
    "ng.FlowLedger(g, residual=np.zeros(1, dtype=np.float64))",
    "ng.spf(g, 0, ng.KspOptions(1, sel))",
    "ng.capacity(g, 0, 1, cfg)",
    "ng.CapacityAware(guard, max_path_cost=guard)",
    'ng.to_units("1.5", 1.0)',
    "ng.FlowPolicyConfig(ng.Global())",
    "ng.Path(g, np.zeros((1, 1), dtype=np.float64))",
    "ng.batch_capacity(g, float_pairs,",
    "cost: int = result.target_cost",
    "first = result.costs[0]",
    "total = flows.sum()",
    "cost_total: int = best.total_cost",
    "weights = best.program.source_weights",
]


def fail(label: str, problems: list[str]) -> int:
    for p in problems:
        print(f"dev/stubs/check_clients.py [{label}]: {p}", file=sys.stderr)
    return 1


def installation(problems: list[str]) -> Path:
    pkg = importlib.import_module("netgraph_core")
    where = Path(os.path.realpath(pkg.__file__ or "")).parent
    roots = {os.path.realpath(sysconfig.get_paths()[k]) for k in ("purelib", "platlib")}
    if not any(str(where).startswith(r + os.sep) for r in roots):
        problems.append(f"netgraph_core imported from {where}, not this site-packages")
    stub = where / "__init__.pyi"
    if not stub.is_file():
        problems.append(f"the installation has no {stub.name}")
    elif stub.read_bytes() != STUB.read_bytes():
        problems.append(f"{stub} differs from the committed {STUB.relative_to(ROOT)}")
    if not (where / "py.typed").is_file():
        problems.append("the installation has no py.typed")
    if (where / "_docs.py").exists():
        problems.append("the installation still carries _docs.py")
    return where


def pyright(work: Path, verbose: bool = False) -> dict:
    """Pyright on the project in ``work``: its verbose log, or its JSON report. Exits with the
    reason unless Pyright ended normally: exit status 0 (no error) or 1 (errors reported), not
    a signal or another status (a fatal error, a bad configuration or command line), and for
    the report, standard output that is one JSON report whose error count agrees with the
    status."""
    cmd = [sys.executable, "-m", "pyright", "--pythonpath", sys.executable]
    cmd += ["--project", str(work), "--verbose" if verbose else "--outputjson"]
    proc = subprocess.run(cmd, capture_output=True, text=True, cwd=work)
    out = proc.stdout + proc.stderr
    (work / ("pyright-verbose.log" if verbose else "pyright.json")).write_text(out)
    if proc.returncode < 0:
        sys.exit(f"pyright was killed by signal {-proc.returncode} on {work}:\n{out}")
    if proc.returncode not in (0, 1):
        sys.exit(f"pyright failed (exit {proc.returncode}) on {work}:\n{out}")
    if verbose:
        return {"log": out}
    try:
        report = json.loads(proc.stdout)
    except json.JSONDecodeError:
        sys.exit(f"pyright gave no JSON (exit {proc.returncode}):\n{out}")
    if not isinstance(report, dict):
        report = {}
    summary = report.get("summary")
    errors = summary.get("errorCount") if isinstance(summary, dict) else None
    if (
        not isinstance(report.get("generalDiagnostics"), list)
        or type(errors) is not int
        or proc.returncode != (errors > 0)
    ):
        sys.exit(
            f"pyright's report is not one of a complete run (exit {proc.returncode}):\n{out}"
        )
    return report


def diagnostics(report: dict) -> dict[str, collections.Counter]:
    got: dict[str, collections.Counter] = collections.defaultdict(collections.Counter)
    for d in report.get("generalDiagnostics", []):
        line = d["range"]["start"]["line"] + 1
        got[Path(d["file"]).name][(line, d.get("rule", d["severity"]))] += 1
    return got


def main() -> int:
    if len(sys.argv) != 3:
        sys.exit("usage: check_clients.py LABEL WORK_DIR")
    label, work = sys.argv[1], Path(sys.argv[2]).resolve()
    problems: list[str] = []
    site = installation(problems)
    shutil.rmtree(work, ignore_errors=True)
    work.mkdir(parents=True)
    names = ["valid_client.py", "invalid_client.py"]
    for name in names:
        shutil.copyfile(CLIENTS / name, work / name)
    version = f"{sys.version_info.major}.{sys.version_info.minor}"
    (work / "pyrightconfig.json").write_text(
        json.dumps(
            {
                "include": names,
                "typeCheckingMode": "standard",
                "pythonVersion": version,
            },
            indent=2,
        )
    )

    log = pyright(work, verbose=True)["log"]
    paths = [ln.strip() for ln in log.splitlines() if ln.startswith("  ")]
    if not any(
        os.path.realpath(p) == str(site.parent) for p in paths if os.path.isdir(p)
    ):
        problems.append(
            f"Pyright's search paths do not hold {site.parent} (pyright-verbose.log)"
        )
    source = os.path.realpath(ROOT / "python")
    if any(os.path.realpath(p) == source for p in paths if os.path.isdir(p)):
        problems.append("Pyright searches the source tree's python/ directory")

    report = pyright(work)
    got = diagnostics(report)
    if got["valid_client.py"]:
        for (line, rule), n in sorted(got["valid_client.py"].items()):
            problems.append(
                f"valid_client.py:{line}: {n} unexpected diagnostic(s) {rule}"
            )

    lines = (work / "invalid_client.py").read_text().splitlines()
    want: collections.Counter = collections.Counter()
    for i, text in enumerate(lines, 1):
        m = EXPECT.search(text)
        if m:
            for rule in re.split(r"\s*,\s*", m.group(1)):
                want[(i, rule)] += 1
    for case in REQUIRED:
        if not any(case in t and EXPECT.search(t) for t in lines):
            problems.append(f"invalid_client.py lost the marked case {case!r}")
    for key in sorted(set(want) | set(got["invalid_client.py"])):
        w, g = want[key], got["invalid_client.py"][key]
        if w != g:
            problems.append(
                f"invalid_client.py:{key[0]}: {key[1]} expected {w}, got {g}"
            )
    details = {
        (d["range"]["start"]["line"] + 1, d.get("rule")): d["message"].splitlines()[0]
        for d in report.get("generalDiagnostics", [])
        if Path(d["file"]).name == "invalid_client.py"
    }

    # The control: one wrong line appended to a copy of the valid client.
    control = work / "control"
    control.mkdir()
    planted = (CLIENTS / "valid_client.py").read_text() + "oops: str = g.num_edges()\n"
    (control / "valid_client.py").write_text(planted)
    (control / "pyrightconfig.json").write_text(
        json.dumps({"include": ["valid_client.py"], "typeCheckingMode": "standard"})
    )
    seen = diagnostics(pyright(control))["valid_client.py"]
    last = len(planted.splitlines())
    if seen != collections.Counter({(last, "reportAssignmentType"): 1}):
        problems.append(f"control: the planted line {last} gave {dict(seen)}")

    run = subprocess.run(
        [sys.executable, "-I", str(work / "valid_client.py")],
        capture_output=True,
        text=True,
    )
    if run.returncode != 0:
        problems.append(f"valid_client.py does not run:\n{run.stderr}")

    if problems:
        return fail(label, problems)
    print(
        f"[{label}] {site}: stub identical to the committed one; valid client: 0 "
        f"diagnostics and it runs (control: a planted line caught); invalid client: "
        f"{sum(want.values())} expected diagnostics on {len({k[0] for k in want})} lines, no other"
    )
    for (line, rule), message in sorted(details.items()):
        print(f"  invalid_client.py:{line}: {rule}: {message}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
