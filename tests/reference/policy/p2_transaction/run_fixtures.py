"""Run every prototype-2 fixture and write ``snapshots/<fixture>.json``.

Usage: ``venv/bin/python run_fixtures.py [fixture ...]``; exit status 1 on any
failed expectation.
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import _p1  # noqa: F401
from fixtures import FIXTURE_MAP, FIXTURES
from recorder import PSnap

HERE = Path(__file__).resolve().parent


def run_one(fx) -> dict:
    s = PSnap(fx.name, fx.family, fx.description, fx.rules)
    fx.run(s)
    data = s.finish()
    old = [k for k, v in FIXTURE_MAP.items() if v == fx.name]
    data["regenerates"] = old
    return data


def main(argv: list[str]) -> int:
    only = set(argv[1:])
    out = HERE / "snapshots"
    out.mkdir(exist_ok=True)
    failed, steps, checks = [], 0, 0
    t0 = time.perf_counter()
    names = {fx.name for fx in FIXTURES}
    missing = sorted(v for v in FIXTURE_MAP.values() if v not in names)
    assert not missing, missing
    for fx in FIXTURES:
        if only and fx.name not in only:
            continue
        data = run_one(fx)
        (out / f"{fx.name}.json").write_text(json.dumps(data, indent=1) + "\n")
        steps += data["counts"]["steps"]
        checks += data["counts"]["checks"]
        print(
            f"{'ok  ' if data['passed'] else 'FAIL'} {fx.name} ({data['counts']['checks']} checks)"
        )
        for f in data["failures"]:
            print(f"     {f}")
        if not data["passed"]:
            failed.append(fx.name)
    n = len(only) if only else len(FIXTURES)
    print(
        f"{n - len(failed)}/{n} fixtures pass ({len(FIXTURE_MAP)} regenerated from the "
        f"grid set), {steps} operations, {checks} checks, {time.perf_counter() - t0:.1f} s"
    )
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
