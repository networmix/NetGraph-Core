"""Small separating checks for clocks, equality and the experiment-only helper."""

import copy
import json
from pathlib import Path
import sys
import threading
import time

from analyze_clocks import complement, union
from mach_clock import Registry
from run_workload import equal
from analyze_sample import classify, parse, verify_ranges
from sensitivity import stats


def main(out):
    import _netgraph_core as native

    checks = []
    assert union([(0, 10), (5, 15)]) == 15
    assert complement([(0, 10), (5, 15)], 0, 20) == 5
    assert union([(0, 10), (2, 5)]) == 10
    checks.append("overlapping/nested wall union and independent complement")
    witnesses = verify_ranges()
    assert (
        classify(
            [
                "??? (in _netgraph_core.cpython-314-darwin.so) (load address 0x1000 + 0xc4f28)"
            ]
        )[0]
        == "solver_search"
    )
    assert (
        classify(
            [
                "??? (in _netgraph_core.cpython-314-darwin.so) (load address 0x1000 + 0xc5588)"
            ]
        )[0]
        == "solver_blocking_flow"
    )
    checks.append({"reference_instruction_witnesses": witnesses})
    nodes = parse(
        "Call graph:\n    7 root1\n    + 7 __psynch_cvwait\n    3 root2\n    + 3 native_leaf\nTotal number in stack\n"
    )
    assert [n.count for n in nodes if n.parent is None] == [7, 3]
    assert (
        classify(["frame (in _netgraph_core.test)", "__psynch_cvwait"])[0]
        == "blocked_wait"
    )
    try:
        parse("Call graph:\n    1 root\n    + 2 impossible\nTotal number in stack\n")
    except AssertionError:
        pass
    else:
        raise AssertionError("malformed sample weights accepted")
    triples = [({"wall_s": 2}, {"wall_s": 3}, {"wall_s": 2}) for _ in range(7)]
    statistic = stats(triples, 100000)
    assert statistic["median_ratio"] == 1.5 and statistic["median_slope_s_per_ms"] == 10
    assert statistic["ratio_ci95"] == (1.5, 1.5)
    checks.append(
        "two sample roots conserved; wait separated; malformed weights refused; ratio/slope derived fixture"
    )
    expected = {"metadata": {"core_extension_sha256": "a"}, "rows": [1, 2, False]}
    equal(expected, copy.deepcopy(expected))
    for mutation in (
        lambda d: d["rows"].__setitem__(0, 3),
        lambda d: d["rows"].__setitem__(2, 0),
        lambda d: d.__delitem__("rows"),
        lambda d: d.__setitem__("extra", 1),
        lambda d: d["rows"].reverse(),
    ):
        changed = copy.deepcopy(expected)
        mutation(changed)
        try:
            equal(expected, changed)
        except AssertionError:
            pass
        else:
            raise AssertionError("poisoned export accepted")
    checks.append("scalar/type/missing/extra/order poison rejected")
    registry = Registry()
    registry.install()
    a = registry.snapshot()
    consumed = native._clocks_native_cpu(40_000_000)
    b = registry.snapshot()
    residual = (b[2] + b[3] - a[2] - a[3]) / 2 - (b[4] - a[4])
    assert consumed >= 40_000_000 and abs(residual - consumed) < 2_000_000, (
        consumed,
        residual,
    )
    checks.append(
        {"joined_native_worker_cpu_ns": consumed, "process_minus_python_ns": residual}
    )
    entered, stop = threading.Event(), threading.Event()
    progress = []

    def contender():
        entered.set()
        while not stop.is_set():
            progress.append(1)
            time.sleep(0)

    thread = threading.Thread(target=contender)
    thread.start()
    entered.wait()
    count = len(progress)
    elapsed = native._clocks_delay(20_000_000)
    assert len(progress) > count and elapsed >= 20_000_000
    stop.set()
    thread.join()
    registry.uninstall()
    checks.append(
        {"GIL_free_sleep_ns": elapsed, "other_Python_progress": len(progress) - count}
    )
    report = {
        "checks": checks,
        "summary": "PASS: overlap/closure controls, five export poison classes, joined native CPU residual and GIL-free sleep",
    }
    with Path(out).open("x") as stream:
        json.dump(report, stream, indent=2)
        stream.write("\n")
    print(json.dumps(report), flush=True)


if __name__ == "__main__":
    main(sys.argv[1])
