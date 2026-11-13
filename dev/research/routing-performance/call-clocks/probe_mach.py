"""Untimed access/precision controls; failure is a blocker, never a zero sample."""

import json
import hashlib
import sys
import threading
import time
from pathlib import Path

from mach_clock import Registry, cpu_ns, own_port


def busy(seconds):
    start = time.thread_time_ns()
    while time.thread_time_ns() - start < seconds * 1e9:
        pass


def main(out):
    registry = Registry()
    registry.install()
    controls = []
    for kind in ("busy", "sleep", "native_busy"):
        ready, done = threading.Event(), threading.Event()
        go, release = threading.Event(), threading.Event()
        result = {}

        def target():
            port = own_port()
            result["port"] = port
            ready.set()
            go.wait()
            start = time.thread_time_ns()
            mach0 = cpu_ns(port)
            if kind == "busy":
                busy(0.04)
            elif kind == "native_busy":
                hashlib.pbkdf2_hmac(
                    "sha256", b"clock-control", b"fixed-salt", 1_000_000
                )
            else:
                time.sleep(0.04)
            result.update(
                self_ns=time.thread_time_ns() - start, mach_ns=cpu_ns(port) - mach0
            )
            done.set()
            release.wait()

        thread = threading.Thread(target=target, name=kind)
        thread.start()
        ready.wait()
        foreign0 = cpu_ns(result["port"])
        go.set()
        reads = []
        while not done.wait(0.001):
            reads.append(cpu_ns(result["port"]))
        foreign1 = cpu_ns(result["port"])
        release.set()
        thread.join()
        assert all(b >= a for a, b in zip(reads, reads[1:])), reads
        error = abs(result["self_ns"] - result["mach_ns"])
        # 1 ms is a predeclared control tolerance, not a per-call error guarantee.
        assert error < 1_000_000, (kind, result)
        foreign_error = abs(foreign1 - foreign0 - result["self_ns"])
        assert foreign_error < 1_000_000, (kind, foreign_error)
        assert (
            result["mach_ns"] > 35_000_000
            if kind != "sleep"
            else result["mach_ns"] < 5_000_000
        )
        controls.append(
            {
                "kind": kind,
                "result": result,
                "foreign_reads": len(reads),
                "error_ns": error,
                "foreign_delta_ns": foreign1 - foreign0,
                "foreign_error_ns": foreign_error,
            }
        )
    snap = registry.snapshot()
    assert all(row["end_mach"] is not None for row in registry.rows[1:])
    registry.uninstall()
    rows = [{k: v for k, v in row.items() if k != "events"} for row in registry.rows]
    report = {
        "controls": controls,
        "snapshot": snap,
        "threads": rows,
        "summary": "PASS: Mach access, busy/sleep self-clock agreement, foreign monotonic reads and exited counters",
    }
    with Path(out).open("x") as stream:
        json.dump(report, stream, indent=2)
        stream.write("\n")
    print(json.dumps(report), flush=True)


if __name__ == "__main__":
    main(sys.argv[1])
