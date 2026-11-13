"""Null-wrapper CPU/clock control with eight registered Python callers."""

import argparse
import json
from pathlib import Path
import threading
import time
from types import SimpleNamespace

from collector import Collector


def run(args):
    collector = Collector(SimpleNamespace(), {})
    collector.registry.install()
    collector.active = True
    barrier = threading.Barrier(args.workers)
    results = []

    def null():
        return None

    wrapped = collector.wrap("calibration.null", null)

    def worker():
        barrier.wait()
        start = time.thread_time_ns()
        for _ in range(args.calls):
            null()
        direct = time.thread_time_ns() - start
        start = time.thread_time_ns()
        for _ in range(args.calls):
            wrapped()
        full = time.thread_time_ns() - start
        results.append(
            {
                "thread": collector.registry.local.row["index"],
                "full_wrap_cpu_ns": (full - direct) / args.calls,
            }
        )

    threads = [threading.Thread(target=worker) for _ in range(args.workers)]
    for thread in threads:
        thread.start()
    for thread in threads:
        thread.join()
    collector.registry.uninstall()
    assert (
        len(results) == args.workers
        and len(collector.events) == args.workers * args.calls * 17
    )
    for row in results:
        row["null_span_cpu_ns"] = (
            sum(
                collector.events[i + 6] - collector.events[i + 5]
                for i in range(0, len(collector.events), 17)
                if collector.events[i + 2] == row["thread"]
            )
            / args.calls
        )
    report = {
        "workers": args.workers,
        "calls_per_thread": args.calls,
        "threads": results,
        "summary": "PASS: registered null-call streams recorded; per-thread CPU overhead, not summed wall correction",
    }
    with args.out.open("x") as stream:
        json.dump(report, stream, indent=2)
        stream.write("\n")
    print(json.dumps(report), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("out", type=Path)
    parser.add_argument("--calls", type=int, default=3000)
    parser.add_argument("--workers", type=int, choices=(1, 8, 28), required=True)
    run(parser.parse_args())
