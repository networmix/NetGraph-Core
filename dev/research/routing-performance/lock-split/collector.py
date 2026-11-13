"""Frozen T4 outer clocks with two native TLS snapshots per event."""

from array import array
import time

from reuse import load

base = load("collector")
install = base.install


class Collector(base.Collector):
    def __init__(self, core, api):
        super().__init__(core, api)
        import _netgraph_core as native

        self.native = native
        self.detached = array("q")

    def measure(self, name, original, args=(), kwargs=None):
        start = self.registry.snapshot()
        before = self.native._split_snapshot()
        assert before[6] == 0, "native clock invalid"
        w0, c0 = time.perf_counter_ns(), time.thread_time_ns()
        try:
            return original(*args, **(kwargs or {}))
        finally:
            c1, w1 = time.thread_time_ns(), time.perf_counter_ns()
            after = self.native._split_snapshot()
            assert after[6] == 0, "native clock invalid"
            stop = self.registry.snapshot()
            with self.registry.lock:
                self.record(name, w0, w1, c0, c1, start, stop)
                self.detached.extend((*before, *after))

    def save(self, out):
        summary = super().save(out)
        with (out / "detached.i64").open("xb") as stream:
            self.detached.tofile(stream)
        summary["detached_width"] = 14
        summary["native_clock_calibration"] = [
            self.native._split_clock_reads(30000) for _ in range(5)
        ]
        summary["native_clock_calibration_calls"] = 30000
        assert self.native._split_snapshot()[6] == 0, "native calibration clock invalid"
        return summary
