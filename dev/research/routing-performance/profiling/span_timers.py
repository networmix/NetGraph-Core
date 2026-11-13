"""All public Core Python boundaries, with exclusive per-thread CPU and wall spans.

CPU outside spans estimates serial GIL work. It also contains implicit native
object destruction and Python/runtime work not reached through a public call.
Wall spans include GIL release/reacquisition and binding conversion.
"""

import threading
import time
from collections import Counter


class SpanTimers:
    def __init__(self):
        self.active = False
        self.local = threading.local()
        self.threads = []
        self.installed = []

    def begin_thread(self):
        record = {
            "name": threading.current_thread().name,
            "ident": threading.get_ident(),
            "events": [],
            "start_wall_ns": time.perf_counter_ns(),
            "start_cpu_ns": time.thread_time_ns(),
        }
        self.local.record = record
        self.threads.append(record)
        return record

    def end_thread(self, record):
        record["end_cpu_ns"] = time.thread_time_ns()
        record["end_wall_ns"] = time.perf_counter_ns()

    def wrap(self, name, original):
        def wrapped(*args, **kwargs):
            if not self.active:
                return original(*args, **kwargs)
            record = self.local.record
            start_wall = time.perf_counter_ns()
            start_cpu = time.thread_time_ns()
            try:
                return original(*args, **kwargs)
            finally:
                stop_cpu = time.thread_time_ns()
                stop_wall = time.perf_counter_ns()
                record["events"].append(
                    (name, start_wall, stop_wall, start_cpu, stop_cpu)
                )

        return wrapped

    def install(self, core, *, properties_only=False):
        for name, value in list(vars(core).items()):
            if (
                type(value).__name__ == "builtin_function_or_method"
                and not properties_only
            ):
                setattr(core, name, self.wrap(name, value))
                self.installed.append(name)
            elif type(value).__name__ == "pybind11_type" and not hasattr(
                value, "__members__"
            ):
                for attr, member in list(vars(value).items()):
                    if isinstance(member, property):
                        wrapped = property(
                            self.wrap(name + "." + attr, member.fget),
                            member.fset,
                            member.fdel,
                            member.__doc__,
                        )
                    elif properties_only or (
                        attr.startswith("_") and attr != "__init__"
                    ):
                        continue
                    elif isinstance(member, staticmethod):
                        wrapped = staticmethod(
                            self.wrap(name + "." + attr, member.__func__)
                        )
                    elif type(member).__name__ == "instancemethod":
                        wrapped = self.wrap(name + "." + attr, getattr(value, attr))
                    else:
                        continue
                    setattr(value, attr, wrapped)
                    self.installed.append(name + "." + attr)

    def calibrate(self, count=30000):
        def null():
            return None

        self.active = True
        record = self.begin_thread()
        wrapped = self.wrap("calibration.null", null)
        w0, c0 = time.perf_counter_ns(), time.thread_time_ns()
        for _ in range(count):
            null()
        direct_wall, direct_cpu = (
            time.perf_counter_ns() - w0,
            time.thread_time_ns() - c0,
        )
        w0, c0 = time.perf_counter_ns(), time.thread_time_ns()
        for _ in range(count):
            wrapped()
        wrapped_wall, wrapped_cpu = (
            time.perf_counter_ns() - w0,
            time.thread_time_ns() - c0,
        )
        inner_wall = sum(e[2] - e[1] for e in record["events"]) / count
        inner_cpu = sum(e[4] - e[3] for e in record["events"]) / count
        self.end_thread(record)
        self.active = False
        self.threads.clear()
        return {
            "calls": count,
            "cost_per_wrap_wall_ns": (wrapped_wall - direct_wall) / count,
            "cost_per_wrap_cpu_ns": (wrapped_cpu - direct_cpu) / count,
            "recorded_null_wall_ns": inner_wall,
            "recorded_null_cpu_ns": inner_cpu,
            "clock_resolution_s": time.get_clock_info("thread_time").resolution,
        }

    def summarize(self):
        rows = []
        for index, record in enumerate(self.threads):
            counts = Counter(e[0] for e in record["events"])
            rows.append(
                {
                    "thread_index": index,
                    "name": record["name"],
                    "cpu_s": (record["end_cpu_ns"] - record["start_cpu_ns"]) / 1e9,
                    "wall_s": (record["end_wall_ns"] - record["start_wall_ns"]) / 1e9,
                    "native_span_wall_s": sum(e[2] - e[1] for e in record["events"])
                    / 1e9,
                    "native_span_cpu_s": sum(e[4] - e[3] for e in record["events"])
                    / 1e9,
                    "calls": len(record["events"]),
                    "call_counts": dict(counts),
                }
            )
        return rows
