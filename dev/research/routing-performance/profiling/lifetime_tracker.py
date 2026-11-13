"""Scratch-consumer lexical scopes retain native wrappers until timed release.

Only Core objects created or returned through Python public boundaries are held.
A state/probe's frame has returned before its holding dictionary is cleared.
Live output wrappers survive through their ordinary references; dead wrappers and
native holders are destroyed inside the release span. No destructor monkeypatch.
"""

import functools
import threading
import time
import weakref
from collections import Counter
from contextlib import contextmanager

from span_timers import SpanTimers

COLLECTOR = None


class LifetimeTimers(SpanTimers):
    def __init__(self, core):
        super().__init__()
        self.core_types = {
            v
            for v in vars(core).values()
            if type(v).__name__ == "pybind11_type" and not hasattr(v, "__members__")
        }
        self.properties = {
            name + "." + field
            for name, kind in vars(core).items()
            if type(kind).__name__ == "pybind11_type" and kind in self.core_types
            for field, member in vars(kind).items()
            if isinstance(member, property)
        }
        self.releases = []

    def begin_thread(self):
        record = super().begin_thread()
        self.local.holds = []
        return record

    def retain(self, value):
        if type(value) in self.core_types:
            hold = self.local.holds[-1] if self.local.holds else None
            if hold is not None:
                hold[id(value)] = value
        elif isinstance(value, (tuple, list)):
            for item in value:
                if type(item) in self.core_types:
                    self.retain(item)

    def wrap(self, name, original):
        timed = super().wrap(name, original)

        def wrapped(*args, **kwargs):
            result = timed(*args, **kwargs)
            if self.active:
                if name.endswith(".__init__"):
                    self.retain(args[0])
                else:
                    self.retain(result)
            return result

        return wrapped

    @contextmanager
    def scope(self, label):
        if not self.active:
            yield
            return
        hold = {}
        self.local.holds.append(hold)
        try:
            yield
        finally:
            self.local.holds.pop()
            refs = [
                (type(value).__name__, weakref.ref(value)) for value in hold.values()
            ]
            record = self.local.record
            w0, c0 = time.perf_counter_ns(), time.thread_time_ns()
            hold.clear()
            c1, w1 = time.thread_time_ns(), time.perf_counter_ns()
            record["events"].append(("CoreObjects.release", w0, w1, c0, c1))
            dead = Counter(kind for kind, reference in refs if reference() is None)
            self.releases.append(
                {
                    "thread": threading.current_thread().name,
                    "scope": label,
                    "held_wrappers": len(refs),
                    "destroyed_wrappers": dict(dead),
                    "wall_ns": w1 - w0,
                    "cpu_ns": c1 - c0,
                }
            )


def lifetime_scope(label):
    def decorate(original):
        @functools.wraps(original)
        def scoped(*args, **kwargs):
            if COLLECTOR is None or not COLLECTOR.active:
                return original(*args, **kwargs)
            with COLLECTOR.scope(label):
                return original(*args, **kwargs)

        return scoped

    return decorate
