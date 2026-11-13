"""T4 thread registry, retaining each thread's own native counters before exit."""

from reuse import load

base = load("mach_clock")


class Registry(base.Registry):
    def _begin(self, thread):
        import _netgraph_core as native

        row = super()._begin(thread)
        row["split_start"] = native._split_snapshot()
        return row

    def end(self, row):
        import _netgraph_core as native

        row["split_stop"] = native._split_snapshot()
        super().end(row)
