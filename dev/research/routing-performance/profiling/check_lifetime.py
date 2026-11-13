"""Separate a last-owner release from a surviving result at a timed scope."""

import json
import weakref
from pathlib import Path

import netgraph_core as core

import lifetime_tracker
from lifetime_tracker import LifetimeTimers

collector = LifetimeTimers(core)
lifetime_tracker.COLLECTOR = collector
collector.install(core)
collector.active = True
record = collector.begin_thread()
with collector.scope("control"):
    doomed = core.OutputOptions()
    doomed_ref = weakref.ref(doomed)
    del doomed
    assert doomed_ref() is not None, (
        "holding scope must postpone the final reference release"
    )
    survivor = core.OutputOptions()
assert doomed_ref() is None
assert survivor.by_source is False
assert collector.releases[0]["destroyed_wrappers"] == {"OutputOptions": 1}
collector.end_thread(record)
collector.active = False
report = {
    "release": collector.releases,
    "threads": collector.summarize(),
    "summary": "last-owner wrapper freed within timed release, survivor stays valid; destruction counted independently by weakrefs",
}
p = Path(__file__).resolve().parents[1] / "lifetime-controls-01.json"
with p.open("x") as stream:
    json.dump(report, stream, indent=2)
print(json.dumps(report), flush=True)
