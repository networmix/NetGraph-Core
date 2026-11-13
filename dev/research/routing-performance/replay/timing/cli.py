"""Fresh NetGraph CLI entry with only before/after GC metadata."""

import gc
import json
from pathlib import Path
import runpy
import sys

destination, *arguments = sys.argv[1:]
before = {
    "counts": gc.get_stats(),
    "enabled": gc.isenabled(),
    "threshold": gc.get_threshold(),
}
sys.argv = ["ngraph", *arguments]
try:
    runpy.run_module("ngraph", run_name="__main__", alter_sys=True)
finally:
    with Path(destination).open("x") as stream:
        json.dump(
            {
                "before": before,
                "after": {
                    "counts": gc.get_stats(),
                    "enabled": gc.isenabled(),
                    "threshold": gc.get_threshold(),
                },
            },
            stream,
        )
