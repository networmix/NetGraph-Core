"""Make prototype 1's library importable (units, apportion, graph, flowcore, search,
kernels, ledger, pinned, certify).

The directory is appended, not prepended, so this prototype's own ``fixtures``,
``run_fixtures`` and ``tests`` modules are never shadowed by prototype 1's.
"""

import sys
from pathlib import Path

P1_DIR = Path(__file__).resolve().parent.parent / "p1_amount"
if str(P1_DIR) not in sys.path:
    sys.path.append(str(P1_DIR))
