"""Whether the extension carries the private hooks (CMake option NETGRAPH_CORE_PRIVATE_HOOKS: ON
for the editable development build, OFF for an installed wheel). Tests that need a hook are
skipped without them, never failed."""

import _netgraph_core
import pytest

HOOKS = hasattr(_netgraph_core, "_max_flows_per_call")

needs_hooks = pytest.mark.skipif(
    not HOOKS, reason="needs the private hooks (NETGRAPH_CORE_PRIVATE_HOOKS)"
)
