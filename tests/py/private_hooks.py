"""Whether the extension carries the private hooks (CMake option NETGRAPH_CORE_PRIVATE_HOOKS: ON
for the editable development build, OFF for an installed wheel). Tests that need a hook are
skipped without them, never failed.

A few scenario tests run in both builds and add checks of the private state observers
(``_state()``: next FlowId, bound and orphaned owners, the requested volume, a flow's members)
under ``if HOOKS:``; without the hooks those lines alone are left out. They are
test_flow_ledger.py's test_lifetime_bound_rows_and_epochs,
test_orphans_from_destruction_reset_and_copy, test_records_fields_and_counters_beyond_64_bits,
test_empty_split_graph_keeps_the_route_cost and test_rejected_requests_leave_the_policy_and_the_ledger,
and test_flow_policy.py's test_lifetimes_keep_alive_and_orphans."""

import _netgraph_core
import pytest

HOOKS = hasattr(_netgraph_core, "_fail_native_allocation")

needs_hooks = pytest.mark.skipif(
    not HOOKS, reason="needs the private hooks (NETGRAPH_CORE_PRIVATE_HOOKS)"
)
