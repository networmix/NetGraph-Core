"""The public surface: the graph on units, Eligibility, the conversions, the route structures,
the selection and placement axes, the ledger, the policy engine, the DagCache and the
evaluators; the extension's private entries, which netgraph_core never re-exports; and the
spellings the interface grooming retired (APPROACH 12), which exist nowhere, without aliases."""

import _netgraph_core
import pytest
from private_hooks import HOOKS

import netgraph_core as ngc

PUBLIC = {
    "Graph",
    "Eligibility",
    "to_units",
    "as_float",
    "PredDAG",
    "Path",
    "EdgeTieBreak",
    "EdgePriority",
    "PathSelection",
    "FlowPlacement",
    "Admission",
    "FlowLedger",
    "FlowRecord",
    "Absolute",
    "Stretch",
    "DiminishingReturns",
    "CostOnly",
    "CapacityAware",
    "Global",
    "Pinned",
    "FlowPolicyConfig",
    "OutputOptions",
    "Termination",
    "PlacementResult",
    "FlowPolicy",
    "MAX_FLOWS_PER_CALL",
    "DagCache",
    "DagCacheStats",
    "CapacityOptions",
    "SensitivityOptions",
    "SensitivityStatus",
    "SensitivityEntry",
    "SensitivityResult",
    "capacity",
    "batch_capacity",
    "sensitivity",
}

# Shipped private entries: the flow core's recorder entry (APPROACH 2.2) and the internal search.
PRIVATE = {"_place_global_units", "_spf"}

# The private hooks (NETGRAPH_CORE_PRIVATE_HOOKS: the development build only): module functions
# and class members.
HOOK_FUNCTIONS = {
    "_max_flows_per_call",
    "_sensitivity_fail_edge",
    "_fail_python_allocation",
    "_round_robin",
    "_highest_averages",
}
HOOK_MEMBERS = {
    "FlowLedger": {"_state"},
    "FlowPolicy": {"_place_traced", "_state", "_owner"},
}


def test_public_surface():
    assert set(ngc.__all__) == PUBLIC | {"__version__"}
    for name in ngc.__all__:
        assert hasattr(ngc, name)
    assert isinstance(ngc.__version__, str)
    assert ngc.MAX_FLOWS_PER_CALL == 10000


def test_extension_surface_is_exactly_the_public_and_private_sets():
    names = {n for n in dir(_netgraph_core) if not n.startswith("__")}
    assert names == PUBLIC | PRIVATE | (HOOK_FUNCTIONS if HOOKS else set())
    for name in PRIVATE | HOOK_FUNCTIONS:
        assert not hasattr(ngc, name)
    # Every private class member is a hook (beside pybind11's own conduit); the hooks are all
    # present or all absent.
    for cls in PUBLIC:
        obj = getattr(_netgraph_core, cls)
        private = {
            n
            for n in dir(obj)
            if n.startswith("_")
            and not n.startswith("__")
            and n != "_pybind11_conduit_v1_"
        }
        want = HOOK_MEMBERS.get(cls, set()) if HOOKS else set()
        assert private == want, cls


# The 1a to 1c spellings renamed before 1d (APPROACH 12) and the old Python keyword: no alias.
RETIRED = {
    "StrictMultiDiGraph",
    "FlowGraph",
    "Outputs",
    "Diminishing",
    "kMaxFlowsPerCall",
}


def test_retired_spellings_are_absent():
    for name in RETIRED:
        assert not hasattr(_netgraph_core, name), name
        assert not hasattr(ngc, name), name
    assert ngc.Graph.__name__ == "Graph"
    assert ngc.FlowLedger.__name__ == "FlowLedger"
    assert ngc.OutputOptions.__name__ == "OutputOptions"
    assert ngc.DiminishingReturns.__name__ == "DiminishingReturns"


def _config():
    return ngc.FlowPolicyConfig(
        routing=ngc.CostOnly(),
        placement=ngc.FlowPlacement.EQUAL_BALANCED,
        admission=ngc.Admission.LOSSLESS,
        selection=ngc.PathSelection(
            ngc.EdgeTieBreak.DETERMINISTIC, ngc.EdgePriority.RESIDUAL
        ),
    )


def test_policy_takes_its_ledger_by_the_ledger_keyword_only():
    g = ngc.Graph(2, [0], [1], [3], [1], 1.0)
    ledger = ngc.FlowLedger(g)
    p = ngc.FlowPolicy(ledger=ledger, sources=0, destinations=1, config=_config())
    assert p.place(2).placed == 2
    assert ledger.records()[0].delivered == 2
    # Positional order unchanged: the ledger first.
    assert ngc.FlowPolicy(ledger, 0, 1, _config()).place(1).placed == 1
    with pytest.raises(TypeError):
        ngc.FlowPolicy(flow_graph=ledger, sources=0, destinations=1, config=_config())
    with pytest.raises(ValueError, match="^ledger must be a FlowLedger$"):
        ngc.FlowPolicy(None, 0, 1, _config())
    assert "ledger: _netgraph_core.FlowLedger" in ngc.FlowPolicy.__init__.__doc__
