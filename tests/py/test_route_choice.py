"""The route choice of ShortestPaths: the single-path order read on the service
view at preparation, one record per source over the route the rule chooses."""

from __future__ import annotations

from helpers import build

import netgraph_core as ngc

TIE = ngc.RouteTieBreak


def test_shortest_paths_single_path():
    # TRI (three parallel cost-1 links of capacities 4, 10, 100 and a cost-2 link of 50) with the
    # service view 4, 8, 5, 50: FEWER_HOPS e0, HIGHER_CAPACITY e2, HIGHER_AVAILABLE_CAPACITY e1,
    # HIGHER_HEADROOM e0 (4/4 = 1 against 4/5 and 1/20).
    g = build(2, [(0, 1, 4, 1), (0, 1, 10, 1), (0, 1, 100, 1), (0, 1, 50, 2)])
    for tie, want in (
        (TIE.FEWER_HOPS, 0),
        (TIE.HIGHER_CAPACITY, 2),
        (TIE.HIGHER_AVAILABLE_CAPACITY, 1),
        (TIE.HIGHER_HEADROOM, 0),
    ):
        ledger = ngc.FlowLedger(g, service_baseline=[4, 8, 5, 50])
        cfg = ngc.FlowPolicyConfig(
            routing=ngc.ShortestPaths(tie_break=tie, multipath=False)
        )
        policy = ngc.FlowPolicy(ledger, 0, 1, cfg)
        assert policy.place(1).call_delivered == 1, tie
        (record,) = policy.records()
        assert record.route.via_edges().tolist() == [want], tie
        assert record.edge_flows == [(want, 1)], tie
