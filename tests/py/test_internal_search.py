"""The private entries of sub-PR 1b: the internal search (``_spf``, APPROACH 6) and the
apportionment (``_round_robin``, ``_highest_averages``, APPROACH 2.1)."""

from __future__ import annotations

import _netgraph_core as ext
import pytest

import netgraph_core as ngc

DET = ngc.EdgeTieBreak.DETERMINISTIC
PREFER = ngc.EdgeTieBreak.PREFER_HIGHER_PRIORITY
CAP = ngc.EdgePriority.CAPACITY


def build(n, edges, q=1.0):
    return ngc.StrictMultiDiGraph(
        n,
        [e[0] for e in edges],
        [e[1] for e in edges],
        [e[2] for e in edges],
        [e[3] for e in edges],
        q,
    )


def ids(g, positions):
    return sorted(int(e) for e, p in enumerate(g.input_index()) if p in positions)


def test_full_search_and_anycast():
    g = build(5, [(0, 1, 1, 1), (1, 2, 1, 1), (0, 3, 1, 2), (3, 2, 1, 5), (1, 4, 1, 5)])
    full = ext._spf(g, 0, "forward", None, ngc.PathSelection(DET, CAP))
    assert full["costs"] == [0, 1, 2, 2, 6]
    assert "target_cost" not in full
    any_ = ext._spf(g, [0], "forward", [2, 4], ngc.PathSelection(DET, CAP))
    assert any_["target_cost"] == 2
    assert sorted(int(e) for e in any_["dag"].via_edges()) == ids(g, {0, 1})
    rev = ext._spf(g, [2, 4], "reverse", [0], ngc.PathSelection(DET, CAP))
    assert rev["dag"] == any_["dag"]
    none = ext._spf(g, 2, "forward", 0, ngc.PathSelection(DET, CAP))
    assert none["target_cost"] is None
    assert list(none["dag"].via_edges()) == []
    unreachable = ext._spf(g, 2, "forward", None, ngc.PathSelection(DET, CAP))
    assert unreachable["costs"] == [None, None, 0, None, None]


def test_staged_single_path_review13_d3():
    # s -> v (1, prio 1), v -> t (3, prio 1), v -> a -> b -> t (1 + 1 + 1, prio 2): the
    # bottleneck of every walk is 1, so the two-hop walk wins under both tie modes.
    g = build(5, [(0, 1, 1, 1), (1, 4, 1, 3), (1, 2, 2, 1), (2, 3, 2, 1), (3, 4, 2, 1)])
    for tie in (DET, PREFER):
        r = ext._spf(g, 0, "forward", 4, ngc.PathSelection(tie, CAP, multipath=False))
        assert r["path"] == ids(g, {0, 1})
        assert r["target_cost"] == 4


def test_search_validation():
    g = build(3, [(0, 1, 2, 1), (1, 2, 2, 1)])
    sel = ngc.PathSelection(DET, CAP)
    with pytest.raises(ValueError, match="requires targets"):
        ext._spf(g, 0, "forward", None, ngc.PathSelection(DET, CAP, multipath=False))
    with pytest.raises(ValueError, match="both endpoint sets"):
        ext._spf(g, [0, 1], "forward", [1], sel)
    with pytest.raises(ValueError, match="empty"):
        ext._spf(g, [], "forward", None, sel)
    with pytest.raises(ValueError, match="residual has length"):
        ext._spf(g, 0, "forward", None, sel, [1])
    with pytest.raises(ValueError, match=r"residual\[1\] = 3"):
        ext._spf(g, 0, "forward", None, sel, [1, 3])
    with pytest.raises(ValueError, match="min_residual"):
        ext._spf(g, 0, "forward", None, sel, None, 0)
    with pytest.raises(ValueError, match="direction"):
        ext._spf(g, 0, "sideways", None, sel)
    gated = ext._spf(g, 0, "forward", None, sel, [2, 1], 2)
    assert gated["costs"] == [0, 1, None]
    # Section 6's order (review-1b R15): the endpoint checks precede the residual's.
    with pytest.raises(ValueError, match="both endpoint sets"):
        ext._spf(g, 0, "forward", 0, sel, [])
    with pytest.raises(ValueError, match="requires targets"):
        ext._spf(
            g, 0, "forward", None, ngc.PathSelection(DET, CAP, multipath=False), []
        )


def test_apportionment_entries():
    assert ext._round_robin(7, 3) == [3, 2, 2]
    assert ext._highest_averages(10, [3, 9, 9]) == [1, 5, 4]
    assert ext._highest_averages(11, [3, 9, 9]) == [1, 5, 5]
    assert ext._highest_averages(4, [1, 3, 3]) == [0, 2, 2]
    big = (1 << 62) - 2
    a = ext._highest_averages(big, [(1 << 60) + 1, (1 << 60) - 1, 7])
    assert sum(a) == big
    with pytest.raises(ValueError, match="outside"):
        ext._highest_averages(1, [0, 1])
    with pytest.raises(ValueError, match="2\\^62"):
        ext._highest_averages(1, [1 << 61, 1 << 61])
