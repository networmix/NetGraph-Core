"""PredecessorDAG, Path, RouteTieBreak, the enums and Eligibility.intersect through the bindings
(docs/semantics.md#3-route-values-preddag-and-path): validation in the documented order on
unnarrowed values, read-only views, structural equality."""

from __future__ import annotations

import gc

import numpy as np
import pytest
from helpers import build, exact

import netgraph_core as ngc

L64 = 1 << 64


# ids by (source, cost, destination, input): 0: 0->1, 1: 0->1, 2: 0->2 (cost 3), 3: 1->2, 4: 2->0
G_EDGES = [(0, 1, 5, 1), (0, 1, 5, 1), (1, 2, 5, 1), (0, 2, 5, 3), (2, 0, 5, 1)]


def test_from_csr_round_trip_and_views():
    g = build(3, G_EDGES)
    dag = ngc.PredecessorDAG.from_csr(g, [0, 0, 2, 4], [0, 0, 0, 1], [0, 1, 2, 3])
    assert dag.graph_uid() == g.uid()
    assert list(dag.parent_offsets()) == [0, 0, 2, 4]
    assert list(dag.parents()) == [0, 0, 0, 1]
    assert list(dag.via_edges()) == [0, 1, 2, 3]
    # parent_offsets() is a new read-only array on each call.
    offsets = dag.parent_offsets()
    assert offsets.dtype == np.int32 and not offsets.flags.writeable
    with pytest.raises(ValueError):
        offsets[0] = 9
    assert dag.parent_offsets() is not offsets
    assert not np.shares_memory(dag.parent_offsets(), offsets)
    # parents() and via_edges() are read-only views that keep the DAG alive.
    parents, view = dag.parents(), dag.via_edges()
    for array in (parents, view):
        assert array.dtype == np.int32 and not array.flags.writeable
        with pytest.raises(ValueError):
            array[0] = 9
    del dag
    gc.collect()
    assert list(view) == [0, 1, 2, 3]  # the views keep the DAG alive
    assert list(parents) == [0, 0, 0, 1]
    same = ngc.PredecessorDAG.from_csr(
        g, np.array([0, 0, 2, 4]), np.array([0, 0, 0, 1]), np.array([0, 1, 2, 3])
    )
    again = ngc.PredecessorDAG.from_csr(g, [0, 0, 2, 4], [0, 0, 0, 1], [0, 1, 2, 3])
    assert same == again and hash(same) == hash(again)
    assert same != ngc.PredecessorDAG.from_csr(g, [0, 0, 0, 0], [], [])
    # Equality includes the graph: the same CSR on a twin graph (the same arrays, another uid) is
    # another value.
    twin = build(3, G_EDGES)
    there = ngc.PredecessorDAG.from_csr(twin, [0, 0, 2, 4], [0, 0, 0, 1], [0, 1, 2, 3])
    assert there.graph_uid() == twin.uid() != g.uid()
    assert list(there.via_edges()) == list(same.via_edges())
    assert there != same and not there == same


@pytest.mark.parametrize(
    ("offsets", "parents", "via", "message"),
    [
        ([0, 0, 0], [], [], "parent_offsets has length 3, expected N + 1 = 4"),
        ([1, 0, 2, 1], [], [], "parent_offsets[0] = 1, expected 0"),
        (
            [L64, 0, 0, 0],
            [],
            [],
            "parent_offsets[0] = 18446744073709551616, expected 0",
        ),
        ([0, 2, 1, 2], [0, 0], [0, 1], "parent_offsets are not monotone at 2 (2 > 1)"),
        (
            [0, L64, 0, 0],
            [],
            [],
            "parent_offsets are not monotone at 2 (18446744073709551616 > 0)",
        ),
        (
            [0, 0, 2, 2],
            [0],
            [0, 1],
            "parent_offsets[N] = 2 but parents has 1 and via_edges has 2 entries",
        ),
        (
            [0, 0, 2, 3],
            [9, 0],
            [0, 1],
            "parent_offsets[N] = 3 but parents has 2 and via_edges has 2 entries",
        ),
        ([0, 0, 2, 2], [3, 0], [0, 1], "parents[0] = 3 out of node range [0, 3)"),
        (
            [0, 0, 2, 2],
            [-L64, 0],
            [0, 1],
            "parents[0] = -18446744073709551616 out of node range [0, 3)",
        ),
        (
            [0, 0, 2, 2],
            [0, 0],
            [0, L64],
            "via_edges[1] = 18446744073709551616 out of edge range [0, 5)",
        ),
        (
            [0, 0, 2, 2],
            [1, 0],
            [0, 1],
            "entry 0 of node 1: via edge 0 joins (0, 1), not (1, 1)",
        ),
        ([0, 0, 2, 2], [0, 0], [1, 1], "duplicate entry (child 1, via 1)"),
        # A cycle (0 -> 1 -> 2 -> 0) before the order check: node 1's entries descend.
        ([0, 1, 3, 4], [2, 0, 0, 1], [4, 1, 0, 3], "the entries contain a cycle"),
        (
            [0, 0, 2, 2],
            [0, 0],
            [1, 0],
            "the entries of node 1 are not in ascending via EdgeId (1 then 0)",
        ),
    ],
)
def test_from_csr_validation_order(offsets, parents, via, message):
    g = build(3, G_EDGES)
    with pytest.raises(ValueError, match=exact(f"PredecessorDAG.from_csr: {message}")):
        ngc.PredecessorDAG.from_csr(g, offsets, parents, via)


def test_from_csr_reports_the_first_of_several_defects():
    # Every defect at once, then each reported one removed: the shape (length, offsets[0],
    # monotone, offsets[N]), then the ranges of every entry (unnarrowed, before any endpoint is
    # read), the endpoints, a duplicate, a cycle, the order within a node. Node 1's entries hold
    # the cycle 0 -> 1 -> 2 -> 0 with node 0's, a duplicate and a descending pair.
    g = build(3, G_EDGES)
    big = 2**70
    steps = [
        ([1, 0, 0], [big], [big], "parent_offsets has length 3, expected N + 1 = 4"),
        ([1, 0, 2, 5], [big], [big], "parent_offsets[0] = 1, expected 0"),
        ([0, 2, 1, 5], [big], [big], "parent_offsets are not monotone at 2 (2 > 1)"),
        (
            [0, 0, 2, 5],
            [big],
            [big],
            "parent_offsets[N] = 5 but parents has 1 and via_edges has 1 entries",
        ),
        (
            [0, 1, 4, 6],
            [2, 0, 0, 0, 0, 1],
            [4, 1, 0, 0, 0, big],
            f"via_edges[5] = {big} out of edge range [0, 5)",
        ),
        (
            [0, 1, 4, 6],
            [2, 0, 0, 0, 0, 1],
            [4, 1, 0, 0, 0, 3],
            "entry 4 of node 2: via edge 0 joins (0, 1), not (0, 2)",
        ),
        (
            [0, 1, 4, 6],
            [2, 0, 0, 0, 0, 1],
            [4, 1, 0, 0, 2, 3],
            "duplicate entry (child 1, via 0)",
        ),
        ([0, 1, 3, 5], [2, 0, 0, 0, 1], [4, 1, 0, 2, 3], "the entries contain a cycle"),
        (
            [0, 0, 2, 4],
            [0, 0, 0, 1],
            [1, 0, 2, 3],
            "the entries of node 1 are not in ascending via EdgeId (1 then 0)",
        ),
    ]
    for offsets, parents, via, message in steps:
        with pytest.raises(ValueError) as info:
            ngc.PredecessorDAG.from_csr(g, offsets, parents, via)
        assert str(info.value) == "PredecessorDAG.from_csr: " + message
    ngc.PredecessorDAG.from_csr(g, [0, 0, 2, 4], [0, 0, 0, 1], [0, 1, 2, 3])


@pytest.mark.parametrize(
    "offsets",
    [
        [0, 2**100, 2**99],
        (0, 2**100, 2**99),
        [0, 2**64 - 1, 2**64 - 2],
        np.array([0, 2**64 - 1, 2**64 - 2], dtype=np.uint64),
        [0, 100, 99],
    ],
)
def test_from_csr_compares_offsets_beyond_int64_exactly(offsets):
    # Two offsets above int64 compare by value, so a decreasing pair is the
    # monotonicity error at its index, before the offsets[N] count check.
    g = build(2, [(0, 1, 1, 1)])
    with pytest.raises(ValueError, match="not monotone at 2"):
        ngc.PredecessorDAG.from_csr(g, offsets, [], [])
    increasing = [0, 2**99, 2**100]
    with pytest.raises(
        ValueError, match=r"offsets\[N\] = 1267650600228229401496703205376"
    ):
        ngc.PredecessorDAG.from_csr(g, increasing, [], [])


class FailingIndex:
    """An integer-like object whose __index__ returns 2**100, except on call `fail_on`, which
    raises `error` (kept to check that the caller receives this very exception)."""

    def __init__(self, fail_on=0, error=None):
        self.calls = 0
        self.fail_on = fail_on
        self.error = error

    def __index__(self):
        self.calls += 1
        if self.calls == self.fail_on:
            raise self.error
        return 2**100


class ConversionFailed(Exception):
    pass


@pytest.mark.parametrize("error_type", [ConversionFailed, MemoryError])
@pytest.mark.parametrize("which", [1, 2])
@pytest.mark.parametrize("fail_on", [1, 2])
def test_from_csr_passes_a_failed_offset_conversion_through(fail_on, which, error_type):
    # Offsets 1 and 2 both overflow int64: they are read once, then converted again to compare
    # them exactly. A conversion that fails, on either read of either offset, raises its own
    # exception to the caller, unchanged and with no other exception attached.
    g = build(2, [(0, 1, 1, 1)])
    offsets = [0, FailingIndex(), FailingIndex()]
    error = error_type("conversion failed")
    offsets[which].fail_on, offsets[which].error = fail_on, error
    with pytest.raises(error_type) as info:
        ngc.PredecessorDAG.from_csr(g, offsets, [], [])
    assert info.value is error
    assert info.value.__cause__ is None and info.value.__context__ is None
    assert offsets[which].calls == fail_on
    # Without a failure the same offsets reach the count check (2**100 is not 0 entries).
    with pytest.raises(
        ValueError, match=r"^PredecessorDAG.from_csr: parent_offsets\[N\] = "
    ):
        ngc.PredecessorDAG.from_csr(g, [0, FailingIndex(), FailingIndex()], [], [])


def test_from_csr_type_and_shape_defects():
    g = build(3, G_EDGES)
    with pytest.raises(
        TypeError, match=exact("parent_offsets must be a 1-D array (got 2 dimensions)")
    ):
        ngc.PredecessorDAG.from_csr(g, np.zeros((2, 2), dtype=np.int64), [], [])
    with pytest.raises(
        TypeError, match=exact("parent_offsets[3] must be an integer (got float)")
    ):
        ngc.PredecessorDAG.from_csr(g, [0, 0, 0, 0.0], [], [])
    with pytest.raises(
        TypeError, match=exact("parent_offsets[3] must be an integer (got bool)")
    ):
        ngc.PredecessorDAG.from_csr(g, [0, 0, 0, True], [], [])


def test_path_cost_hops_and_to_dag():
    g = build(3, G_EDGES)
    p = ngc.Path(g, [[1, 0], [3]])
    assert (p.source(), p.destination(), p.cost()) == (0, 2, 2)
    assert p.hops() == [[1, 0], [3]]  # the caller's member order is kept
    dag = p.to_dag()
    assert list(dag.parent_offsets()) == [0, 0, 2, 3]
    assert list(dag.via_edges()) == [0, 1, 3]
    assert p.graph_uid() == g.uid()
    assert ngc.Path(g, [[1, 0], [3]]) == p
    assert ngc.Path(g, [[0, 1], [3]]) != p
    assert hash(ngc.Path(g, [[1, 0], [3]])) == hash(p)
    # Equality includes the graph: the same hops on a twin graph (the same arrays, another uid).
    twin = build(3, G_EDGES)
    there = ngc.Path(twin, [[1, 0], [3]])
    assert there.graph_uid() == twin.uid() != g.uid()
    assert there.hops() == p.hops() and there.cost() == p.cost()
    assert there != p and not there == p


@pytest.mark.parametrize(
    ("hops", "message"),
    [
        ([[], [7]], "hop 1 member 7 out of edge range [0, 5)"),
        ([[L64]], "hop 0 member 18446744073709551616 out of edge range [0, 5)"),
        ([[0, 0], []], "duplicate member 0"),
        ([], "the hop list is empty"),
        ([[0], []], "hop 1 is empty"),
        (
            [[0, 2]],
            "the members of hop 0 do not share (u, v): edge 2 joins (0, 2), edge 0 "
            "joins (0, 1)",
        ),
        ([[0], [2]], "hops 0 and 1 are not contiguous (1 then 0)"),
        ([[2], [4]], "node 0 is repeated"),
    ],
)
def test_path_validation_order(hops, message):
    g = build(3, G_EDGES)
    with pytest.raises(ValueError, match=exact(f"Path: {message}")):
        ngc.Path(g, hops)


def test_path_reports_the_first_of_several_defects():
    # Several defects at once, then each reported one removed; the later defect sits in an
    # earlier hop where it can. Hop 2 of the fourth step mixes 0 -> 2 and 0 -> 1, hops 0 and 1
    # do not meet and node 0 repeats.
    g = build(3, G_EDGES)
    big = 2**70
    steps = [
        ([[0, 0, 2], [], [big]], f"hop 2 member {big} out of edge range [0, 5)"),
        ([[0, 0, 2], [], [4]], "duplicate member 0"),
        ([[0, 2], [], [4]], "hop 1 is empty"),
        (
            [[1], [4], [2, 0]],
            "the members of hop 2 do not share (u, v): edge 0 joins (0, 1), edge 2 "
            "joins (0, 2)",
        ),
        ([[1], [4], [2]], "hops 0 and 1 are not contiguous (1 then 2)"),
        ([[2], [4], [0]], "node 0 is repeated"),
    ]
    for hops, message in steps:
        with pytest.raises(ValueError) as info:
            ngc.Path(g, hops)
        assert str(info.value) == "Path: " + message


def test_route_tie_break_and_enums():
    # RouteTieBreak has exactly four values, FEWER_HOPS first (the default everywhere); the
    # route choice is fields of the models and the search options, not an object.
    assert list(ngc.RouteTieBreak.__members__) == [
        "FEWER_HOPS",
        "HIGHER_CAPACITY",
        "HIGHER_AVAILABLE_CAPACITY",
        "HIGHER_HEADROOM",
    ]
    assert ngc.ShortestPaths().tie_break == ngc.RouteTieBreak.FEWER_HOPS
    assert {m.name for m in ngc.RouteAllocation.__members__.values()} == {
        "OPTIMIZED",
        "WEIGHTED",
        "EQUAL",
    }
    assert set(ngc.Admission.__members__) == {"LOSSLESS", "LOSSY"}
    assert set(ngc.SourceAllocation.__members__) == {"EQUAL", "OPTIMIZED", "WEIGHTED"}
    for retired in ("PathSelection", "EdgeTieBreak", "EdgePriority"):
        assert retired not in dir(ngc), retired


def test_eligibility_intersect_returns_the_operand_object():
    g = build(3, G_EDGES)
    everything = g.eligibility()
    a = g.eligibility(excluded_edges=[0])
    b = g.eligibility(excluded_edges=[1])
    ab = g.eligibility(excluded_edges=[0, 1])
    assert ngc.Eligibility.intersect(everything, a) is a
    assert ngc.Eligibility.intersect(a, everything) is a
    assert ngc.Eligibility.intersect(a, a) is a
    # Two equal values built separately (distinct uids): the first operand wins the tie.
    a2 = g.eligibility(excluded_edges=[0])
    assert a2.uid() != a.uid()
    assert ngc.Eligibility.intersect(a, a2) is a
    assert ngc.Eligibility.intersect(a2, a) is a2
    assert ngc.Eligibility.intersect(a, ab) is ab
    x = ngc.Eligibility.intersect(a, b)
    y = ngc.Eligibility.intersect(a, b)
    assert x.uid() not in (0, a.uid(), b.uid(), y.uid())
    assert list(x.edge_allowed()) == list(ab.edge_allowed())
    other = build(3, G_EDGES)
    with pytest.raises(ValueError, match="different graphs"):
        ngc.Eligibility.intersect(a, other.eligibility())
