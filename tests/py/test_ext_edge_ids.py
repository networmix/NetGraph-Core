"""Tests for external edge ID preservation in StrictMultiDiGraph."""

import numpy as np

import netgraph_core


def test_ext_edge_ids_without_reverse():
    """Test that external edge IDs are preserved through reordering."""
    # Create edges in non-sorted order to trigger reordering
    num_nodes = 4
    src = np.array([2, 0, 1, 0], dtype=np.int32)
    dst = np.array([3, 1, 2, 2], dtype=np.int32)
    capacity = np.array([10.0, 20.0, 15.0, 25.0], dtype=np.float64)
    cost = np.array([5, 1, 3, 1], dtype=np.int64)
    ext_edge_ids = np.array([1000, 2000, 3000, 4000], dtype=np.int64)

    g = netgraph_core.StrictMultiDiGraph.from_arrays(
        num_nodes=num_nodes,
        src=src,
        dst=dst,
        capacity=capacity,
        cost=cost,
        ext_edge_ids=ext_edge_ids,
    )

    assert g.num_nodes() == 4
    assert g.num_edges() == 4

    # Get reordered data
    src_view = g.edge_src_view()
    dst_view = g.edge_dst_view()
    cost_view = g.cost_view()
    ext_ids_view = g.ext_edge_ids_view()

    # Core sorts by (cost, src, dst), so expected order:
    # (0→1, cost=1, ext_id=2000)
    # (0→2, cost=1, ext_id=4000)
    # (1→2, cost=3, ext_id=3000)
    # (2→3, cost=5, ext_id=1000)
    expected_order = [
        (0, 1, 1, 2000),
        (0, 2, 1, 4000),
        (1, 2, 3, 3000),
        (2, 3, 5, 1000),
    ]

    for i, (exp_src, exp_dst, exp_cost, exp_ext_id) in enumerate(expected_order):
        assert src_view[i] == exp_src
        assert dst_view[i] == exp_dst
        assert cost_view[i] == exp_cost
        assert ext_ids_view[i] == exp_ext_id


def test_ext_edge_ids_with_manual_reverse():
    """Test ext_edge_ids when reverse edges are manually added with distinct IDs."""
    num_nodes = 3

    # Manually add forward and reverse edges with distinct ext_edge_ids
    src = np.array([0, 1, 1, 2], dtype=np.int32)  # fwd0, rev0, fwd1, rev1
    dst = np.array([1, 0, 2, 1], dtype=np.int32)
    capacity = np.array([10.0, 10.0, 20.0, 20.0], dtype=np.float64)
    cost = np.array([1, 1, 2, 2], dtype=np.int64)
    ext_edge_ids = np.array([0, 1, 2, 3], dtype=np.int64)  # Distinct for fwd/rev

    g = netgraph_core.StrictMultiDiGraph.from_arrays(
        num_nodes=num_nodes,
        src=src,
        dst=dst,
        capacity=capacity,
        cost=cost,
        ext_edge_ids=ext_edge_ids,
    )

    assert g.num_edges() == 4

    ext_ids_view = g.ext_edge_ids_view()
    src_view = g.edge_src_view()
    dst_view = g.edge_dst_view()

    # Verify all ext_edge_ids are distinct
    assert len(set(ext_ids_view)) == 4, "ext_edge_ids should all be distinct"

    # Verify we can distinguish forward from reverse
    # In sorted order by cost, src, dst:
    # (0→1, cost=1, ext_id=0), (1→0, cost=1, ext_id=1), (1→2, cost=2, ext_id=2), (2→1, cost=2, ext_id=3)
    expected = [(0, 1, 0), (1, 0, 1), (1, 2, 2), (2, 1, 3)]
    for i, (exp_src, exp_dst, exp_ext_id) in enumerate(expected):
        assert src_view[i] == exp_src
        assert dst_view[i] == exp_dst
        assert ext_ids_view[i] == exp_ext_id


def test_ext_edge_ids_encoding_scheme():
    """Test a typical encoding scheme: (linkIndex << 1) | directionBit."""
    num_nodes = 3

    # Simulate two links with forward/reverse encoding
    # Link 0: fwd=(0<<1)|0=0, rev=(0<<1)|1=1
    # Link 1: fwd=(1<<1)|0=2, rev=(1<<1)|1=3
    src = np.array([0, 1, 1, 2], dtype=np.int32)
    dst = np.array([1, 0, 2, 1], dtype=np.int32)
    capacity = np.array([10.0, 10.0, 20.0, 20.0], dtype=np.float64)
    cost = np.array([1, 1, 1, 1], dtype=np.int64)
    ext_edge_ids = np.array([0, 1, 2, 3], dtype=np.int64)

    g = netgraph_core.StrictMultiDiGraph.from_arrays(
        num_nodes=num_nodes,
        src=src,
        dst=dst,
        capacity=capacity,
        cost=cost,
        ext_edge_ids=ext_edge_ids,
    )

    ext_ids_view = g.ext_edge_ids_view()
    src_view = g.edge_src_view()
    dst_view = g.edge_dst_view()

    # Decode ext_edge_ids back to (linkIndex, directionBit)
    for i in range(len(ext_ids_view)):
        ext_id = int(ext_ids_view[i])
        link_idx = ext_id >> 1
        dir_bit = ext_id & 1
        direction = "rev" if dir_bit else "fwd"

        # Verify encoding is consistent
        assert link_idx in (0, 1), f"Invalid link_idx {link_idx}"
        assert direction in ("fwd", "rev"), f"Invalid direction {direction}"

        # Verify forward edges have dir_bit=0, reverse have dir_bit=1
        if src_view[i] < dst_view[i]:  # Assuming forward means src < dst
            assert dir_bit == 0, f"Forward edge {i} should have dir_bit=0"
        else:
            assert dir_bit == 1, f"Reverse edge {i} should have dir_bit=1"


if __name__ == "__main__":
    test_ext_edge_ids_without_reverse()
    test_ext_edge_ids_with_manual_reverse()
    test_ext_edge_ids_encoding_scheme()
    print("All ext_edge_ids tests passed!")
