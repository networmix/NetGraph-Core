#!/usr/bin/env python3
"""Benchmark for measuring netgraph_core performance and profiling overhead.

Runs SPF and flow placement operations on a mesh graph topology.
Useful for:
  - Measuring profiling infrastructure overhead (compare with/without NGRAPH_CORE_PROFILE=1)
  - Regression testing after code changes
  - Comparing performance across versions

Usage:
    # Default run (400-node mesh, 500 SPF + 100 flow ops)
    python dev/benchmark_profiling_overhead.py

    # Quick smoke test
    python dev/benchmark_profiling_overhead.py --quick

    # Custom parameters
    python dev/benchmark_profiling_overhead.py --mesh-size 30 --spf-iters 1000 --flow-iters 200

    # Compare profiling overhead
    python dev/benchmark_profiling_overhead.py                      # baseline
    NGRAPH_CORE_PROFILE=1 python dev/benchmark_profiling_overhead.py  # with profiling

    # JSON output for scripting
    python dev/benchmark_profiling_overhead.py --json
"""

import argparse
import json
import statistics
import sys
import time

import numpy as np

import netgraph_core as ngc


def create_mesh_graph(n: int) -> ngc.StrictMultiDiGraph:
    """Create an n x n mesh graph with bidirectional edges.

    Args:
        n: Grid dimension (creates n*n nodes)

    Returns:
        Graph with 4-connected mesh topology, unit costs, 100.0 capacity per edge.
    """
    num_nodes = n * n
    edges = []

    for row in range(n):
        for col in range(n):
            node = row * n + col
            if col < n - 1:
                right = row * n + (col + 1)
                edges.append((node, right, 1, 100.0))
                edges.append((right, node, 1, 100.0))
            if row < n - 1:
                down = (row + 1) * n + col
                edges.append((node, down, 1, 100.0))
                edges.append((down, node, 1, 100.0))

    src = np.array([e[0] for e in edges], dtype=np.int32)
    dst = np.array([e[1] for e in edges], dtype=np.int32)
    cost = np.array([e[2] for e in edges], dtype=np.int64)
    cap = np.array([e[3] for e in edges], dtype=np.float64)
    ext_edge_ids = np.arange(len(edges), dtype=np.int64)

    return ngc.StrictMultiDiGraph.from_arrays(
        num_nodes, src, dst, cap, cost, ext_edge_ids
    )


def benchmark_spf(
    graph_handle, algs: ngc.Algorithms, num_nodes: int, iterations: int
) -> list[float]:
    """Run SPF iterations and return per-call times in ms."""
    times = []
    edge_sel = ngc.EdgeSelection(
        multi_edge=False,
        require_capacity=False,
        tie_break=ngc.EdgeTieBreak.DETERMINISTIC,
    )

    for i in range(iterations):
        src = i % num_nodes
        start = time.perf_counter()
        algs.spf(graph_handle, src, selection=edge_sel)
        times.append((time.perf_counter() - start) * 1000)

    return times


def benchmark_flow_placement(
    g: ngc.StrictMultiDiGraph, graph_handle, algs: ngc.Algorithms, iterations: int
) -> list[float]:
    """Run flow placement iterations and return per-call times in ms."""
    num_nodes = g.num_nodes()
    times = []

    for i in range(iterations):
        fg = ngc.FlowGraph(g)
        cfg = ngc.FlowPolicyConfig()
        cfg.flow_placement = ngc.FlowPlacement.EQUAL_BALANCED
        cfg.max_flow_count = 16
        cfg.min_flow_count = 16
        cfg.selection = ngc.EdgeSelection(
            multi_edge=False,
            require_capacity=True,
            tie_break=ngc.EdgeTieBreak.PREFER_HIGHER_RESIDUAL,
        )
        cfg.require_capacity = True
        cfg.multipath = False
        cfg.reoptimize_flows_on_each_placement = True
        policy = ngc.FlowPolicy(algs, graph_handle, cfg)

        src = i % num_nodes
        dst = (num_nodes - 1 - i) % num_nodes
        if src == dst:
            dst = (dst + 1) % num_nodes

        start = time.perf_counter()
        policy.place_demand(fg, src, dst, 0, 1000.0)
        times.append((time.perf_counter() - start) * 1000)

    return times


def compute_stats(name: str, times: list[float]) -> dict:
    """Compute statistics for a benchmark run."""
    return {
        "name": name,
        "iterations": len(times),
        "mean_ms": statistics.mean(times),
        "median_ms": statistics.median(times),
        "stdev_ms": statistics.stdev(times) if len(times) > 1 else 0,
        "min_ms": min(times),
        "max_ms": max(times),
        "total_ms": sum(times),
    }


def print_stats(stats: dict) -> None:
    """Print benchmark statistics."""
    print(f"\n{stats['name']}:")
    print(f"  Iterations: {stats['iterations']}")
    print(f"  Mean:   {stats['mean_ms']:.3f} ms")
    print(f"  Median: {stats['median_ms']:.3f} ms")
    print(f"  Stdev:  {stats['stdev_ms']:.3f} ms")
    print(f"  Min:    {stats['min_ms']:.3f} ms")
    print(f"  Max:    {stats['max_ms']:.3f} ms")
    print(f"  Total:  {stats['total_ms']:.1f} ms")


def main():
    parser = argparse.ArgumentParser(
        description="Benchmark netgraph_core performance and profiling overhead.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s                     # default benchmark
  %(prog)s --quick             # fast smoke test
  %(prog)s --mesh-size 30      # larger graph
  %(prog)s --json              # machine-readable output
  NGRAPH_CORE_PROFILE=1 %(prog)s  # measure with profiling enabled
""",
    )
    parser.add_argument(
        "--mesh-size",
        type=int,
        default=20,
        metavar="N",
        help="mesh grid size (N x N nodes, default: 20)",
    )
    parser.add_argument(
        "--spf-iters",
        type=int,
        default=500,
        metavar="N",
        help="number of SPF iterations (default: 500)",
    )
    parser.add_argument(
        "--flow-iters",
        type=int,
        default=100,
        metavar="N",
        help="number of flow placement iterations (default: 100)",
    )
    parser.add_argument(
        "--quick",
        action="store_true",
        help="quick run with minimal iterations (mesh=10, spf=50, flow=10)",
    )
    parser.add_argument(
        "--json",
        action="store_true",
        help="output results as JSON (for scripting)",
    )
    args = parser.parse_args()

    if args.quick:
        args.mesh_size = 10
        args.spf_iters = 50
        args.flow_iters = 10

    profiling_enabled = getattr(ngc, "profiling_enabled", lambda: None)()

    if not args.json:
        print("=" * 60)
        print("NetGraph-Core Benchmark")
        print("=" * 60)
        print(f"\nVersion: {ngc.__version__}")
        if profiling_enabled is not None:
            print(f"Profiling: {'ON' if profiling_enabled else 'OFF'}")
        print(
            f"\nGraph: {args.mesh_size}x{args.mesh_size} mesh ({args.mesh_size**2} nodes)"
        )
        print(f"SPF iterations: {args.spf_iters}")
        print(f"Flow iterations: {args.flow_iters}")
        print("\nCreating graph...", end=" ", flush=True)

    g = create_mesh_graph(args.mesh_size)
    num_nodes = g.num_nodes()
    num_edges = g.num_edges()

    if not args.json:
        print(f"done ({num_nodes} nodes, {num_edges} edges)")

    backend = ngc.Backend.cpu()
    algs = ngc.Algorithms(backend)
    graph_handle = algs.build_graph(g)

    # Warmup (results discarded, profiling stats reset after)
    if not args.json:
        print("\nWarmup...", end=" ", flush=True)
        sys.stdout.flush()
    _ = benchmark_spf(graph_handle, algs, num_nodes, min(10, args.spf_iters))
    _ = benchmark_flow_placement(g, graph_handle, algs, min(5, args.flow_iters))
    if hasattr(ngc, "profiling_reset"):
        ngc.profiling_reset()
    if not args.json:
        print("done", flush=True)
        print("\n" + "-" * 40)
        print("Running benchmarks...", flush=True)
        sys.stdout.flush()

    spf_times = benchmark_spf(graph_handle, algs, num_nodes, args.spf_iters)
    spf_stats = compute_stats("Shortest Paths (SPF)", spf_times)

    flow_times = benchmark_flow_placement(g, graph_handle, algs, args.flow_iters)
    flow_stats = compute_stats("Flow Placement", flow_times)

    total_ops = args.spf_iters + args.flow_iters
    total_time_ms = spf_stats["total_ms"] + flow_stats["total_ms"]
    throughput = total_ops / (total_time_ms / 1000) if total_time_ms > 0 else 0

    if args.json:
        result = {
            "version": ngc.__version__,
            "profiling_enabled": profiling_enabled,
            "config": {
                "mesh_size": args.mesh_size,
                "num_nodes": num_nodes,
                "num_edges": num_edges,
                "spf_iterations": args.spf_iters,
                "flow_iterations": args.flow_iters,
            },
            "spf": {k: v for k, v in spf_stats.items() if k != "name"},
            "flow": {k: v for k, v in flow_stats.items() if k != "name"},
            "summary": {
                "total_ops": total_ops,
                "total_time_ms": total_time_ms,
                "throughput_ops_per_sec": throughput,
            },
        }
        print(json.dumps(result, indent=2))
    else:
        print_stats(spf_stats)
        print_stats(flow_stats)

        print("\n" + "=" * 60)
        print("SUMMARY")
        print("=" * 60)
        print(f"Total operations: {total_ops}")
        print(f"Total time: {total_time_ms:.1f} ms ({total_time_ms / 1000:.2f} s)")
        print(f"Throughput: {throughput:.1f} ops/sec")

        if profiling_enabled:
            print("\n" + "-" * 40)
            print("Profiling stats:", flush=True)
            sys.stdout.flush()
            ngc.profiling_dump()
            sys.stderr.flush()


if __name__ == "__main__":
    main()
