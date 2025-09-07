/* FlowRecord — immutable identity + current path DAG and placed amount. */
#pragma once

#include <cstdint>
#include <utility>

#include "netgraph/core/shortest_paths.hpp"
#include "netgraph/core/types.hpp"

namespace netgraph::core {

// FlowRecord is a logical unit of routed demand. It stores identity, current
// path (as a PredDAG), path cost, and the volume placed so far. It does not
// mutate any graph state; placement is performed via FlowGraph.
struct FlowRecord {
  FlowIndex index {};
  NodeId src { -1 };
  NodeId dst { -1 };
  PredDAG dag {};
  Cost cost { 0 };
  double placed_flow { 0.0 };

  FlowRecord() noexcept = default;
  FlowRecord(FlowIndex idx, NodeId src, NodeId dst, PredDAG d, Cost c)
      : index(idx), src(src), dst(dst), dag(std::move(d)), cost(c), placed_flow(0.0) {}
};

} // namespace netgraph::core
