/* Core type aliases and helper structs. */
#pragma once

#include <cstdint>
#include <functional>

namespace netgraph::core {

// Node and edge identifiers are signed 32-bit integers.
using NodeId = std::int32_t;
using EdgeId = std::int32_t;
using Cost   = std::int64_t;  // Path cost (64-bit for large accumulations)
using Cap    = double;         // Edge capacity
using Flow = double;           // Flow amount (same unit as capacity)

using FlowClass = std::int32_t;  // Flow priority/class bucket
using FlowId = std::int64_t;     // Unique flow identifier

// FlowIndex uniquely identifies a flow: (src, dst, class, id).
// Used as a key in unordered_map.
struct FlowIndex {
  NodeId src;
  NodeId dst;
  FlowClass flowClass;  // Priority bucket (e.g., for QoS classes)
  FlowId    flowId;     // Unique ID within a FlowPolicy
  friend bool operator==(const FlowIndex& a, const FlowIndex& b) noexcept {
    return a.src==b.src && a.dst==b.dst && a.flowClass==b.flowClass && a.flowId==b.flowId;
  }
};

// Hash function for FlowIndex (enables use in std::unordered_map).
struct FlowIndexHash {
  std::size_t operator()(const FlowIndex& k) const noexcept {
    std::size_t h = 0;
    auto combine = [&h](std::size_t v) {
      // Hash combine formula
      h ^= v + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
    };
    combine(std::hash<NodeId>{}(k.src));
    combine(std::hash<NodeId>{}(k.dst));
    combine(std::hash<FlowClass>{}(k.flowClass));
    combine(std::hash<FlowId>{}(k.flowId));
    return h;
  }
};

// Flow placement strategy for distributing demand across multiple paths.
//
// IMPORTANT SEMANTICS:
//
// - EqualBalanced models *ECMP admission on a fixed SPF DAG* ("one-shot ECMP").
//   We compute a single global scale so that no edge is oversubscribed under
//   equal per-edge splits within each (u->v) group, place once, and stop.
//   If you keep injecting after the first bottleneck saturates, you would
//   change the split set (or oversubscribe and drop), which is intentionally
//   out of scope for EqualBalanced.
//
// - Proportional may be used iteratively (e.g., for max-flow).
enum class FlowPlacement {
  Proportional = 1,    // Distribute flow proportionally to residual capacity (like ECMP with weights)
  EqualBalanced = 2    // Split equally per parallel edge on a fixed DAG (single-pass ECMP admission)
};

// Tie-breaking rule when multiple equal-cost edges exist between the same (u,v) pair.
enum class EdgeTieBreak {
  Deterministic = 1,         // Use edge order from graph construction (reproducible)
  PreferHigherResidual = 2   // Prefer edge with more available capacity
};

// Edge selection policy for shortest path algorithms.
struct EdgeSelection {
  // multi_edge: if true, keep all equal-cost parallel edges per (u,v) pair.
  //             if false, select one edge per (u,v) using tie_break rule.
  bool multi_edge { true };
  // require_capacity: if true, only consider edges with residual capacity > kMinCap.
  bool require_capacity { false };
  EdgeTieBreak tie_break { EdgeTieBreak::Deterministic };
};

} // namespace netgraph::core
