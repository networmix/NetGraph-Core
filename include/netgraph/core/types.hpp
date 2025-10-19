/* Core type aliases and helper structs.
 *
 * For Python developers:
 * - NodeId/EdgeId: int32 (matches np.int32)
 * - Cost: int64 (matches np.int64)
 * - Cap/Flow: double (matches np.float64)
 * - std::span<T>: lightweight view over contiguous arrays (like memoryview, no copy)
 * - std::optional<T>: nullable value (like T | None)
 */
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
// Used as a key in unordered_map (like Python dict with custom __hash__).
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
      // Hash combine formula (similar to Python's hash tuple)
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
enum class FlowPlacement {
  Proportional = 1,    // Distribute flow proportionally to residual capacity (like ECMP with weights)
  EqualBalanced = 2    // Split flow equally across paths (strict equal-cost multipath)
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
