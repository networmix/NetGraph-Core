/* Core type aliases and small helper structs used across the library. */
#pragma once

#include <cstdint>
#include <functional>

namespace netgraph::core {

// Use signed 32-bit ids for node and edge identifiers at API boundaries and internally.
using NodeId = std::int32_t;
using EdgeId = std::int32_t;
using Cost   = std::int64_t;
using Cap    = double;
// Semantic alias for flow amounts (same unit as capacity)
using Flow = double;

// Identity of a flow: endpoints + class (priority bucket) + per-policy unique id
struct FlowIndex {
  NodeId src;
  NodeId dst;
  std::int32_t flowClass; // small priority bucket
  std::int64_t flowId;    // per-policy unique id
  friend bool operator==(const FlowIndex& a, const FlowIndex& b) noexcept {
    return a.src==b.src && a.dst==b.dst && a.flowClass==b.flowClass && a.flowId==b.flowId;
  }
};

struct FlowIndexHash {
  std::size_t operator()(const FlowIndex& k) const noexcept {
    std::size_t h = 0;
    auto combine = [&h](std::size_t v) {
      // Standard hash combine (boost-like)
      h ^= v + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
    };
    combine(std::hash<NodeId>{}(k.src));
    combine(std::hash<NodeId>{}(k.dst));
    combine(std::hash<std::int32_t>{}(k.flowClass));
    combine(std::hash<std::int64_t>{}(k.flowId));
    return h;
  }
};

enum class FlowPlacement {
  Proportional = 1,
  EqualBalanced = 2
};

enum class EdgeTieBreak { Deterministic = 1, PreferHigherResidual = 2 };

struct EdgeSelection {
  // When true, keep all equal-cost parallel edges per (u,v) adjacency.
  // When false, pick a single edge per (u,v) according to tie_break.
  bool multi_edge { true };
  bool require_capacity { false };
  EdgeTieBreak tie_break { EdgeTieBreak::Deterministic };
};

} // namespace netgraph::core
