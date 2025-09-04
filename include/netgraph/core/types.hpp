/* Core type aliases and small helper structs used across the library. */
#pragma once

#include <cstdint>

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
    // simple mix
    std::size_t h = static_cast<std::size_t>(static_cast<std::uint64_t>(k.src) * 1469598103934665603ULL);
    h ^= static_cast<std::size_t>(static_cast<std::uint32_t>(k.dst)) + 0x9e3779b97f4a7c15ULL + (h<<6) + (h>>2);
    h ^= static_cast<std::size_t>(static_cast<std::uint32_t>(k.flowClass)) + 0x9e3779b97f4a7c15ULL + (h<<6) + (h>>2);
    h ^= static_cast<std::size_t>(static_cast<std::uint64_t>(k.flowId)) + 0x9e3779b97f4a7c15ULL + (h<<6) + (h>>2);
    return h;
  }
};

enum class FlowPlacement {
  Proportional = 1,
  EqualBalanced = 2
};

enum class EdgeTieBreak { Deterministic = 1, PreferHigherResidual = 2 };

struct EdgeSelection {
  bool multipath { true };
  bool require_capacity { false };
  EdgeTieBreak tie_break { EdgeTieBreak::Deterministic };
};

} // namespace netgraph::core
