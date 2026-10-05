/* A reusable NodeId -> local index map for work on a small part of a large graph (internal). */
#pragma once

#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

#include "netgraph/core/types.hpp"

#if defined(__GNUC__) || defined(__clang__)
#define NETGRAPH_CORE_ALWAYS_INLINE __attribute__((always_inline)) inline
#else
#define NETGRAPH_CORE_ALWAYS_INLINE inline
#endif

namespace netgraph::core::detail {

// Local indices 0.. in first-touch order, backed by dense arrays of the graph's node count
// that are reset through the touched entries only, so a call costs O(touched) after the first
// use. One instance per thread (`local()`); not reentrant: a caller holding it must not call
// another routine that uses it.
//
// reset() sizes both arrays for the graph, so add() never allocates: it is a few instructions,
// noexcept and inlined into the callers' per-edge loops. An out-of-line add() (it grew its list
// with push_back) saved and restored callee-saved registers on every call, and in some process
// layouts those stack reloads stalled: the flow core's subset solve then ran about 1.45 times
// slower for a whole process (the bimodal `kernel.proportional.grid` blocks; NOTES
// "Bimodality").
class NodeMap {
public:
  void reset(std::size_t num_nodes) {
    for (std::size_t i = 0; i < count_; ++i) slot_[static_cast<std::size_t>(nodes_[i])] = -1;
    count_ = 0;
    // Each array is checked against its own size: a failed growth leaves no array short of a
    // size the other one promises.
    if (slot_.size() < num_nodes) slot_.resize(num_nodes, -1);
    if (nodes_.size() < num_nodes) nodes_.resize(num_nodes);
  }
  NETGRAPH_CORE_ALWAYS_INLINE std::int32_t add(NodeId v) noexcept {
    auto& s = slot_[static_cast<std::size_t>(v)];
    if (s < 0) {
      s = static_cast<std::int32_t>(count_);
      nodes_[count_++] = v;
    }
    return s;
  }
  [[nodiscard]] std::int32_t of(NodeId v) const noexcept { return slot_[static_cast<std::size_t>(v)]; }
  [[nodiscard]] std::size_t size() const noexcept { return count_; }
  [[nodiscard]] std::span<const NodeId> nodes() const noexcept { return {nodes_.data(), count_}; }

  static NodeMap& local() {
    thread_local NodeMap map;
    return map;
  }

private:
  std::vector<std::int32_t> slot_ {};
  std::vector<NodeId> nodes_ {};
  std::size_t count_ {0};
};

} // namespace netgraph::core::detail
