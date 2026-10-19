/* Per-edge amount arrays (internal): one check and one pair of messages for an array of
   units indexed by EdgeId (a search's residual, the ledger's baseline, the flow core's view).
   The bindings check unnarrowed Python integers with their own labels before calling in. */
#pragma once

#include <cstddef>
#include <span>
#include <stdexcept>
#include <string>

#include "netgraph/core/graph.hpp"
#include "netgraph/core/types.hpp"

namespace netgraph::core::detail {

// std::invalid_argument unless `amounts` has one entry per edge of graph, each in
// [0, the edge's capacity]: the length first, then the first entry outside, in EdgeId order. `label`
// names the array in the message ("residual", "view").
inline void check_edge_amounts(const Graph& graph, std::span<const Amount> amounts, const char* label) {
  const auto capacity = graph.edge_capacity_units();
  if (amounts.size() != capacity.size()) {
    throw std::invalid_argument(std::string(label) + " has length " + std::to_string(amounts.size()) + ", expected " +
                                std::to_string(capacity.size()) + " (one entry per edge)");
  }
  for (std::size_t e = 0; e < amounts.size(); ++e) {
    if (amounts[e] < 0 || amounts[e] > capacity[e]) {
      throw std::invalid_argument(std::string(label) + "[" + std::to_string(e) + "] = " + std::to_string(amounts[e]) +
                                  " is outside [0, edge_capacity_units[" + std::to_string(e) + "] = " +
                                  std::to_string(capacity[e]) + "]");
    }
  }
}

} // namespace netgraph::core::detail
