// The nonzero entries of a member-aligned array, preserving member order and owning the result.
#pragma once

#include <algorithm>
#include <cstddef>
#include <span>
#include <utility>
#include <vector>

#include "netgraph/core/types.hpp"

namespace netgraph::core::detail {

inline std::vector<std::pair<EdgeId, Amount>> nonzero_pairs(std::span<const EdgeId> members, std::span<const Amount> amounts) {
  std::vector<std::pair<EdgeId, Amount>> out;
  out.reserve(static_cast<std::size_t>(std::count_if(amounts.begin(), amounts.end(), [](Amount x) { return x != 0; })));
  for (std::size_t m = 0; m < members.size(); ++m) {
    if (amounts[m] != 0) out.emplace_back(members[m], amounts[m]);
  }
  return out;
}

} // namespace netgraph::core::detail
