#pragma once

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <span>
#include <utility>
#include <vector>

namespace netgraph::core::detail {

// Iterative Tarjan over a validated CSR graph: offsets has n + 1 entries, heads are
// vertex indices, and both vertices and arc positions fit int32. usable(a) selects
// arcs without copying the adjacency. Labels are assigned in component completion
// order; the return value is their count. The predicate is inlined at each call site.
template <class Usable>
std::int32_t strong_components(std::span<const std::int32_t> offsets,
                               std::span<const std::int32_t> heads, Usable usable,
                               std::vector<std::int32_t>& labels) {
  const auto n = offsets.size() - 1;
  labels.assign(n, -1);
  std::vector<std::int32_t> index(n, -1);
  std::vector<std::int32_t> low(n, 0);
  std::vector<std::uint8_t> on(n, 0);
  std::vector<std::int32_t> stack;
  std::vector<std::pair<std::int32_t, std::int32_t>> work;  // (vertex, next adjacency position)
  std::int32_t counter = 0;
  std::int32_t num_components = 0;
  for (std::size_t s = 0; s < n; ++s) {
    if (index[s] >= 0) continue;
    work.emplace_back(static_cast<std::int32_t>(s), offsets[s]);
    index[s] = low[s] = counter++;
    stack.push_back(static_cast<std::int32_t>(s));
    on[s] = 1;
    while (!work.empty()) {
      auto& [v, k] = work.back();
      const auto vi = static_cast<std::size_t>(v);
      if (k < offsets[vi + 1]) {
        const auto a = static_cast<std::size_t>(k++);
        if (!usable(a)) continue;
        const auto w = static_cast<std::size_t>(heads[a]);
        if (index[w] < 0) {
          index[w] = low[w] = counter++;
          stack.push_back(static_cast<std::int32_t>(w));
          on[w] = 1;
          work.emplace_back(static_cast<std::int32_t>(w), offsets[w]);
        } else if (on[w] != 0) {
          low[vi] = std::min(low[vi], index[w]);
        }
        continue;
      }
      if (low[vi] == index[vi]) {
        while (true) {
          const auto w = static_cast<std::size_t>(stack.back());
          stack.pop_back();
          on[w] = 0;
          labels[w] = num_components;
          if (w == vi) break;
        }
        ++num_components;
      }
      const auto done = v;
      work.pop_back();
      if (!work.empty()) {
        const auto u = static_cast<std::size_t>(work.back().first);
        low[u] = std::min(low[u], low[static_cast<std::size_t>(done)]);
      }
    }
  }
  return num_components;
}

} // namespace netgraph::core::detail
