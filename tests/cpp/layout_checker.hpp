// An independent checker of the graph's stored layout (the EdgeId order and the CSR arrays that
// graph.hpp describes).
//
// It reads the arrays as spans and the caller's input separately, and derives every expectation
// from the layout's definition, never from the constructor's helpers:
//   1. sizes first, before any element is read: N >= 0, E >= 0, the input arrays of length E,
//      every per-edge array and the input map of length E, both offset arrays of length N + 1,
//      the reverse slots of length E, the internal reverse-source copy of length E;
//   2. then ranges, before any value is used as an index: offsets in [0, E], [0] = 0, [N] = E,
//      non-decreasing; endpoints in [0, N); the map and the slots in [0, E) and permutations;
//      the copied sources in [0, N);
//   3. then the order: the ids of row u have source u and strictly increasing (cost, destination, input);
//      the slots of v hold edges with destination v in strictly increasing (cost, source, input);
//      each copied source equals edge_sources()[incoming_edges()[slot]]; every property is the input
//      array permuted by the map; max_edge_cost and edge_cost_gcd are the full scans (0 without edges).
// A failing step stops the check (later steps would read through unchecked values). The result
// is the list of defect codes of the failing step, empty when the layout is correct.
#pragma once

#include <cstdint>
#include <numeric>
#include <span>
#include <string>
#include <tuple>
#include <vector>

#include "netgraph/core/types.hpp"

namespace netgraph::core::test {

// The caller's input, in input coordinates.
struct LayoutInput {
  std::int64_t num_nodes {0};
  std::span<const NodeId> src {};
  std::span<const NodeId> dst {};
  std::span<const Amount> capacity {};
  std::span<const Cost> cost {};
};

// What a graph stores (or claims to store).
struct LayoutSpans {
  std::int64_t num_nodes {0};
  std::int64_t num_edges {0};
  std::span<const Amount> capacity {};
  std::span<const Cost> cost {};
  std::span<const NodeId> src {};
  std::span<const NodeId> dst {};
  std::span<const std::int32_t> input_positions {};
  std::span<const std::int32_t> outgoing_offsets {};
  std::span<const std::int32_t> incoming_offsets {};
  std::span<const EdgeId> incoming_edges {};
  std::span<const NodeId> incoming_edge_sources {};  // the internal copy of the in-edges' sources
  Cost max_edge_cost {0};
  Cost edge_cost_gcd {0};
};

namespace layout_checker_detail {

inline void offsets_ok(std::span<const std::int32_t> off, std::int64_t n, std::int64_t m, const char* name,
                       std::vector<std::string>& bad) {
  bool range = true;
  bool monotone = true;
  for (std::size_t k = 0; k < off.size(); ++k) {
    if (off[k] < 0 || off[k] > m) range = false;
    if (k > 0 && off[k] < off[k - 1]) monotone = false;
  }
  if (!range) bad.push_back(std::string(name) + ":range");
  if (off[0] != 0) bad.push_back(std::string(name) + ":start");
  if (off[static_cast<std::size_t>(n)] != m) bad.push_back(std::string(name) + ":end");
  if (!monotone) bad.push_back(std::string(name) + ":monotone");
}

template <typename T>
bool permutation_of_range(std::span<const T> a, std::int64_t m) {
  std::vector<std::uint8_t> seen(static_cast<std::size_t>(m), 0);
  for (const T x : a) {
    if (x < 0 || x >= m) return false;
    if (seen[static_cast<std::size_t>(x)] != 0) return false;
    seen[static_cast<std::size_t>(x)] = 1;
  }
  return true;
}

} // namespace layout_checker_detail

inline std::vector<std::string> check_layout(const LayoutSpans& g, const LayoutInput& in) {
  using namespace layout_checker_detail;
  std::vector<std::string> bad;
  const std::int64_t n = g.num_nodes;
  const std::int64_t m = g.num_edges;
  // 1. Sizes (no element is read).
  if (n < 0 || n != in.num_nodes) bad.emplace_back("size:num_nodes");
  if (m < 0 || static_cast<std::size_t>(m) != in.src.size()) bad.emplace_back("size:num_edges");
  if (!bad.empty()) return bad;
  const auto me = static_cast<std::size_t>(m);
  const auto np1 = static_cast<std::size_t>(n) + 1;
  auto size_is = [&](std::size_t got, std::size_t want, const char* name) {
    if (got != want) bad.emplace_back(std::string("size:") + name);
  };
  size_is(in.dst.size(), me, "input.dst");
  size_is(in.capacity.size(), me, "input.capacity");
  size_is(in.cost.size(), me, "input.cost");
  size_is(g.capacity.size(), me, "capacity");
  size_is(g.cost.size(), me, "cost");
  size_is(g.src.size(), me, "src");
  size_is(g.dst.size(), me, "dst");
  size_is(g.input_positions.size(), me, "input_index");
  size_is(g.outgoing_offsets.size(), np1, "row_offsets");
  size_is(g.incoming_offsets.size(), np1, "in_row_offsets");
  size_is(g.incoming_edges.size(), me, "in_adj_edge_index");
  size_is(g.incoming_edge_sources.size(), me, "in_col_indices");
  if (!bad.empty()) return bad;

  // 2. Ranges and permutations (values are read, none is used as an index yet).
  offsets_ok(g.outgoing_offsets, n, m, "row_offsets", bad);
  offsets_ok(g.incoming_offsets, n, m, "in_row_offsets", bad);
  for (std::size_t e = 0; e < me; ++e) {
    if (g.src[e] < 0 || g.src[e] >= n) {
      bad.emplace_back("range:src");
      break;
    }
  }
  for (std::size_t e = 0; e < me; ++e) {
    if (g.dst[e] < 0 || g.dst[e] >= n) {
      bad.emplace_back("range:dst");
      break;
    }
  }
  if (!permutation_of_range(g.input_positions, m)) bad.emplace_back("map:permutation");
  if (!permutation_of_range(g.incoming_edges, m)) bad.emplace_back("slot:permutation");
  for (const NodeId x : g.incoming_edge_sources) {
    if (x < 0 || x >= n) {
      bad.emplace_back("range:in_col_indices");
      break;
    }
  }
  if (!bad.empty()) return bad;

  // 3. Properties, order, the copy, the scans.
  for (std::size_t e = 0; e < me; ++e) {
    const auto i = static_cast<std::size_t>(g.input_positions[e]);
    if (g.capacity[e] != in.capacity[i]) {
      bad.emplace_back("property:capacity");
      break;
    }
  }
  for (std::size_t e = 0; e < me; ++e) {
    if (g.cost[e] != in.cost[static_cast<std::size_t>(g.input_positions[e])]) {
      bad.emplace_back("property:cost");
      break;
    }
  }
  for (std::size_t e = 0; e < me; ++e) {
    if (g.src[e] != in.src[static_cast<std::size_t>(g.input_positions[e])]) {
      bad.emplace_back("property:src");
      break;
    }
  }
  for (std::size_t e = 0; e < me; ++e) {
    if (g.dst[e] != in.dst[static_cast<std::size_t>(g.input_positions[e])]) {
      bad.emplace_back("property:dst");
      break;
    }
  }
  // Forward rows: the ids of row u, in (cost, destination, input) strictly increasing, source u.
  bool fwd_member = true;
  bool fwd_order = true;
  for (std::size_t u = 0; u + 1 < np1; ++u) {
    for (auto e = static_cast<std::size_t>(g.outgoing_offsets[u]); e < static_cast<std::size_t>(g.outgoing_offsets[u + 1]); ++e) {
      if (g.src[e] != static_cast<NodeId>(u)) fwd_member = false;
      if (e > static_cast<std::size_t>(g.outgoing_offsets[u]) &&
          !(std::tuple(g.cost[e - 1], g.dst[e - 1], g.input_positions[e - 1]) <
            std::tuple(g.cost[e], g.dst[e], g.input_positions[e]))) {
        fwd_order = false;
      }
    }
  }
  if (!fwd_member) bad.emplace_back("forward:membership");
  if (!fwd_order) bad.emplace_back("forward:order");
  // Reverse rows: the slots of v, edges with destination v in (cost, source, input) strictly increasing.
  bool rev_member = true;
  bool rev_order = true;
  for (std::size_t v = 0; v + 1 < np1; ++v) {
    const auto lo = static_cast<std::size_t>(g.incoming_offsets[v]);
    const auto hi = static_cast<std::size_t>(g.incoming_offsets[v + 1]);
    for (std::size_t s = lo; s < hi; ++s) {
      const auto e = static_cast<std::size_t>(g.incoming_edges[s]);
      if (g.dst[e] != static_cast<NodeId>(v)) rev_member = false;
      if (s > lo) {
        const auto p = static_cast<std::size_t>(g.incoming_edges[s - 1]);
        if (!(std::tuple(g.cost[p], g.src[p], g.input_positions[p]) < std::tuple(g.cost[e], g.src[e], g.input_positions[e]))) {
          rev_order = false;
        }
      }
    }
  }
  if (!rev_member) bad.emplace_back("reverse:membership");
  if (!rev_order) bad.emplace_back("reverse:order");
  for (std::size_t s = 0; s < me; ++s) {
    if (g.incoming_edge_sources[s] != g.src[static_cast<std::size_t>(g.incoming_edges[s])]) {
      bad.emplace_back("copy:source");
      break;
    }
  }
  // The scans, from the input (0 on an empty graph).
  Cost mx = 0;
  Cost gcd = 0;
  for (const Cost c : in.cost) {
    mx = c > mx ? c : mx;
    gcd = std::gcd(gcd, c);
  }
  if (g.max_edge_cost != mx) bad.emplace_back("scan:max_cost");
  if (g.edge_cost_gcd != gcd) bad.emplace_back("scan:cost_gcd");
  return bad;
}

} // namespace netgraph::core::test
