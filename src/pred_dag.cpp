/*
  PredDAG and Path (APPROACH 5): immutable, graph-bound route structures validated once at
  construction, in the documented order, before any traversal.
*/
#include "netgraph/core/pred_dag.hpp"

#include <algorithm>
#include <stdexcept>
#include <string>

#include "netgraph/core/detail/pred_dag_builder.hpp"

namespace netgraph::core {

namespace {

std::string idx(std::size_t i) { return std::to_string(i); }

} // namespace

PredDAG PredDAG::from_csr(const StrictMultiDiGraph& g, std::span<const std::int32_t> offsets,
                          std::span<const NodeId> parents, std::span<const EdgeId> via) {
  const auto n = static_cast<std::size_t>(g.num_nodes());
  const auto m = static_cast<std::int64_t>(g.num_edges());
  // Shape: offsets of length N + 1, from 0, monotone, ending at len(parents) == len(via).
  if (offsets.size() != n + 1) {
    throw std::invalid_argument("PredDAG.from_csr: offsets has length " + idx(offsets.size()) + ", expected N + 1 = " +
                                idx(n + 1));
  }
  if (offsets[0] != 0) {
    throw std::invalid_argument("PredDAG.from_csr: offsets[0] = " + std::to_string(offsets[0]) + ", expected 0");
  }
  for (std::size_t v = 0; v < n; ++v) {
    if (offsets[v + 1] < offsets[v]) {
      throw std::invalid_argument("PredDAG.from_csr: offsets are not monotone at " + idx(v + 1) + " (" +
                                  std::to_string(offsets[v]) + " > " + std::to_string(offsets[v + 1]) + ")");
    }
  }
  const auto total = static_cast<std::size_t>(offsets[n]);
  if (parents.size() != total || via.size() != total) {
    throw std::invalid_argument("PredDAG.from_csr: offsets[N] = " + idx(total) + " but parents has " +
                                idx(parents.size()) + " and via has " + idx(via.size()) + " entries");
  }
  // Entries: ids in range, then each via joining (parent, child).
  for (std::size_t k = 0; k < total; ++k) {
    if (parents[k] < 0 || static_cast<std::size_t>(parents[k]) >= n) {
      throw std::invalid_argument("PredDAG.from_csr: parents[" + idx(k) + "] = " + std::to_string(parents[k]) +
                                  " out of node range [0, " + idx(n) + ")");
    }
    if (via[k] < 0 || via[k] >= m) {
      throw std::invalid_argument("PredDAG.from_csr: via[" + idx(k) + "] = " + std::to_string(via[k]) +
                                  " out of edge range [0, " + std::to_string(m) + ")");
    }
  }
  const auto edge_src = g.edge_src();
  const auto edge_dst = g.edge_dst();
  for (std::size_t v = 0; v < n; ++v) {
    for (auto k = static_cast<std::size_t>(offsets[v]); k < static_cast<std::size_t>(offsets[v + 1]); ++k) {
      const auto e = static_cast<std::size_t>(via[k]);
      if (edge_src[e] != parents[k] || static_cast<std::size_t>(edge_dst[e]) != v) {
        throw std::invalid_argument("PredDAG.from_csr: entry " + idx(k) + " of node " + idx(v) + ": via edge " +
                                    idx(e) + " joins (" + std::to_string(edge_src[e]) + ", " +
                                    std::to_string(edge_dst[e]) + "), not (" + std::to_string(parents[k]) + ", " +
                                    idx(v) + ")");
      }
    }
  }
  // No duplicate (child, via): an edge has one child, so a duplicate is a repeated via.
  {
    std::vector<std::uint8_t> seen(static_cast<std::size_t>(m), 0);
    for (std::size_t k = 0; k < total; ++k) {
      const auto e = static_cast<std::size_t>(via[k]);
      if (seen[e] != 0) {
        throw std::invalid_argument("PredDAG.from_csr: duplicate entry (child " + std::to_string(edge_dst[e]) +
                                    ", via " + idx(e) + ")");
      }
      seen[e] = 1;
    }
  }
  // Acyclic: Kahn over the entries (a self-loop is a cycle).
  {
    std::vector<std::int32_t> out_count(n + 1, 0);
    std::vector<std::int32_t> indeg(n, 0);
    for (std::size_t v = 0; v < n; ++v) {
      indeg[v] = offsets[v + 1] - offsets[v];
      for (auto k = static_cast<std::size_t>(offsets[v]); k < static_cast<std::size_t>(offsets[v + 1]); ++k) {
        ++out_count[static_cast<std::size_t>(parents[k]) + 1];
      }
    }
    for (std::size_t v = 0; v < n; ++v) out_count[v + 1] += out_count[v];
    std::vector<std::int32_t> children(total);
    std::vector<std::int32_t> pos(out_count.begin(), out_count.end() - 1);
    for (std::size_t v = 0; v < n; ++v) {
      for (auto k = static_cast<std::size_t>(offsets[v]); k < static_cast<std::size_t>(offsets[v + 1]); ++k) {
        children[static_cast<std::size_t>(pos[static_cast<std::size_t>(parents[k])]++)] = static_cast<std::int32_t>(v);
      }
    }
    std::vector<std::int32_t> ready;
    for (std::size_t v = 0; v < n; ++v) {
      if (indeg[v] == 0) ready.push_back(static_cast<std::int32_t>(v));
    }
    std::size_t done = 0;
    while (!ready.empty()) {
      const auto u = static_cast<std::size_t>(ready.back());
      ready.pop_back();
      ++done;
      for (auto k = out_count[u]; k < out_count[u + 1]; ++k) {
        const auto c = static_cast<std::size_t>(children[static_cast<std::size_t>(k)]);
        if (--indeg[c] == 0) ready.push_back(static_cast<std::int32_t>(c));
      }
    }
    if (done != n) throw std::invalid_argument("PredDAG.from_csr: the entries contain a cycle");
  }
  // Canonical form: within a node, ascending via EdgeId.
  for (std::size_t v = 0; v < n; ++v) {
    for (auto k = static_cast<std::size_t>(offsets[v]) + 1; k < static_cast<std::size_t>(offsets[v + 1]); ++k) {
      if (via[k] <= via[k - 1]) {
        throw std::invalid_argument("PredDAG.from_csr: the entries of node " + idx(v) +
                                    " are not in ascending via EdgeId (" + std::to_string(via[k - 1]) + " then " +
                                    std::to_string(via[k]) + ")");
      }
    }
  }
  PredDAG out;
  out.graph_uid_ = g.uid();
  out.offsets_.assign(offsets.begin(), offsets.end());
  out.parents_.assign(parents.begin(), parents.end());
  out.via_.assign(via.begin(), via.end());
  out.edges_ = out.via_;
  std::sort(out.edges_.begin(), out.edges_.end());
  return out;
}

namespace detail {

PredDAG PredDAGBuilder::from_edges(const StrictMultiDiGraph& g, std::span<const EdgeId> edges) {
  const auto n = static_cast<std::size_t>(g.num_nodes());
  const auto edge_src = g.edge_src();
  const auto edge_dst = g.edge_dst();
  PredDAG out;
  out.graph_uid_ = g.uid();
  out.offsets_.assign(n + 1, 0);
  std::vector<EdgeId> sorted(edges.begin(), edges.end());
  std::sort(sorted.begin(), sorted.end());
  for (EdgeId e : sorted) ++out.offsets_[static_cast<std::size_t>(edge_dst[static_cast<std::size_t>(e)]) + 1];
  for (std::size_t v = 0; v < n; ++v) out.offsets_[v + 1] += out.offsets_[v];
  out.parents_.resize(sorted.size());
  out.via_.resize(sorted.size());
  std::vector<std::int32_t> pos(out.offsets_.begin(), out.offsets_.end() - 1);
  for (EdgeId e : sorted) {  // ascending: each node's entries stay ascending
    const auto ei = static_cast<std::size_t>(e);
    const auto k = static_cast<std::size_t>(pos[static_cast<std::size_t>(edge_dst[ei])]++);
    out.parents_[k] = edge_src[ei];
    out.via_[k] = e;
  }
  out.edges_ = std::move(sorted);
  return out;
}

PredDAG PredDAGBuilder::from_children(std::uint64_t graph_uid, const std::vector<std::vector<EdgeId>>& per_child,
                                      const std::vector<NodeId>& parent_of) {
  PredDAG out;
  out.graph_uid_ = graph_uid;
  out.offsets_.assign(per_child.size() + 1, 0);
  for (std::size_t v = 0; v < per_child.size(); ++v) {
    out.offsets_[v + 1] = out.offsets_[v] + static_cast<std::int32_t>(per_child[v].size());
    for (EdgeId e : per_child[v]) {
      out.parents_.push_back(parent_of[v]);
      out.via_.push_back(e);
    }
  }
  out.edges_ = out.via_;
  std::sort(out.edges_.begin(), out.edges_.end());
  return out;
}

std::span<const EdgeId> PredDAGBuilder::ascending(const PredDAG& dag) noexcept { return dag.edges_; }

} // namespace detail

Path::Path(const StrictMultiDiGraph& g, std::vector<std::vector<EdgeId>> hops) {
  const auto m = static_cast<std::int64_t>(g.num_edges());
  // Every member in range before any endpoint is read.
  for (std::size_t h = 0; h < hops.size(); ++h) {
    for (EdgeId e : hops[h]) {
      if (e < 0 || e >= m) {
        throw std::invalid_argument("Path: hop " + idx(h) + " member " + std::to_string(e) +
                                    " out of edge range [0, " + std::to_string(m) + ")");
      }
    }
  }
  {
    std::vector<std::uint8_t> seen(static_cast<std::size_t>(m), 0);
    for (const auto& hop : hops) {
      for (EdgeId e : hop) {
        if (seen[static_cast<std::size_t>(e)] != 0) {
          throw std::invalid_argument("Path: duplicate member " + std::to_string(e));
        }
        seen[static_cast<std::size_t>(e)] = 1;
      }
    }
  }
  if (hops.empty()) throw std::invalid_argument("Path: the hop list is empty");
  for (std::size_t h = 0; h < hops.size(); ++h) {
    if (hops[h].empty()) throw std::invalid_argument("Path: hop " + idx(h) + " is empty");
  }
  const auto edge_src = g.edge_src();
  const auto edge_dst = g.edge_dst();
  const auto cost = g.cost();
  for (std::size_t h = 0; h < hops.size(); ++h) {
    const auto first = static_cast<std::size_t>(hops[h].front());
    for (EdgeId e : hops[h]) {
      const auto ei = static_cast<std::size_t>(e);
      if (edge_src[ei] != edge_src[first] || edge_dst[ei] != edge_dst[first]) {
        throw std::invalid_argument("Path: the members of hop " + idx(h) + " do not share (u, v): edge " +
                                    std::to_string(e) + " joins (" + std::to_string(edge_src[ei]) + ", " +
                                    std::to_string(edge_dst[ei]) + "), edge " + idx(first) + " joins (" +
                                    std::to_string(edge_src[first]) + ", " + std::to_string(edge_dst[first]) + ")");
      }
    }
  }
  for (std::size_t h = 1; h < hops.size(); ++h) {
    const NodeId prev = edge_dst[static_cast<std::size_t>(hops[h - 1].front())];
    const NodeId next = edge_src[static_cast<std::size_t>(hops[h].front())];
    if (prev != next) {
      throw std::invalid_argument("Path: hops " + idx(h - 1) + " and " + idx(h) + " are not contiguous (" +
                                  std::to_string(prev) + " then " + std::to_string(next) + ")");
    }
  }
  std::vector<NodeId> nodes;
  nodes.reserve(hops.size() + 1);
  nodes.push_back(edge_src[static_cast<std::size_t>(hops.front().front())]);
  Cost total = 0;
  for (const auto& hop : hops) {
    nodes.push_back(edge_dst[static_cast<std::size_t>(hop.front())]);
    Cost best = cost[static_cast<std::size_t>(hop.front())];
    for (EdgeId e : hop) best = std::min(best, cost[static_cast<std::size_t>(e)]);
    total += best;
  }
  {
    std::vector<NodeId> sorted(nodes);
    std::sort(sorted.begin(), sorted.end());
    const auto dup = std::adjacent_find(sorted.begin(), sorted.end());
    if (dup != sorted.end()) throw std::invalid_argument("Path: node " + std::to_string(*dup) + " is repeated");
  }
  graph_uid_ = g.uid();
  num_nodes_ = g.num_nodes();
  hops_ = std::move(hops);
  nodes_ = std::move(nodes);
  cost_ = total;
}

PredDAG Path::to_dag() const {
  // Built here (not through the builder) because a Path keeps no graph reference: the parent of
  // every member of hop h is nodes_[h], its child nodes_[h + 1]; no node repeats, so each node
  // holds the members of at most one hop.
  std::vector<std::vector<EdgeId>> per_child(static_cast<std::size_t>(num_nodes_));
  std::vector<NodeId> parent_of(static_cast<std::size_t>(num_nodes_), -1);
  for (std::size_t h = 0; h < hops_.size(); ++h) {
    const auto child = static_cast<std::size_t>(nodes_[h + 1]);
    per_child[child] = hops_[h];
    std::sort(per_child[child].begin(), per_child[child].end());
    parent_of[child] = nodes_[h];
  }
  return detail::PredDAGBuilder::from_children(graph_uid_, per_child, parent_of);
}

} // namespace netgraph::core
