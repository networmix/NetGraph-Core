/*
  FlowGraph — authoritative flow ledger layered over FlowState.

  Tracks per-flow edge deltas to support exact removal and path inspection,
  while delegating residual and aggregate flow management to FlowState.
*/
#include "netgraph/core/flow_graph.hpp"
#include "netgraph/core/constants.hpp"

#include <algorithm>

namespace netgraph::core {

FlowGraph::FlowGraph(const StrictMultiDiGraph& g)
  : g_(&g), fs_(g) {}

Flow FlowGraph::place(const FlowIndex& idx, NodeId src, NodeId dst,
                      const PredDAG& dag, Flow amount,
                      FlowPlacement placement, bool shortest_path) {
  if (amount <= 0.0) return 0.0;
  auto& bucket = ledger_[idx];
  bucket.clear();
  Flow placed = fs_.place_on_dag(src, dst, dag, amount, placement, shortest_path, &bucket);
  // Drop tiny entries and coalesce if needed
  if (!bucket.empty()) {
    std::vector<std::pair<EdgeId, Flow>> compact;
    compact.reserve(bucket.size());
    std::sort(bucket.begin(), bucket.end(), [](auto const& a, auto const& b){ return a.first < b.first; });
    for (std::size_t i=0;i<bucket.size();) {
      EdgeId e = bucket[i].first; double sum = 0.0; std::size_t j=i;
      while (j<bucket.size() && bucket[j].first==e) { sum += bucket[j].second; ++j; }
      if (sum >= kMinFlow) compact.emplace_back(e, static_cast<Flow>(sum));
      i=j;
    }
    bucket.swap(compact);
  }
  if (bucket.empty() && placed < kMinFlow) {
    ledger_.erase(idx);
    return 0.0;
  }
  return placed;
}

void FlowGraph::remove(const FlowIndex& idx) {
  auto it = ledger_.find(idx);
  if (it == ledger_.end()) return;
  auto& deltas = it->second;
  if (!deltas.empty()) {
    // revert
    std::vector<std::pair<EdgeId, Flow>> neg(deltas.begin(), deltas.end());
    fs_.apply_deltas(neg, /*add=*/false);
  }
  ledger_.erase(it);
}

void FlowGraph::remove_by_class(std::int32_t flowClass) {
  std::vector<FlowIndex> to_rm;
  to_rm.reserve(ledger_.size());
  for (auto const& kv : ledger_) if (kv.first.flowClass == flowClass) to_rm.push_back(kv.first);
  for (auto const& idx : to_rm) remove(idx);
}

void FlowGraph::reset() {
  fs_.reset();
  ledger_.clear();
}

std::vector<std::pair<EdgeId, Flow>> FlowGraph::get_flow_edges(const FlowIndex& idx) const {
  auto it = ledger_.find(idx);
  if (it == ledger_.end()) return {};
  return it->second;
}

std::vector<EdgeId> FlowGraph::get_flow_path(const FlowIndex& idx) const {
  auto it = ledger_.find(idx);
  if (it == ledger_.end()) return {};
  const auto& deltas = it->second;
  // Build a single-path check by constructing adjacency with positive flow
  std::unordered_map<NodeId, std::vector<std::pair<NodeId, EdgeId>>> adj;
  const auto& row = g_->row_offsets_view();
  const auto& col = g_->col_indices_view();
  const auto& aei = g_->adj_edge_index_view();
  // Reverse map: EdgeId -> (u,v). Build once
  std::vector<NodeId> src_of(static_cast<std::size_t>(g_->num_edges()));
  std::vector<NodeId> dst_of(static_cast<std::size_t>(g_->num_edges()));
  for (std::int32_t u = 0; u < g_->num_nodes(); ++u) {
    auto s = static_cast<std::size_t>(row[static_cast<std::size_t>(u)]);
    auto e = static_cast<std::size_t>(row[static_cast<std::size_t>(u)+1]);
    for (std::size_t j = s; j < e; ++j) {
      auto v = static_cast<std::int32_t>(col[j]);
      auto eid = static_cast<std::size_t>(aei[j]);
      src_of[eid] = u; dst_of[eid] = v;
    }
  }
  for (auto const& pr : deltas) {
    if (pr.second < kMinFlow) continue;
    auto eid = static_cast<std::size_t>(pr.first);
    NodeId u = src_of[eid]; NodeId v = dst_of[eid];
    adj[u].emplace_back(v, pr.first);
  }
  // Try to find a path: each internal node must have out-degree==1
  // Find a start node: one with out-degree==1 and not an in-edge target
  std::unordered_map<NodeId, int> indeg;
  for (auto const& kv : adj) for (auto const& pr : kv.second) indeg[pr.first]++;
  NodeId start = -1;
  for (auto const& kv : adj) {
    if (kv.second.size() == 1) {
      if (indeg.find(kv.first) == indeg.end()) { start = kv.first; break; }
    }
  }
  if (start < 0) return {};
  std::vector<EdgeId> path;
  NodeId cur = start;
  std::size_t guard = 0;
  while (adj.count(cur)) {
    if (adj[cur].size() != 1) return {};
    auto [nxt, eid] = adj[cur][0];
    path.push_back(eid);
    cur = nxt;
    if (++guard > static_cast<std::size_t>(g_->num_edges())) return {};
  }
  return path;
}

} // namespace netgraph::core
