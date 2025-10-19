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

  // Get or create ledger entry for this flow. The ledger tracks per-edge deltas.
  // Python developers: ledger_ is like a defaultdict(list).
  auto& bucket = ledger_[idx];
  bucket.clear();

  // Delegate placement to FlowState, which returns the actual placed flow and
  // populates bucket with per-edge deltas (EdgeId, Flow) pairs.
  Flow placed = fs_.place_on_dag(src, dst, dag, amount, placement, shortest_path, &bucket);

  // Coalesce and filter: drop tiny entries and merge duplicate EdgeIds.
  if (!bucket.empty()) {
    std::vector<std::pair<EdgeId, Flow>> compact;
    compact.reserve(bucket.size());
    // Sort by EdgeId to group duplicates together.
    std::sort(bucket.begin(), bucket.end(), [](auto const& a, auto const& b){ return a.first < b.first; });
    // Merge consecutive entries with the same EdgeId.
    for (std::size_t i=0;i<bucket.size();) {
      EdgeId e = bucket[i].first; double sum = 0.0; std::size_t j=i;
      while (j<bucket.size() && bucket[j].first==e) { sum += bucket[j].second; ++j; }
      if (sum >= kMinFlow) compact.emplace_back(e, static_cast<Flow>(sum));
      i=j;
    }
    bucket.swap(compact);  // replace bucket with compacted version
  }

  // Clean up if no meaningful flow was placed.
  if (bucket.empty() && placed < kMinFlow) {
    ledger_.erase(idx);
    return 0.0;
  }
  return placed;
}

void FlowGraph::remove(const FlowIndex& idx) {
  auto it = ledger_.find(idx);
  if (it == ledger_.end()) return;  // flow not found
  auto& deltas = it->second;

  // Revert this flow's deltas from the FlowState by subtracting them.
  // This restores residual capacity and removes the flow's contribution.
  if (!deltas.empty()) {
    std::vector<std::pair<EdgeId, Flow>> neg(deltas.begin(), deltas.end());
    fs_.apply_deltas(neg, /*add=*/false);  // false = subtract
  }
  ledger_.erase(it);
}

void FlowGraph::remove_by_class(FlowClass flowClass) {
  std::vector<FlowIndex> to_rm;
  to_rm.reserve(ledger_.size());
  for (auto const& kv : ledger_) if (kv.first.flowClass == flowClass) to_rm.push_back(kv.first);
  for (auto const& idx : to_rm) remove(idx);
}

void FlowGraph::reset() noexcept {
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

  // Reconstruct a simple path from the flow's edge deltas.
  // This only succeeds if the flow forms a single path (not a DAG).

  // Step 1: Build adjacency map from edges with positive flow.
  std::unordered_map<NodeId, std::vector<std::pair<NodeId, EdgeId>>> adj;
  const auto& row = g_->row_offsets_view();
  const auto& col = g_->col_indices_view();
  const auto& aei = g_->adj_edge_index_view();

  // Build reverse map: EdgeId -> (src, dst).
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

  // Build adjacency list from deltas.
  for (auto const& pr : deltas) {
    if (pr.second < kMinFlow) continue;
    auto eid = static_cast<std::size_t>(pr.first);
    NodeId u = src_of[eid]; NodeId v = dst_of[eid];
    adj[u].emplace_back(v, pr.first);
  }

  // Step 2: Find starting node (out-degree == 1, in-degree == 0).
  std::unordered_map<NodeId, int> indeg;
  for (auto const& kv : adj)
    for (auto const& pr : kv.second)
      indeg[pr.first]++;

  NodeId start = -1;
  for (auto const& kv : adj) {
    // Candidate: has exactly one outgoing edge and zero incoming edges.
    if (kv.second.size() == 1) {
      if (indeg.find(kv.first) == indeg.end()) {
        start = kv.first;
        break;
      }
    }
  }
  if (start < 0) return {};  // no clear start node

  // Step 3: Walk the path from start, collecting edges.
  std::vector<EdgeId> path;
  NodeId cur = start;
  std::size_t guard = 0;  // prevent infinite loops
  while (adj.count(cur)) {
    if (adj[cur].size() != 1) return {};  // branching detected, not a simple path
    auto [nxt, eid] = adj[cur][0];
    path.push_back(eid);
    cur = nxt;
    if (++guard > static_cast<std::size_t>(g_->num_edges())) return {};  // cycle detected
  }
  return path;
}

} // namespace netgraph::core
