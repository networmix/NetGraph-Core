/*
  FlowGraph — authoritative flow ledger layered over FlowState.

  Tracks per-flow edge allocations to support exact removal and path inspection,
  while delegating residual and aggregate flow management to FlowState.
*/
#include "netgraph/core/flow_graph.hpp"
#include "netgraph/core/constants.hpp"

#include <algorithm>

namespace netgraph::core {

FlowGraph::FlowGraph(const StrictMultiDiGraph& g)
  : g_(&g), fs_(g) {
}

Flow FlowGraph::place(const FlowIndex& idx, NodeId src, NodeId dst,
                      const PredDAG& dag, Flow amount,
                      FlowPlacement placement) {
  if (amount <= 0.0) return 0.0;

  // Get or create ledger entry for this flow. The ledger tracks per-edge
  // cumulative amounts contributed by this flow (not just last placement).
  // We append fresh deltas to the existing bucket, then coalesce by EdgeId.
  auto& bucket = ledger_[idx];

  // Delegate placement to FlowState, which returns the actual placed flow and
  // populates bucket with per-edge allocations (EdgeId, Flow) pairs.
  Flow placed = fs_.place_on_dag(src, dst, dag, amount, placement, &bucket);

  // Coalesce and filter: merge duplicate EdgeIds. Keep any positive totals
  // (do not drop sub-kMinFlow amounts) to preserve exact reversibility.
  if (!bucket.empty()) {
    std::vector<std::pair<EdgeId, Flow>> compact;
    compact.reserve(bucket.size());
    // Sort by EdgeId to group duplicates together.
    std::sort(bucket.begin(), bucket.end(), [](auto const& a, auto const& b){ return a.first < b.first; });
    // Merge consecutive entries with the same EdgeId.
    for (std::size_t i=0;i<bucket.size();) {
      EdgeId e = bucket[i].first; double sum = 0.0; std::size_t j=i;
      while (j<bucket.size() && bucket[j].first==e) { sum += bucket[j].second; ++j; }
      if (sum > 0.0) compact.emplace_back(e, static_cast<Flow>(sum));
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
  const auto& deltas = it->second;
  // Revert this flow's allocations from the FlowState by subtracting them.
  if (!deltas.empty()) {
    fs_.apply_deltas(deltas, /*add=*/false);  // false = subtract
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

  // Reconstruct a simple path from the flow's edge allocations.
  // This only succeeds if the flow forms a single path (not a DAG).

  // Step 1: Build adjacency map from edges with positive flow.
  std::unordered_map<NodeId, std::vector<std::pair<NodeId, EdgeId>>> adj;
  // Build adjacency list from allocations using cached src/dst.
  for (auto const& pr : deltas) {
    if (pr.second < kEpsilon) continue;
    auto eid = static_cast<std::size_t>(pr.first);
    NodeId u = g_->edge_src_view()[eid]; NodeId v = g_->edge_dst_view()[eid];
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
