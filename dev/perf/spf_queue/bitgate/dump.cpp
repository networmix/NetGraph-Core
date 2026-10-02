// Old-versus-new output bit-gate. Built once against the previous sources and
// once against the working tree; the SHA-256 of stdout must match. Exercises
// forward/reverse SPF, max-flow, batch max-flow, sensitivity, FlowPolicy
// placement sequences and KSP through the public C++ API only.
#include "netgraph/core/algorithms.hpp"
#include "netgraph/core/backend.hpp"
#include "netgraph/core/constants.hpp"
#include "netgraph/core/flow_graph.hpp"
#include "netgraph/core/flow_policy.hpp"
#include "netgraph/core/k_shortest_paths.hpp"
#include "netgraph/core/max_flow.hpp"
#include "netgraph/core/options.hpp"
#include "netgraph/core/shortest_paths.hpp"
#include "netgraph/core/strict_multidigraph.hpp"

#include <cstdio>
#include <memory>
#include <random>
#include <vector>

using namespace netgraph::core;

struct Case {
  StrictMultiDiGraph g;
  std::unique_ptr<bool[]> nm, em;
  std::vector<Cap> residual;
  NodeId src = 0, dst = 1;
  std::span<const bool> nmask() const { return {nm.get(), size_t(g.num_nodes())}; }
  std::span<const bool> emask() const { return {em.get(), size_t(g.num_edges())}; }
};

static Case make(unsigned seed, int n, int degree, Cost max_cost, int zero_pct, bool duplex, bool pseudo) {
  std::mt19937 rng(seed);
  std::vector<std::int32_t> s, d; std::vector<Cap> cap; std::vector<Cost> cost;
  auto add = [&](int u, int v, Cost c, double k) { s.push_back(u); d.push_back(v); cost.push_back(c); cap.push_back(k); };
  for (int u = 0; u < n; ++u) for (int k = 0; k < degree; ++k) {
    int v = k == 0 ? (u + 1) % n : int(rng() % n);
    Cost c = (int(rng() % 100) < zero_pct) ? 0 : 1 + Cost(rng() % unsigned(max_cost));
    double w = 1.0 + double(rng() % 100) / 7.0;
    add(u, v, c, w); if (duplex) add(v, u, c, w);
  }
  Case c;
  NodeId src = 0, dst = n / 2;
  int total = n;
  if (pseudo) {
    int k = std::max(1, n / 8);
    for (int u = 0; u < k; ++u) add(n, u, 0, 1e9);
    for (int u = n - k; u < n; ++u) add(u, n + 1, 0, 1e9);
    src = n; dst = n + 1; total = n + 2;
  }
  c.g = StrictMultiDiGraph::from_arrays(total, s, d, cap, cost);
  c.src = src; c.dst = dst;
  int N = c.g.num_nodes(), E = c.g.num_edges();
  c.nm = std::make_unique<bool[]>(N); c.em = std::make_unique<bool[]>(E);
  for (int i = 0; i < N; ++i) c.nm[i] = rng() % 29 != 0 || i == src || i == dst;
  for (int i = 0; i < E; ++i) c.em[i] = rng() % 13 != 0;
  c.residual.assign(c.g.capacity_view().begin(), c.g.capacity_view().end());
  for (int i = 0; i < E; ++i) if (rng() % 11 == 0) c.residual[i] = (rng() % 2) ? 0.0 : kMinCap / 2;
  return c;
}

static void dump_res(const std::pair<std::vector<Cost>, PredDAG>& r) {
  for (auto v : r.first) std::printf("%lld,", (long long)v);
  std::printf("|");
  for (auto v : r.second.parent_offsets) std::printf("%d,", v);
  std::printf("|");
  for (auto v : r.second.parents) std::printf("%d,", v);
  std::printf("|");
  for (auto v : r.second.via_edges) std::printf("%d,", v);
  std::printf("\n");
}

static void dump_summary(double total, const FlowSummary& sm) {
  std::printf("flow %.17g cut:", total);
  for (auto e : sm.min_cut.edges) std::printf("%d,", e);
  std::printf(" costs:");
  for (size_t i = 0; i < sm.costs.size(); ++i) std::printf("%lld=%.17g,", (long long)sm.costs[i], sm.flows[i]);
  std::printf(" ef:");
  for (auto f : sm.edge_flows) std::printf("%.17g,", f);
  std::printf("\n");
}

int main() {
  auto be = make_cpu_backend();
  auto algs = std::make_shared<Algorithms>(be);
  std::vector<Case> cases;
  cases.push_back(make(1, 60, 4, 31, 0, true, false));
  cases.push_back(make(2, 80, 3, 31, 20, true, true));
  cases.push_back(make(3, 120, 5, 1, 0, true, true));
  cases.push_back(make(4, 90, 4, 20, 10, false, false));
  cases.push_back(make(5, 200, 3, 1000, 5, true, true));
  cases.push_back(make(6, 40, 6, 4, 30, true, false));
  for (auto& c : cases) {
    const int n = c.g.num_nodes();
    std::printf("# case n=%d e=%d\n", n, c.g.num_edges());
    // Forward and reverse SPF over the option matrix.
    for (NodeId s : {c.src, NodeId(3), NodeId(n - 1)}) for (int mp = 0; mp < 2; ++mp) for (int me = 0; me < 2; ++me)
    for (int tb = 0; tb < 2; ++tb) for (int rc = 0; rc < 2; ++rc) for (int ur = 0; ur < 2; ++ur) for (int um = 0; um < 2; ++um) {
      EdgeSelection sel; sel.multi_edge = me; sel.require_capacity = rc;
      sel.tie_break = tb ? EdgeTieBreak::PreferHigherResidual : EdgeTieBreak::Deterministic;
      auto res = ur ? std::span<const Cap>(c.residual) : std::span<const Cap>{};
      auto nm = um ? c.nmask() : std::span<const bool>{};
      auto em = um ? c.emask() : std::span<const bool>{};
      for (std::optional<NodeId> dst : {std::optional<NodeId>{}, std::optional<NodeId>{c.dst}, std::optional<NodeId>{NodeId((s + 1) % n)}})
        dump_res(shortest_paths(c.g, s, dst, mp, sel, res, nm, em));
      const auto row = c.g.row_offsets_view(); const auto aei = c.g.adj_edge_index_view();
      std::vector<EdgeId> fan(aei.begin() + row[c.src], aei.begin() + row[c.src + 1]);
      dump_res(shortest_paths_to(c.g, s, mp, sel, res, nm, em, {}));
      if (c.src >= n - 2) dump_res(shortest_paths_to(c.g, c.dst, mp, sel, res, nm, em, fan));   // pseudo source fan-out
    }
    // Max-flow matrix.
    for (int pl = 1; pl <= 3; ++pl) for (int sp = 0; sp < 2; ++sp) for (int rc = 0; rc < 2; ++rc) for (int um = 0; um < 2; ++um) {
      auto nm = um ? c.nmask() : std::span<const bool>{};
      auto em = um ? c.emask() : std::span<const bool>{};
      auto [total, sm] = calc_max_flow(c.g, c.src, c.dst, static_cast<FlowPlacement>(pl), sp, rc, true, true, true, nm, em);
      dump_summary(total, sm);
    }
    auto gh = algs->build_graph(c.g);
    MaxFlowOptions mo; mo.with_edge_flows = true;
    std::vector<std::pair<NodeId, NodeId>> pairs;
    for (int i = 0; i < 12; ++i) pairs.emplace_back(NodeId((i * 7) % n), NodeId((i * 11 + 5) % n));
    for (auto& sm : algs->batch_max_flow(gh, pairs, mo)) dump_summary(sm.total_flow, sm);
    for (auto& [e, f] : algs->sensitivity_analysis(gh, c.src, c.dst, mo)) std::printf("sens %d %.17g\n", e, f);
    // FlowPolicy sequences.
    for (int pl = 1; pl <= 3; ++pl) for (int mp = 0; mp < 2; ++mp) for (int rc = 0; rc < 2; ++rc) {
      FlowGraph fg(c.g);
      ExecutionContext ctx(algs, gh);
      FlowPolicyConfig cfg;
      cfg.flow_placement = static_cast<FlowPlacement>(pl);
      cfg.multipath = mp; cfg.require_capacity = rc;
      cfg.min_flow_count = mp ? 1 : 4; cfg.max_flow_count = mp ? std::optional<int>{} : std::optional<int>{4};
      cfg.node_mask = c.nmask(); cfg.edge_mask = c.emask();
      FlowPolicy pol(ctx, cfg);
      auto r1 = pol.place_demand(fg, c.src, c.dst, 0, 25.0);
      auto r2 = pol.place_demand(fg, c.src, c.dst, 0, 40.0);
      auto r3 = pol.rebalance_demand(fg, c.src, c.dst, 0, 5.0);
      std::printf("policy %d %d %d: %.17g %.17g %.17g %.17g %.17g %.17g placed=%.17g flows=%zu\n", pl, mp, rc,
                  r1.first, r1.second, r2.first, r2.second, r3.first, r3.second, pol.placed_demand(), pol.flows().size());
      std::vector<std::pair<long long, std::pair<long long, double>>> recs;
      for (auto& kv : pol.flows()) recs.push_back({kv.first.flowId, {(long long)kv.second.cost, kv.second.placed_flow}});
      std::sort(recs.begin(), recs.end());
      for (auto& r : recs) std::printf("  flow %lld cost %lld placed %.17g\n", r.first, r.second.first, r.second.second);
      for (auto f : fg.edge_flow_view()) std::printf("%.17g,", f);
      std::printf("\n");
      pol.remove_demand(fg);
      double left = 0; for (auto f : fg.edge_flow_view()) left += f;
      std::printf("after remove %.17g\n", left);
    }
    // KSP.
    for (int k : {1, 3, 6}) for (int uq = 0; uq < 2; ++uq) {
      auto paths = k_shortest_paths(c.g, c.src, c.dst, k, std::optional<double>{1.5}, uq, c.nmask(), c.emask());
      std::printf("ksp k=%d unique=%d n=%zu\n", k, uq, paths.size());
      for (auto& p : paths) dump_res(p);
    }
  }
  return 0;
}
