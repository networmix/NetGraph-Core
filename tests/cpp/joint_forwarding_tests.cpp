// The loss rule over several demands at once (detail::lossy of ForwardingDemand values: the demand
// level above sources and entries, A+ and the cycle rule across demands, per-demand black holes),
// natively: random tunnel demands against an independent reference evaluator (cycles, black holes,
// serial pairs and round counts included); one pooled demand against the pooled kernel it replaces
// in a list; random pooled and tunnel demands without a cyclic group against an independent
// synchronous sweep; hand-derived cases (the demand level above tunnels and sources, a pooled demand
// sharing a link, pooled roots blocking serial pairs, the source-labelled round state of a pooled
// cycle, link sums beyond 2^63 and 2^64). The references call none of the library's apportionment
// helpers (program_reference's round robin and highest averages).
#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cstdint>
#include <functional>
#include <map>
#include <limits>
#include <memory>
#include <optional>
#include <random>
#include <set>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "netgraph/core/detail/forwarding.hpp"
#include "netgraph/core/detail/private_hooks.hpp"
#include "netgraph/core/detail/program_data.hpp"
#include "netgraph/core/detail/route_kernels.hpp"
#include "netgraph/core/flow_ledger.hpp"
#include "netgraph/core/flow_policy.hpp"
#include "netgraph/core/forwarding_program.hpp"
#include "test_utils.hpp"

using namespace netgraph::core;
using namespace netgraph::core::test;

namespace {

namespace ref = program_reference;

using Amounts = std::vector<Amount>;
using Pairs = std::vector<std::pair<EdgeId, Amount>>;
using detail::PreparedRoute;
using detail::ProgramFlow;

const std::string kPeriodTwo = "lossy forwarding iteration entered a period-2 cycle";
const std::string kRoundLimit = "lossy forwarding iteration did not reach a fixed point within the round limit";

// ---------------------------------------------------------------- the reference of simple paths

// One component of the reference: one tunnel of one demand with its own offer, its route
// (EdgeIds), its label (demand, source, tunnel) and the links its demand black-holes.
struct Component {
  Amount offer {0};
  std::vector<EdgeId> route {};
  std::tuple<int, int, int> label {};
  std::set<EdgeId> dark {};
};

struct RefResult {
  std::string error {};  // empty: solved
  std::int64_t rounds {0};  // every round of every cyclic group
  std::vector<Amounts> arrival {};
  std::vector<Amounts> carried {};
  Amounts delivered {};
};

// The loss rule at a link: demands by H, then sources, then tunnels, over the arrivals of one link's members.
Amounts hierarchical(Amount total, const Amounts& arrivals, const std::vector<std::tuple<int, int, int>>& labels) {
  Amounts out(arrivals.size(), 0);
  std::function<void(const std::vector<std::size_t>&, Amount, int)> split = [&](const std::vector<std::size_t>& ids, Amount amount,
                                                                              int depth) {
    if (depth == 3) {
      out[ids.front()] = amount;
      return;
    }
    std::map<int, std::vector<std::size_t>> groups;
    for (auto i : ids) {
      const int key = depth == 0 ? std::get<0>(labels[i]) : depth == 1 ? std::get<1>(labels[i]) : std::get<2>(labels[i]);
      groups[key].push_back(i);
    }
    Amounts w;
    for (const auto& [k, g] : groups) {
      Amount s = 0;
      for (auto i : g) s += arrivals[i];
      w.push_back(s);
    }
    auto q = ref::highest_averages(amount, w);
    if (q.empty()) q.assign(w.size(), 0);
    std::size_t j = 0;
    for (const auto& [k, g] : groups) split(g, q[j++], depth + 1);
  };
  std::vector<std::size_t> all(arrivals.size());
  for (std::size_t i = 0; i < all.size(); ++i) all[i] = i;
  split(all, total, 0);
  return out;
}

// A reference evaluator of the loss rule for simple-path components, written apart from the
// library (its own SCCs, splits and rounds): A+ where arrivals are determined, serial normalization, bounded synchronous
// rounds from zero in every cyclic group, the period-two and round-limit errors.
RefResult reference(const std::vector<Component>& comps, const Amounts& cap) {
  const std::size_t n = comps.size();
  auto black = [&](std::size_t i, std::size_t j) { return comps[i].dark.count(comps[i].route[j]) != 0; };
  // Dead: can carry nothing (a zero offer, or behind its own black hole); no part in the loads,
  // the waits or serial normalization.
  auto dead = [&](std::size_t i, std::size_t j) {
    if (comps[i].offer == 0) return true;
    for (std::size_t k = 0; k < j; ++k) {
      if (black(i, k)) return true;
    }
    return false;
  };
  std::map<EdgeId, std::vector<std::pair<std::size_t, std::size_t>>> occ;
  for (std::size_t i = 0; i < n; ++i) {
    for (std::size_t j = 0; j < comps[i].route.size(); ++j) {
      if (!black(i, j) && !dead(i, j)) occ[comps[i].route[j]].emplace_back(i, j);
    }
  }
  std::set<EdgeId> binding;
  for (const auto& [e, members] : occ) {
    Amount upper = 0;
    for (const auto& [i, j] : members) upper += comps[i].offer;
    if (upper > cap[static_cast<std::size_t>(e)]) binding.insert(e);
  }
  for (const auto& [e, members] : occ) {
    std::vector<std::pair<std::size_t, std::size_t>> nexts;
    bool all = true;
    for (const auto& [i, j] : members) {
      if (j + 1 >= comps[i].route.size()) all = false;
      nexts.emplace_back(i, j + 1);
    }
    if (!all) continue;
    std::set<EdgeId> g;
    for (const auto& [i, j] : nexts) g.insert(comps[i].route[j]);
    if (g.size() != 1) continue;
    const EdgeId f = *g.begin();
    auto at = occ.find(f);
    auto theirs = at == occ.end() ? std::vector<std::pair<std::size_t, std::size_t>> {} : at->second;
    std::sort(nexts.begin(), nexts.end());
    std::sort(theirs.begin(), theirs.end());
    if (f != e && cap[static_cast<std::size_t>(f)] == cap[static_cast<std::size_t>(e)] && nexts == theirs) binding.erase(f);
  }
  // Nodes: routers (0, i, j) and services (1, e, 0).
  using Node = std::tuple<int, std::int64_t, std::int64_t>;
  std::map<Node, std::set<Node>> graph;
  auto arc = [&](const Node& a, const Node& b) {
    graph[a].insert(b);
    graph[b];
  };
  for (std::size_t i = 0; i < n; ++i) {
    if (comps[i].offer == 0) continue;
    for (std::size_t j = 0; j < comps[i].route.size(); ++j) {
      const EdgeId e = comps[i].route[j];
      const Node r {0, static_cast<std::int64_t>(i), static_cast<std::int64_t>(j)};
      const Node r2 {0, static_cast<std::int64_t>(i), static_cast<std::int64_t>(j + 1)};
      if (binding.count(e) != 0 && !black(i, j) && !dead(i, j)) {
        arc(r, Node {1, e, 0});
        arc(Node {1, e, 0}, r2);
      } else {
        arc(r, r2);
      }
    }
  }
  // Strongly connected components (recursive Tarjan; the graphs are small).
  std::map<Node, int> index;
  std::map<Node, int> low;
  std::vector<Node> stack;
  std::set<Node> on;
  std::vector<std::vector<Node>> groups;
  std::function<void(const Node&)> strong = [&](const Node& v) {
    index[v] = low[v] = static_cast<int>(index.size());
    stack.push_back(v);
    on.insert(v);
    for (const auto& w : graph[v]) {
      if (index.count(w) == 0) {
        strong(w);
        low[v] = std::min(low[v], low[w]);
      } else if (on.count(w) != 0) {
        low[v] = std::min(low[v], index[w]);
      }
    }
    if (low[v] == index[v]) {
      std::vector<Node> group;
      while (true) {
        const Node w = stack.back();
        stack.pop_back();
        on.erase(w);
        group.push_back(w);
        if (w == v) break;
      }
      groups.push_back(group);
    }
  };
  for (const auto& [v, _] : graph) {
    if (index.count(v) == 0) strong(v);
  }
  std::map<Node, std::size_t> where;
  for (std::size_t g = 0; g < groups.size(); ++g) {
    for (const auto& v : groups[g]) where[v] = g;
  }
  std::vector<std::set<std::size_t>> preds(groups.size());
  for (const auto& [a, bs] : graph) {
    for (const auto& b : bs) {
      if (where[a] != where[b]) preds[where[b]].insert(where[a]);
    }
  }
  RefResult out;
  out.arrival.resize(n);
  out.carried.resize(n);
  for (std::size_t i = 0; i < n; ++i) {
    out.arrival[i].assign(comps[i].route.size(), 0);
    out.carried[i].assign(comps[i].route.size(), 0);
  }
  auto router = [&](std::size_t i, std::size_t j) {
    if (j >= comps[i].route.size()) return;
    out.arrival[i][j] = j == 0 ? comps[i].offer : out.carried[i][j - 1];
    if (black(i, j)) {
      out.carried[i][j] = 0;
    } else if (binding.count(comps[i].route[j]) == 0) {
      out.carried[i][j] = out.arrival[i][j];
    }
  };
  auto service = [&](EdgeId e) {
    const auto& members = occ[e];
    Amounts w;
    std::vector<std::tuple<int, int, int>> labels;
    Amount sum = 0;
    for (const auto& [i, j] : members) {
      w.push_back(out.arrival[i][j]);
      labels.push_back(comps[i].label);
      sum += out.arrival[i][j];
    }
    return hierarchical(std::min(cap[static_cast<std::size_t>(e)], sum), w, labels);
  };
  std::set<std::size_t> done;
  while (done.size() < groups.size()) {
    // Any ready group: the order of independent groups changes no value and no round count.
    std::size_t g = 0;
    while (done.count(g) != 0 || !std::includes(done.begin(), done.end(), preds[g].begin(), preds[g].end())) ++g;
    done.insert(g);
    const auto& nodes = groups[g];
    if (nodes.size() == 1 && graph[nodes.front()].count(nodes.front()) == 0) {
      const auto& v = nodes.front();
      if (std::get<0>(v) == 0) {
        router(static_cast<std::size_t>(std::get<1>(v)), static_cast<std::size_t>(std::get<2>(v)));
      } else {
        const EdgeId e = static_cast<EdgeId>(std::get<1>(v));
        const auto c = service(e);
        for (std::size_t k = 0; k < c.size(); ++k) out.carried[occ[e][k].first][occ[e][k].second] = c[k];
      }
      continue;
    }
    std::vector<EdgeId> es;
    std::vector<std::pair<std::size_t, std::size_t>> inner;
    for (const auto& v : nodes) {
      if (std::get<0>(v) == 1) {
        es.push_back(static_cast<EdgeId>(std::get<1>(v)));
      } else {
        inner.emplace_back(static_cast<std::size_t>(std::get<1>(v)), static_cast<std::size_t>(std::get<2>(v)));
      }
    }
    std::sort(es.begin(), es.end());
    std::sort(inner.begin(), inner.end());
    std::vector<std::pair<std::size_t, std::size_t>> coords;
    for (EdgeId e : es) {
      for (const auto& c : occ[e]) coords.push_back(c);
    }
    Amounts previous(coords.size(), 0);
    std::optional<Amounts> before;
    for (const auto& [i, j] : coords) out.carried[i][j] = 0;
    bool solved = false;
    const auto limit = 64 * static_cast<std::int64_t>(es.size());
    for (std::int64_t step = 1; step <= limit; ++step) {
      ++out.rounds;
      for (const auto& [i, j] : inner) router(i, j);
      Amounts next;
      for (EdgeId e : es) {
        const auto c = service(e);
        next.insert(next.end(), c.begin(), c.end());
      }
      for (std::size_t k = 0; k < coords.size(); ++k) out.carried[coords[k].first][coords[k].second] = next[k];
      if (next == previous) {
        solved = true;
        break;
      }
      if (before.has_value() && next == *before) {
        out.error = kPeriodTwo;
        return out;
      }
      before = previous;
      previous = next;
    }
    if (!solved) {
      out.error = kRoundLimit;
      return out;
    }
  }
  out.delivered.resize(n);
  for (std::size_t i = 0; i < n; ++i) out.delivered[i] = out.carried[i].empty() ? comps[i].offer : out.carried[i].back();
  return out;
}

// The reference's own controls: D-crossing and D-C2 (two fixed points, period two
// at round 3), a dark arrival dropped where it arrives.
TEST(JointForwarding, ReferenceControls) {
  // D-crossing: A on e, f (8); B on e (4); C on f (4); V(e) = 9, V(f) = 7.
  const auto r = reference({{8, {0, 1}, {0, 0, 0}, {}}, {4, {0}, {1, 0, 0}, {}}, {4, {1}, {2, 0, 0}, {}}}, {9, 7});
  EXPECT_EQ(r.error, "");
  EXPECT_EQ(r.delivered, (Amounts {4, 3, 3}));
  const auto c2 = reference({{1, {0, 1}, {0, 0, 0}, {}}, {2, {1, 0}, {1, 0, 0}, {}}}, {1, 2});
  EXPECT_EQ(c2.error, kPeriodTwo);
  EXPECT_EQ(c2.rounds, 3);
  const auto d = reference({{5, {0, 1}, {0, 0, 0}, {0}}, {3, {1}, {1, 0, 0}, {}}}, {10, 4});
  EXPECT_EQ(d.delivered, (Amounts {0, 3}));
  EXPECT_EQ(d.arrival[0], (Amounts {5, 0}));
}

// ---------------------------------------------------------------- the library's side

// A demand of tunnels for the library: per head its injection, per tunnel (head index, EdgeIds,
// weight), its black-holed links.
struct TunnelDemand {
  std::vector<NodeId> heads {};
  Amounts injections {};
  std::vector<std::tuple<std::size_t, std::vector<EdgeId>, Amount>> tunnels {};
  std::vector<EdgeId> dark {};
};

struct Built {
  std::vector<std::shared_ptr<const ForwardingProgram>> programs;
  std::vector<std::shared_ptr<const Eligibility>> black;
  std::vector<std::unique_ptr<detail::ForwardingView>> views;
};

Built build(const Graph& g, const std::vector<TunnelDemand>& demands, NodeId dst) {
  Built b;
  for (const auto& d : demands) {
    std::vector<RouteDefinition> routes;
    for (const auto& [h, es, w] : d.tunnels) {
      std::vector<std::vector<EdgeId>> hops;
      for (EdgeId e : es) hops.push_back({e});
      routes.push_back(RouteDefinition {d.heads[h], std::make_shared<const Path>(g, hops), w});
    }
    b.programs.push_back(ForwardingProgram::from_routes(g, d.heads, std::vector<NodeId> {dst}, routes));
    b.black.push_back(d.dark.empty() ? nullptr : std::make_shared<const Eligibility>(g.eligibility({}, d.dark)));
    b.views.push_back(std::make_unique<detail::ForwardingView>(
        detail::make_forwarding_view(g, b.programs.back(), nullptr, std::vector<NodeId> {dst}, b.black.back().get())));
  }
  return b;
}

detail::AvailableCapacityFn capacities(const Graph& g) {
  return [&g](std::span<const EdgeId> edges, std::span<Amount> out) {
    for (std::size_t i = 0; i < edges.size(); ++i) out[i] = g.edge_capacity_units()[static_cast<std::size_t>(edges[i])];
  };
}

std::vector<ProgramFlow> joint(const Graph& g, const Built& b, const std::vector<Amounts>& injections) {
  std::vector<detail::ForwardingDemand> demands;
  for (std::size_t d = 0; d < b.views.size(); ++d) demands.push_back({b.views[d].get(), injections[d]});
  return detail::lossy(demands, capacities(g));
}

// Random tunnel demands (one to three, each with one or two heads and one or two tunnels per head)
// on random small graphs with cycles, some links black-holed for one demand, against the
// reference: every hop's arrival and carried amount, each tunnel's delivery, the errors, and the
// rounds of every cyclic group (hooks). Each tunnel ends at the one destination.
TEST(JointForwarding, TunnelDemandsMatchTheReference) {
  std::mt19937 rng(20261020);
  int checked = 0;
  int cyclic = 0;      // trials whose evaluation iterated a group
  int errors = 0;      // trials that raised a cycle-rule error (rare: the deterministic cases are below)
  int dark = 0;        // trials with a black-holed link on some forwarded tunnel
  int crossings = 0;   // trials where some link carries two demands
  for (int trial = 0; trial < 3000; ++trial) {
    const std::int32_t inner = 5;
    const NodeId dst = inner;
    std::vector<InEdge> rows;
    for (NodeId u = 0; u < inner; ++u) {
      for (NodeId w = 0; w < inner; ++w) {
        if (u != w && rng() % 2 == 0) rows.push_back({u, w, static_cast<Amount>(rng() % 4), 1});
      }
      rows.push_back({u, dst, static_cast<Amount>(rng() % 7), 1});
    }
    const auto g = shared(make_graph(inner + 1, rows));
    // A random simple walk from h over inner nodes (1 to 3 inner edges), then its exit to dst.
    auto walk = [&](NodeId h) {
      std::vector<EdgeId> path;
      std::set<NodeId> seen {h};
      NodeId u = h;
      const int len = static_cast<int>(rng() % 5);
      for (int k = 0; k < len; ++k) {
        std::vector<EdgeId> out;
        for (EdgeId e = 0; e < g->num_edges(); ++e) {
          const auto ei = static_cast<std::size_t>(e);
          if (g->edge_sources()[ei] == u && g->edge_destinations()[ei] != dst && seen.count(g->edge_destinations()[ei]) == 0) out.push_back(e);
        }
        if (out.empty()) break;
        const EdgeId e = out[rng() % out.size()];
        path.push_back(e);
        u = g->edge_destinations()[static_cast<std::size_t>(e)];
        seen.insert(u);
      }
      for (EdgeId e = 0; e < g->num_edges(); ++e) {
        if (g->edge_sources()[static_cast<std::size_t>(e)] == u && g->edge_destinations()[static_cast<std::size_t>(e)] == dst) {
          path.push_back(e);
          break;
        }
      }
      return path;
    };
    std::vector<TunnelDemand> demands(1 + rng() % 4);
    std::vector<Component> comps;
    for (std::size_t d = 0; d < demands.size(); ++d) {
      auto& dm = demands[d];
      std::set<NodeId> heads;
      const int nh = 1 + static_cast<int>(rng() % 2);
      while (static_cast<int>(heads.size()) < nh) heads.insert(static_cast<NodeId>(rng() % inner));
      dm.heads.assign(heads.begin(), heads.end());
      for (std::size_t h = 0; h < dm.heads.size(); ++h) {
        const int nt = 1 + static_cast<int>(rng() % 2);
        for (int t = 0; t < nt; ++t) dm.tunnels.emplace_back(h, walk(dm.heads[h]), 1 + static_cast<Amount>(rng() % 3));
        dm.injections.push_back(static_cast<Amount>(rng() % 13));
      }
      if (rng() % 3 == 0) {
        std::set<EdgeId> black;
        for (const auto& [h, es, w] : dm.tunnels) {
          if (rng() % 2 == 0) black.insert(es[rng() % es.size()]);
        }
        dm.dark.assign(black.begin(), black.end());
      }
      // The components: per head, its tunnels' highest-averages amounts of its injection.
      for (std::size_t h = 0; h < dm.heads.size(); ++h) {
        Amounts w;
        std::vector<std::size_t> ids;
        for (std::size_t t = 0; t < dm.tunnels.size(); ++t) {
          if (std::get<0>(dm.tunnels[t]) == h) {
            ids.push_back(t);
            w.push_back(std::get<2>(dm.tunnels[t]));
          }
        }
        const auto split = ref::highest_averages(dm.injections[h], w);
        for (std::size_t k = 0; k < ids.size(); ++k) {
          const auto& [hh, es, ww] = dm.tunnels[ids[k]];
          comps.push_back({split[k], es, {static_cast<int>(d), static_cast<int>(h), static_cast<int>(ids[k])},
                           std::set<EdgeId>(dm.dark.begin(), dm.dark.end())});
        }
      }
    }
    Amounts cap(g->edge_capacity_units().begin(), g->edge_capacity_units().end());
    const auto want = reference(comps, cap);
    const auto b = build(*g, demands, dst);
    std::vector<Amounts> injections;
    for (const auto& dm : demands) injections.push_back(dm.injections);
    const std::string where = "trial " + std::to_string(trial);
    std::vector<ProgramFlow> got;
    std::string error;
    try {
      got = joint(*g, b, injections);
    } catch (const std::invalid_argument& e) {
      error = e.what();
    }
    ASSERT_EQ(error, want.error) << where;
#if NETGRAPH_CORE_PRIVATE_HOOKS
    EXPECT_EQ(detail::thread_lossy_counters().rounds, want.rounds) << where;
    cyclic += detail::thread_lossy_counters().groups > 0 ? 1 : 0;
#endif
    errors += error.empty() ? 0 : 1;
    if (!error.empty()) {
      ++checked;
      continue;
    }
    // Per component (in demand, then head, then tunnel order of the components list): its entry.
    std::size_t c = 0;
    std::map<EdgeId, std::set<int>> users;
    for (std::size_t d = 0; d < demands.size(); ++d) {
      const auto& dm = demands[d];
      for (std::size_t h = 0; h < dm.heads.size(); ++h) {
        for (std::size_t t = 0; t < dm.tunnels.size(); ++t) {
          if (std::get<0>(dm.tunnels[t]) != h) continue;
          const auto& f = got[d].entries[t];
          const auto& es = std::get<1>(dm.tunnels[t]);
          Pairs flows;
          Pairs drops;
          for (std::size_t j = 0; j < es.size(); ++j) {
            if (want.carried[c][j] != 0) {
              flows.emplace_back(es[j], want.carried[c][j]);
              users[es[j]].insert(static_cast<int>(d));
            }
            if (want.arrival[c][j] != want.carried[c][j]) drops.emplace_back(es[j], want.arrival[c][j] - want.carried[c][j]);
            dark += comps[c].dark.count(es[j]) != 0 && comps[c].offer > 0 ? 1 : 0;
          }
          std::sort(flows.begin(), flows.end());
          std::sort(drops.begin(), drops.end());
          EXPECT_EQ(f.assigned, comps[c].offer) << where << " component " << c;
          EXPECT_EQ(f.delivered, want.delivered[c]) << where << " component " << c;
          EXPECT_EQ(f.edge_flows, flows) << where << " component " << c;
          EXPECT_EQ(f.edge_drops, drops) << where << " component " << c;
          ++c;
        }
      }
    }
    crossings += std::any_of(users.begin(), users.end(), [](const auto& u) { return u.second.size() > 1; }) ? 1 : 0;
    ++checked;
    if (::testing::Test::HasFailure()) return;
  }
  RecordProperty("checked", checked);
  RecordProperty("cyclic", cyclic);
  RecordProperty("errors", errors);
  RecordProperty("dark", dark);
  RecordProperty("crossings", crossings);
  EXPECT_EQ(checked, 3000);
#if NETGRAPH_CORE_PRIVATE_HOOKS
  EXPECT_GT(cyclic, 20);  // not vacuous: groups iterated
#endif
  EXPECT_GT(dark, 100);
  EXPECT_GT(crossings, 300);
}

// The cycle rule's errors across demands, against the reference: D-C2 on the two-barrier cycle
// (A: X then Y, offer 1; B: Y then X, offer 2; V = 1, 2): period two at round 3; three barriers
// that never settle (two demands of 1 and 2): the round limit after 64 * 3 rounds. Nothing returns.
TEST(JointForwarding, CycleErrorsMatchTheReference) {
  // CYC: e0 0->1 (X), e1 1->2, e2 1->5, e3 2->3 (Y), e4 3->0, e5 3->5, e6 4->0, e7 4->2.
  const auto cyc = shared(make_graph(6, {{0, 1, 1, 1}, {1, 2, 1000, 1}, {1, 5, 1000, 1}, {2, 3, 2, 1}, {3, 0, 1000, 1},
                                         {3, 5, 1000, 1}, {4, 0, 1000, 1}, {4, 2, 1000, 1}}));
  ASSERT_EQ(id_of_input(*cyc), (std::vector<EdgeId> {0, 1, 2, 3, 4, 5, 6, 7}));
  const std::vector<EdgeId> t0 = {6, 0, 1, 3, 5};
  const std::vector<EdgeId> t1 = {7, 3, 4, 0, 2};
  const Amounts cap {1, 1000, 1000, 2, 1000, 1000, 1000, 1000};
  const auto want = reference({{1, t0, {0, 0, 0}, {}}, {2, t1, {1, 0, 0}, {}}}, cap);
  EXPECT_EQ(want.error, kPeriodTwo);
  EXPECT_EQ(want.rounds, 3);
  const auto b = build(*cyc, {{{4}, {1}, {{0, t0, 1}}, {}}, {{4}, {2}, {{0, t1, 1}}, {}}}, 5);
  EXPECT_EQ(error_of([&] { (void)joint(*cyc, b, {{1}, {2}}); }), kPeriodTwo);
#if NETGRAPH_CORE_PRIVATE_HOOKS
  EXPECT_EQ(detail::thread_lossy_counters().rounds, 3);
#endif
  const auto lim = shared(make_graph(8, {{0, 1, 1, 1}, {1, 4, 1000, 1}, {2, 3, 2, 1}, {3, 0, 1000, 1}, {3, 7, 1000, 1},
                                         {4, 5, 1, 1}, {5, 2, 1000, 1}, {5, 7, 1000, 1}, {6, 0, 1000, 1}, {6, 2, 1000, 1}}));
  ASSERT_EQ(id_of_input(*lim), (std::vector<EdgeId> {0, 1, 2, 3, 4, 5, 6, 7, 8, 9}));
  const std::vector<EdgeId> p0 = {8, 0, 1, 5, 6, 2, 4};
  const std::vector<EdgeId> p1 = {9, 2, 3, 0, 1, 5, 7};
  const Amounts lcap {1, 1000, 2, 1000, 1000, 1, 1000, 1000, 1000, 1000};
  const auto limit = reference({{1, p0, {0, 0, 0}, {}}, {2, p1, {1, 0, 0}, {}}}, lcap);
  EXPECT_EQ(limit.error, kRoundLimit);
  EXPECT_EQ(limit.rounds, 192);
  const auto lb = build(*lim, {{{6}, {1}, {{0, p0, 1}}, {}}, {{6}, {2}, {{0, p1, 1}}, {}}}, 7);
  EXPECT_EQ(error_of([&] { (void)joint(*lim, lb, {{1}, {2}}); }), kRoundLimit);
#if NETGRAPH_CORE_PRIVATE_HOOKS
  EXPECT_EQ(detail::thread_lossy_counters().rounds, 192);
#endif
}

// ---------------------------------------------------------------- hand-derived cases

// The demand level above tunnels: A has two unit tunnels on link 0 (capacity 1), B one tunnel; offers 2 and 2.
// e0: H(1, (2, 2)) = (1, 0) by demand (A first), within A H(1, (1, 1)) = (1, 0): A's tunnel 0 gets
// the unit (flat over tunnels H(1, (1, 1, 2)) would give it to B).
TEST(JointForwarding, DemandLevelAboveTunnels) {
  const auto g = shared(make_graph(2, {{0, 1, 1, 1}}));
  const auto b = build(*g, {{{0}, {2}, {{0, {0}, 1}, {0, {0}, 1}}, {}}, {{0}, {2}, {{0, {0}, 1}}, {}}}, 1);
  const auto got = joint(*g, b, {{2}, {2}});
  EXPECT_EQ(got[0].entries[0].delivered, 1);
  EXPECT_EQ(got[0].entries[1].delivered, 0);
  EXPECT_EQ(got[1].entries[0].delivered, 0);
  EXPECT_EQ(got[1].entries[0].edge_drops, (Pairs {{0, 2}}));
  // The demand order is the given order: B first takes the unit.
  const Amounts two {2};
  const std::vector<detail::ForwardingDemand> swapped {{b.views[1].get(), two}, {b.views[0].get(), two}};
  const auto other = detail::lossy(swapped, capacities(*g));
  EXPECT_EQ(other[0].entries[0].delivered, 1);
  EXPECT_EQ(other[1].delivered, 0);
}

// The demand level above sources: A's heads 1 and 2 (one unit each) against B's head 5 (2 units) on link 2
// (capacity 1): H(1, (2, 2)) = (1, 0), within A H(1, (1, 1)) = (1, 0): head 1 delivers.
TEST(JointForwarding, DemandLevelAboveSources) {
  const auto g = shared(make_graph(6, {{1, 3, 100, 1}, {2, 3, 100, 1}, {3, 4, 1, 1}, {5, 3, 100, 1}}));
  ASSERT_EQ(id_of_input(*g), (std::vector<EdgeId> {0, 1, 2, 3}));
  const auto b = build(*g, {{{1, 2}, {1, 1}, {{0, {0, 2}, 1}, {1, {1, 2}, 1}}, {}}, {{5}, {2}, {{0, {3, 2}, 1}}, {}}}, 4);
  const auto got = joint(*g, b, {{1, 1}, {2}});
  EXPECT_EQ(got[0].entries[0].delivered, 1);
  EXPECT_EQ(got[0].entries[1].delivered, 0);
  EXPECT_EQ(got[1].delivered, 0);
}

// Link sums beyond 2^63: w = 2^62 - 1 per demand on one link; two, three and five
// demands: (5, 5) of 10, (4, 3, 3), (2, 2, 2, 2, 2) and with capacity 11 (3, 2, 2, 2, 2) (the
// extra unit to the first demand). Each drops the rest.
TEST(JointForwarding, LinkSumsBeyond64Bits) {
  const Amount w = kL62 - 1;
  for (const auto& [n, cap, want] : std::vector<std::tuple<int, Amount, Amounts>> {
           {2, 10, {5, 5}}, {3, 10, {4, 3, 3}}, {5, 10, {2, 2, 2, 2, 2}}, {5, 11, {3, 2, 2, 2, 2}}}) {
    const auto g = shared(make_graph(2, {{0, 1, cap, 1}}));
    std::vector<TunnelDemand> demands(static_cast<std::size_t>(n), TunnelDemand {{0}, {w}, {{0, {0}, 1}}, {}});
    const auto b = build(*g, demands, 1);
    const auto got = joint(*g, b, std::vector<Amounts>(static_cast<std::size_t>(n), Amounts {w}));
    for (int d = 0; d < n; ++d) {
      EXPECT_EQ(got[static_cast<std::size_t>(d)].delivered, want[static_cast<std::size_t>(d)]) << n << " demands, demand " << d;
      EXPECT_EQ(got[static_cast<std::size_t>(d)].entries[0].edge_drops, (Pairs {{0, w - want[static_cast<std::size_t>(d)]}}));
    }
  }
}

// ---------------------------------------------------------------- pooled demands

// A ShortestPaths program of `sources` to `dst` (pooled: one split context), prepared on a fresh
// ledger of g.
std::shared_ptr<const ForwardingProgram> pooled_program(const std::shared_ptr<const Graph>& g, const std::vector<NodeId>& sources,
                                                        NodeId dst, RouteAllocation allocation = RouteAllocation::Equal) {
  FlowLedger ledger(g);
  FlowPolicy p(ledger, sources, std::vector<NodeId> {dst}, FlowPolicyConfig {.routing = shortest_paths(allocation), .admission = Admission::Lossy});
  (void)p.replace(0);
  return p.program();
}

bool same_flow(const ProgramFlow& a, const ProgramFlow& b, const std::string& where) {
  EXPECT_EQ(a.admitted, b.admitted) << where;
  EXPECT_EQ(a.offered, b.offered) << where;
  EXPECT_EQ(a.source_drops, b.source_drops) << where;
  EXPECT_EQ(a.delivered, b.delivered) << where;
  EXPECT_EQ(a.offered_by_source, b.offered_by_source) << where;
  EXPECT_EQ(a.source_drops_by_source, b.source_drops_by_source) << where;
  EXPECT_EQ(a.entries.size(), b.entries.size()) << where;
  for (std::size_t e = 0; e < std::min(a.entries.size(), b.entries.size()); ++e) {
    EXPECT_EQ(a.entries[e].assigned, b.entries[e].assigned) << where << " entry " << e;
    EXPECT_EQ(a.entries[e].delivered, b.entries[e].delivered) << where << " entry " << e;
    EXPECT_EQ(a.entries[e].source_drops, b.entries[e].source_drops) << where << " entry " << e;
    EXPECT_EQ(a.entries[e].edge_flows, b.entries[e].edge_flows) << where << " entry " << e;
    EXPECT_EQ(a.entries[e].edge_drops, b.entries[e].edge_drops) << where << " entry " << e;
  }
  return !::testing::Test::HasFailure();
}

// One pooled demand: ShortestPaths EQUAL and WEIGHTED programs and from_dag programs of
// three or four sources on random DAGs, with a random mask and random black holes, forwarded alone
// on the joint evaluator, equal the pooled kernel's forwarding field by field (per head: injection,
// source loss, delivery, carried and dropped per edge).
TEST(JointForwarding, OnePooledDemandEqualsTheKernel) {
  std::mt19937 rng(5101);
  int checked = 0;
  int attributed = 0;  // trials where some head both delivers and drops
  int black = 0;
  for (int trial = 0; trial < 400; ++trial) {
    const std::int32_t n = 7;
    const NodeId dst = n - 1;
    std::vector<InEdge> rows;
    for (NodeId u = 0; u < n; ++u) {
      for (NodeId w = u + 1; w < n; ++w) {
        if (rng() % 3 != 0) rows.push_back({u, w, static_cast<Amount>(rng() % 6), 1 + static_cast<Cost>(rng() % 2)});
      }
    }
    const auto g = shared(make_graph(n, rows));
    std::vector<NodeId> sources = {0, 1, 2};
    if (rng() % 2 == 0) sources.push_back(3);
    std::shared_ptr<const ForwardingProgram> prog;
    const int kind = static_cast<int>(rng() % 3);
    if (kind < 2) {
      prog = pooled_program(g, sources, dst, kind == 0 ? RouteAllocation::Equal : RouteAllocation::Weighted);
    } else {
      try {
        prog = ForwardingProgram::from_dag(*g, sources, std::vector<NodeId> {dst}, all_edges(*g), RouteOptions {});
      } catch (const std::invalid_argument&) {
        continue;
      }
    }
    if (!prog || !detail::ProgramAccess::data(*prog).pools_source_inflows) continue;
    std::shared_ptr<const Eligibility> mask;
    if (rng() % 3 == 0) {
      const std::vector<EdgeId> masked {static_cast<EdgeId>(rng() % g->num_edges())};
      mask = std::make_shared<const Eligibility>(g->eligibility({}, masked));
    }
    std::shared_ptr<const Eligibility> dark;
    if (rng() % 2 == 0) {
      std::vector<EdgeId> es {static_cast<EdgeId>(rng() % g->num_edges())};
      if (rng() % 2 == 0) es.push_back(static_cast<EdgeId>(rng() % g->num_edges()));
      std::vector<NodeId> nodes;
      if (rng() % 4 == 0) nodes.push_back(sources[rng() % sources.size()]);
      dark = std::make_shared<const Eligibility>(g->eligibility(nodes, es));
    }
    const auto view = detail::make_forwarding_view(*g, prog, mask.get(), vec(prog->destinations()), dark.get());
    if (!view.pooled) continue;
    Amounts inject(view.num_sources());
    for (auto& x : inject) x = static_cast<Amount>(rng() % 9);
    const auto free = capacities(*g);
    const auto kernel = detail::lossy(view, free, inject);
    const detail::ForwardingDemand one {&view, inject};
    const auto got = detail::lossy(std::span<const detail::ForwardingDemand>(&one, 1), free);
    ASSERT_EQ(got.size(), 1U);
    if (!same_flow(got.front(), kernel, "trial " + std::to_string(trial))) return;
    for (const auto& f : kernel.entries) attributed += f.delivered > 0 && !f.edge_drops.empty() ? 1 : 0;
    black += dark ? 1 : 0;
    ++checked;
  }
  RecordProperty("checked", checked);
  EXPECT_GT(checked, 200);
  EXPECT_GT(attributed, 50);
  EXPECT_GT(black, 80);
}

// An independent sweep of pooled and tunnel demands on one view: from zero carried, every sweep
// recomputes every arrival from the previous carried amounts (a tunnel's hop from its previous hop;
// a pooled router splits its heads' inflows' total by round robin or highest averages over its
// members' frozen weights and attributes each member's amount to the heads by highest averages over
// their remaining inflows, member by member) and every link's carried amounts by the loss rule over its
// non-dark arrivals (labels demand, head, entry); a dark arrival carries 0. It stops when a sweep
// changes nothing. When the dependencies are acyclic it reaches A+'s answer (depth by depth). It
// reads the plans' structure only; the splits and the rule are its own.
struct SweepPooled {
  int demand {0};
  const PreparedRoute* plan {nullptr};
  std::vector<std::int32_t> root_positions;  // per head (-1: no member)
  Amounts shares;                           // per head
  std::set<EdgeId> dark;
};
struct SweepTunnel {
  int demand {0};
  int head {0};
  int index {0};
  Amount share {0};
  std::vector<EdgeId> route;
  std::set<EdgeId> dark;
};
struct SweepResult {
  bool settled {false};
  std::vector<std::vector<std::vector<Amount>>> pooled_carried;  // per pooled demand, member, head
  std::vector<std::vector<std::vector<Amount>>> pooled_arrival;
  std::vector<Amounts> tunnel_carried;
  std::vector<Amounts> tunnel_arrival;
  std::vector<Amounts> pooled_delivered;  // per pooled demand, head
};

SweepResult sweep(const std::vector<SweepPooled>& pooled, const std::vector<SweepTunnel>& tunnels, const Amounts& cap) {
  SweepResult r;
  r.pooled_carried.resize(pooled.size());
  r.pooled_arrival.resize(pooled.size());
  for (std::size_t q = 0; q < pooled.size(); ++q) {
    const auto M = pooled[q].plan->forwarding_edges.size();
    r.pooled_carried[q].assign(M, Amounts(pooled[q].shares.size(), 0));
    r.pooled_arrival[q] = r.pooled_carried[q];
  }
  r.tunnel_carried.resize(tunnels.size());
  r.tunnel_arrival.resize(tunnels.size());
  for (std::size_t t = 0; t < tunnels.size(); ++t) {
    r.tunnel_carried[t].assign(tunnels[t].route.size(), 0);
    r.tunnel_arrival[t] = r.tunnel_carried[t];
  }
  for (int round = 0; round < 400; ++round) {
    // Arrivals from the previous carried amounts.
    for (std::size_t q = 0; q < pooled.size(); ++q) {
      const auto& d = pooled[q];
      const auto& plan = *d.plan;
      const std::size_t P = plan.topological_order.size();
      const std::size_t k = d.shares.size();
      std::vector<Amounts> in(P, Amounts(k, 0));
      for (std::size_t h = 0; h < k; ++h) {
        if (d.root_positions[h] >= 0) in[static_cast<std::size_t>(d.root_positions[h])][h] += d.shares[h];
      }
      for (std::size_t m = 0; m < plan.forwarding_edges.size(); ++m) {
        for (std::size_t h = 0; h < k; ++h) in[static_cast<std::size_t>(plan.edge_destination_positions[m])][h] += r.pooled_carried[q][m][h];
      }
      for (std::size_t pos = 0; pos < P; ++pos) {
        if (plan.is_destination[pos] != 0) continue;
        const auto b = static_cast<std::size_t>(plan.outgoing_begin[pos]);
        const auto e = static_cast<std::size_t>(plan.outgoing_end[pos]);
        Amount total = 0;
        for (Amount x : in[pos]) total += x;
        Amounts part;
        if (plan.allocation == RouteAllocation::Weighted) {
          part = ref::highest_averages(total, Amounts(plan.weights.begin() + static_cast<std::ptrdiff_t>(b), plan.weights.begin() + static_cast<std::ptrdiff_t>(e)));
        } else {
          part = ref::round_robin(total, e - b);
        }
        if (part.empty()) part.assign(e - b, 0);
        Amounts remaining = in[pos];
        for (std::size_t m = b; m < e; ++m) {
          auto row = ref::highest_averages(part[m - b], remaining);
          if (row.empty()) row.assign(k, 0);
          for (std::size_t h = 0; h < k; ++h) {
            remaining[h] -= row[h];
            r.pooled_arrival[q][m][h] = row[h];
          }
        }
      }
    }
    for (std::size_t t = 0; t < tunnels.size(); ++t) {
      for (std::size_t j = 0; j < tunnels[t].route.size(); ++j) r.tunnel_arrival[t][j] = j == 0 ? tunnels[t].share : r.tunnel_carried[t][j - 1];
    }
    // Every link's rule over its non-dark arrivals, in (demand, head, entry) order.
    struct Member {
      std::tuple<int, int, int> label;
      Amount* carried;
      Amount arrival;
    };
    std::map<EdgeId, std::vector<Member>> links;
    auto next_pc = r.pooled_carried;
    auto next_tc = r.tunnel_carried;
    for (std::size_t q = 0; q < pooled.size(); ++q) {
      const auto& plan = *pooled[q].plan;
      for (std::size_t m = 0; m < plan.forwarding_edges.size(); ++m) {
        const EdgeId e = plan.forwarding_edges[m];
        for (std::size_t h = 0; h < pooled[q].shares.size(); ++h) {
          if (pooled[q].dark.count(e) != 0) {
            next_pc[q][m][h] = 0;
            continue;
          }
          links[e].push_back({{pooled[q].demand, static_cast<int>(h), 0}, &next_pc[q][m][h], r.pooled_arrival[q][m][h]});
        }
      }
    }
    for (std::size_t t = 0; t < tunnels.size(); ++t) {
      for (std::size_t j = 0; j < tunnels[t].route.size(); ++j) {
        const EdgeId e = tunnels[t].route[j];
        if (tunnels[t].dark.count(e) != 0) {
          next_tc[t][j] = 0;
          continue;
        }
        links[e].push_back({{tunnels[t].demand, tunnels[t].head, tunnels[t].index}, &next_tc[t][j], r.tunnel_arrival[t][j]});
      }
    }
    for (auto& [e, members] : links) {
      std::stable_sort(members.begin(), members.end(), [](const Member& a, const Member& b) { return a.label < b.label; });
      Amounts w;
      std::vector<std::tuple<int, int, int>> labels;
      Amount sum = 0;
      for (const auto& x : members) {
        w.push_back(x.arrival);
        labels.push_back(x.label);
        sum += x.arrival;
      }
      const auto c = hierarchical(std::min(cap[static_cast<std::size_t>(e)], sum), w, labels);
      for (std::size_t i = 0; i < members.size(); ++i) *members[i].carried = c[i];
    }
    const bool same = next_pc == r.pooled_carried && next_tc == r.tunnel_carried && round > 0;
    r.pooled_carried = std::move(next_pc);
    r.tunnel_carried = std::move(next_tc);
    if (same) {
      r.settled = true;
      break;
    }
  }
  r.pooled_delivered.resize(pooled.size());
  for (std::size_t q = 0; q < pooled.size(); ++q) {
    const auto& plan = *pooled[q].plan;
    r.pooled_delivered[q].assign(pooled[q].shares.size(), 0);
    for (std::size_t m = 0; m < plan.forwarding_edges.size(); ++m) {
      if (plan.is_destination[static_cast<std::size_t>(plan.edge_destination_positions[m])] == 0) continue;
      for (std::size_t h = 0; h < pooled[q].shares.size(); ++h) r.pooled_delivered[q][h] += r.pooled_carried[q][m][h];
    }
  }
  return r;
}

// Random pooled demands (ShortestPaths EQUAL or WEIGHTED, two or three sources) beside random
// tunnel demands on random DAGs, each demand with its own black holes, in random creation orders,
// against the sweep wherever the evaluation needs no rounds: every head's and tunnel's arrivals,
// carried amounts and deliveries.
TEST(JointForwarding, PooledAndTunnelDemandsMatchASweep) {
  std::mt19937 rng(77);
  int checked = 0;
  int shared_pooled = 0;  // trials where a pooled member carries another demand too
  int attributed = 0;     // trials where a pooled head drops while another delivers
  for (int trial = 0; trial < 600; ++trial) {
    const std::int32_t n = 7;
    const NodeId dst = n - 1;
    std::vector<InEdge> rows;
    for (NodeId u = 0; u < n; ++u) {
      for (NodeId w = u + 1; w < n; ++w) {
        if (rng() % 3 != 0) rows.push_back({u, w, static_cast<Amount>(rng() % 5), 1 + static_cast<Cost>(rng() % 2)});
      }
    }
    const auto g = shared(make_graph(n, rows));
    const Amounts cap(g->edge_capacity_units().begin(), g->edge_capacity_units().end());
    // Each demand: pooled (kind 0) or tunnels (kind 1).
    const int D = 2 + static_cast<int>(rng() % 3);
    std::vector<std::shared_ptr<const ForwardingProgram>> programs;
    std::vector<std::shared_ptr<const Eligibility>> black;
    std::vector<std::unique_ptr<detail::ForwardingView>> views;
    std::vector<Amounts> injections;
    std::vector<SweepPooled> sp;
    std::vector<SweepTunnel> st;
    std::vector<int> kind;
    bool skip = false;
    for (int d = 0; d < D && !skip; ++d) {
      std::vector<EdgeId> dark_edges;
      if (rng() % 3 == 0) dark_edges.push_back(static_cast<EdgeId>(rng() % g->num_edges()));
      black.push_back(dark_edges.empty() ? nullptr : std::make_shared<const Eligibility>(g->eligibility({}, dark_edges)));
      const std::set<EdgeId> dark_set(dark_edges.begin(), dark_edges.end());
      if (rng() % 2 == 0) {
        std::vector<NodeId> sources = {0, 1};
        if (rng() % 2 == 0) sources.push_back(2);
        programs.push_back(pooled_program(g, sources, dst, rng() % 2 == 0 ? RouteAllocation::Equal : RouteAllocation::Weighted));
        views.push_back(std::make_unique<detail::ForwardingView>(
            detail::make_forwarding_view(*g, programs.back(), nullptr, std::vector<NodeId> {dst}, black.back().get())));
        if (!views.back()->pooled) {
          skip = true;
          break;
        }
        Amounts inj(views.back()->num_sources());
        for (auto& x : inj) x = static_cast<Amount>(rng() % 9);
        SweepPooled q;
        q.demand = d;
        q.plan = views.back()->prepared_route(0).get();
        for (std::size_t h = 0; h < inj.size(); ++h) {
          const auto e = static_cast<std::size_t>(views.back()->entries_for_source(h).front());
          const auto rr = views.back()->context_source_index(e);
          q.root_positions.push_back(rr < 0 ? -1 : q.plan->source_positions[static_cast<std::size_t>(rr)]);
        }
        q.shares = inj;
        q.dark = dark_set;
        sp.push_back(q);
        injections.push_back(inj);
        kind.push_back(0);
        continue;
      }
      // One or two tunnels from a random head to dst along random paths of the DAG.
      const NodeId head = static_cast<NodeId>(rng() % 3);
      std::vector<RouteDefinition> routes;
      const int nt = 1 + static_cast<int>(rng() % 2);
      Amounts weights;
      std::vector<std::vector<EdgeId>> paths;
      for (int t = 0; t < nt; ++t) {
        std::vector<EdgeId> path;
        NodeId u = head;
        while (u != dst) {
          std::vector<EdgeId> out;
          for (EdgeId e = 0; e < g->num_edges(); ++e) {
            if (g->edge_sources()[static_cast<std::size_t>(e)] == u) out.push_back(e);
          }
          if (out.empty()) break;
          const EdgeId e = out[rng() % out.size()];
          path.push_back(e);
          u = g->edge_destinations()[static_cast<std::size_t>(e)];
        }
        if (u != dst || path.empty()) continue;
        std::vector<std::vector<EdgeId>> hops;
        for (EdgeId e : path) hops.push_back({e});
        weights.push_back(1 + static_cast<Amount>(rng() % 3));
        routes.push_back(RouteDefinition {head, std::make_shared<const Path>(*g, hops), weights.back()});
        paths.push_back(path);
      }
      if (routes.empty()) {
        skip = true;
        break;
      }
      programs.push_back(ForwardingProgram::from_routes(*g, std::vector<NodeId> {head}, std::vector<NodeId> {dst}, routes));
      views.push_back(std::make_unique<detail::ForwardingView>(
          detail::make_forwarding_view(*g, programs.back(), nullptr, std::vector<NodeId> {dst}, black.back().get())));
      const Amount x = static_cast<Amount>(rng() % 9);
      injections.push_back({x});
      const auto split = ref::highest_averages(x, weights);
      for (std::size_t t = 0; t < paths.size(); ++t) st.push_back({d, 0, static_cast<int>(t), split[t], paths[t], dark_set});
      kind.push_back(1);
    }
    if (skip) continue;
    std::vector<detail::ForwardingDemand> demands;
    for (std::size_t d = 0; d < views.size(); ++d) demands.push_back({views[d].get(), injections[d]});
    const auto got = detail::lossy(demands, capacities(*g));
#if NETGRAPH_CORE_PRIVATE_HOOKS
    if (detail::thread_lossy_counters().groups != 0) continue;
#endif
    const auto want = sweep(sp, st, cap);
    ASSERT_TRUE(want.settled) << "trial " << trial;
    const std::string where = "trial " + std::to_string(trial);
    std::size_t qi = 0;
    std::size_t ti = 0;
    std::map<EdgeId, std::set<int>> users;
    for (std::size_t d = 0; d < views.size(); ++d) {
      if (kind[d] == 0) {
        const auto& plan = *sp[qi].plan;
        bool drops = false;
        bool delivers = false;
        for (std::size_t h = 0; h < sp[qi].shares.size(); ++h) {
          const auto& f = got[d].entries[static_cast<std::size_t>(views[d]->entries_for_source(h).front())];
          Pairs flows;
          Pairs dr;
          for (std::size_t m = 0; m < plan.forwarding_edges.size(); ++m) {
            const Amount c = want.pooled_carried[qi][m][h];
            const Amount a = want.pooled_arrival[qi][m][h];
            if (c != 0) {
              flows.emplace_back(plan.forwarding_edges[m], c);
              users[plan.forwarding_edges[m]].insert(static_cast<int>(d));
            }
            if (a != c) dr.emplace_back(plan.forwarding_edges[m], a - c);
          }
          EXPECT_EQ(f.edge_flows, flows) << where << " demand " << d << " head " << h;
          EXPECT_EQ(f.edge_drops, dr) << where << " demand " << d << " head " << h;
          EXPECT_EQ(f.delivered, want.pooled_delivered[qi][h]) << where << " demand " << d << " head " << h;
          drops = drops || !dr.empty();
          delivers = delivers || f.delivered > 0;
        }
        attributed += drops && delivers ? 1 : 0;
        ++qi;
        continue;
      }
      for (std::size_t t = 0; ti < st.size() && st[ti].demand == static_cast<int>(d); ++t, ++ti) {
        const auto& f = got[d].entries[t];
        Pairs flows;
        Pairs dr;
        for (std::size_t j = 0; j < st[ti].route.size(); ++j) {
          if (want.tunnel_carried[ti][j] != 0) {
            flows.emplace_back(st[ti].route[j], want.tunnel_carried[ti][j]);
            users[st[ti].route[j]].insert(static_cast<int>(d));
          }
          if (want.tunnel_arrival[ti][j] != want.tunnel_carried[ti][j]) dr.emplace_back(st[ti].route[j], want.tunnel_arrival[ti][j] - want.tunnel_carried[ti][j]);
        }
        std::sort(flows.begin(), flows.end());
        std::sort(dr.begin(), dr.end());
        EXPECT_EQ(f.edge_flows, flows) << where << " demand " << d << " tunnel " << t;
        EXPECT_EQ(f.edge_drops, dr) << where << " demand " << d << " tunnel " << t;
        EXPECT_EQ(f.delivered, want.tunnel_carried[ti].empty() ? st[ti].share : want.tunnel_carried[ti].back()) << where;
      }
    }
    for (const auto& q : sp) {
      for (EdgeId e : q.plan->forwarding_edges) {
        if (users[e].size() > 1) {
          ++shared_pooled;
          break;
        }
      }
    }
    ++checked;
    if (::testing::Test::HasFailure()) return;
  }
  RecordProperty("checked", checked);
  RecordProperty("shared_pooled", shared_pooled);
  RecordProperty("attributed", attributed);
  EXPECT_GT(checked, 250);
  EXPECT_GT(shared_pooled, 100);
  EXPECT_GT(attributed, 50);
}

// A pooled demand shares a link with a tunnel. G7: e0 0->2, e1 1->2, e2 2->3 (1), e3 2->3, e4
// 4->2 (all others wide); A pooled sources {0, 1} -> 3 (RR(2, 2) = (1, 1)), B tunnel [e4, e2] (2).
// Node 2 splits A's 2 as (1, 1) over e2, e3; e2: A 1, B 2, capacity 1: by demand H(1, (1, 2)) =
// (0, 1); e3 carries A's 1. Attribution at node 2: e2's unit to source 0 (H(1, (1, 1))), e3's to
// source 1.
TEST(JointForwarding, PooledDemandSharesALink) {
  const auto g = shared(make_graph(5, {{0, 2, 100, 1}, {1, 2, 100, 1}, {2, 3, 1, 1}, {2, 3, 100, 1}, {4, 2, 100, 1}}));
  ASSERT_EQ(id_of_input(*g), (std::vector<EdgeId> {0, 1, 2, 3, 4}));
  const auto a = pooled_program(g, {0, 1}, 3);
  ASSERT_TRUE(detail::ProgramAccess::data(*a).pools_source_inflows);
  const auto va = detail::make_forwarding_view(*g, a, nullptr, std::vector<NodeId> {3}, nullptr);
  const auto bp = build(*g, {{{4}, {2}, {{0, {4, 2}, 1}}, {}}}, 3);
  const Amounts ia {1, 1};
  const Amounts ib {2};
  const std::vector<detail::ForwardingDemand> demands {{&va, ia}, {bp.views[0].get(), ib}};
  const auto got = detail::lossy(demands, capacities(*g));
  EXPECT_EQ(got[0].delivered, 1);
  EXPECT_EQ(got[0].entries[0].delivered, 0);
  EXPECT_EQ(got[0].entries[0].edge_flows, (Pairs {{0, 1}}));
  EXPECT_EQ(got[0].entries[0].edge_drops, (Pairs {{2, 1}}));
  EXPECT_EQ(got[0].entries[1].delivered, 1);
  EXPECT_EQ(got[0].entries[1].edge_flows, (Pairs {{1, 1}, {3, 1}}));
  EXPECT_TRUE(got[0].entries[1].edge_drops.empty());
  EXPECT_EQ(got[1].delivered, 1);
  EXPECT_EQ(got[1].entries[0].edge_drops, (Pairs {{2, 1}}));
}

// Every root is an injection. Chain 0->1->2->3 of capacity 2, a pooled demand of sources {0,
// 1, 2} with shares (2, 2, 2): arrivals (2, 4, 4), carries (2, 2, 2); neither (e0, e1) nor (e1, e2)
// is a transparent pair (nodes 1 and 2 inject). Per source: e1 carries H(2, (2, 2)) = (1, 1), e2
// H(2, (1, 1, 2)) = (1, 0, 1): delivered (1, 0, 1). Alone and beside a demand on a disjoint link,
// equal to the kernel.
TEST(JointForwarding, PooledRootsBlockSerialPairs) {
  const auto g = shared(make_graph(6, {{0, 1, 2, 1}, {1, 2, 2, 1}, {2, 3, 2, 1}, {4, 5, 1, 1}}));
  const auto a = pooled_program(g, {0, 1, 2}, 3);
  const auto va = detail::make_forwarding_view(*g, a, nullptr, std::vector<NodeId> {3}, nullptr);
  const Amounts ia {2, 2, 2};
  const auto kernel = detail::lossy(va, capacities(*g), ia);
  EXPECT_EQ(kernel.delivered, 2);
  const detail::ForwardingDemand alone {&va, ia};
  const auto one = detail::lossy(std::span<const detail::ForwardingDemand>(&alone, 1), capacities(*g));
  same_flow(one.front(), kernel, "alone");
  EXPECT_EQ(one.front().entries[0].delivered, 1);
  EXPECT_EQ(one.front().entries[1].delivered, 0);
  EXPECT_EQ(one.front().entries[2].delivered, 1);
  const auto other = build(*g, {{{4}, {3}, {{0, {3}, 1}}, {}}}, 5);
  const Amounts io {3};
  const std::vector<detail::ForwardingDemand> both {{&va, ia}, {other.views[0].get(), io}};
  const auto two = detail::lossy(both, capacities(*g));
  same_flow(two.front(), kernel, "beside another demand");
  EXPECT_EQ(two[1].delivered, 1);
}

// A head whose share is 0 injects nothing and carries nothing (its occurrences are dead), and a dark
// member's traffic goes no further. Chain 0->1->2->3, a pooled demand of sources {0, 1, 2}, four
// routers. Capacities (1, 1, 1), shares (2, 0, 0): only e0 can drop and e1, e2 continue it (nodes 1
// and 2 inject nothing), so one barrier: five events (hooks: the event limit 5 is met). Capacities
// (2, 1, 2), shares (0, 2, 0): e1's first occurrence is source 0's, which carries nothing, yet e1
// waits for source 1's arrival: delivered 1, source 1 drops 1 on e1. Capacities (1, 1, 1), shares
// (2, 0, 0), e0 black-holed: source 0's traffic dies on e0, so e1 and e2 can carry nothing: no
// barrier, four events. Each equals the kernel, and the program of its one live source alone
// (the same plan without the heads that carry nothing): that source's amounts and (hooks) the
// links that wait, the groups and the dependency arcs.
TEST(JointForwarding, ZeroShareHeadsOfAPooledDemand) {
  struct Restore {
    ~Restore() {
#if NETGRAPH_CORE_PRIVATE_HOOKS
      detail::thread_forwarding_limits() = {};
#endif
    }
  } restore;
  struct Case {
    Amount v1;
    Amounts shares;
    bool black;
    std::int64_t events;
    Amounts delivered;
  };
  for (const auto& c : {Case {1, {2, 0, 0}, false, 5, {1, 0, 0}}, Case {1, {0, 2, 0}, false, 5, {0, 1, 0}},
                        Case {1, {2, 0, 0}, true, 4, {0, 0, 0}}}) {
    const Amount v0 = c.shares[1] != 0 ? 2 : 1;
    const auto g = shared(make_graph(4, {{0, 1, v0, 1}, {1, 2, c.v1, 1}, {2, 3, v0, 1}}));
    const auto a = pooled_program(g, {0, 1, 2}, 3);
    const auto black = c.black ? std::make_shared<const Eligibility>(g->eligibility({}, std::vector<EdgeId> {0})) : nullptr;
    const auto va = detail::make_forwarding_view(*g, a, nullptr, std::vector<NodeId> {3}, black.get());
    const auto kernel = detail::lossy(va, capacities(*g), c.shares);
#if NETGRAPH_CORE_PRIVATE_HOOKS
    detail::thread_forwarding_limits()[static_cast<std::size_t>(detail::ForwardingCounter::Event)] = c.events;
#endif
    const detail::ForwardingDemand alone {&va, c.shares};
    const auto one = detail::lossy(std::span<const detail::ForwardingDemand>(&alone, 1), capacities(*g));
#if NETGRAPH_CORE_PRIVATE_HOOKS
    detail::thread_forwarding_limits() = {};
#endif
    const std::string where = "shares " + std::to_string(c.shares[0]) + std::to_string(c.shares[1]) + (c.black ? " black" : "");
    same_flow(one.front(), kernel, where);
    for (std::size_t h = 0; h < 3; ++h) EXPECT_EQ(one.front().entries[h].delivered, c.delivered[h]) << where << " source " << h;
#if NETGRAPH_CORE_PRIVATE_HOOKS
    const auto structure = std::make_tuple(detail::thread_lossy_counters().groups, detail::thread_lossy_counters().barriers,
                                           detail::thread_lossy_counters().dependencies);
#endif
    const std::size_t live = c.shares[0] != 0 ? 0 : 1;
    const auto b = pooled_program(g, {static_cast<NodeId>(live)}, 3);
    const auto vb = detail::make_forwarding_view(*g, b, nullptr, std::vector<NodeId> {3}, black.get());
    const Amounts ib {c.shares[live]};
    const detail::ForwardingDemand only {&vb, ib};
    const auto lone = detail::lossy(std::span<const detail::ForwardingDemand>(&only, 1), capacities(*g));
    EXPECT_EQ(one.front().entries[live].delivered, lone.front().entries[0].delivered) << where;
    EXPECT_EQ(one.front().entries[live].edge_flows, lone.front().entries[0].edge_flows) << where;
    EXPECT_EQ(one.front().entries[live].edge_drops, lone.front().entries[0].edge_drops) << where;
#if NETGRAPH_CORE_PRIVATE_HOOKS
    EXPECT_EQ(structure, std::make_tuple(detail::thread_lossy_counters().groups, detail::thread_lossy_counters().barriers,
                                         detail::thread_lossy_counters().dependencies))
        << where;
#endif
  }
}

// Serial normalization reads only what can carry something: a router's merges and branches are
// those of its live members. Each case against the same program without its dead member, pooled
// (from_dag) and not (the DAG route beside a tunnel of weight 0), on links of capacity 0 that drop
// everything (so every link with traffic can drop, and a link that passes its predecessor's traffic
// through is no barrier): the delivery and (hooks) the links that wait, the groups and the arcs.
// - Dead branch: e0 0->1, e1 1->2, e2 1->3 (capacity 9); one unit, round robin sends it to e1, so
//   e2 carries nothing and e0, e1 stay a pair: one barrier, e0.
// - Dead merge in one split route: e0 0->1 and e1 0->2 (capacity 9), e2 1->3, e3 2->3, e4 3->4; one
//   unit takes e0 and e2, so e1 and e3 carry nothing and e2, e4 stay a pair: one barrier, e2.
// - Dead merge of a pooled source: e0 0->2 (capacity 9), e1 1->2, e2 2->3; sources 0 and 1. Source 0's
//   share is 0, or its 1 dies on e0, its own black hole (its router arc stays: one more arc); e1, e2
//   stay a pair: one barrier, e1.
TEST(JointForwarding, DeadBranchesAndMergesKeepPairs) {
  struct Seen {
    Amount delivered {0};
    std::vector<EdgeId> barriers;
    std::int64_t groups {-1};
    std::int64_t dependencies {-1};
  };
  auto run = [](const std::shared_ptr<const Graph>& g, const std::shared_ptr<const ForwardingProgram>& prog, const Amounts& inject,
                std::vector<NodeId> dests, const std::vector<EdgeId>& black) {
    const auto holes = black.empty() ? nullptr : std::make_shared<const Eligibility>(g->eligibility({}, black));
    const auto view = detail::make_forwarding_view(*g, prog, nullptr, dests, holes.get());
    const detail::ForwardingDemand demand {&view, inject};
    const auto out = detail::lossy(std::span<const detail::ForwardingDemand>(&demand, 1), capacities(*g));
    Seen seen {out.front().delivered, {}, -1, -1};
#if NETGRAPH_CORE_PRIVATE_HOOKS
    seen.barriers = detail::thread_lossy_counters().barriers;
    seen.groups = detail::thread_lossy_counters().groups;
    seen.dependencies = detail::thread_lossy_counters().dependencies;
#endif
    return seen;
  };
  // A program of one DAG route from `source`, pooled or beside a weight-0 tunnel over `path`.
  auto program = [](const std::shared_ptr<const Graph>& g, NodeId source, const std::vector<NodeId>& dests, const std::vector<EdgeId>& members,
                    bool pooled, const std::vector<EdgeId>& path) {
    const NodeId s[1] = {source};
    if (pooled) return ForwardingProgram::from_dag(*g, s, dests, dag_of(*g, members), RouteOptions {});
    std::vector<std::vector<EdgeId>> hops;
    for (EdgeId e : path) hops.push_back({e});
    const std::vector<RouteDefinition> defs {RouteDefinition {source, dag_of(*g, members), 1},
                                             RouteDefinition {source, std::make_shared<const Path>(*g, hops), 0}};
    return ForwardingProgram::from_routes(*g, s, dests, defs);
  };
  auto expect_same = [](const Seen& with, const Seen& without, std::int64_t extra_arcs, const std::string& where) {
    EXPECT_EQ(with.delivered, without.delivered) << where;
    EXPECT_EQ(with.barriers, without.barriers) << where;
    EXPECT_EQ(with.groups, without.groups) << where;
#if NETGRAPH_CORE_PRIVATE_HOOKS
    EXPECT_EQ(with.dependencies, without.dependencies + extra_arcs) << where;
#else
    (void)extra_arcs;
#endif
  };
  for (bool pooled : {true, false}) {
    const std::string kind = pooled ? "pooled" : "not pooled";
    {
      const auto g = shared(make_graph(4, {{0, 1, 0, 1}, {1, 2, 0, 1}, {1, 3, 9, 1}}));
      const auto id = id_of_input(*g);
      const auto with = run(g, program(g, 0, {2, 3}, {id[0], id[1], id[2]}, pooled, {id[0], id[1]}), {1}, {2, 3}, {});
      const auto without = run(g, program(g, 0, {2}, {id[0], id[1]}, pooled, {id[0], id[1]}), {1}, {2}, {});
      expect_same(with, without, 0, "dead branch, " + kind);
#if NETGRAPH_CORE_PRIVATE_HOOKS
      EXPECT_EQ(with.barriers, (std::vector<EdgeId> {id[0]})) << kind;
#endif
    }
    {
      const auto g = shared(make_graph(5, {{0, 1, 9, 1}, {0, 2, 9, 1}, {1, 3, 0, 1}, {2, 3, 0, 1}, {3, 4, 0, 1}}));
      const auto id = id_of_input(*g);
      const std::vector<EdgeId> path {id[0], id[2], id[4]};
      const auto with = run(g, program(g, 0, {4}, {id[0], id[1], id[2], id[3], id[4]}, pooled, path), {1}, {4}, {});
      const auto without = run(g, program(g, 0, {4}, path, pooled, path), {1}, {4}, {});
      expect_same(with, without, 0, "dead merge in one split route, " + kind);
#if NETGRAPH_CORE_PRIVATE_HOOKS
      EXPECT_EQ(with.barriers, (std::vector<EdgeId> {id[2]})) << kind;
#endif
    }
  }
  const auto g = shared(make_graph(4, {{0, 2, 9, 1}, {1, 2, 0, 1}, {2, 3, 0, 1}}));
  const auto id = id_of_input(*g);
  const NodeId both[2] = {0, 1};
  const NodeId live[1] = {1};
  const NodeId d[1] = {3};
  const auto two = ForwardingProgram::from_dag(*g, both, d, dag_of(*g, {id[0], id[1], id[2]}), RouteOptions {});
  const auto without = run(g, ForwardingProgram::from_dag(*g, live, d, dag_of(*g, {id[1], id[2]}), RouteOptions {}), {1}, {3}, {});
  expect_same(run(g, two, {0, 1}, {3}, {}), without, 0, "dead merge of a zero-share source");
  expect_same(run(g, two, {1, 1}, {3}, {id[0]}), without, 1, "dead merge of a black-holed source");
#if NETGRAPH_CORE_PRIVATE_HOOKS
  EXPECT_EQ(without.barriers, (std::vector<EdgeId> {id[1]}));
#endif
}

// Adding an element that can carry nothing changes nothing, on every path of the rule: random small
// graphs with cycles (six inner nodes, each with an exit to the destination),
// two to four demands that are tunnels, split routes (a DAG route beside a tunnel of weight 0) or
// pooled contexts (from_dag over two sources), each over one or two random walks of up to four inner
// links, and one added element:
// - a tunnel whose offer is 0;
// - a tunnel behind its own black hole (against the tunnel cut after its black hole);
// - a pooled source whose share is 0, its walk merging into its demand's DAG;
// - such a source behind its own black hole (its black-holed first link keeps one router arc);
// - a branch that round robin leaves empty (one unit at its root, the branch on the higher EdgeId),
//   its walk merging back or reaching the destination, pooled and not;
// - a source of share 0 inside a demand: a tunnel or a DAG route of a new head, or a pooled source
//   on the demand's DAG (its root inside a chain).
// One demand in a quarter of the draws; in half of those the element joins that demand (the
// scalar path).
// With and without the element: every demand's delivery and carried amounts per link, the drops of
// every demand whose own traffic is unchanged, the error, and (hooks) the rounds, the groups, the
// links that wait and the dependency arcs.
TEST(JointForwarding, DeadElementsChangeNothingOnEveryPath) {
  std::mt19937 rng(20261021);
  const NodeId inner = 6;
  const NodeId dst = inner;
  struct Part {
    std::shared_ptr<const ForwardingProgram> prog;
    std::shared_ptr<const Eligibility> black;
    Amounts inject;
  };
  struct Seen {
    std::string error;
    Amounts delivered;
    std::vector<std::map<EdgeId, Amount>> carried;
    std::vector<std::map<EdgeId, Amount>> drops;
    std::int64_t rounds {-1};
    std::int64_t groups {-1};
    std::vector<EdgeId> barriers;
    std::int64_t dependencies {-1};
  };
  std::array<int, 7> kinds {};
  int checked = 0;
  int skipped = 0;
  int cyclic = 0;
  int errors = 0;
  for (int trial = 0; trial < 8000; ++trial) {
    std::vector<InEdge> rows;
    for (NodeId u = 0; u < inner; ++u) {
      for (NodeId w = 0; w < inner; ++w) {
        if (u != w && rng() % 2 == 0) rows.push_back({u, w, static_cast<Amount>(rng() % 4), 1});
      }
      rows.push_back({u, dst, static_cast<Amount>(1 + rng() % 5), 1});
    }
    const auto g = shared(make_graph(inner + 1, rows));
    std::vector<std::vector<EdgeId>> outs(static_cast<std::size_t>(inner));
    std::vector<EdgeId> exit(static_cast<std::size_t>(inner), -1);
    for (EdgeId e = 0; e < g->num_edges(); ++e) {
      const auto u = g->edge_sources()[static_cast<std::size_t>(e)];
      const auto w = g->edge_destinations()[static_cast<std::size_t>(e)];
      if (w == dst) {
        exit[static_cast<std::size_t>(u)] = e;
      } else {
        outs[static_cast<std::size_t>(u)].push_back(e);
      }
    }
    auto head_of = [&](EdgeId e) { return g->edge_destinations()[static_cast<std::size_t>(e)]; };
    // A random simple walk from u (0 to 3 inner links, avoiding `avoid`), ending at the destination
    // or, when it can, at a node of `merge` (taken half of the time it is offered).
    auto walk = [&](NodeId u, std::set<NodeId> avoid, const std::set<NodeId>& merge) {
      std::vector<EdgeId> path;
      avoid.insert(u);
      const int len = static_cast<int>(rng() % 5);
      for (int k = 0; k <= len; ++k) {
        std::vector<EdgeId> into;
        std::vector<EdgeId> free;
        for (EdgeId e : outs[static_cast<std::size_t>(u)]) {
          const auto w = head_of(e);
          if (merge.count(w) != 0) into.push_back(e);
          if (avoid.count(w) == 0 && merge.count(w) == 0) free.push_back(e);
        }
        if (!into.empty() && rng() % 2 == 0) {
          path.push_back(into[rng() % into.size()]);
          return path;
        }
        if (k == len || free.empty()) break;
        const EdgeId e = free[rng() % free.size()];
        path.push_back(e);
        u = head_of(e);
        avoid.insert(u);
      }
      path.push_back(exit[static_cast<std::size_t>(u)]);
      return path;
    };
    auto nodes_of = [&](const std::vector<EdgeId>& es) {
      std::set<NodeId> out;
      for (EdgeId e : es) {
        out.insert(g->edge_sources()[static_cast<std::size_t>(e)]);
        out.insert(head_of(e));
      }
      out.erase(dst);
      return out;
    };
    auto path_of = [&](const std::vector<EdgeId>& es) {
      std::vector<std::vector<EdgeId>> hops;
      for (EdgeId e : es) hops.push_back({e});
      return std::make_shared<const Path>(*g, hops);
    };
    // The union of two member sets; a union with a cycle is no route (thrown as invalid_argument: the
    // trial is skipped).
    auto unite = [&](std::vector<EdgeId> a, const std::vector<EdgeId>& b) {
      a.insert(a.end(), b.begin(), b.end());
      std::sort(a.begin(), a.end());
      a.erase(std::unique(a.begin(), a.end()), a.end());
      std::vector<int> indegree(static_cast<std::size_t>(inner) + 1, 0);
      for (EdgeId e : a) ++indegree[static_cast<std::size_t>(head_of(e))];
      std::vector<NodeId> ready;
      for (NodeId v = 0; v <= inner; ++v) {
        if (indegree[static_cast<std::size_t>(v)] == 0) ready.push_back(v);
      }
      std::size_t seen = 0;
      while (!ready.empty()) {
        const NodeId v = ready.back();
        ready.pop_back();
        ++seen;
        for (EdgeId e : a) {
          if (g->edge_sources()[static_cast<std::size_t>(e)] == v && --indegree[static_cast<std::size_t>(head_of(e))] == 0) ready.push_back(head_of(e));
        }
      }
      if (seen != static_cast<std::size_t>(inner) + 1) throw std::invalid_argument("cyclic union");
      return a;
    };
    const NodeId d[1] = {dst};
    std::vector<Part> with;
    std::vector<Part> without;
    int pooled_at = -1;
    std::vector<NodeId> pooled_sources;
    std::vector<EdgeId> pooled_members;
    std::vector<std::vector<RouteDefinition>> defs_of;  // per demand built from routes: its definitions
    try {
      const int n = 1 + static_cast<int>(rng() % 4);
      for (int i = 0; i < n; ++i) {
        const auto kind = rng() % 3;
        const auto h = static_cast<NodeId>(rng() % inner);
        if (kind == 0) {
          std::vector<RouteDefinition> defs {RouteDefinition {h, path_of(walk(h, {}, {})), 1 + static_cast<Amount>(rng() % 3)}};
          if (rng() % 2 == 0) defs.push_back(RouteDefinition {h, path_of(walk(h, {}, {})), 1 + static_cast<Amount>(rng() % 3)});
          const NodeId s[1] = {h};
          with.push_back({ForwardingProgram::from_routes(*g, s, d, defs), nullptr, {static_cast<Amount>(rng() % 6)}});
          defs_of.push_back(defs);
        } else if (kind == 1) {
          const auto w1 = walk(h, {}, {});
          const auto members = rng() % 2 == 0 ? unite(w1, walk(h, {}, {})) : w1;
          const NodeId s[1] = {h};
          const std::vector<RouteDefinition> defs {RouteDefinition {h, dag_of(*g, members), 1}, RouteDefinition {h, path_of(w1), 0}};
          with.push_back({ForwardingProgram::from_routes(*g, s, d, defs), nullptr, {static_cast<Amount>(rng() % 6)}});
          defs_of.push_back(defs);
        } else {
          NodeId h2 = static_cast<NodeId>(rng() % inner);
          if (h2 == h) h2 = (h + 1) % inner;
          const auto members = unite(walk(h, {}, {}), walk(h2, {}, {}));
          std::vector<NodeId> s {std::min(h, h2), std::max(h, h2)};
          with.push_back({ForwardingProgram::from_dag(*g, s, d, dag_of(*g, members), RouteOptions {}), nullptr,
                          {static_cast<Amount>(rng() % 6), static_cast<Amount>(rng() % 6)}});
          defs_of.emplace_back();
          if (pooled_at < 0) {
            pooled_at = i;
            pooled_sources = s;
            pooled_members = members;
          }
        }
      }
      without = with;
      auto kind = static_cast<std::size_t>(rng() % 7);
      if (n == 1 && kind != 6 && rng() % 2 == 0) kind = pooled_at >= 0 ? 2 + rng() % 2 : 6;  // the element joins the one demand
      if ((kind == 2 || kind == 3) && pooled_at < 0) kind = 0;
      const auto h = static_cast<NodeId>(rng() % inner);
      const NodeId s[1] = {h};
      if (kind == 0) {
        with.push_back({ForwardingProgram::from_routes(*g, s, d, std::vector<RouteDefinition> {RouteDefinition {h, path_of(walk(h, {}, {})), 1}}),
                        nullptr, {0}});
      } else if (kind == 1) {
        const auto w = walk(h, {}, {});
        const auto p = rng() % w.size();
        const auto black = std::make_shared<const Eligibility>(g->eligibility({}, std::vector<EdgeId> {w[p]}));
        const Amount z = 1 + static_cast<Amount>(rng() % 3);
        with.push_back({ForwardingProgram::from_routes(*g, s, d, std::vector<RouteDefinition> {RouteDefinition {h, path_of(w), 1}}), black, {z}});
        const std::vector<EdgeId> cut(w.begin(), w.begin() + static_cast<std::ptrdiff_t>(p) + 1);
        const NodeId end[1] = {head_of(cut.back())};
        without.push_back({ForwardingProgram::from_routes(*g, s, end, std::vector<RouteDefinition> {RouteDefinition {h, path_of(cut), 1}}), black, {z}});
      } else if (kind == 2 || kind == 3) {
        const auto live = nodes_of(pooled_members);
        std::vector<NodeId> free;
        for (NodeId v = 0; v < inner; ++v) {
          if (live.count(v) == 0) free.push_back(v);
        }
        if (free.empty()) {
          ++skipped;
          continue;
        }
        const NodeId h0 = free[rng() % free.size()];
        const auto w = walk(h0, live, live);
        std::vector<NodeId> sources = pooled_sources;
        sources.push_back(h0);
        std::sort(sources.begin(), sources.end());
        auto& part = with[static_cast<std::size_t>(pooled_at)];
        const auto old_inject = part.inject;
        part.prog = ForwardingProgram::from_dag(*g, sources, d, dag_of(*g, unite(pooled_members, w)), RouteOptions {});
        part.inject.clear();
        for (NodeId v : part.prog->sources()) {
          if (v == h0) {
            part.inject.push_back(kind == 2 ? 0 : 1 + static_cast<Amount>(rng() % 3));
          } else {
            part.inject.push_back(old_inject[v == pooled_sources[0] ? 0 : 1]);
          }
        }
        if (kind == 3) part.black = std::make_shared<const Eligibility>(g->eligibility({}, std::vector<EdgeId> {w.front()}));
      } else if (kind == 6) {
        // A source of share 0 inside demand i: a new head's tunnel or DAG route, or a pooled source
        // on the demand's DAG.
        const auto i = static_cast<std::size_t>(rng() % static_cast<std::size_t>(n));
        auto& part = with[i];
        std::vector<NodeId> sources(part.prog->sources().begin(), part.prog->sources().end());
        NodeId h0 = -1;
        if (defs_of[i].empty()) {
          std::vector<NodeId> inside;
          for (NodeId v : nodes_of(i == static_cast<std::size_t>(pooled_at) ? pooled_members : std::vector<EdgeId> {})) {
            if (std::find(sources.begin(), sources.end(), v) == sources.end()) inside.push_back(v);
          }
          if (inside.empty()) {
            ++skipped;
            continue;
          }
          h0 = inside[rng() % inside.size()];
          sources.push_back(h0);
          std::sort(sources.begin(), sources.end());
          part.prog = ForwardingProgram::from_dag(*g, sources, d, dag_of(*g, pooled_members), RouteOptions {});
        } else {
          h0 = static_cast<NodeId>(rng() % inner);
          if (std::find(sources.begin(), sources.end(), h0) != sources.end()) {
            ++skipped;
            continue;
          }
          auto defs = defs_of[i];
          const auto w = walk(h0, {}, {});
          if (rng() % 2 == 0) {
            defs.push_back(RouteDefinition {h0, path_of(w), 1});
          } else {
            defs.push_back(RouteDefinition {h0, dag_of(*g, w), 1});
          }
          sources.push_back(h0);
          std::sort(sources.begin(), sources.end());
          part.prog = ForwardingProgram::from_routes(*g, sources, d, defs);
        }
        const auto old = without[i].prog->sources();
        const auto old_inject = part.inject;
        part.inject.clear();
        for (NodeId v : part.prog->sources()) {
          const auto at = std::find(old.begin(), old.end(), v);
          part.inject.push_back(v == h0 ? 0 : old_inject[static_cast<std::size_t>(at - old.begin())]);
        }
      } else {
        // The live walk W from h; the branch on a higher EdgeId out of h than W's first link.
        const auto w = walk(h, {}, {});
        std::vector<EdgeId> later;
        for (EdgeId e : outs[static_cast<std::size_t>(h)]) {
          if (e > w.front() && std::find(w.begin(), w.end(), e) == w.end()) later.push_back(e);
        }
        if (exit[static_cast<std::size_t>(h)] > w.front() && exit[static_cast<std::size_t>(h)] != w.front()) later.push_back(exit[static_cast<std::size_t>(h)]);
        if (later.empty()) {
          ++skipped;
          continue;
        }
        const EdgeId branch = later[rng() % later.size()];
        auto live = nodes_of(w);
        std::vector<EdgeId> dead {branch};
        const auto b = head_of(branch);
        if (b != dst && live.count(b) == 0) {
          auto avoid = live;
          const auto rest = walk(b, avoid, live);
          dead.insert(dead.end(), rest.begin(), rest.end());
        }
        auto program = [&](const std::vector<EdgeId>& members) {
          if (kind == 4) return ForwardingProgram::from_dag(*g, s, d, dag_of(*g, members), RouteOptions {});
          return ForwardingProgram::from_routes(*g, s, d,
                                                std::vector<RouteDefinition> {RouteDefinition {h, dag_of(*g, members), 1},
                                                                              RouteDefinition {h, path_of(w), 0}});
        };
        with.push_back({program(unite(w, dead)), nullptr, {1}});
        without.push_back({program(w), nullptr, {1}});
      }
      ++kinds[kind];
      auto evaluate = [&](const std::vector<Part>& parts) {
        std::vector<std::unique_ptr<detail::ForwardingView>> views;
        std::vector<detail::ForwardingDemand> demands;
        for (const auto& p : parts) {
          views.push_back(std::make_unique<detail::ForwardingView>(detail::make_forwarding_view(
              *g, p.prog, nullptr, std::vector<NodeId>(p.prog->destinations().begin(), p.prog->destinations().end()), p.black.get())));
          demands.push_back({views.back().get(), p.inject});
        }
        Seen seen;
        std::vector<ProgramFlow> got;
        try {
          got = detail::lossy(demands, capacities(*g));
        } catch (const std::invalid_argument& e) {
          seen.error = e.what();
        }
#if NETGRAPH_CORE_PRIVATE_HOOKS
        seen.rounds = detail::thread_lossy_counters().rounds;
        seen.groups = detail::thread_lossy_counters().groups;
        seen.barriers = detail::thread_lossy_counters().barriers;
        seen.dependencies = detail::thread_lossy_counters().dependencies;
#endif
        for (const auto& f : got) {
          seen.delivered.push_back(f.delivered);
          std::map<EdgeId, Amount> c;
          std::map<EdgeId, Amount> dr;
          for (const auto& entry : f.entries) {
            for (const auto& [e, x] : entry.edge_flows) c[e] += x;
            for (const auto& [e, x] : entry.edge_drops) dr[e] += x;
          }
          seen.carried.push_back(c);
          seen.drops.push_back(dr);
        }
        return seen;
      };
      const auto a = evaluate(with);
      const auto b = evaluate(without);
      const std::string where = "trial " + std::to_string(trial) + ", kind " + std::to_string(kind);
      EXPECT_EQ(a.error, b.error) << where;
      EXPECT_EQ(a.rounds, b.rounds) << where;
      EXPECT_EQ(a.groups, b.groups) << where;
      EXPECT_EQ(a.barriers, b.barriers) << where;
#if NETGRAPH_CORE_PRIVATE_HOOKS
      EXPECT_EQ(a.dependencies, b.dependencies + (kind == 3 ? 1 : 0)) << where;
#endif
      if (a.error.empty() && b.error.empty()) {
        // Demands present in both runs: the live ones, and the added tunnel cut or branch removed.
        const std::size_t both = std::min(a.delivered.size(), b.delivered.size());
        for (std::size_t i = 0; i < both; ++i) {
          EXPECT_EQ(a.delivered[i], b.delivered[i]) << where << ", demand " << i;
          EXPECT_EQ(a.carried[i], b.carried[i]) << where << ", demand " << i;
          if (!(kind == 3 && static_cast<int>(i) == pooled_at)) EXPECT_EQ(a.drops[i], b.drops[i]) << where << ", demand " << i;
        }
      }
      errors += a.error.empty() ? 0 : 1;
      cyclic += a.groups > 0 ? 1 : 0;
      ++checked;
    } catch (const std::invalid_argument&) {
      ++skipped;  // a random member set the builders refuse (a cycle, no route to the destination)
    }
    if (::testing::Test::HasFailure()) return;
  }
  RecordProperty("checked", checked);
  RecordProperty("skipped", skipped);
  RecordProperty("cyclic", cyclic);
  RecordProperty("errors", errors);
  EXPECT_GT(checked, 3500);
  for (std::size_t k = 0; k < kinds.size(); ++k) EXPECT_GT(kinds[k], 300) << "kind " << k;
#if NETGRAPH_CORE_PRIVATE_HOOKS
  EXPECT_GT(cyclic, 50);  // not vacuous: groups iterated
#endif
}

// The source-labelled round state. E = e0 0->1 (1), F = e1 2->3 (1); wide e2 4->0, e3 1->2, e4
// 5->2, e5 3->6, e6 7->2, e7 3->0, e8 1->8. A pooled sources {4, 5} -> 6 (offer 2: one each), B a
// tunnel [e6, e1, e7, e0, e8] from 7 to 8 (offer 1). One cyclic group {E, F}; the state (E: A4, B;
// F: A4, A5, B) settles at round 3 (x1 ((1, 0), (0, 1, 0)), x2 = x3 ((1, 0), (1, 0, 0))), while
// per-demand aggregates would repeat at round 2. A delivers 1 (source 4), source 5 drops 1 on F; B
// drops 1 on F.
TEST(JointForwarding, PooledCycleKeepsSourceState) {
  const auto g = shared(make_graph(9, {{0, 1, 1, 1},
                                       {2, 3, 1, 1},
                                       {4, 0, 100, 1},
                                       {1, 2, 100, 1},
                                       {5, 2, 100, 1},
                                       {3, 6, 100, 1},
                                       {7, 2, 100, 1},
                                       {3, 0, 100, 1},
                                       {1, 8, 100, 1}}));
  const auto id = id_of_input(*g);
  auto ids = [&](std::vector<std::int32_t> in) {
    std::vector<EdgeId> out;
    for (auto i : in) out.push_back(id[static_cast<std::size_t>(i)]);
    return out;
  };
  const auto a = pooled_program(g, {4, 5}, 6);
  ASSERT_TRUE(detail::ProgramAccess::data(*a).pools_source_inflows);
  const auto va = detail::make_forwarding_view(*g, a, nullptr, std::vector<NodeId> {6}, nullptr);
  const auto bp = build(*g, {{{7}, {1}, {{0, ids({6, 1, 7, 0, 8}), 1}}, {}}}, 8);
  const Amounts ia {1, 1};
  const Amounts ib {1};
  const std::vector<detail::ForwardingDemand> demands {{&va, ia}, {bp.views[0].get(), ib}};
  const auto got = detail::lossy(demands, capacities(*g));
#if NETGRAPH_CORE_PRIVATE_HOOKS
  EXPECT_EQ(detail::thread_lossy_counters().groups, 1);
  EXPECT_EQ(detail::thread_lossy_counters().rounds, 3);
#endif
  EXPECT_EQ(got[0].delivered, 1);
  EXPECT_EQ(got[0].entries[0].delivered, 1);
  EXPECT_EQ(got[0].entries[1].delivered, 0);
  EXPECT_EQ(got[0].entries[1].edge_drops, (Pairs {{id[1], 1}}));
  EXPECT_EQ(got[1].delivered, 0);
  EXPECT_EQ(got[1].entries[0].edge_drops, (Pairs {{id[1], 1}}));
}

} // namespace


// Source labels also reach destinations, which have no outgoing occurrence. A long chain
// followed by a disconnected live component therefore crosses the occurrence index's width.
TEST(JointForwarding, PooledHeadLabelOffsetsIncludeDestinationsWithoutNarrowing) {
  constexpr std::int64_t chain_heads = 65535;
  constexpr std::int64_t entries = chain_heads + 1;
  constexpr std::int64_t occurrences = chain_heads * (chain_heads + 1) / 2 + 1;
  constexpr std::int64_t routers = chain_heads + 3;
  EXPECT_NO_THROW(detail::check_forwarding_count(detail::ForwardingCounter::Entry, entries));
  EXPECT_NO_THROW(detail::check_forwarding_count(detail::ForwardingCounter::Occurrence, occurrences));
  EXPECT_NO_THROW(detail::check_forwarding_count(detail::ForwardingCounter::Router, routers));
  using Offset = detail::HeadLabelOffset;
  // Each label is a root entry or is reached over an incoming occurrence. The cumulative
  // offset must hold their sum, although each separate preflight counter fits in int32.
  EXPECT_GE(static_cast<std::uint64_t>(std::numeric_limits<Offset>::max()),
            2 * static_cast<std::uint64_t>(detail::kForwardingLimit));
  std::int64_t exact_labels = 0;
  Offset offset = 0;
  for (std::int64_t labels = 1; labels <= chain_heads; ++labels) {
    exact_labels += labels;
    offset = static_cast<Offset>(exact_labels);
    ASSERT_EQ(static_cast<std::int64_t>(offset), exact_labels);
  }
  exact_labels += chain_heads;  // destination labels have no outgoing occurrence
  offset = static_cast<Offset>(exact_labels);
  EXPECT_EQ(exact_labels, 2147516415LL);
  EXPECT_EQ(static_cast<std::int64_t>(offset), exact_labels);
  EXPECT_GT(static_cast<std::uint64_t>(offset), static_cast<std::uint64_t>(detail::kForwardingLimit));
  const Offset disconnected_begin = offset;
  offset = static_cast<Offset>(++exact_labels);
  EXPECT_EQ(offset - disconnected_begin, Offset {1});
  offset = static_cast<Offset>(++exact_labels);  // its destination
  EXPECT_EQ(static_cast<std::int64_t>(offset), 2147516417LL);
  // The bound includes destination labels; it is not the occurrence-only bound.
  EXPECT_LE(static_cast<std::int64_t>(offset), entries + occurrences);
  EXPECT_EQ(occurrences, 2147450881LL);
}
