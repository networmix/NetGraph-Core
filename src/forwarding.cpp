/*
  The forwarding engine of programs (forwarding.hpp): views, the injection, the lossless bound
  and the loss rule.

  A pooled program (one split context) runs on its plan's kernels: allocate and the T_max search
  for the bound, forward() with its attribution for the losses. Every other program (tunnels and
  supplied split routes, each split route a one-root context) runs on the general network of one
  call: the union of its live entries' links, each link with its occurrences (a tunnel's hop, a
  split route's member), and per split route its routers. Each occurrence belongs to one source
  and one entry, so the loss rule at a link is two levels of highest averages.

  The evaluation (A+): the no-drop loads U0 first; serial normalization marks the second link
  of every transparent pair (equal free capacity, every occurrence of the first continuing on the
  second and every occurrence of the second coming from the first, with nothing injected, absorbed,
  merged or split between them) as a pass-through; every other link with U0 > V is a barrier. The
  events are the barriers and the routers; an event waits for its predecessors (a router for its
  in-members' carried amounts, a barrier for all its arrivals, a tunnel's next barrier for its
  previous one). Tarjan's algorithm groups them into strongly connected components (SCCs), which
  run in dependency order, the smallest key first among the ready ones ((0, smallest barrier
  EdgeId) for a group with a barrier, else (1, its first event)). A cyclic group of b barriers is
  iterated from zero carried: one round replays the group's private forwarding from the
  previous carried state x, then computes every new carry with the loss rule into a separate state;
  equal states solve it, a state equal to the one two rounds before raises the period-two error, and
  64 b rounds without either raise the round-limit error. Once every group is solved, a call that
  iterated checks the whole result independently (every router's split, every link's rule and
  capacity, every continuation, the conservation of every entry) before anything is returned.

  Amounts are int64 units below 2^62: a call's total bounds every arrival and load. The per-thread
  scratch of the general network is rebuilt by every call before it is read.
*/
#include "netgraph/core/detail/forwarding.hpp"

#include <algorithm>
#include <limits>
#include <queue>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>

#include "netgraph/core/detail/int128.hpp"
#include "netgraph/core/detail/route_builders.hpp"

namespace netgraph::core::detail {

#if NETGRAPH_CORE_PRIVATE_HOOKS
LossyCounters& thread_lossy_counters() noexcept {
  thread_local LossyCounters c;
  return c;
}

std::array<std::int64_t, 5>& thread_forwarding_limits() noexcept {
  thread_local std::array<std::int64_t, 5> limits {};
  return limits;
}
#endif

void check_forwarding_count(ForwardingCounter counter, std::int64_t n) {
  std::int64_t limit = kForwardingLimit;
#if NETGRAPH_CORE_PRIVATE_HOOKS
  if (const auto lower = thread_forwarding_limits()[static_cast<std::size_t>(counter)]; lower > 0) limit = lower;
#endif
  if (n <= limit) return;
  static constexpr const char* kNames[] = {"entry", "occurrence", "router", "event", "dependency"};
  throw std::invalid_argument(std::string("forwarding ") + kNames[static_cast<std::size_t>(counter)] + " limit exceeded: " +
                              std::to_string(n) + " > " + std::to_string(limit));
}

namespace {

// ================================================================ the view

bool has_black_holes(const Eligibility* dark) noexcept { return dark != nullptr && dark->uid() != 0; }

// The roster slots that reach an allowed destination over their basis under the mask: one
// reverse search from the allowed destinations over every allowed edge (topology), or per slot a
// forward search over its own allowed basis edges.
std::vector<std::uint8_t> replay_roster(const Graph& graph, const ProgramData& p, const std::vector<std::uint8_t>& edge_allowed,
                                        const std::vector<std::uint8_t>& node_allowed,
                                        const std::vector<std::uint8_t>& is_dst) {
  const std::size_t k = p.reachable_sources.size();
  std::vector<std::uint8_t> keep(k, 0);
  const auto edge_sources = graph.edge_sources();
  const auto edge_destinations = graph.edge_destinations();
  if (p.uses_topology_basis) {
    const auto in_row = graph.incoming_offsets();
    const auto in_adj = graph.incoming_edges();
    std::vector<std::uint8_t> seen(static_cast<std::size_t>(graph.num_nodes()), 0);
    std::vector<NodeId> stack;
    for (std::size_t v = 0; v < is_dst.size(); ++v) {
      if (is_dst[v] != 0) {
        seen[v] = 1;
        stack.push_back(static_cast<NodeId>(v));
      }
    }
    while (!stack.empty()) {
      const auto v = static_cast<std::size_t>(stack.back());
      stack.pop_back();
      for (auto i = in_row[v]; i < in_row[v + 1]; ++i) {
        const auto e = static_cast<std::size_t>(in_adj[static_cast<std::size_t>(i)]);
        if (edge_allowed[e] == 0) continue;
        const auto u = static_cast<std::size_t>(edge_sources[e]);
        // A path stops at the first destination: a destination is not passed through.
        if (seen[u] != 0 || is_dst[u] != 0) continue;
        seen[u] = 1;
        stack.push_back(static_cast<NodeId>(u));
      }
    }
    for (std::size_t s = 0; s < k; ++s) {
      const auto h = static_cast<std::size_t>(p.reachable_sources[s]);
      keep[s] = node_allowed[h] != 0 && seen[h] != 0 ? 1 : 0;
    }
    return keep;
  }
  std::vector<std::uint8_t> seen(static_cast<std::size_t>(graph.num_nodes()), 0);
  std::vector<NodeId> touched;
  std::vector<NodeId> stack;
  std::vector<std::pair<std::span<const EdgeId>, bool>> parts;
  const auto row = graph.outgoing_offsets();
  for (std::size_t s = 0; s < k; ++s) {
    const NodeId h = p.reachable_sources[s];
    if (node_allowed[static_cast<std::size_t>(h)] == 0) continue;
    parts.clear();
    p.basis_parts(s, parts);
    // A forward search over the slot's allowed basis edges (an ascending part's edges of one tail
    // are consecutive; a tunnel's path is scanned whole).
    touched.assign(1, h);
    stack.assign(1, h);
    seen[static_cast<std::size_t>(h)] = 1;
    bool reached = false;
    auto visit = [&](EdgeId x) {
      const auto e = static_cast<std::size_t>(x);
      if (edge_allowed[e] == 0) return;
      const auto v = static_cast<std::size_t>(edge_destinations[e]);
      if (is_dst[v] != 0) {
        reached = true;
        return;
      }
      if (seen[v] != 0) return;
      seen[v] = 1;
      touched.push_back(static_cast<NodeId>(v));
      stack.push_back(static_cast<NodeId>(v));
    };
    while (!stack.empty() && !reached) {
      const NodeId u = stack.back();
      stack.pop_back();
      const auto lo = row[static_cast<std::size_t>(u)];
      const auto hi = row[static_cast<std::size_t>(u) + 1];
      for (const auto& [part, ascending] : parts) {
        if (ascending) {
          for (auto it = std::lower_bound(part.begin(), part.end(), lo); it != part.end() && *it < hi && !reached; ++it) visit(*it);
        } else {
          for (EdgeId x : part) {
            if (x >= lo && x < hi && !reached) visit(x);
          }
        }
      }
    }
    for (NodeId v : touched) seen[static_cast<std::size_t>(v)] = 0;
    keep[s] = reached ? 1 : 0;
  }
  return keep;
}

// An edge set under the mask, pruned to the walks from `source` to the allowed destinations, as a
// record route with its cost (an empty DAG and kUnreachable when nothing remains).
std::pair<std::shared_ptr<const PredecessorDAG>, Cost> masked_route(const Graph& graph, std::span<const EdgeId> edges, NodeId source,
                                                             std::span<const NodeId> destinations, const Eligibility* eligibility) {
  const PredecessorDAG dag = PredecessorDAGBuilder::of_edges(graph, edges);
  auto kept = prune(graph, dag, source, destinations, eligibility);
  if (kept.empty()) return {std::make_shared<const PredecessorDAG>(PredecessorDAGBuilder::empty(graph)), kUnreachable};
  auto out = std::make_shared<const PredecessorDAG>(PredecessorDAGBuilder::of_edges(graph, kept));
  return {out, min_route_cost(graph, *out, source, destinations)};
}

void finish_view(ForwardingView& v) {
  const ProgramData& p = *v.data;
  // The surviving entries per view head, ascending index.
  v.source_entry_offsets_.assign(v.source_indices_.size() + 1, 0);
  v.source_entry_positions_.clear();
  for (std::size_t h = 0; h < v.source_indices_.size(); ++h) {
    for (std::int32_t e : p.entries_for_source(static_cast<std::size_t>(v.source_indices_[h]))) {
      if (v.status_[static_cast<std::size_t>(e)] != TunnelStatus::Failed) v.source_entry_positions_.push_back(e);
    }
    v.source_entry_offsets_[h + 1] = static_cast<std::int32_t>(v.source_entry_positions_.size());
  }
  // Black holes: the black-holed heads, and every up tunnel whose path reaches one.
  v.blackholed_sources_.clear();
  if (v.black_holes != nullptr) {
    v.blackholed_sources_.assign(v.source_indices_.size(), 0);
    for (std::size_t h = 0; h < v.source_indices_.size(); ++h) v.blackholed_sources_[h] = v.black_holes->allows_node(v.source(h)) ? 0 : 1;
    for (std::size_t t = 0; t < p.tunnels.size(); ++t) {
      const auto e = static_cast<std::size_t>(p.tunnel_entry_positions[t]);
      if (v.status_[e] != TunnelStatus::Up) continue;
      // An edge is black-holed when the Eligibility disallows it (itself or an endpoint).
      for (EdgeId x : path_of(p.tunnels[t])) {
        if (!v.black_holes->allows_edge(x)) {
          v.status_[e] = TunnelStatus::Blackholed;
          break;
        }
      }
    }
  }
}

// A borrowed view's state written out (the program's own: every head, every entry up), for a
// view that then changes it (black holes).
void materialize(ForwardingView& v) {
  if (!v.borrowed) return;
  const ProgramData& p = *v.data;
  const std::size_t n = p.entries.size();
  v.source_indices_.resize(p.reachable_sources.size());
  for (std::size_t s = 0; s < p.reachable_sources.size(); ++s) v.source_indices_[s] = static_cast<std::int32_t>(s);
  v.source_weights_ = p.source_weight_values;
  v.status_.resize(n);
  v.routes_.resize(n);
  v.min_route_costs_.resize(n);
  v.context_source_indices_.resize(n);
  for (std::size_t e = 0; e < n; ++e) {
    v.status_[e] = v.status(e);
    v.routes_[e] = p.entries[e].route;
    v.min_route_costs_[e] = p.entries[e].min_route_cost;
    v.context_source_indices_[e] = v.context_source_index(e);
  }
  v.prepared_routes_.clear();
  for (const auto& c : p.split_contexts) v.prepared_routes_.push_back(c.prepared_route);
  v.borrowed = false;
}

} // namespace

ForwardingView make_forwarding_view(const Graph& graph, std::shared_ptr<const ForwardingProgram> program, const Eligibility* mask,
                         std::span<const NodeId> destinations, const Eligibility* dark, bool records) {
  ForwardingView v;
  v.data = &ProgramAccess::data(*program);
  v.program = std::move(program);
  v.black_holes = has_black_holes(dark) ? dark : nullptr;
  const ProgramData& p = *v.data;
  const std::size_t n = p.entries.size();
  v.pooled = p.pools_source_inflows;
  if (mask == nullptr || mask->uid() == 0) {
    // The program's own state: its roster, its routes and plans, borrowed from the program (a
    // tunnel's record route is made here; black holes change the statuses).
    v.borrowed = true;
    if (v.black_holes == nullptr && !(records && !p.tunnels.empty())) return v;
    materialize(v);
    for (std::size_t e = 0; records && e < n; ++e) {
      if (p.entries[e].tunnel_position < 0) continue;
      const auto path = path_of(p.tunnels[static_cast<std::size_t>(p.entries[e].tunnel_position)]);
      v.routes_[e] = std::make_shared<const PredecessorDAG>(path.empty() ? PredecessorDAGBuilder::empty(graph) : PredecessorDAGBuilder::of_edges(graph, path));
    }
    finish_view(v);
    return v;
  }
  v.status_.assign(n, TunnelStatus::Up);
  v.routes_.resize(n);
  v.min_route_costs_.assign(n, kUnreachable);
  v.context_source_indices_.assign(n, -1);
  // The replay: the mask and the program's own exclusions together.
  const auto num_nodes = static_cast<std::size_t>(graph.num_nodes());
  const auto num_edges = static_cast<std::size_t>(graph.num_edges());
  std::vector<std::uint8_t> node_allowed(num_nodes, 0);
  std::vector<std::uint8_t> edge_allowed(num_edges, 0);
  for (std::size_t x = 0; x < num_nodes; ++x) node_allowed[x] = mask->allows_node(static_cast<NodeId>(x)) ? 1 : 0;
  for (NodeId x : p.excluded_nodes) node_allowed[static_cast<std::size_t>(x)] = 0;
  for (std::size_t e = 0; e < num_edges; ++e) edge_allowed[e] = mask->allows_edge(static_cast<EdgeId>(e)) ? 1 : 0;
  for (EdgeId e : p.excluded_edges) edge_allowed[static_cast<std::size_t>(e)] = 0;
  std::vector<NodeId> live_dests;
  std::vector<std::uint8_t> is_dst(num_nodes, 0);
  for (NodeId d : destinations) {
    if (node_allowed[static_cast<std::size_t>(d)] == 0) continue;
    live_dests.push_back(d);
    is_dst[static_cast<std::size_t>(d)] = 1;
  }
  // The program's exclusions as an Eligibility of the replay (null when the program has none).
  std::optional<Eligibility> both;
  const Eligibility* replay = mask;
  if (!p.excluded_nodes.empty() || !p.excluded_edges.empty()) {
    both = Eligibility::intersect(*mask, graph.eligibility(p.excluded_nodes, p.excluded_edges));
    replay = &*both;
  }
  const auto keep = replay_roster(graph, p, edge_allowed, node_allowed, is_dst);
  for (std::size_t s = 0; s < p.reachable_sources.size(); ++s) {
    if (keep[s] == 0) continue;
    v.source_indices_.push_back(static_cast<std::int32_t>(s));
    if (!p.source_weight_values.empty()) v.source_weights_.push_back(p.source_weight_values[s]);
  }
  std::vector<std::int32_t> view_head(p.reachable_sources.size(), -1);
  for (std::size_t h = 0; h < v.source_indices_.size(); ++h) view_head[static_cast<std::size_t>(v.source_indices_[h])] = static_cast<std::int32_t>(h);
  // Entries: a tunnel fails when it touches an excluded node or edge (an edge is excluded with its
  // endpoints); a split route keeps its allowed edges on the walks from its head to an allowed
  // destination and fails when none remains. Every entry of a head outside the replay roster fails.
  // A tunnel that setup left down stays down with its status (replay never brings one up).
  for (std::size_t e = 0; e < n; ++e) {
    const auto& entry = p.entries[e];
    bool failed = view_head[static_cast<std::size_t>(entry.roster_index)] < 0;
    if (entry.tunnel_position >= 0) {
      const auto& t = p.tunnels[static_cast<std::size_t>(entry.tunnel_position)];
      if (t.status != TunnelStatus::Up) {
        v.status_[e] = t.status;
        if (records) v.routes_[e] = std::make_shared<const PredecessorDAG>(PredecessorDAGBuilder::empty(graph));
        continue;
      }
      for (EdgeId x : path_of(t)) failed = failed || edge_allowed[static_cast<std::size_t>(x)] == 0;
      if (!failed) {
        if (records) v.routes_[e] = std::make_shared<const PredecessorDAG>(PredecessorDAGBuilder::of_edges(graph, path_of(t)));
        v.min_route_costs_[e] = entry.min_route_cost;
      }
    } else if (!failed) {
      // A split route's route under the mask (its edges before weight pruning).
      const auto edges = PredecessorDAGBuilder::ascending(*entry.route);
      auto [route, cost] = masked_route(graph, edges, entry.head, live_dests, replay);
      failed = cost == kUnreachable;
      v.routes_[e] = std::move(route);
      v.min_route_costs_[e] = cost;
    }
    if (failed) {
      v.status_[e] = TunnelStatus::Failed;
      v.routes_[e] = records ? std::make_shared<const PredecessorDAG>(PredecessorDAGBuilder::empty(graph)) : nullptr;
      v.min_route_costs_[e] = kUnreachable;
    }
  }
  // The split contexts restricted to the mask, rooted at their surviving entries' heads.
  v.prepared_routes_.reserve(p.split_contexts.size());
  std::vector<NodeId> roots;
  std::vector<std::int32_t> entries;
  for (const auto& c : p.split_contexts) {
    roots.clear();
    entries.clear();
    for (std::size_t r = 0; r < c.source_entry_positions.size(); ++r) {
      const auto e = static_cast<std::size_t>(c.source_entry_positions[r]);
      if (v.status_[e] == TunnelStatus::Failed) continue;
      roots.push_back(c.prepared_route->sources[r]);
      entries.push_back(c.source_entry_positions[r]);
    }
    v.prepared_routes_.push_back(std::make_shared<const PreparedRoute>(restrict_prepared_route(graph, *c.prepared_route, roots, live_dests, replay)));
    for (std::size_t r = 0; r < entries.size(); ++r) v.context_source_indices_[static_cast<std::size_t>(entries[r])] = static_cast<std::int32_t>(r);
  }
  finish_view(v);
  return v;
}

ForwardingView with_black_holes(const ForwardingView& view, const Eligibility* dark) {
  ForwardingView v = view;
  v.black_holes = has_black_holes(dark) ? dark : nullptr;
  if (v.borrowed && v.black_holes == nullptr) return v;
  materialize(v);
  // The statuses again from the replay's (a black hole never fails a tunnel; a mask does).
  for (std::size_t t = 0; t < v.data->tunnels.size(); ++t) {
    const auto e = static_cast<std::size_t>(v.data->tunnel_entry_positions[t]);
    if (v.status_[e] == TunnelStatus::Blackholed) v.status_[e] = TunnelStatus::Up;
  }
  finish_view(v);
  return v;
}

// ================================================================ the injection

void source_shares(const ForwardingView& view, Amount x, std::span<Amount> out) {
  std::fill(out.begin(), out.end(), 0);
  if (out.empty()) return;
  if (view.source_weights().empty()) {
    (void)round_robin(x, out);
  } else {
    highest_averages(x, view.source_weights(), out);
  }
}

void source_offer_increments(const ForwardingView& view, const Counter128& X, Amount v, std::span<Amount> out) {
  const std::size_t k = out.size();
  std::fill(out.begin(), out.end(), 0);
  if (k == 0) return;
  const auto weights = view.source_weights();
  if (weights.empty()) {
    round_robin_increments(X, v, out);
    return;
  }
  // H(X + c W) = H(X) + c w (the divisor rule's thresholds shift by c), so only r = X mod W and
  // r + v enter: H(X + v) - H(X) = floor((r + v) / W) w + H((r + v) mod W) - H(r). Captured
  // weights are a controller's solved supplies: W is its delivery, below the graph's capacity sum
  // (< 2^62), so W, r + v and c w fit 64 bits.
  std::uint64_t W = 0;
  for (Amount w : weights) W += static_cast<std::uint64_t>(w);
  if (W == 0) return;
  const std::uint64_t r = mod_u128(U128 {X.hi, X.lo}, W);
  const std::uint64_t t = r + static_cast<std::uint64_t>(v);
  const auto c = static_cast<Amount>(t / W);
  std::vector<Amount> a(k);
  std::vector<Amount> b(k);
  highest_averages(static_cast<Amount>(r), weights, a);
  highest_averages(static_cast<Amount>(t % W), weights, b);
  for (std::size_t i = 0; i < k; ++i) out[i] = c * weights[i] + b[i] - a[i];
}

namespace {

// Each head's share over its surviving entries (highest averages by weight) into `share` (per
// program entry); what a head cannot place (black-holed, no surviving entry of positive weight)
// goes to `source_drops_by_source`.
void entry_shares(const ForwardingView& view, std::span<const Amount> heads, std::vector<Amount>& share,
                  std::vector<Amount>& source_drops_by_source, std::vector<Amount>& w, std::vector<Amount>& part) {
  const ProgramData& p = *view.data;
  share.assign(p.entries.size(), 0);
  source_drops_by_source.assign(heads.size(), 0);
  for (std::size_t h = 0; h < heads.size(); ++h) {
    const Amount x = heads[h];
    if (x == 0) continue;
    if (view.is_blackholed_source(h)) {
      source_drops_by_source[h] = x;
      continue;
    }
    const auto entries = view.entries_for_source(h);
    w.resize(entries.size());
    bool positive = false;  // a weight sum may exceed 2^63: only whether one weight is positive
    for (std::size_t i = 0; i < entries.size(); ++i) {
      w[i] = p.entries[static_cast<std::size_t>(entries[i])].weight;
      positive = positive || w[i] > 0;
    }
    if (!positive) {
      source_drops_by_source[h] = x;
      continue;
    }
    part.resize(entries.size());
    highest_averages(x, w, part);
    for (std::size_t i = 0; i < entries.size(); ++i) share[static_cast<std::size_t>(entries[i])] = part[i];
  }
}

// The positive entries of a member-aligned array as (EdgeId, amount) pairs.
std::vector<std::pair<EdgeId, Amount>> positive_pairs(std::span<const EdgeId> members, std::span<const Amount> a) {
  std::vector<std::pair<EdgeId, Amount>> out;
  out.reserve(static_cast<std::size_t>(std::count_if(a.begin(), a.end(), [](Amount x) { return x != 0; })));
  for (std::size_t m = 0; m < members.size(); ++m) {
    if (a[m] != 0) out.emplace_back(members[m], a[m]);
  }
  return out;
}

void sort_pairs(std::vector<std::pair<EdgeId, Amount>>& v) {
  std::sort(v.begin(), v.end(), [](const auto& a, const auto& b) { return a.first < b.first; });
}

// A forwarding's outcome before anything is sent; the heads' own losses only for Lossy (no
// Lossless forwarding loses anything).
ProgramFlow empty_flow(const ForwardingView& view, bool lossy) {
  ProgramFlow out;
  out.entries.resize(view.data->entries.size());
  out.offered_by_source.assign(view.num_sources(), 0);
  if (lossy) out.source_drops_by_source.assign(view.num_sources(), 0);
  return out;
}

// ================================================================ the pooled form

// The free capacity of the pooled plan's members with the black holes at 0.
void pooled_view(const PreparedRoute& plan, const ForwardingView& view, const AvailableCapacityFn& free, std::vector<Amount>& out) {
  out.resize(plan.forwarding_edges.size());
  free(plan.forwarding_edges, out);
  if (view.black_holes == nullptr) return;
  for (std::size_t m = 0; m < plan.forwarding_edges.size(); ++m) {
    if (!view.black_holes->allows_edge(plan.forwarding_edges[m])) out[m] = 0;
  }
}

// The pooled form's per-root amounts as entry flows (root r carries entry source_entry_positions[r]).
void pooled_entries(const ForwardingView& view, const PreparedRoute& plan, const Attribution* attr,
                    std::span<const Amount> alloc, std::span<const Amount> drops, ProgramFlow& out) {
  for (std::size_t h = 0; h < view.num_sources(); ++h) {
    const auto e = static_cast<std::size_t>(view.entries_for_source(h).front());
    auto& f = out.entries[e];
    if (attr == nullptr) {
      // One root holds the whole forwarding (the attribution's degenerate form).
      f.edge_flows = positive_pairs(plan.forwarding_edges, alloc);
      if (!drops.empty()) f.edge_drops = positive_pairs(plan.forwarding_edges, drops);
      continue;
    }
    for (const auto& x : attr->entries[h]) {
      const EdgeId edge = plan.forwarding_edges[static_cast<std::size_t>(x.edge_position)];
      if (x.carried != 0) f.edge_flows.emplace_back(edge, x.carried);
      if (x.drops != 0) f.edge_drops.emplace_back(edge, x.drops);
    }
  }
}

// The first black-holed head's index (the equal injection gives head i a unit from total i + 1),
// or the number of heads.
Amount black_hole_offer_cap(const ForwardingView& view) {
  for (std::size_t h = 0; h < view.num_sources(); ++h) {
    if (view.is_blackholed_source(h)) return static_cast<Amount>(h);
  }
  return std::numeric_limits<Amount>::max();
}

ProgramFlow pooled_lossless(const ForwardingView& view, const AvailableCapacityFn& free, Amount lo, std::optional<Amount> hi) {
  const PreparedRoute& plan = *view.prepared_route(0);
  ProgramFlow out = empty_flow(view, false);
  const std::size_t k = view.num_sources();
  // A black-holed head may receive nothing: the total stays at most its index.
  const Amount cap = std::min(hi.value_or(std::numeric_limits<Amount>::max()), black_hole_offer_cap(view));
  std::vector<Amount> V;
  pooled_view(plan, view, free, V);
  std::vector<Amount> own(plan.forwarding_edges.size(), 0);
  if (lo != 0 && !allocate(plan, lo, own)) throw std::logic_error("forwarding: the admitted total does not fit");
  Amount total = lo;
  if (cap > lo && !plan.forwarding_edges.empty()) {
    const std::optional<Amount> request =
        cap == std::numeric_limits<Amount>::max() ? std::nullopt : std::optional<Amount>(cap - lo);
    auto step = lossless_step(plan, lo, V, own, request);
    total = step.admitted;
    own = std::move(step.edge_flows);
  }
  out.admitted = total;
  source_shares(view, total, out.offered_by_source);
  out.offered = total;
  out.delivered = total;
  for (std::size_t h = 0; h < k; ++h) {
    const auto e = static_cast<std::size_t>(view.entries_for_source(h).front());
    out.entries[e].assigned = out.offered_by_source[h];
    out.entries[e].delivered = out.offered_by_source[h];
  }
  if (k > 1 && total > 0) {
    // Several roots: forward the total again to attribute it.
    Attribution attr;
    const auto fw = forward(plan, Injection {out.offered_by_source, false}, {}, &attr);
    if (fw.edge_flows != own) throw std::logic_error("forwarding: the attributed loads differ from the admitted ones");
    pooled_entries(view, plan, &attr, {}, {}, out);
  } else if (k == 1) {
    pooled_entries(view, plan, nullptr, own, {}, out);
  }
  return out;
}

ProgramFlow pooled_lossy(const ForwardingView& view, const AvailableCapacityFn& free, std::span<const Amount> injections) {
  const PreparedRoute& plan = *view.prepared_route(0);
  ProgramFlow out = empty_flow(view, true);
  const std::size_t k = view.num_sources();
  std::vector<Amount> inject(injections.begin(), injections.end());
  for (std::size_t h = 0; h < k; ++h) {
    out.offered_by_source[h] = injections[h];
    out.offered += injections[h];
    if (view.is_blackholed_source(h)) {
      // A black-holed head loses its own share at the source.
      out.source_drops_by_source[h] = injections[h];
      out.source_drops += injections[h];
      inject[h] = 0;
    }
  }
  out.admitted = out.offered;
  std::vector<Amount> V;
  pooled_view(plan, view, free, V);
  Attribution attr;
  const auto fw = forward(plan, Injection {inject, false}, V, k > 1 ? &attr : nullptr);
  out.source_drops += fw.source_drops;
  out.delivered = fw.delivered;
  for (std::size_t h = 0; h < k; ++h) {
    const auto e = static_cast<std::size_t>(view.entries_for_source(h).front());
    auto& f = out.entries[e];
    f.assigned = inject[h];
    if (k > 1) {
      f.delivered = attr.delivered[h];
      f.source_drops = attr.source_drops_by_source[h];
    } else {
      f.delivered = fw.delivered;
      f.source_drops = fw.source_drops;
    }
  }
  pooled_entries(view, plan, k > 1 ? &attr : nullptr, fw.edge_flows, fw.edge_drops, out);
  return out;
}

// ================================================================ the general network

// The live entries of one call (surviving, positive weight, a usable route) and the union of their
// links. A tunnel's occurrences are its hops in path order, a split route's its plan's members.
struct ForwardingNetwork {
  struct Item {
    std::int32_t entry_position {0};
    std::int32_t source_index {0};      // view head
    std::int32_t tunnel_position {-1};
    const PreparedRoute* prepared_route {nullptr};
    std::int32_t source_position {-1};
    std::int32_t first_occurrence {0};      // its first occurrence
    std::int32_t num_occurrences {0};
    std::int32_t first_router_position {0};      // its first router (split routes: one per plan position)
    // A split route's entry whose root is not in its plan: its share has no route.
    [[nodiscard]] bool routeless() const noexcept { return tunnel_position < 0 && source_position < 0; }
    [[nodiscard]] std::int32_t num_routers() const noexcept {
      return prepared_route == nullptr ? 0 : static_cast<std::int32_t>(prepared_route->topological_order.size());
    }
  };
  std::vector<Item> entry_occurrences;
  std::vector<std::int32_t> entry_occurrence_positions;  // per program entry: its item, or -1 (or -2: a routeless entry)
  std::vector<EdgeId> edges;          // ascending
  std::vector<Amount> available_capacity;              // per link
  std::vector<std::int32_t> occurrence_edge_positions; // per occurrence
  std::vector<std::int32_t> occurrence_entry_positions;
  std::vector<std::int32_t> occurrence_route_positions;  // tunnel: hop; split route: member
  std::vector<std::int32_t> edge_occurrence_offsets;  // per link + 1 into occurrences_by_edge
  std::vector<std::int32_t> occurrences_by_edge;     // occurrences, by (head NodeId, entry index)
  std::int32_t num_router_positions {0};
  std::vector<std::int32_t> router_entry_positions;    // per router position: its item
  std::vector<std::int32_t> incoming_edge_positions;   // per router position: its one in-member, -1 none, -2 several
};

std::int32_t edge_position(const ForwardingNetwork& net, EdgeId e) {
  return static_cast<std::int32_t>(std::lower_bound(net.edges.begin(), net.edges.end(), e) - net.edges.begin());
}

// Every entry of the forwarding with a positive weight, by (head, entry index), to `visit` as its
// item before its positions in the network are set (first_occurrence, first_router_position): a
// tunnel with its hops, or a split route with its root's position in the plan and its members
// (none when routeless).
template <class Visit>
void for_each_forwarded_entry(const ForwardingView& view, Visit&& visit) {
  const ProgramData& p = *view.data;
  for (std::size_t h = 0; h < view.num_sources(); ++h) {
    for (const std::int32_t ei : view.entries_for_source(h)) {
      const auto e = static_cast<std::size_t>(ei);
      const auto& entry = p.entries[e];
      if (entry.weight == 0) continue;
      ForwardingNetwork::Item it;
      it.entry_position = ei;
      it.source_index = static_cast<std::int32_t>(h);
      if (entry.tunnel_position >= 0) {
        it.tunnel_position = entry.tunnel_position;
        it.num_occurrences = static_cast<std::int32_t>(path_of(p.tunnels[static_cast<std::size_t>(entry.tunnel_position)]).size());
      } else {
        it.prepared_route = view.prepared_route(static_cast<std::size_t>(entry.split_context_position)).get();
        const auto r = view.context_source_index(e);
        it.source_position = r < 0 ? -1 : it.prepared_route->source_positions[static_cast<std::size_t>(r)];
        if (it.source_position >= 0) it.num_occurrences = static_cast<std::int32_t>(it.prepared_route->forwarding_edges.size());
      }
      visit(it);
    }
  }
}

ForwardingNetwork build_forwarding_network(const ForwardingView& view, const AvailableCapacityFn& free) {
  const ProgramData& p = *view.data;
  // The network's counts, exactly, after validation and before any of its arrays is sized: the
  // forwarded entries (positive weight, with a route), their occurrences and the split routes'
  // routers. A refusal changes nothing.
  std::int64_t num_entries = 0;
  std::int64_t num_occurrences = 0;
  std::int64_t num_routers = 0;
  for_each_forwarded_entry(view, [&](const ForwardingNetwork::Item& it) {
    if (it.routeless()) return;
    ++num_entries;
    num_occurrences += it.num_occurrences;
    num_routers += it.num_routers();
  });
  check_forwarding_count(ForwardingCounter::Entry, num_entries);
  check_forwarding_count(ForwardingCounter::Occurrence, num_occurrences);
  check_forwarding_count(ForwardingCounter::Router, num_routers);
  ForwardingNetwork net;
  net.entry_occurrence_positions.assign(p.entries.size(), -1);
  for_each_forwarded_entry(view, [&](ForwardingNetwork::Item it) {
    const auto e = static_cast<std::size_t>(it.entry_position);
    if (it.routeless()) {
      net.entry_occurrence_positions[e] = -2;  // its share has no route: lost at the source
      return;
    }
    if (it.prepared_route != nullptr) {
      it.first_router_position = net.num_router_positions;
      net.num_router_positions += it.num_routers();
    }
    net.entry_occurrence_positions[e] = static_cast<std::int32_t>(net.entry_occurrences.size());
    net.entry_occurrences.push_back(it);
  });
  std::int32_t nocc = 0;
  for (auto& it : net.entry_occurrences) {
    it.first_occurrence = nocc;
    nocc += it.num_occurrences;
    if (it.tunnel_position >= 0) {
      const auto edges = path_of(p.tunnels[static_cast<std::size_t>(it.tunnel_position)]);
      net.edges.insert(net.edges.end(), edges.begin(), edges.end());
    } else {
      net.edges.insert(net.edges.end(), it.prepared_route->forwarding_edges.begin(), it.prepared_route->forwarding_edges.end());
    }
  }
  std::sort(net.edges.begin(), net.edges.end());
  net.edges.erase(std::unique(net.edges.begin(), net.edges.end()), net.edges.end());
  net.available_capacity.resize(net.edges.size());
  free(net.edges, net.available_capacity);
  for (std::size_t l = 0; view.black_holes != nullptr && l < net.edges.size(); ++l) {
    if (!view.black_holes->allows_edge(net.edges[l])) net.available_capacity[l] = 0;
  }
  const auto n = static_cast<std::size_t>(nocc);
  net.occurrence_edge_positions.resize(n);
  net.occurrence_entry_positions.resize(n);
  net.occurrence_route_positions.resize(n);
  for (std::size_t i = 0; i < net.entry_occurrences.size(); ++i) {
    const auto& it = net.entry_occurrences[i];
    for (std::int32_t j = 0; j < it.num_occurrences; ++j) {
      const auto o = static_cast<std::size_t>(it.first_occurrence + j);
      const EdgeId e = it.tunnel_position >= 0 ? path_of(p.tunnels[static_cast<std::size_t>(it.tunnel_position)])[static_cast<std::size_t>(j)]
                                      : it.prepared_route->forwarding_edges[static_cast<std::size_t>(j)];
      net.occurrence_edge_positions[o] = edge_position(net, e);
      net.occurrence_entry_positions[o] = static_cast<std::int32_t>(i);
      net.occurrence_route_positions[o] = j;
    }
  }
  net.router_entry_positions.resize(static_cast<std::size_t>(net.num_router_positions));
  net.incoming_edge_positions.assign(static_cast<std::size_t>(net.num_router_positions), -1);
  for (std::size_t i = 0; i < net.entry_occurrences.size(); ++i) {
    const auto& it = net.entry_occurrences[i];
    if (it.tunnel_position >= 0) continue;
    for (std::size_t pos = 0; pos < it.prepared_route->topological_order.size(); ++pos) net.router_entry_positions[static_cast<std::size_t>(it.first_router_position) + pos] = static_cast<std::int32_t>(i);
    for (std::int32_t m = 0; m < it.num_occurrences; ++m) {
      auto& slot = net.incoming_edge_positions[static_cast<std::size_t>(it.first_router_position + it.prepared_route->edge_destination_positions[static_cast<std::size_t>(m)])];
      slot = slot == -1 ? m : -2;
    }
  }
  net.edge_occurrence_offsets.assign(net.edges.size() + 1, 0);
  for (auto l : net.occurrence_edge_positions) ++net.edge_occurrence_offsets[static_cast<std::size_t>(l) + 1];
  for (std::size_t l = 0; l < net.edges.size(); ++l) net.edge_occurrence_offsets[l + 1] += net.edge_occurrence_offsets[l];
  net.occurrences_by_edge.resize(n);
  std::vector<std::int32_t> fill(net.edge_occurrence_offsets.begin(), net.edge_occurrence_offsets.end() - 1);
  // Occurrences are numbered by item, items by (head, entry index): each link's list comes out in
  // the loss rule's order (sources ascending, then entries ascending).
  for (std::size_t o = 0; o < n; ++o) net.occurrences_by_edge[static_cast<std::size_t>(fill[static_cast<std::size_t>(net.occurrence_edge_positions[o])]++)] = static_cast<std::int32_t>(o);
  return net;
}

// The no-drop loads of the entries' shares per link (and per occurrence when `occ` is given).
// A split route allocates its share over its plan from its root.
void no_drop(const ForwardingNetwork& net, std::span<const Amount> share, std::vector<Amount>& load, std::vector<Amount>* occ,
             std::vector<Amount>& scratch) {
  load.assign(net.edges.size(), 0);
  if (occ != nullptr) occ->assign(net.occurrence_edge_positions.size(), 0);
  for (const auto& it : net.entry_occurrences) {
    const Amount s = share[static_cast<std::size_t>(it.entry_position)];
    if (s == 0) continue;
    if (it.tunnel_position >= 0) {
      for (std::int32_t j = 0; j < it.num_occurrences; ++j) {
        const auto o = static_cast<std::size_t>(it.first_occurrence + j);
        load[static_cast<std::size_t>(net.occurrence_edge_positions[o])] += s;
        if (occ != nullptr) (*occ)[o] = s;
      }
      continue;
    }
    scratch.resize(it.prepared_route->forwarding_edges.size());
    (void)allocate(*it.prepared_route, s, scratch);
    for (std::int32_t j = 0; j < it.num_occurrences; ++j) {
      const auto o = static_cast<std::size_t>(it.first_occurrence + j);
      load[static_cast<std::size_t>(net.occurrence_edge_positions[o])] += scratch[static_cast<std::size_t>(j)];
      if (occ != nullptr) (*occ)[o] = scratch[static_cast<std::size_t>(j)];
    }
  }
}

// The router split of a split route's plan position p: round robin or highest averages over
// its out-members.
void split_at(const PreparedRoute& plan, std::size_t pos, Amount amount, std::span<Amount> out) {
  const auto b = static_cast<std::size_t>(plan.outgoing_begin[pos]);
  const auto e = static_cast<std::size_t>(plan.outgoing_end[pos]);
  std::fill(out.begin(), out.end(), 0);
  if (amount == 0 || e == b) return;
  if (plan.allocation == RouteAllocation::Weighted) {
    highest_averages(amount, std::span<const Amount>(plan.weights.data() + b, e - b), out);
  } else {
    (void)round_robin(amount, out);
  }
}

// The loss rule at one link: T = min(V, sum of the arrivals); highest averages over the sources'
// arrivals, then within a source over its entries'. `occs` lists the link's occurrences by
// (source, entry); writes their carried amounts.
void loss_rule(const ForwardingNetwork& net, std::span<const std::int32_t> occs, Amount V, std::span<const Amount> arrival,
               std::span<Amount> carried, std::vector<Amount>& wa, std::vector<Amount>& wb) {
  Amount sum = 0;
  for (auto o : occs) sum += arrival[static_cast<std::size_t>(o)];
  if (sum <= V) {
    for (auto o : occs) carried[static_cast<std::size_t>(o)] = arrival[static_cast<std::size_t>(o)];
    return;
  }
  // The source groups: consecutive runs of one head.
  wa.clear();
  std::vector<std::size_t> run;  // the first index of each source's run
  for (std::size_t i = 0; i < occs.size(); ++i) {
    const auto head = net.entry_occurrences[static_cast<std::size_t>(net.occurrence_entry_positions[static_cast<std::size_t>(occs[i])])].source_index;
    if (i == 0 || head != net.entry_occurrences[static_cast<std::size_t>(net.occurrence_entry_positions[static_cast<std::size_t>(occs[i - 1])])].source_index) {
      run.push_back(i);
      wa.push_back(0);
    }
    wa.back() += arrival[static_cast<std::size_t>(occs[i])];
  }
  run.push_back(occs.size());
  wb.resize(wa.size());
  highest_averages(V, wa, wb);
  std::vector<Amount> w;
  std::vector<Amount> c;
  for (std::size_t s = 0; s + 1 < run.size(); ++s) {
    const std::size_t b = run[s];
    const std::size_t e = run[s + 1];
    w.resize(e - b);
    c.resize(e - b);
    for (std::size_t i = b; i < e; ++i) w[i - b] = arrival[static_cast<std::size_t>(occs[i])];
    highest_averages(wb[s], w, c);
    for (std::size_t i = b; i < e; ++i) carried[static_cast<std::size_t>(occs[i])] = c[i - b];
  }
}

// ---------------------------------------------------------------- the general lossless bound

ProgramFlow general_lossless(const ForwardingView& view, const AvailableCapacityFn& free, Amount lo, std::optional<Amount> hi) {
  const ForwardingNetwork net = build_forwarding_network(view, free);
  const std::size_t k = view.num_sources();
  std::vector<Amount> heads(k);
  std::vector<Amount> share;
  std::vector<Amount> source_drops_by_source;
  std::vector<Amount> w;
  std::vector<Amount> part;
  std::vector<Amount> load;
  std::vector<Amount> scratch;
  auto fits = [&](Amount t) {
    source_shares(view, t, heads);
    // Units no head receives (every captured weight 0) have no route.
    Amount given = 0;
    for (Amount x : heads) given += x;
    if (given != t) return false;
    entry_shares(view, heads, share, source_drops_by_source, w, part);
    for (Amount x : source_drops_by_source) {
      if (x != 0) return false;
    }
    for (std::size_t e = 0; e < share.size(); ++e) {
      if (share[e] != 0 && net.entry_occurrence_positions[e] == -2) return false;  // a share with no route
    }
    no_drop(net, share, load, nullptr, scratch);
    for (std::size_t l = 0; l < load.size(); ++l) {
      if (load[l] > net.available_capacity[l]) return false;
    }
    return true;
  };
  // An upper bound: every unit leaves on an entry's first link, which carries it without a drop;
  // with no head, nothing enters.
  Amount top = 0;
  {
    std::vector<std::uint8_t> first(net.edges.size(), 0);
    for (const auto& it : net.entry_occurrences) {
      if (it.tunnel_position >= 0) {
        first[static_cast<std::size_t>(net.occurrence_edge_positions[static_cast<std::size_t>(it.first_occurrence)])] = 1;
      } else {
        const auto rp = static_cast<std::size_t>(it.source_position);
        for (auto m = it.prepared_route->outgoing_begin[rp]; m < it.prepared_route->outgoing_end[rp]; ++m) {
          first[static_cast<std::size_t>(net.occurrence_edge_positions[static_cast<std::size_t>(it.first_occurrence + m)])] = 1;
        }
      }
    }
    for (std::size_t l = 0; l < first.size(); ++l) {
      if (first[l] != 0) top += net.available_capacity[l];
    }
  }
  Amount upper = hi.has_value() ? std::min(*hi, std::max(top, lo)) : std::max(top, lo);
  Amount best = lo;
  if (upper > lo) {
    if (fits(upper)) {
      best = upper;
    } else {
      // The largest fitting total in [lo, upper - 1]: lo fits and upper does not (a prefix).
      Amount a = lo;
      Amount b = upper;
      while (b - a > 1) {
        const Amount mid = a + (b - a) / 2;
        if (fits(mid)) {
          a = mid;
        } else {
          b = mid;
        }
      }
      best = a;
    }
  }
  // The loads of the admitted total.
  ProgramFlow out = empty_flow(view, false);
  out.admitted = best;
  source_shares(view, best, out.offered_by_source);
  entry_shares(view, out.offered_by_source, share, source_drops_by_source, w, part);
  std::vector<Amount> occ;
  no_drop(net, share, load, &occ, scratch);
  for (std::size_t l = 0; l < load.size(); ++l) {
    if (load[l] > net.available_capacity[l]) throw std::logic_error("forwarding: the admitted total does not fit");
  }
  for (Amount x : out.offered_by_source) out.offered += x;
  out.delivered = out.offered;
  for (const auto& it : net.entry_occurrences) {
    auto& f = out.entries[static_cast<std::size_t>(it.entry_position)];
    f.assigned = share[static_cast<std::size_t>(it.entry_position)];
    f.delivered = f.assigned;
    for (std::int32_t j = 0; j < it.num_occurrences; ++j) {
      const auto o = static_cast<std::size_t>(it.first_occurrence + j);
      if (occ[o] != 0) f.edge_flows.emplace_back(net.edges[static_cast<std::size_t>(net.occurrence_edge_positions[o])], occ[o]);
    }
    if (it.tunnel_position >= 0) sort_pairs(f.edge_flows);
  }
  return out;
}

// ---------------------------------------------------------------- the general loss rule (A+)

class LossyForwardingEvaluator {
public:
  LossyForwardingEvaluator(const ForwardingView& view, const ForwardingNetwork& net) : view_(view), net_(net) {}

  ProgramFlow run(std::span<const Amount> injections);

private:
  // Event ids: barriers [0, nb), routers [nb, nb + num_pos).
  [[nodiscard]] std::int32_t router(const ForwardingNetwork::Item& it, std::int32_t pos) const noexcept {
    return num_barriers_ + it.first_router_position + pos;
  }
  [[nodiscard]] std::int32_t event_of_occ(std::int32_t o) const noexcept {
    return barrier_[static_cast<std::size_t>(net_.occurrence_edge_positions[static_cast<std::size_t>(o)])];
  }
  // A push into event `ev` is suppressed while a group iterates and ev lies outside it.
  [[nodiscard]] bool blocked(std::int32_t ev) const noexcept { return iterating_ && scc_[static_cast<std::size_t>(ev)] != current_component_; }

  void normalize_serial_links(std::span<const Amount> U0);
  void build_events();
  void find_strong_components();
  void carry(std::int32_t o, Amount c);
  void split(std::int32_t ev);
  void evaluate_barrier(std::int32_t ev);
  void iterate(std::int32_t scc);
  void verify(std::span<const Amount> share) const;

  const ForwardingView& view_;
  const ForwardingNetwork& net_;
  std::vector<std::int32_t> occurrence_source_positions_;         // per occurrence of a split route: its member's tail position
  std::vector<std::uint8_t> transparent_;  // per link: the second link of a transparent pair
  std::vector<std::int32_t> barrier_;      // per link: its event, or -1 (a pass-through)
  std::vector<std::int32_t> barrier_link_;
  std::int32_t num_barriers_ {0};
  std::int32_t num_events_ {0};
  std::vector<std::int32_t> dependency_offsets_;
  std::vector<std::int32_t> dependencies_;
  std::vector<std::int32_t> scc_;
  std::int32_t num_scc_ {0};
  bool iterating_ {false};
  std::int32_t current_component_ {-1};
  std::vector<Amount> arrival_;
  std::vector<Amount> carried_;
  std::vector<Amount> inflow_;    // per router position
  std::vector<Amount> delivered_; // per item
  std::vector<Amount> outgoing_shares_;
  std::vector<Amount> source_arrivals_;
  std::vector<Amount> source_carries_;
  bool iterated_ {false};
};

void LossyForwardingEvaluator::normalize_serial_links(std::span<const Amount> U0) {
  const std::size_t L = net_.edges.size();
  transparent_.assign(L, 0);
  // Per link: the one link every occurrence continues on, and the one every occurrence comes from
  // (-1: none yet; -2: they differ, or an occurrence is injected, absorbed, merged or split there).
  std::vector<std::int32_t> next(L, -1);
  std::vector<std::int32_t> prev(L, -1);
  auto merge = [](std::int32_t& slot, std::int32_t v) {
    if (slot == -1) {
      slot = v;
    } else if (slot != v) {
      slot = -2;
    }
  };
  for (std::size_t o = 0; o < net_.occurrence_edge_positions.size(); ++o) {
    const auto l = static_cast<std::size_t>(net_.occurrence_edge_positions[o]);
    const auto& it = net_.entry_occurrences[static_cast<std::size_t>(net_.occurrence_entry_positions[o])];
    const std::int32_t j = net_.occurrence_route_positions[o];
    std::int32_t succ = -2;
    std::int32_t pred = -2;
    if (it.tunnel_position >= 0) {
      if (j + 1 < it.num_occurrences) succ = net_.occurrence_edge_positions[static_cast<std::size_t>(it.first_occurrence + j + 1)];
      if (j > 0) pred = net_.occurrence_edge_positions[static_cast<std::size_t>(it.first_occurrence + j - 1)];
    } else {
      const auto& plan = *it.prepared_route;
      // Continues: its head is no destination and no root, holds this member alone as its
      // in-member and has one out-member.
      const auto q = static_cast<std::size_t>(plan.edge_destination_positions[static_cast<std::size_t>(j)]);
      if (plan.is_destination[q] == 0 && static_cast<std::int32_t>(q) != it.source_position &&
          plan.outgoing_end[q] - plan.outgoing_begin[q] == 1 && net_.incoming_edge_positions[static_cast<std::size_t>(it.first_router_position) + q] == j) {
        succ = net_.occurrence_edge_positions[static_cast<std::size_t>(it.first_occurrence + plan.outgoing_begin[q])];
      }
      // Comes from: its tail is no root, has one in-member and this one out-member.
      const auto t = static_cast<std::size_t>(occurrence_source_positions_[static_cast<std::size_t>(it.first_occurrence + j)]);
      const auto in = net_.incoming_edge_positions[static_cast<std::size_t>(it.first_router_position) + t];
      if (static_cast<std::int32_t>(t) != it.source_position && in >= 0 && plan.outgoing_end[t] - plan.outgoing_begin[t] == 1) {
        pred = net_.occurrence_edge_positions[static_cast<std::size_t>(it.first_occurrence + in)];
      }
    }
    merge(next[l], succ);
    merge(prev[l], pred);
  }
  for (std::size_t l = 0; l < L; ++l) {
    const auto n = next[l];
    if (n < 0) continue;
    const auto ni = static_cast<std::size_t>(n);
    if (prev[ni] != static_cast<std::int32_t>(l) || net_.available_capacity[l] != net_.available_capacity[ni]) continue;
    transparent_[ni] = 1;
  }
  barrier_.assign(L, -1);
  barrier_link_.clear();
  for (std::size_t l = 0; l < L; ++l) {
    if (transparent_[l] == 0 && U0[l] > net_.available_capacity[l]) {
      barrier_[l] = static_cast<std::int32_t>(barrier_link_.size());
      barrier_link_.push_back(static_cast<std::int32_t>(l));
    }
  }
  num_barriers_ = static_cast<std::int32_t>(barrier_link_.size());
  // The events are known only after the no-drop loads: checked before the event structures are
  // sized (the network is built; a refusal changes nothing).
  check_forwarding_count(ForwardingCounter::Event, std::int64_t {num_barriers_} + net_.num_router_positions);
  num_events_ = num_barriers_ + net_.num_router_positions;
}

void LossyForwardingEvaluator::build_events() {
  // Every dependency arc, to `emit`: an event before the next one on its tunnel; on a split route
  // a member's tail router before its event (a member that can drop) and that event before the
  // head router, or the tail before the head.
  auto arcs = [&](auto&& emit) {
    for (const auto& it : net_.entry_occurrences) {
      if (it.tunnel_position >= 0) {
        std::int32_t last = -1;
        for (std::int32_t j = 0; j < it.num_occurrences; ++j) {
          const auto ev = event_of_occ(it.first_occurrence + j);
          if (ev < 0) continue;
          if (last >= 0) emit(last, ev);
          last = ev;
        }
        continue;
      }
      const auto& plan = *it.prepared_route;
      for (std::size_t pos = 0; pos < plan.topological_order.size(); ++pos) {
        for (auto m = plan.outgoing_begin[pos]; m < plan.outgoing_end[pos]; ++m) {
          const auto o = it.first_occurrence + m;
          const auto from = router(it, static_cast<std::int32_t>(pos));
          const auto to = router(it, plan.edge_destination_positions[static_cast<std::size_t>(m)]);
          const auto ev = event_of_occ(o);
          if (ev >= 0) {
            emit(from, ev);
            emit(ev, to);
          } else {
            emit(from, to);
          }
        }
      }
    }
  };
  std::int64_t count = 0;
  arcs([&](std::int32_t, std::int32_t) { ++count; });
  check_forwarding_count(ForwardingCounter::Dependency, count);
  std::vector<std::pair<std::int32_t, std::int32_t>> edges;
  edges.reserve(static_cast<std::size_t>(count));
  arcs([&](std::int32_t from, std::int32_t to) { edges.emplace_back(from, to); });
  dependency_offsets_.assign(static_cast<std::size_t>(num_events_) + 1, 0);
  for (const auto& [a, b] : edges) ++dependency_offsets_[static_cast<std::size_t>(a) + 1];
  for (std::size_t v = 0; v < static_cast<std::size_t>(num_events_); ++v) dependency_offsets_[v + 1] += dependency_offsets_[v];
  dependencies_.resize(edges.size());
  std::vector<std::int32_t> fill(dependency_offsets_.begin(), dependency_offsets_.end() - 1);
  for (const auto& [a, b] : edges) dependencies_[static_cast<std::size_t>(fill[static_cast<std::size_t>(a)]++)] = b;
}

// Tarjan's strongly connected components, iterative (no recursion depth).
void LossyForwardingEvaluator::find_strong_components() {
  const auto n = static_cast<std::size_t>(num_events_);
  scc_.assign(n, -1);
  std::vector<std::int32_t> index(n, -1);
  std::vector<std::int32_t> low(n, 0);
  std::vector<std::uint8_t> on(n, 0);
  std::vector<std::int32_t> stack;
  std::vector<std::pair<std::int32_t, std::int32_t>> work;  // (event, next adjacency position)
  std::int32_t counter = 0;
  num_scc_ = 0;
  for (std::size_t s = 0; s < n; ++s) {
    if (index[s] >= 0) continue;
    work.emplace_back(static_cast<std::int32_t>(s), dependency_offsets_[s]);
    index[s] = low[s] = counter++;
    stack.push_back(static_cast<std::int32_t>(s));
    on[s] = 1;
    while (!work.empty()) {
      auto& [v, k] = work.back();
      const auto vi = static_cast<std::size_t>(v);
      if (k < dependency_offsets_[vi + 1]) {
        const auto w = static_cast<std::size_t>(dependencies_[static_cast<std::size_t>(k)]);
        ++k;
        if (index[w] < 0) {
          index[w] = low[w] = counter++;
          stack.push_back(static_cast<std::int32_t>(w));
          on[w] = 1;
          work.emplace_back(static_cast<std::int32_t>(w), dependency_offsets_[w]);
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
          scc_[w] = num_scc_;
          if (w == vi) break;
        }
        ++num_scc_;
      }
      const auto done = v;
      work.pop_back();
      if (!work.empty()) {
        const auto u = static_cast<std::size_t>(work.back().first);
        low[u] = std::min(low[u], low[static_cast<std::size_t>(done)]);
      }
    }
  }
}

// Occurrence o now carries c: forward it to what it feeds (a split route's head router, a
// tunnel's next hops up to its next barrier or its end).
void LossyForwardingEvaluator::carry(std::int32_t o, Amount c) {
  carried_[static_cast<std::size_t>(o)] = c;
  const auto& it = net_.entry_occurrences[static_cast<std::size_t>(net_.occurrence_entry_positions[static_cast<std::size_t>(o)])];
  const std::int32_t j = net_.occurrence_route_positions[static_cast<std::size_t>(o)];
  if (it.tunnel_position < 0) {
    const auto q = it.prepared_route->edge_destination_positions[static_cast<std::size_t>(j)];
    const auto ev = router(it, q);
    if (!blocked(ev)) inflow_[static_cast<std::size_t>(it.first_router_position + q)] += c;
    return;
  }
  for (std::int32_t h = j + 1;; ++h) {
    if (h == it.num_occurrences) {
      // The tunnel's end is no event of a group.
      if (!iterating_) delivered_[static_cast<std::size_t>(net_.occurrence_entry_positions[static_cast<std::size_t>(o)])] = c;
      return;
    }
    const auto o2 = it.first_occurrence + h;
    const auto ev = event_of_occ(o2);
    if (ev >= 0) {
      if (!blocked(ev)) arrival_[static_cast<std::size_t>(o2)] = c;
      return;
    }
    arrival_[static_cast<std::size_t>(o2)] = c;
    carried_[static_cast<std::size_t>(o2)] = c;
  }
}

// A router splits what it holds once over its out-members; a pass-through member carries its
// whole arrival on at once.
void LossyForwardingEvaluator::split(std::int32_t ev) {
  const auto pos_global = ev - num_barriers_;
  const auto& it = net_.entry_occurrences[static_cast<std::size_t>(net_.router_entry_positions[static_cast<std::size_t>(pos_global)])];
  const auto& plan = *it.prepared_route;
  const auto pos = static_cast<std::size_t>(pos_global - it.first_router_position);
  if (plan.is_destination[pos] != 0) return;
  const auto b = static_cast<std::size_t>(plan.outgoing_begin[pos]);
  const auto e = static_cast<std::size_t>(plan.outgoing_end[pos]);
  outgoing_shares_.resize(e - b);
  split_at(plan, pos, inflow_[static_cast<std::size_t>(pos_global)], outgoing_shares_);
  for (std::size_t m = b; m < e; ++m) {
    const auto o = it.first_occurrence + static_cast<std::int32_t>(m);
    const Amount x = outgoing_shares_[m - b];
    const auto bev = event_of_occ(o);
    if (bev >= 0) {
      if (!blocked(bev)) arrival_[static_cast<std::size_t>(o)] = x;
      continue;
    }
    arrival_[static_cast<std::size_t>(o)] = x;
    carry(o, x);
  }
}

void LossyForwardingEvaluator::evaluate_barrier(std::int32_t ev) {
  const auto l = static_cast<std::size_t>(barrier_link_[static_cast<std::size_t>(ev)]);
  const auto occs = std::span<const std::int32_t>(net_.occurrences_by_edge).subspan(
      static_cast<std::size_t>(net_.edge_occurrence_offsets[l]), static_cast<std::size_t>(net_.edge_occurrence_offsets[l + 1] - net_.edge_occurrence_offsets[l]));
  loss_rule(net_, occs, net_.available_capacity[l], arrival_, carried_, source_arrivals_, source_carries_);
  for (auto o : occs) carry(o, carried_[static_cast<std::size_t>(o)]);
}

// A cyclic group: rounds from zero carried at its barriers.
void LossyForwardingEvaluator::iterate(std::int32_t scc) {
  iterated_ = true;
  current_component_ = scc;
  std::vector<std::int32_t> barriers;
  std::vector<std::int32_t> routers;
  for (std::int32_t ev = 0; ev < num_events_; ++ev) {
    if (scc_[static_cast<std::size_t>(ev)] != scc) continue;
    (ev < num_barriers_ ? barriers : routers).push_back(ev);
  }
  // The group's routers in an order that respects router -> router edges (barriers cut).
  {
    std::vector<std::int32_t> indeg(static_cast<std::size_t>(num_events_), 0);
    for (auto r : routers) {
      for (auto k = dependency_offsets_[static_cast<std::size_t>(r)]; k < dependency_offsets_[static_cast<std::size_t>(r) + 1]; ++k) {
        const auto w = dependencies_[static_cast<std::size_t>(k)];
        if (w >= num_barriers_ && scc_[static_cast<std::size_t>(w)] == scc) ++indeg[static_cast<std::size_t>(w)];
      }
    }
    std::vector<std::int32_t> order;
    std::vector<std::int32_t> ready;
    for (auto r : routers) {
      if (indeg[static_cast<std::size_t>(r)] == 0) ready.push_back(r);
    }
    while (!ready.empty()) {
      const auto r = ready.back();
      ready.pop_back();
      order.push_back(r);
      for (auto k = dependency_offsets_[static_cast<std::size_t>(r)]; k < dependency_offsets_[static_cast<std::size_t>(r) + 1]; ++k) {
        const auto w = dependencies_[static_cast<std::size_t>(k)];
        if (w >= num_barriers_ && scc_[static_cast<std::size_t>(w)] == scc && --indeg[static_cast<std::size_t>(w)] == 0) ready.push_back(w);
      }
    }
    routers.swap(order);
  }
  // The state x: every carried amount at the group's barriers, by EdgeId then source then entry
  // (the barriers ascend by EdgeId; a link's occurrences are in that order).
  std::vector<std::int32_t> occs;
  for (auto b : barriers) {
    const auto l = static_cast<std::size_t>(barrier_link_[static_cast<std::size_t>(b)]);
    for (auto i = net_.edge_occurrence_offsets[l]; i < net_.edge_occurrence_offsets[l + 1]; ++i) occs.push_back(net_.occurrences_by_edge[static_cast<std::size_t>(i)]);
  }
  std::vector<Amount> inflow0(routers.size());
  for (std::size_t i = 0; i < routers.size(); ++i) inflow0[i] = inflow_[static_cast<std::size_t>(routers[i] - num_barriers_)];
  std::vector<Amount> x(occs.size(), 0);
  std::vector<Amount> next(occs.size(), 0);
  std::vector<Amount> before(occs.size(), 0);  // x(r - 2)
  const auto limit = static_cast<std::int64_t>(64) * static_cast<std::int64_t>(barriers.size());
  auto replay = [&]() {
    for (std::size_t i = 0; i < routers.size(); ++i) inflow_[static_cast<std::size_t>(routers[i] - num_barriers_)] = inflow0[i];
    for (std::size_t i = 0; i < occs.size(); ++i) carry(occs[i], x[i]);
    for (auto r : routers) split(r);
  };
  iterating_ = true;
  bool solved = false;
  for (std::int64_t round = 1; round <= limit; ++round) {
#if NETGRAPH_CORE_PRIVATE_HOOKS
    ++thread_lossy_counters().rounds;
#endif
    replay();
    // Every new carry into a separate state.
    for (auto b : barriers) {
      const auto l = static_cast<std::size_t>(barrier_link_[static_cast<std::size_t>(b)]);
      const auto span = std::span<const std::int32_t>(net_.occurrences_by_edge).subspan(
          static_cast<std::size_t>(net_.edge_occurrence_offsets[l]), static_cast<std::size_t>(net_.edge_occurrence_offsets[l + 1] - net_.edge_occurrence_offsets[l]));
      loss_rule(net_, span, net_.available_capacity[l], arrival_, carried_, source_arrivals_, source_carries_);
    }
    for (std::size_t i = 0; i < occs.size(); ++i) next[i] = carried_[static_cast<std::size_t>(occs[i])];
    if (next == x) {
      solved = true;
      break;
    }
    if (round >= 2 && next == before) {
      iterating_ = false;
      throw std::invalid_argument("lossy forwarding iteration entered a period-2 cycle");
    }
    before.swap(x);
    x.swap(next);
  }
  if (!solved) {
    iterating_ = false;
    throw std::invalid_argument("lossy forwarding iteration did not reach a fixed point within the round limit");
  }
  // The solved state's forwarding, now also into the group's successors.
  iterating_ = false;
  current_component_ = -1;
  for (std::size_t i = 0; i < routers.size(); ++i) inflow_[static_cast<std::size_t>(routers[i] - num_barriers_)] = inflow0[i];
  for (std::size_t i = 0; i < occs.size(); ++i) carry(occs[i], x[i]);
  for (auto r : routers) split(r);
}

// The independent check of an iterated result (docs/semantics.md, "Losses inside one program"):
// every entry's injection, every router's split of what reaches it, every link's loss rule and
// capacity, every continuation, and each entry's conservation (delivered + drops = its share).
void LossyForwardingEvaluator::verify(std::span<const Amount> share) const {
  const auto fail = [] { throw std::runtime_error("lossy forwarding fixed-point verification failed"); };
  // Continuations and splits.
  std::vector<Amount> in(static_cast<std::size_t>(net_.num_router_positions), 0);
  std::vector<Amount> part;
  for (std::size_t i = 0; i < net_.entry_occurrences.size(); ++i) {
    const auto& it = net_.entry_occurrences[i];
    const Amount s = share[static_cast<std::size_t>(it.entry_position)];
    if (it.tunnel_position >= 0) {
      Amount prev = s;
      Amount dropped = 0;
      for (std::int32_t j = 0; j < it.num_occurrences; ++j) {
        const auto o = static_cast<std::size_t>(it.first_occurrence + j);
        if (arrival_[o] != prev || carried_[o] > arrival_[o] || carried_[o] < 0) fail();
        dropped += arrival_[o] - carried_[o];
        prev = carried_[o];
      }
      if (delivered_[i] != prev || delivered_[i] + dropped != s) fail();
      continue;
    }
    const auto& plan = *it.prepared_route;
    in[static_cast<std::size_t>(it.first_router_position + it.source_position)] += s;
    Amount dropped = 0;
    Amount reached = 0;
    for (std::size_t pos = 0; pos < plan.topological_order.size(); ++pos) {
      const Amount amount = in[static_cast<std::size_t>(it.first_router_position) + pos];
      if (plan.is_destination[pos] != 0) {
        reached += amount;
        continue;
      }
      const auto b = static_cast<std::size_t>(plan.outgoing_begin[pos]);
      const auto e = static_cast<std::size_t>(plan.outgoing_end[pos]);
      part.resize(e - b);
      split_at(plan, pos, amount, part);
      for (std::size_t m = b; m < e; ++m) {
        const auto o = static_cast<std::size_t>(it.first_occurrence) + m;
        if (arrival_[o] != part[m - b] || carried_[o] > arrival_[o] || carried_[o] < 0) fail();
        dropped += arrival_[o] - carried_[o];
        in[static_cast<std::size_t>(it.first_router_position + plan.edge_destination_positions[m])] += carried_[o];
      }
    }
    if (delivered_[i] != reached || reached + dropped != s) fail();
  }
  // Every link: its loss rule over the arrivals (a transparent link passes everything on) and
  // its capacity.
  std::vector<Amount> expect(carried_.size(), 0);
  std::vector<Amount> wa;
  std::vector<Amount> wb;
  for (std::size_t l = 0; l < net_.edges.size(); ++l) {
    const auto occs = std::span<const std::int32_t>(net_.occurrences_by_edge).subspan(
        static_cast<std::size_t>(net_.edge_occurrence_offsets[l]), static_cast<std::size_t>(net_.edge_occurrence_offsets[l + 1] - net_.edge_occurrence_offsets[l]));
    Amount sum = 0;
    for (auto o : occs) sum += carried_[static_cast<std::size_t>(o)];
    if (sum > net_.available_capacity[l]) fail();
    if (transparent_[l] != 0) {
      for (auto o : occs) {
        if (carried_[static_cast<std::size_t>(o)] != arrival_[static_cast<std::size_t>(o)]) fail();
      }
      continue;
    }
    loss_rule(net_, occs, net_.available_capacity[l], arrival_, expect, wa, wb);
    for (auto o : occs) {
      if (expect[static_cast<std::size_t>(o)] != carried_[static_cast<std::size_t>(o)]) fail();
    }
  }
}

ProgramFlow LossyForwardingEvaluator::run(std::span<const Amount> injections) {
  const std::size_t k = view_.num_sources();
  ProgramFlow out = empty_flow(view_, true);
  for (std::size_t h = 0; h < k; ++h) {
    out.offered_by_source[h] = injections[h];
    out.offered += injections[h];
  }
  out.admitted = out.offered;
  std::vector<Amount> share;
  std::vector<Amount> source_drops_by_source;
  std::vector<Amount> w;
  std::vector<Amount> part;
  entry_shares(view_, injections, share, source_drops_by_source, w, part);
  for (std::size_t h = 0; h < k; ++h) {
    out.source_drops_by_source[h] = source_drops_by_source[h];
    out.source_drops += source_drops_by_source[h];
  }
  for (std::size_t e = 0; e < share.size(); ++e) {
    out.entries[e].assigned = share[e];
    if (share[e] != 0 && net_.entry_occurrence_positions[e] == -2) {
      out.entries[e].source_drops = share[e];
      out.source_drops += share[e];
    }
  }
  // A+: the no-drop loads first, then the barriers.
  std::vector<Amount> U0;
  std::vector<Amount> scratch;
  no_drop(net_, share, U0, nullptr, scratch);
  occurrence_source_positions_.assign(net_.occurrence_edge_positions.size(), -1);
  for (const auto& it : net_.entry_occurrences) {
    if (it.tunnel_position >= 0) continue;
    for (std::size_t pos = 0; pos < it.prepared_route->topological_order.size(); ++pos) {
      for (auto m = it.prepared_route->outgoing_begin[pos]; m < it.prepared_route->outgoing_end[pos]; ++m) occurrence_source_positions_[static_cast<std::size_t>(it.first_occurrence + m)] = static_cast<std::int32_t>(pos);
    }
  }
  normalize_serial_links(U0);
  const auto n = net_.occurrence_edge_positions.size();
  arrival_.assign(n, 0);
  carried_.assign(n, 0);
  inflow_.assign(static_cast<std::size_t>(net_.num_router_positions), 0);
  delivered_.assign(net_.entry_occurrences.size(), 0);
  build_events();
  find_strong_components();
  // The injections: a tunnel's share enters its first hop, a split route's its root router.
  for (std::size_t i = 0; i < net_.entry_occurrences.size(); ++i) {
    const auto& it = net_.entry_occurrences[i];
    const Amount s = share[static_cast<std::size_t>(it.entry_position)];
    if (it.tunnel_position >= 0) {
      const auto o = it.first_occurrence;
      arrival_[static_cast<std::size_t>(o)] = s;
      if (event_of_occ(o) < 0) carry(o, s);
    } else {
      inflow_[static_cast<std::size_t>(it.first_router_position + it.source_position)] += s;
    }
  }
  // The groups in dependency order, the smallest key first among the ready ones.
  const auto ns = static_cast<std::size_t>(num_scc_);
  std::vector<std::vector<std::int32_t>> members(ns);
  for (std::int32_t ev = 0; ev < num_events_; ++ev) members[static_cast<std::size_t>(scc_[static_cast<std::size_t>(ev)])].push_back(ev);
  std::vector<std::int32_t> indeg(ns, 0);
  std::vector<std::uint8_t> cyclic(ns, 0);
  for (std::int32_t ev = 0; ev < num_events_; ++ev) {
    const auto a = scc_[static_cast<std::size_t>(ev)];
    for (auto kk = dependency_offsets_[static_cast<std::size_t>(ev)]; kk < dependency_offsets_[static_cast<std::size_t>(ev) + 1]; ++kk) {
      const auto b = scc_[static_cast<std::size_t>(dependencies_[static_cast<std::size_t>(kk)])];
      if (a == b) {
        cyclic[static_cast<std::size_t>(a)] = 1;
      } else {
        ++indeg[static_cast<std::size_t>(b)];
      }
    }
  }
  using Key = std::tuple<int, std::int64_t, std::int32_t>;
  auto key = [&](std::int32_t s) -> Key {
    const auto& m = members[static_cast<std::size_t>(s)];
    std::int64_t edge = std::numeric_limits<std::int64_t>::max();
    for (auto ev : m) {
      if (ev < num_barriers_) edge = std::min<std::int64_t>(edge, net_.edges[static_cast<std::size_t>(barrier_link_[static_cast<std::size_t>(ev)])]);
    }
    if (edge != std::numeric_limits<std::int64_t>::max()) return {0, edge, s};
    return {1, m.front(), s};
  };
  std::priority_queue<Key, std::vector<Key>, std::greater<>> ready;
  for (std::size_t s = 0; s < ns; ++s) {
    if (indeg[s] == 0) ready.push(key(static_cast<std::int32_t>(s)));
  }
  while (!ready.empty()) {
    const auto s = std::get<2>(ready.top());
    ready.pop();
    const auto& m = members[static_cast<std::size_t>(s)];
    if (cyclic[static_cast<std::size_t>(s)] != 0) {
#if NETGRAPH_CORE_PRIVATE_HOOKS
      ++thread_lossy_counters().groups;
#endif
      iterate(s);
    } else if (m.front() < num_barriers_) {
      evaluate_barrier(m.front());
    } else {
      split(m.front());
    }
    for (auto ev : m) {
      for (auto kk = dependency_offsets_[static_cast<std::size_t>(ev)]; kk < dependency_offsets_[static_cast<std::size_t>(ev) + 1]; ++kk) {
        const auto b = scc_[static_cast<std::size_t>(dependencies_[static_cast<std::size_t>(kk)])];
        if (b != s && --indeg[static_cast<std::size_t>(b)] == 0) ready.push(key(b));
      }
    }
  }
  // Split routes deliver what reaches their destinations.
  for (std::size_t i = 0; i < net_.entry_occurrences.size(); ++i) {
    const auto& it = net_.entry_occurrences[i];
    if (it.tunnel_position >= 0) continue;
    Amount d = 0;
    for (std::size_t pos = 0; pos < it.prepared_route->topological_order.size(); ++pos) {
      if (it.prepared_route->is_destination[pos] != 0) d += inflow_[static_cast<std::size_t>(it.first_router_position) + pos];
    }
    delivered_[i] = d;
  }
  if (iterated_) {
#if NETGRAPH_CORE_PRIVATE_HOOKS
    if (thread_lossy_counters().corrupt && !carried_.empty()) {
      // The planted fault of the check's test: one carried amount changed after the rounds.
      for (auto& c : carried_) {
        if (c > 0) {
          --c;
          break;
        }
      }
    }
#endif
    verify(share);
  }
  for (std::size_t i = 0; i < net_.entry_occurrences.size(); ++i) {
    const auto& it = net_.entry_occurrences[i];
    auto& f = out.entries[static_cast<std::size_t>(it.entry_position)];
    f.delivered = delivered_[i];
    out.delivered += delivered_[i];
    for (std::int32_t j = 0; j < it.num_occurrences; ++j) {
      const auto o = static_cast<std::size_t>(it.first_occurrence + j);
      const EdgeId e = net_.edges[static_cast<std::size_t>(net_.occurrence_edge_positions[o])];
      if (carried_[o] != 0) f.edge_flows.emplace_back(e, carried_[o]);
      if (arrival_[o] != carried_[o]) f.edge_drops.emplace_back(e, arrival_[o] - carried_[o]);
    }
    if (it.tunnel_position >= 0) {
      sort_pairs(f.edge_flows);
      sort_pairs(f.edge_drops);
    }
  }
  return out;
}

} // namespace

ProgramFlow lossless(const ForwardingView& view, const AvailableCapacityFn& free, Amount lo, std::optional<Amount> hi) {
  if (view.pooled) return pooled_lossless(view, free, lo, hi);
  return general_lossless(view, free, lo, hi);
}

ProgramFlow lossy(const ForwardingView& view, const AvailableCapacityFn& free, std::span<const Amount> injections) {
  if (injections.size() != view.num_sources()) throw std::logic_error("forwarding: the injections are not head-aligned");
#if NETGRAPH_CORE_PRIVATE_HOOKS
  thread_lossy_counters().groups = 0;
  thread_lossy_counters().rounds = 0;
#endif
  if (view.pooled) return pooled_lossy(view, free, injections);
  const ForwardingNetwork net = build_forwarding_network(view, free);
  LossyForwardingEvaluator run(view, net);
  return run.run(injections);
}

} // namespace netgraph::core::detail
