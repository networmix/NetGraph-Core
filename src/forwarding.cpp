/*
  The forwarding engine of programs (forwarding.hpp): views, the injection, the lossless bound
  and the loss rule.

  A pooled program (one split context) called alone runs on its plan's kernels: allocate and the
  T_max search for the bound, forward() with its attribution for the losses. Every other program
  (tunnels and supplied split routes, each split route a one-root context), and every demand of a
  SIMULTANEOUS list (lossy() over ForwardingDemand values, pooled ones included), runs on the
  general network of one call: the union of the demands' live entries' links, each link with its
  occurrences (a tunnel's hop, a split route's member; a pooled context's member once per source
  whose traffic reaches it), and per split route its routers (a pooled router keeps one inflow per
  source and attributes each member's amount to them by highest averages over their remaining
  inflows). Each occurrence belongs to one demand, one source and one entry, so the loss rule
  at a link is three levels of highest averages: demands (in creation order), sources, entries. A
  demand's black holes make its occurrences there dark: what reaches one is dropped at once; a
  dark occurrence shares, loads and waits on no link (Lossless reads a black-holed link as serving
  0 instead).

  The evaluation (A+): the loads U0 first, every demand forwarded with no link shared and its
  traffic stopped at its own black holes; an occurrence that carries 0 there (a share of 0, a member
  whose whole-unit share is 0, traffic behind its own black hole) is dead and takes no part in what
  follows. Serial normalization marks the second link of every transparent pair (equal free
  capacity, every live occurrence of the first continuing on the second and every live occurrence
  of the second coming from the first, with nothing injected, absorbed, merged or split between
  them, counting only the routers' live members) as a pass-through; every other link with U0 > V
  is a barrier. The events are the barriers and the routers; an event waits for its predecessors
  (a router for its in-members' carried amounts, a barrier for all its live arrivals, a tunnel's
  next barrier for its previous one). Tarjan's algorithm groups them into strongly connected
  components (SCCs), which run in dependency order, the smallest key first among the ready ones
  ((0, smallest barrier EdgeId) for a group with a barrier, else (1, its first event)). A cyclic
  group of b barriers is iterated from zero carried: one round replays the group's private
  forwarding from the previous carried state x, then computes every new carry with the loss rule
  into a separate state; equal states solve it, a state equal to the one two rounds before raises
  the period-two error, and 64 b rounds without either raise the round-limit error. Once every
  group is solved, a call that iterated checks the whole result independently (every router's
  split, every link's rule and capacity, every continuation, the conservation of every entry)
  before anything is returned.

  Amounts are int64 units below 2^62: a demand's offer bounds its every arrival and load. A link's
  sum over several demands is only compared with a free capacity below 2^62, so it saturates at
  2^62. The per-thread scratch of the general network is rebuilt by every call before it is read.
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
#include "netgraph/core/detail/sparse.hpp"
#include "netgraph/core/detail/route_builders.hpp"
#include "netgraph/core/detail/strong_components.hpp"

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
      f.edge_flows = nonzero_pairs(plan.forwarding_edges, alloc);
      if (!drops.empty()) f.edge_drops = nonzero_pairs(plan.forwarding_edges, drops);
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

// A link sum over several demands can pass 2^63 (each demand's arrivals at a link stay below its
// offer, below 2^62); every such sum is only compared with a free capacity below 2^62, so it
// saturates at 2^62, exact for every comparison.
constexpr Amount kSumCap = Amount {1} << 62;
Amount saturating_add(Amount a, Amount b) noexcept { return b >= kSumCap - a ? kSumCap : a + b; }

// The live entries of one forwarding of one or more demands (each surviving, of positive weight,
// with a usable route) and the union of their links. An item is a tunnel (its occurrences are its
// hops in path order), a one-root split route (its plan's members), or a pooled context (one item
// for the whole context: an occurrence per member and per head whose traffic can reach that
// member, by member then head, and per router position that is no destination a slot per head
// reaching it, which holds that head's inflow there). Items are ordered by (demand, head, entry
// index), so each link's occurrences come out in the loss rule's order. An occurrence on a link
// its demand's black holes exclude is dark (LOSSY only): it belongs to no link's list, so it never
// shares, loads or blocks that link; what reaches it is dropped there. LOSSLESS keeps no dark
// occurrence: a black-holed link serves 0 (the lossless bound).
struct ForwardingNetwork {
  struct Item {
    std::int32_t demand {0};
    std::int32_t entry_position {-1};  // its program entry (a pooled context: -1, one entry per head)
    std::int32_t source_index {0};     // view head (a pooled context: its first)
    std::int32_t tunnel_position {-1};
    const PreparedRoute* prepared_route {nullptr};
    std::int32_t source_position {-1};  // a one-root split route: its root's plan position
    bool pooled {false};
    std::int32_t first_occurrence {0};      // its first occurrence
    std::int32_t num_occurrences {0};
    std::int32_t first_router_position {0};  // its first router (split routes: one per plan position)
    std::size_t first_member {0};    // pooled: its members' entries in member_occurrence_offsets (members + 1)
    std::size_t first_position {0};  // pooled: its positions' entries in slot_offsets (positions + 1)
    // A one-root split route's entry whose root is not in its plan: its share has no route.
    [[nodiscard]] bool routeless() const noexcept { return tunnel_position < 0 && !pooled && source_position < 0; }
    [[nodiscard]] std::int32_t num_routers() const noexcept {
      return prepared_route == nullptr ? 0 : static_cast<std::int32_t>(prepared_route->topological_order.size());
    }
  };
  std::vector<const ForwardingView*> views;  // per demand
  std::vector<Item> entry_occurrences;
  // Per demand, per program entry: its item, or -1 (not forwarded), or -2 (a share with no route).
  std::vector<std::vector<std::int32_t>> entry_occurrence_positions;
  std::vector<EdgeId> edges;          // ascending
  std::vector<Amount> available_capacity;              // per link
  std::vector<std::int32_t> occurrence_edge_positions; // per occurrence
  std::vector<std::int32_t> occurrence_entry_positions;  // per occurrence: its item
  std::vector<std::int32_t> occurrence_route_positions;  // tunnel: hop; split route: member
  std::vector<std::int32_t> occurrence_heads;            // per occurrence: its view head
  std::vector<std::uint8_t> occurrence_dark;             // per occurrence: dropped where it arrives (empty: none can be)
  bool any_dark {false};                                 // some occurrence is dark
  std::int64_t max_link_occurrences {0};                 // the most sharing occurrences of one link
  std::vector<std::int32_t> edge_occurrence_offsets;  // per link + 1 into occurrences_by_edge
  std::vector<std::int32_t> occurrences_by_edge;     // the sharing occurrences, by (demand, head, entry)
  std::int32_t num_router_positions {0};
  std::vector<std::int32_t> router_entry_positions;    // per router position: its item
  // Pooled contexts: per member its first occurrence (members + 1 per item); per position its
  // first slot (positions + 1 per item; a destination has none); per slot its head; per
  // occurrence of a pooled item the slot of its member's head position for its head (-1 at a
  // destination; -1 for other items).
  std::vector<std::int32_t> member_occurrence_offsets;
  std::vector<std::int32_t> slot_offsets;
  std::vector<std::int32_t> slot_heads;
  std::vector<std::int32_t> occurrence_slots;  // empty without a pooled context
  [[nodiscard]] bool dark(std::size_t o) const noexcept { return any_dark && occurrence_dark[o] != 0; }
  std::vector<std::int32_t> member_tail_positions;  // pooled: per member (aligned with member_occurrence_offsets) its tail position
};

std::int32_t edge_position(const ForwardingNetwork& net, EdgeId e) {
  return static_cast<std::int32_t>(std::lower_bound(net.edges.begin(), net.edges.end(), e) - net.edges.begin());
}

// A pooled view's heads as roots of its plan: per view head, its root's plan position (-1: no
// member, its share has no route).
std::vector<std::int32_t> pooled_root_positions(const ForwardingView& view) {
  const PreparedRoute& plan = *view.prepared_route(0);
  std::vector<std::int32_t> out(view.num_sources(), -1);
  for (std::size_t h = 0; h < out.size(); ++h) {
    const auto e = static_cast<std::size_t>(view.entries_for_source(h).front());
    const auto r = view.context_source_index(e);
    if (r >= 0) out[h] = plan.source_positions[static_cast<std::size_t>(r)];
  }
  return out;
}

// Every entry of a non-pooled forwarding with a positive weight, by (head, entry index), to
// `visit` as its item before its positions in the network are set (first_occurrence,
// first_router_position): a tunnel with its hops, or a split route with its root's position in
// the plan and its members (none when routeless).
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

// A pooled plan's per-member tails (each member's tail position).
std::vector<std::int32_t> member_tails(const PreparedRoute& plan) {
  std::vector<std::int32_t> tail(plan.forwarding_edges.size(), -1);
  for (std::size_t pos = 0; pos < plan.topological_order.size(); ++pos) {
    for (auto m = plan.outgoing_begin[pos]; m < plan.outgoing_end[pos]; ++m) tail[static_cast<std::size_t>(m)] = static_cast<std::int32_t>(pos);
  }
  return tail;
}

// A pooled context's counts before anything of the network is allocated: its forwarded entries
// (the heads with a route) and its occurrences (per such head, the members its traffic can reach:
// the members whose tail it reaches in the plan). Exact, in 64 bits.
std::pair<std::int64_t, std::int64_t> pooled_counts(const ForwardingView& view) {
  const PreparedRoute& plan = *view.prepared_route(0);
  const auto roots = pooled_root_positions(view);
  const std::size_t P = plan.topological_order.size();
  std::vector<std::int32_t> seen(P, -1);
  std::vector<std::int32_t> stack;
  std::int64_t entries = 0;
  std::int64_t occurrences = 0;
  for (std::size_t h = 0; h < roots.size(); ++h) {
    if (roots[h] < 0) continue;
    ++entries;
    stack.assign(1, roots[h]);
    seen[static_cast<std::size_t>(roots[h])] = static_cast<std::int32_t>(h);
    while (!stack.empty()) {
      const auto pos = static_cast<std::size_t>(stack.back());
      stack.pop_back();
      occurrences += plan.outgoing_end[pos] - plan.outgoing_begin[pos];
      for (auto m = plan.outgoing_begin[pos]; m < plan.outgoing_end[pos]; ++m) {
        const auto q = static_cast<std::size_t>(plan.edge_destination_positions[static_cast<std::size_t>(m)]);
        if (seen[q] == static_cast<std::int32_t>(h)) continue;
        seen[q] = static_cast<std::int32_t>(h);
        stack.push_back(static_cast<std::int32_t>(q));
      }
    }
  }
  return {entries, occurrences};
}

// The network of `views` (one per demand, in the loss rule's demand order) on `free`. `lossy`: black
// holes make dark occurrences; otherwise a black-holed link serves 0 (the lossless bound).
ForwardingNetwork build_forwarding_network(std::span<const ForwardingView* const> views, const AvailableCapacityFn& free,
                                           bool lossy) {
  // The network's counts, exactly, after validation and before any of its arrays is sized: the
  // forwarded entries (positive weight, with a route; a pooled context's heads with a route), their
  // occurrences and the split routes' routers (a pooled context's once). A refusal changes nothing.
  std::int64_t num_entries = 0;
  std::int64_t num_occurrences = 0;
  std::int64_t num_routers = 0;
  for (const ForwardingView* view : views) {
    if (view->pooled) {
      const auto [e, o] = pooled_counts(*view);
      num_entries += e;
      num_occurrences += o;
      num_routers += static_cast<std::int64_t>(view->prepared_route(0)->topological_order.size());
      continue;
    }
    for_each_forwarded_entry(*view, [&](const ForwardingNetwork::Item& it) {
      if (it.routeless()) return;
      ++num_entries;
      num_occurrences += it.num_occurrences;
      num_routers += it.num_routers();
    });
  }
  check_forwarding_count(ForwardingCounter::Entry, num_entries);
  check_forwarding_count(ForwardingCounter::Occurrence, num_occurrences);
  check_forwarding_count(ForwardingCounter::Router, num_routers);
  ForwardingNetwork net;
  net.views.assign(views.begin(), views.end());
  net.entry_occurrence_positions.resize(views.size());
  for (std::size_t d = 0; d < views.size(); ++d) {
    const ForwardingView& view = *views[d];
    auto& positions = net.entry_occurrence_positions[d];
    positions.assign(view.data->entries.size(), -1);
    if (view.pooled) {
      ForwardingNetwork::Item it;
      it.demand = static_cast<std::int32_t>(d);
      it.pooled = true;
      it.prepared_route = view.prepared_route(0).get();
      it.first_router_position = net.num_router_positions;
      net.num_router_positions += it.num_routers();
      const auto roots = pooled_root_positions(view);
      bool any = false;
      for (std::size_t h = 0; h < roots.size(); ++h) {
        const auto e = static_cast<std::size_t>(view.entries_for_source(h).front());
        if (roots[h] < 0) {
          positions[e] = -2;  // its share has no route: lost at the source
          continue;
        }
        if (!any) it.source_index = static_cast<std::int32_t>(h);
        any = true;
        positions[e] = static_cast<std::int32_t>(net.entry_occurrences.size());
      }
      net.entry_occurrences.push_back(it);
      continue;
    }
    for_each_forwarded_entry(view, [&](ForwardingNetwork::Item it) {
      it.demand = static_cast<std::int32_t>(d);
      const auto e = static_cast<std::size_t>(it.entry_position);
      if (it.routeless()) {
        positions[e] = -2;  // its share has no route: lost at the source
        return;
      }
      if (it.prepared_route != nullptr) {
        it.first_router_position = net.num_router_positions;
        net.num_router_positions += it.num_routers();
      }
      positions[e] = static_cast<std::int32_t>(net.entry_occurrences.size());
      net.entry_occurrences.push_back(it);
    });
  }
  // Pooled contexts: the heads reaching each position (in topological order: its roots and the
  // heads of its in-members' tails), the members' occurrences and the positions' slots.
  std::vector<HeadLabelOffset> heads_offsets;  // per position of the current item + 1
  std::vector<std::int32_t> heads;          // the heads per position, ascending
  std::vector<std::int32_t> stamp;
  std::vector<std::int32_t> in_offsets;
  std::vector<std::int32_t> in_members;
  std::int32_t nocc = 0;
  for (auto& it : net.entry_occurrences) {
    it.first_occurrence = nocc;
    if (!it.pooled) {
      nocc += it.num_occurrences;
      continue;
    }
    const ForwardingView& view = *net.views[static_cast<std::size_t>(it.demand)];
    const PreparedRoute& plan = *it.prepared_route;
    const std::size_t P = plan.topological_order.size();
    const std::size_t M = plan.forwarding_edges.size();
    const auto roots = pooled_root_positions(view);
    const auto tails = member_tails(plan);
    in_offsets.assign(P + 1, 0);
    for (std::size_t m = 0; m < M; ++m) ++in_offsets[static_cast<std::size_t>(plan.edge_destination_positions[m]) + 1];
    for (std::size_t q = 0; q < P; ++q) in_offsets[q + 1] += in_offsets[q];
    in_members.resize(M);
    {
      std::vector<std::int32_t> fill(in_offsets.begin(), in_offsets.end() - 1);
      for (std::size_t m = 0; m < M; ++m) in_members[static_cast<std::size_t>(fill[static_cast<std::size_t>(plan.edge_destination_positions[m])]++)] = static_cast<std::int32_t>(m);
    }
    std::vector<std::vector<std::int32_t>> roots_at(P);
    for (std::size_t h = 0; h < roots.size(); ++h) {
      if (roots[h] >= 0) roots_at[static_cast<std::size_t>(roots[h])].push_back(static_cast<std::int32_t>(h));
    }
    stamp.assign(roots.size(), -1);
    heads_offsets.assign(P + 1, 0);
    heads.clear();
    for (std::size_t q = 0; q < P; ++q) {
      const auto first = heads.size();
      auto add = [&](std::int32_t h) {
        if (stamp[static_cast<std::size_t>(h)] == static_cast<std::int32_t>(q)) return;
        stamp[static_cast<std::size_t>(h)] = static_cast<std::int32_t>(q);
        heads.push_back(h);
      };
      for (auto h : roots_at[q]) add(h);
      for (auto i = in_offsets[q]; i < in_offsets[q + 1]; ++i) {
        const auto t = static_cast<std::size_t>(tails[static_cast<std::size_t>(in_members[static_cast<std::size_t>(i)])]);
        for (auto k = heads_offsets[t]; k < heads_offsets[t + 1]; ++k) add(heads[static_cast<std::size_t>(k)]);
      }
      std::sort(heads.begin() + static_cast<std::ptrdiff_t>(first), heads.end());
      heads_offsets[q + 1] = heads.size();
    }
    it.first_member = net.member_occurrence_offsets.size();
    std::int32_t count = 0;
    for (std::size_t m = 0; m < M; ++m) {
      net.member_occurrence_offsets.push_back(nocc + count);
      net.member_tail_positions.push_back(tails[m]);
      const auto t = static_cast<std::size_t>(tails[m]);
      // This range contributes to the preflight-checked occurrence total, unlike the
      // cumulative label offset (which also includes labels at destinations).
      count += static_cast<std::int32_t>(heads_offsets[t + 1] - heads_offsets[t]);
    }
    net.member_occurrence_offsets.push_back(nocc + count);
    net.member_tail_positions.push_back(-1);
    it.num_occurrences = count;
    it.first_position = net.slot_offsets.size();
    for (std::size_t q = 0; q < P; ++q) {
      net.slot_offsets.push_back(static_cast<std::int32_t>(net.slot_heads.size()));
      if (plan.is_destination[q] != 0) continue;
      for (auto k = heads_offsets[q]; k < heads_offsets[q + 1]; ++k) net.slot_heads.push_back(heads[static_cast<std::size_t>(k)]);
    }
    net.slot_offsets.push_back(static_cast<std::int32_t>(net.slot_heads.size()));
    nocc += count;
  }
  if (nocc != num_occurrences) throw std::logic_error("forwarding: the network's occurrences differ from their count");
  for (const auto& it : net.entry_occurrences) {
    if (it.tunnel_position >= 0) {
      const auto edges = path_of(net.views[static_cast<std::size_t>(it.demand)]->data->tunnels[static_cast<std::size_t>(it.tunnel_position)]);
      net.edges.insert(net.edges.end(), edges.begin(), edges.end());
    } else if (it.num_occurrences != 0) {
      net.edges.insert(net.edges.end(), it.prepared_route->forwarding_edges.begin(), it.prepared_route->forwarding_edges.end());
    }
  }
  std::sort(net.edges.begin(), net.edges.end());
  net.edges.erase(std::unique(net.edges.begin(), net.edges.end()), net.edges.end());
  net.available_capacity.resize(net.edges.size());
  free(net.edges, net.available_capacity);
  if (!lossy) {
    // The lossless bound reads a black-holed link as serving 0 (one view: lossless has no list).
    for (const ForwardingView* view : views) {
      for (std::size_t l = 0; view->black_holes != nullptr && l < net.edges.size(); ++l) {
        if (!view->black_holes->allows_edge(net.edges[l])) net.available_capacity[l] = 0;
      }
    }
  }
  const auto n = static_cast<std::size_t>(nocc);
  net.occurrence_edge_positions.resize(n);
  net.occurrence_entry_positions.resize(n);
  net.occurrence_route_positions.resize(n);
  net.occurrence_heads.resize(n);
  // The dark flags only where a black hole can apply (LOSSY), the slots only for pooled contexts.
  bool can_be_dark = false;
  bool any_pooled = false;
  for (const ForwardingView* v : net.views) {
    can_be_dark = can_be_dark || (lossy && v->black_holes != nullptr);
    any_pooled = any_pooled || v->pooled;
  }
  if (can_be_dark) net.occurrence_dark.assign(n, 0);
  if (any_pooled) net.occurrence_slots.assign(n, -1);
  for (std::size_t i = 0; i < net.entry_occurrences.size(); ++i) {
    const auto& it = net.entry_occurrences[i];
    const ForwardingView& view = *net.views[static_cast<std::size_t>(it.demand)];
    const Eligibility* dark = lossy ? view.black_holes : nullptr;
    auto set = [&](std::size_t o, EdgeId e, std::int32_t j, std::int32_t head) {
      net.occurrence_edge_positions[o] = edge_position(net, e);
      net.occurrence_entry_positions[o] = static_cast<std::int32_t>(i);
      net.occurrence_route_positions[o] = j;
      net.occurrence_heads[o] = head;
      if (dark != nullptr && !dark->allows_edge(e)) {
        net.occurrence_dark[o] = 1;
        net.any_dark = true;
      }
    };
    if (!it.pooled) {
      for (std::int32_t j = 0; j < it.num_occurrences; ++j) {
        const auto o = static_cast<std::size_t>(it.first_occurrence + j);
        const EdgeId e = it.tunnel_position >= 0 ? path_of(view.data->tunnels[static_cast<std::size_t>(it.tunnel_position)])[static_cast<std::size_t>(j)]
                                                 : it.prepared_route->forwarding_edges[static_cast<std::size_t>(j)];
        set(o, e, j, it.source_index);
      }
      continue;
    }
    const PreparedRoute& plan = *it.prepared_route;
    const std::size_t M = plan.forwarding_edges.size();
    for (std::size_t m = 0; m < M; ++m) {
      const auto t = static_cast<std::size_t>(net.member_tail_positions[it.first_member + m]);
      const auto q = static_cast<std::size_t>(plan.edge_destination_positions[m]);
      const auto first = static_cast<std::size_t>(net.member_occurrence_offsets[it.first_member + m]);
      const auto tail_slots = static_cast<std::size_t>(net.slot_offsets[it.first_position + t]);
      const auto num = static_cast<std::size_t>(net.member_occurrence_offsets[it.first_member + m + 1]) - first;
      const auto head_slots = std::span<const std::int32_t>(net.slot_heads).subspan(
          static_cast<std::size_t>(net.slot_offsets[it.first_position + q]),
          static_cast<std::size_t>(net.slot_offsets[it.first_position + q + 1] - net.slot_offsets[it.first_position + q]));
      for (std::size_t k = 0; k < num; ++k) {
        const auto o = first + k;
        const auto head = net.slot_heads[tail_slots + k];  // the heads of the tail, ascending
        set(o, plan.forwarding_edges[m], static_cast<std::int32_t>(m), head);
        if (plan.is_destination[q] == 0) {
          const auto at = std::lower_bound(head_slots.begin(), head_slots.end(), head);
          net.occurrence_slots[o] = net.slot_offsets[it.first_position + q] + static_cast<std::int32_t>(at - head_slots.begin());
        }
      }
    }
  }
  net.router_entry_positions.resize(static_cast<std::size_t>(net.num_router_positions));
  for (std::size_t i = 0; i < net.entry_occurrences.size(); ++i) {
    const auto& it = net.entry_occurrences[i];
    if (it.tunnel_position >= 0) continue;
    const auto& plan = *it.prepared_route;
    for (std::size_t pos = 0; pos < plan.topological_order.size(); ++pos) net.router_entry_positions[static_cast<std::size_t>(it.first_router_position) + pos] = static_cast<std::int32_t>(i);
  }
  net.edge_occurrence_offsets.assign(net.edges.size() + 1, 0);
  for (std::size_t o = 0; o < n; ++o) {
    if (!net.dark(o)) ++net.edge_occurrence_offsets[static_cast<std::size_t>(net.occurrence_edge_positions[o]) + 1];
  }
  for (std::size_t l = 0; l < net.edges.size(); ++l) {
    net.max_link_occurrences = std::max<std::int64_t>(net.max_link_occurrences, net.edge_occurrence_offsets[l + 1]);
    net.edge_occurrence_offsets[l + 1] += net.edge_occurrence_offsets[l];
  }
  net.occurrences_by_edge.resize(static_cast<std::size_t>(net.edge_occurrence_offsets.back()));
  std::vector<std::int32_t> fill(net.edge_occurrence_offsets.begin(), net.edge_occurrence_offsets.end() - 1);
  // Occurrences are numbered by item, items by (demand, head, entry index), a pooled item's by
  // member then head: each link's list comes out in the loss rule's order (demands, then heads
  // ascending, then entries ascending).
  for (std::size_t o = 0; o < n; ++o) {
    if (net.dark(o)) continue;
    net.occurrences_by_edge[static_cast<std::size_t>(fill[static_cast<std::size_t>(net.occurrence_edge_positions[o])]++)] = static_cast<std::int32_t>(o);
  }
  return net;
}

// One no-drop pass; Checked: dark occurrences stop their traffic and link sums saturate, else plain
// sums (no_drop takes the plain pass only where nothing is dark and no sum can reach kSumCap, where
// both give the same loads).
template <bool Checked>
void no_drop_pass(const ForwardingNetwork& net, std::span<const std::vector<Amount>> shares, std::vector<Amount>& load,
                  std::vector<Amount>* occ, std::vector<Amount>& scratch) {
  auto add = [&](std::size_t o, Amount x) {
    if constexpr (Checked) {
      if (net.dark(o)) return;
      auto& y = load[static_cast<std::size_t>(net.occurrence_edge_positions[o])];
      y = saturating_add(y, x);
    } else {
      load[static_cast<std::size_t>(net.occurrence_edge_positions[o])] += x;
    }
  };
  // A split route's plan with a dark member: forwarded with those members' view at 0, so what
  // they receive goes no further and the routers split only what arrives. Per member its arrival
  // (in scratch); false without a dark member (the forward is not needed).
  std::vector<Amount> view_m;
  auto dark_forward = [&](const ForwardingNetwork::Item& it, std::span<const Amount> inject) {
    if constexpr (Checked) {
      if (!net.any_dark) return false;
      const PreparedRoute& plan = *it.prepared_route;
      const std::size_t M = plan.forwarding_edges.size();
      view_m.assign(M, std::numeric_limits<Amount>::max());
      bool any = false;
      for (std::size_t m = 0; m < M; ++m) {
        const auto o = it.pooled ? static_cast<std::size_t>(net.member_occurrence_offsets[it.first_member + m])
                                 : static_cast<std::size_t>(it.first_occurrence) + m;
        if (it.pooled && o == static_cast<std::size_t>(net.member_occurrence_offsets[it.first_member + m + 1])) continue;  // no head reaches it
        if (net.dark(o)) {
          view_m[m] = 0;
          any = true;
        }
      }
      if (!any) return false;
      const auto fw = forward(plan, Injection {inject, false}, view_m);
      scratch.resize(M);
      for (std::size_t m = 0; m < M; ++m) scratch[m] = fw.edge_flows[m] + fw.edge_drops[m];
      return true;
    } else {
      (void)it;
      (void)inject;
      return false;
    }
  };
  std::vector<Amount> inject;
  std::vector<std::int32_t> head_root;   // a pooled context: per view head, its root's position (-1: no share)
  std::vector<std::uint8_t> slot_live;   // a pooled context: per slot, its head's traffic can reach it
  for (const auto& it : net.entry_occurrences) {
    const auto& share = shares[static_cast<std::size_t>(it.demand)];
    if (it.pooled) {
      const ForwardingView& view = *net.views[static_cast<std::size_t>(it.demand)];
      const PreparedRoute& plan = *it.prepared_route;
      inject.assign(plan.sources.size(), 0);
      head_root.assign(view.num_sources(), -1);
      bool any = false;
      for (std::size_t h = 0; h < view.num_sources(); ++h) {
        const auto e = static_cast<std::size_t>(view.entries_for_source(h).front());
        const auto r = view.context_source_index(e);
        if (r < 0 || plan.source_positions[static_cast<std::size_t>(r)] < 0) continue;
        inject[static_cast<std::size_t>(r)] = share[e];
        if (share[e] != 0) head_root[h] = plan.source_positions[static_cast<std::size_t>(r)];
        any = any || share[e] != 0;
      }
      if (!any) continue;
      if (!dark_forward(it, inject)) scratch = forward(plan, Injection {inject, false}, {}).edge_flows;
      for (std::size_t m = 0; m < plan.forwarding_edges.size(); ++m) {
        const auto o = static_cast<std::size_t>(net.member_occurrence_offsets[it.first_member + m]);
        const auto end = static_cast<std::size_t>(net.member_occurrence_offsets[it.first_member + m + 1]);
        if (scratch[m] == 0 || o == end) continue;
        add(o, scratch[m]);
      }
      if (occ == nullptr) continue;
      // Per head: its occurrence on a member can carry something when the member's amount is
      // positive and the head's traffic reaches the member's tail, from its own root (a positive
      // share) or over a member it can carry something on that is not dark (what reaches a dark
      // member goes no further). Positions run in topological order; a member's occurrences and its
      // tail's slots are the heads reaching that tail, both ascending.
      const auto first_slot = net.slot_offsets[it.first_position];
      slot_live.assign(static_cast<std::size_t>(net.slot_offsets[it.first_position + plan.topological_order.size()] - first_slot), 0);
      for (std::size_t h = 0; h < head_root.size(); ++h) {
        if (head_root[h] < 0 || plan.is_destination[static_cast<std::size_t>(head_root[h])] != 0) continue;
        const auto q = it.first_position + static_cast<std::size_t>(head_root[h]);
        const auto b = net.slot_heads.begin();
        const auto at = std::lower_bound(b + net.slot_offsets[q], b + net.slot_offsets[q + 1], static_cast<std::int32_t>(h));
        slot_live[static_cast<std::size_t>(at - b - first_slot)] = 1;
      }
      for (std::size_t pos = 0; pos < plan.topological_order.size(); ++pos) {
        if (plan.is_destination[pos] != 0) continue;
        const auto tail_slot = static_cast<std::size_t>(net.slot_offsets[it.first_position + pos] - first_slot);
        for (auto m = static_cast<std::size_t>(plan.outgoing_begin[pos]); m < static_cast<std::size_t>(plan.outgoing_end[pos]); ++m) {
          if (scratch[m] == 0) continue;
          const auto o = static_cast<std::size_t>(net.member_occurrence_offsets[it.first_member + m]);
          const auto end = static_cast<std::size_t>(net.member_occurrence_offsets[it.first_member + m + 1]);
          for (auto x = o; x < end; ++x) {
            if (slot_live[tail_slot + (x - o)] == 0) continue;
            (*occ)[x] = scratch[m];
            const auto head_slot = net.occurrence_slots[x];
            if (head_slot >= 0 && !net.dark(x)) slot_live[static_cast<std::size_t>(head_slot - first_slot)] = 1;
          }
        }
      }
      continue;
    }
    const Amount s = share[static_cast<std::size_t>(it.entry_position)];
    if (s == 0) continue;
    if (it.tunnel_position >= 0) {
      for (std::int32_t j = 0; j < it.num_occurrences; ++j) {
        const auto o = static_cast<std::size_t>(it.first_occurrence + j);
        if (occ != nullptr) (*occ)[o] = s;
        add(o, s);
        // What reaches its own black hole goes no further: the later hops carry nothing.
        if (Checked && net.dark(o)) break;
      }
      continue;
    }
    const PreparedRoute& plan = *it.prepared_route;
    bool forwarded = false;
    if (Checked && net.any_dark) {
      inject.assign(plan.sources.size(), 0);
      for (std::size_t r = 0; r < plan.sources.size(); ++r) {
        if (plan.source_positions[r] == it.source_position) inject[r] = s;
      }
      forwarded = dark_forward(it, inject);
    }
    if (!forwarded) {
      scratch.resize(plan.forwarding_edges.size());
      (void)allocate(plan, s, scratch);
    }
    for (std::int32_t j = 0; j < it.num_occurrences; ++j) {
      const auto o = static_cast<std::size_t>(it.first_occurrence + j);
      if (occ != nullptr) (*occ)[o] = scratch[static_cast<std::size_t>(j)];
      add(o, scratch[static_cast<std::size_t>(j)]);
    }
  }
}

// The no-drop loads per link of every demand's shares (`shares[d]`: per program entry of demand
// d), each demand's traffic dropped only at its own black holes, and per occurrence when `occ` is
// given its possible amount (a pooled context's occurrence: its member's flow where its head's
// traffic can reach the member, else 0). A split route allocates its share over its plan from its root; a pooled
// context forwards its heads' shares over its plan without a limit; a tunnel's hops after a dark
// one, and a dark member's onward traffic, carry nothing. A dark occurrence (its arrival as its
// amount) loads no link; link sums saturate (kSumCap). An occurrence whose amount is 0 can carry
// nothing in any evaluation (splits never give less for more): it is dead, and takes no part in
// the waits or in serial normalization.
void no_drop(const ForwardingNetwork& net, std::span<const std::vector<Amount>> shares, std::vector<Amount>& load,
             std::vector<Amount>* occ, std::vector<Amount>& scratch) {
  load.assign(net.edges.size(), 0);
  if (occ != nullptr) occ->assign(net.occurrence_edge_positions.size(), 0);
  // Every occurrence carries at most the sum of all shares (a pooled member at most its demand's),
  // so a link's sum stays below total * max_link_occurrences: below kSumCap, plain sums are exact.
  Amount total = 0;
  for (const auto& share : shares) {
    for (Amount x : share) total = saturating_add(total, x);
  }
  const bool plain = !net.any_dark && (net.max_link_occurrences == 0 || total < kSumCap / net.max_link_occurrences);
  if (plain) {
    no_drop_pass<false>(net, shares, load, occ, scratch);
  } else {
    no_drop_pass<true>(net, shares, load, occ, scratch);
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

// The loss rule's scratch (per call).
struct LossScratch {
  std::vector<Amount> wd;  // per demand run: its arrivals
  std::vector<Amount> td;  // per demand run: its share of T
  std::vector<std::size_t> demand_run;
  std::vector<Amount> wa;
  std::vector<Amount> wb;
  std::vector<std::size_t> run;
  std::vector<Amount> w;
  std::vector<Amount> c;
};

// Within one demand's run occs[b, e): highest averages of T over its sources' arrivals, then
// within a source over its entries'.
void source_levels(const ForwardingNetwork& net, std::span<const std::int32_t> occs, std::size_t b, std::size_t e, Amount T,
                   std::span<const Amount> arrival, std::span<Amount> carried, LossScratch& s) {
  // One positive arrival takes T whole (highest averages over one positive weight).
  if (e - b == 1 && arrival[static_cast<std::size_t>(occs[b])] > 0) {
    carried[static_cast<std::size_t>(occs[b])] = T;
    return;
  }
  s.wa.clear();
  s.run.clear();  // the first index of each source's run
  for (std::size_t i = b; i < e; ++i) {
    const auto head = net.occurrence_heads[static_cast<std::size_t>(occs[i])];
    if (i == b || head != net.occurrence_heads[static_cast<std::size_t>(occs[i - 1])]) {
      s.run.push_back(i);
      s.wa.push_back(0);
    }
    s.wa.back() += arrival[static_cast<std::size_t>(occs[i])];
  }
  s.run.push_back(e);
  s.wb.resize(s.wa.size());
  highest_averages(T, s.wa, s.wb);
  for (std::size_t k = 0; k + 1 < s.run.size(); ++k) {
    const std::size_t rb = s.run[k];
    const std::size_t re = s.run[k + 1];
    if (re - rb == 1 && arrival[static_cast<std::size_t>(occs[rb])] > 0) {
      carried[static_cast<std::size_t>(occs[rb])] = s.wb[k];
      continue;
    }
    s.w.resize(re - rb);
    s.c.resize(re - rb);
    for (std::size_t i = rb; i < re; ++i) s.w[i - rb] = arrival[static_cast<std::size_t>(occs[i])];
    highest_averages(s.wb[k], s.w, s.c);
    for (std::size_t i = rb; i < re; ++i) carried[static_cast<std::size_t>(occs[i])] = s.c[i - rb];
  }
}

// The loss rule at one link: T = min(V, sum of the arrivals); highest averages over the
// demands' arrivals (one demand takes T whole), then within a demand over its sources', then
// within a source over its entries'. `occs` lists the link's sharing occurrences by (demand,
// source, entry); writes their carried amounts.
void loss_rule(const ForwardingNetwork& net, std::span<const std::int32_t> occs, Amount V, std::span<const Amount> arrival,
               std::span<Amount> carried, LossScratch& s) {
  Amount sum = 0;
  for (auto o : occs) sum = saturating_add(sum, arrival[static_cast<std::size_t>(o)]);
  if (sum <= V) {
    for (auto o : occs) carried[static_cast<std::size_t>(o)] = arrival[static_cast<std::size_t>(o)];
    return;
  }
  // The demand runs: consecutive runs of one demand (each demand's arrivals sum below 2^62).
  s.wd.clear();
  s.demand_run.clear();
  for (std::size_t i = 0; i < occs.size(); ++i) {
    const auto d = net.entry_occurrences[static_cast<std::size_t>(net.occurrence_entry_positions[static_cast<std::size_t>(occs[i])])].demand;
    if (i == 0 || d != net.entry_occurrences[static_cast<std::size_t>(net.occurrence_entry_positions[static_cast<std::size_t>(occs[i - 1])])].demand) {
      s.demand_run.push_back(i);
      s.wd.push_back(0);
    }
    s.wd.back() += arrival[static_cast<std::size_t>(occs[i])];
  }
  s.demand_run.push_back(occs.size());
  if (s.wd.size() == 1) {
    source_levels(net, occs, 0, occs.size(), V, arrival, carried, s);
    return;
  }
  s.td.resize(s.wd.size());
  highest_averages(V, s.wd, s.td);
  // source_levels uses none of the demand level's scratch.
  for (std::size_t k = 0; k + 1 < s.demand_run.size(); ++k) {
    source_levels(net, occs, s.demand_run[k], s.demand_run[k + 1], s.td[k], arrival, carried, s);
  }
}

// ---------------------------------------------------------------- the general lossless bound

ProgramFlow general_lossless(const ForwardingView& view, const AvailableCapacityFn& free, Amount lo, std::optional<Amount> hi) {
  const ForwardingView* const views[] = {&view};
  const ForwardingNetwork net = build_forwarding_network(views, free, false);
  const auto& positions = net.entry_occurrence_positions.front();
  const std::size_t k = view.num_sources();
  std::vector<Amount> heads(k);
  std::vector<std::vector<Amount>> shares(1);
  auto& share = shares.front();
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
      if (share[e] != 0 && positions[e] == -2) return false;  // a share with no route
    }
    no_drop(net, shares, load, nullptr, scratch);
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
  no_drop(net, shares, load, &occ, scratch);
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

// The loss rule over every demand of one network at once: A+ where arrivals are
// determined, the cyclic groups iterated, an independent check of an iterated result. One demand
// of a non-pooled program is the scalar call's evaluation.
class LossyForwardingEvaluator {
public:
  explicit LossyForwardingEvaluator(const ForwardingNetwork& net) : net_(net) {}

  // Per demand its per-head injections; one ProgramFlow per demand.
  std::vector<ProgramFlow> run(std::span<const std::span<const Amount>> injections);

private:
  // Event ids: barriers [0, nb), routers [nb, nb + num_pos).
  [[nodiscard]] std::int32_t router(const ForwardingNetwork::Item& it, std::int32_t pos) const noexcept {
    return num_barriers_ + it.first_router_position + pos;
  }
  [[nodiscard]] bool dark(std::int32_t o) const noexcept { return net_.dark(static_cast<std::size_t>(o)); }
  // A dark occurrence is no event (what reaches it is dropped at once).
  [[nodiscard]] std::int32_t event_of_occ(std::int32_t o) const noexcept {
    if (dark(o)) return -1;
    return barrier_[static_cast<std::size_t>(net_.occurrence_edge_positions[static_cast<std::size_t>(o)])];
  }
  // A push into event `ev` is suppressed while a group iterates and ev lies outside it.
  [[nodiscard]] bool blocked(std::int32_t ev) const noexcept { return iterating_ && scc_[static_cast<std::size_t>(ev)] != current_component_; }
  [[nodiscard]] const ForwardingNetwork::Item& item_of(std::int32_t o) const noexcept {
    return net_.entry_occurrences[static_cast<std::size_t>(net_.occurrence_entry_positions[static_cast<std::size_t>(o)])];
  }
  // A pooled item's slots at position pos (a destination: none).
  [[nodiscard]] std::pair<std::int32_t, std::int32_t> slots(const ForwardingNetwork::Item& it, std::size_t pos) const noexcept {
    return {net_.slot_offsets[it.first_position + pos], net_.slot_offsets[it.first_position + pos + 1]};
  }
  [[nodiscard]] std::int32_t member_occurrence(const ForwardingNetwork::Item& it, std::size_t m) const noexcept {
    return net_.member_occurrence_offsets[it.first_member + m];
  }

  void normalize_serial_links(std::span<const Amount> U0);
  void build_events();
  void find_strong_components();
  void reach(std::int32_t o, Amount x);
  void carry(std::int32_t o, Amount c);
  void split(std::int32_t ev);
  void evaluate_barrier(std::int32_t ev);
  void iterate(std::int32_t scc);
  void verify() const;

  const ForwardingNetwork& net_;
  std::vector<std::vector<Amount>> shares_;  // per demand, per program entry
  std::vector<std::int32_t> occurrence_source_positions_;  // per occurrence of a split route: its member's tail position
  std::vector<std::uint8_t> root_position_;  // per router position: a root of its item's heads
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
  std::vector<Amount> possible_;  // per occurrence: its possible amount (0: dead)
  std::vector<Amount> arrival_;
  std::vector<Amount> carried_;
  std::vector<Amount> inflow_;       // per router position (one-root split routes)
  std::vector<Amount> slot_inflow_;  // per slot (pooled contexts)
  std::vector<Amount> delivered_;    // per item (tunnels and one-root split routes)
  std::vector<std::vector<Amount>> head_delivered_;  // per demand: a pooled context's per head
  std::vector<Amount> outgoing_shares_;
  std::vector<Amount> remaining_;
  std::vector<Amount> row_;
  LossScratch loss_;
  bool iterated_ {false};
};

void LossyForwardingEvaluator::normalize_serial_links(std::span<const Amount> U0) {
  const std::size_t L = net_.edges.size();
  transparent_.assign(L, 0);
  // Per link: the one link every sharing occurrence continues on, and the one every sharing
  // occurrence comes from (-1: none yet; -2: they differ, or an occurrence is injected, absorbed,
  // merged or split there, or its neighbour on its route is dark). Only what can carry something
  // counts: a dead occurrence counts for neither link, and a router's merges and branches are
  // those of its live members.
  std::vector<std::int32_t> next(L, -1);
  std::vector<std::int32_t> prev(L, -1);
  auto merge = [](std::int32_t& slot, std::int32_t v) {
    if (slot == -1) {
      slot = v;
    } else if (slot != v) {
      slot = -2;
    }
  };
  // The link of occurrence o, or -2 when it is dark.
  auto link_of = [&](std::int32_t o) { return dark(o) ? -2 : net_.occurrence_edge_positions[static_cast<std::size_t>(o)]; };
  // Per router position, over the live members only (an occurrence of the member can carry
  // something): its one in-member (-1 none, -2 several) and its out-members (their count, and the
  // one when there is one). A dark member takes its share out of its tail but brings nothing into
  // its head (what reaches it is dropped there).
  const auto R = static_cast<std::size_t>(net_.num_router_positions);
  std::vector<std::int32_t> live_in(R, -1);
  std::vector<std::int32_t> live_out(R, -1);
  std::vector<std::int32_t> live_outs(R, 0);
  for (const auto& it : net_.entry_occurrences) {
    if (it.tunnel_position >= 0 || it.num_occurrences == 0) continue;
    const auto& plan = *it.prepared_route;
    const auto base = static_cast<std::size_t>(it.first_router_position);
    for (std::size_t pos = 0; pos < plan.topological_order.size(); ++pos) {
      for (auto m = plan.outgoing_begin[pos]; m < plan.outgoing_end[pos]; ++m) {
        const auto first = it.pooled ? member_occurrence(it, static_cast<std::size_t>(m)) : it.first_occurrence + m;
        const auto end = it.pooled ? member_occurrence(it, static_cast<std::size_t>(m) + 1) : first + 1;
        bool live = false;
        for (auto x = first; x < end && !live; ++x) live = possible_[static_cast<std::size_t>(x)] != 0;
        if (!live) continue;
        ++live_outs[base + pos];
        live_out[base + pos] = m;
        if (dark(first)) continue;
        auto& slot = live_in[base + static_cast<std::size_t>(plan.edge_destination_positions[static_cast<std::size_t>(m)])];
        slot = slot == -1 ? m : -2;
      }
    }
  }
  for (std::size_t o = 0; o < net_.occurrence_edge_positions.size(); ++o) {
    if (net_.dark(o) || possible_[o] == 0) continue;
    const auto l = static_cast<std::size_t>(net_.occurrence_edge_positions[o]);
    const auto& it = net_.entry_occurrences[static_cast<std::size_t>(net_.occurrence_entry_positions[o])];
    const std::int32_t j = net_.occurrence_route_positions[o];
    std::int32_t succ = -2;
    std::int32_t pred = -2;
    if (it.tunnel_position >= 0) {
      if (j + 1 < it.num_occurrences) succ = link_of(static_cast<std::int32_t>(o) + 1);
      if (j > 0) pred = link_of(static_cast<std::int32_t>(o) - 1);
    } else {
      const auto& plan = *it.prepared_route;
      const auto base = static_cast<std::size_t>(it.first_router_position);
      // The occurrence of member m on this occurrence's route (a pooled context: the same head's).
      auto on_member = [&](std::int32_t m) -> std::int32_t {
        if (!it.pooled) return it.first_occurrence + m;
        // A member's occurrences ascend by head.
        const auto b = net_.occurrence_heads.begin();
        const auto first = b + member_occurrence(it, static_cast<std::size_t>(m));
        const auto last = b + member_occurrence(it, static_cast<std::size_t>(m) + 1);
        const auto at = std::lower_bound(first, last, net_.occurrence_heads[o]);
        return at != last && *at == net_.occurrence_heads[o] ? static_cast<std::int32_t>(at - b) : -1;
      };
      // Continues: its head is no destination and no root, holds this member alone as its live
      // in-member and has one live out-member.
      const auto q = static_cast<std::size_t>(plan.edge_destination_positions[static_cast<std::size_t>(j)]);
      if (plan.is_destination[q] == 0 && root_position_[base + q] == 0 && live_outs[base + q] == 1 && live_in[base + q] == j) {
        const auto x = on_member(live_out[base + q]);
        if (x >= 0) succ = link_of(x);
      }
      // Comes from: its tail is no root, has one live in-member and this one live out-member.
      const auto t = static_cast<std::size_t>(occurrence_source_positions_[o]);
      const auto in = live_in[base + t];
      if (root_position_[base + t] == 0 && in >= 0 && live_outs[base + t] == 1) {
        const auto x = on_member(in);
        if (x >= 0) pred = link_of(x);
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
#if NETGRAPH_CORE_PRIVATE_HOOKS
  thread_lossy_counters().barriers.clear();
  for (auto l : barrier_link_) thread_lossy_counters().barriers.push_back(net_.edges[static_cast<std::size_t>(l)]);
#endif
  // The events are known only after the no-drop loads: checked before the event structures are
  // sized (the network is built; a refusal changes nothing).
  check_forwarding_count(ForwardingCounter::Event, std::int64_t {num_barriers_} + net_.num_router_positions);
  num_events_ = num_barriers_ + net_.num_router_positions;
}

void LossyForwardingEvaluator::build_events() {
  // Every dependency arc, to `emit`: an event before the next one on its tunnel; on a split route
  // (a pooled context: per member, whatever its heads) a member's tail router before its event (a
  // member that can drop) and that event before the head router, or the tail before the head.
  // Only through an occurrence that can carry something: a dead occurrence (a tunnel's hop whose
  // share is 0 or that lies behind its own black hole, a member no head can carry anything on)
  // adds no arc, so it makes no link wait and joins no group.
  auto arcs = [&](auto&& emit) {
    for (const auto& it : net_.entry_occurrences) {
      if (it.tunnel_position >= 0) {
        std::int32_t last = -1;
        for (std::int32_t j = 0; j < it.num_occurrences; ++j) {
          if (possible_[static_cast<std::size_t>(it.first_occurrence + j)] == 0) break;  // the rest are dead too
          const auto ev = event_of_occ(it.first_occurrence + j);
          if (ev < 0) continue;
          if (last >= 0) emit(last, ev);
          last = ev;
        }
        continue;
      }
      if (it.num_occurrences == 0) continue;
      const auto& plan = *it.prepared_route;
      for (std::size_t pos = 0; pos < plan.topological_order.size(); ++pos) {
        for (auto m = plan.outgoing_begin[pos]; m < plan.outgoing_end[pos]; ++m) {
          std::int32_t o = it.first_occurrence + m;
          if (it.pooled) {
            // Its first occurrence that can carry something (none: no head reaches it, or none
            // can carry anything on it).
            const auto end = member_occurrence(it, static_cast<std::size_t>(m) + 1);
            o = member_occurrence(it, static_cast<std::size_t>(m));
            while (o < end && possible_[static_cast<std::size_t>(o)] == 0) ++o;
            if (o == end) continue;
          }
          if (possible_[static_cast<std::size_t>(o)] == 0) continue;
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
#if NETGRAPH_CORE_PRIVATE_HOOKS
  ++thread_lossy_counters().event_structures;
  thread_lossy_counters().dependencies = count;
#endif
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

void LossyForwardingEvaluator::find_strong_components() {
  num_scc_ = strong_components(dependency_offsets_, dependencies_,
                               [](std::size_t) { return true; }, scc_);
}

// Amount x reaches occurrence o (a router's split, or a tunnel's injection): a dark occurrence
// drops it there and passes nothing on; a barrier's waits for its link's rule; a pass-through
// carries it on at once.
void LossyForwardingEvaluator::reach(std::int32_t o, Amount x) {
  const auto oi = static_cast<std::size_t>(o);
  if (dark(o)) {
    arrival_[oi] = x;
    carry(o, 0);
    return;
  }
  const auto ev = event_of_occ(o);
  if (ev >= 0) {
    if (!blocked(ev)) arrival_[oi] = x;
    return;
  }
  arrival_[oi] = x;
  carry(o, x);
}

// Occurrence o now carries c: forward it to what it feeds (a split route's head router, a pooled
// context's slot of its head there, a tunnel's next hops up to its next barrier or its end; a dark
// hop drops what reaches it).
void LossyForwardingEvaluator::carry(std::int32_t o, Amount c) {
  carried_[static_cast<std::size_t>(o)] = c;
  const auto& it = item_of(o);
  const std::int32_t j = net_.occurrence_route_positions[static_cast<std::size_t>(o)];
  if (it.tunnel_position < 0) {
    const auto q = it.prepared_route->edge_destination_positions[static_cast<std::size_t>(j)];
    const auto ev = router(it, q);
    if (it.pooled) {
      // A pooled destination's deliveries are read from the carried amounts at the end.
      if (it.prepared_route->is_destination[static_cast<std::size_t>(q)] != 0) return;
      if (!blocked(ev)) slot_inflow_[static_cast<std::size_t>(net_.occurrence_slots[static_cast<std::size_t>(o)])] += c;
      return;
    }
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
    const auto o2i = static_cast<std::size_t>(o2);
    if (dark(o2)) {
      arrival_[o2i] = c;
      carried_[o2i] = 0;
      c = 0;
      continue;
    }
    const auto ev = event_of_occ(o2);
    if (ev >= 0) {
      if (!blocked(ev)) arrival_[o2i] = c;
      return;
    }
    arrival_[o2i] = c;
    carried_[o2i] = c;
  }
}

// A router splits what it holds once over its out-members; a pass-through member carries its
// whole arrival on at once. A pooled router splits its heads' inflows' total and attributes each
// member's amount to the heads by highest averages over their remaining inflows, member by member
// in member order (the pooled kernel's attribution).
void LossyForwardingEvaluator::split(std::int32_t ev) {
  const auto pos_global = ev - num_barriers_;
  const auto& it = net_.entry_occurrences[static_cast<std::size_t>(net_.router_entry_positions[static_cast<std::size_t>(pos_global)])];
  const auto& plan = *it.prepared_route;
  const auto pos = static_cast<std::size_t>(pos_global - it.first_router_position);
  if (plan.is_destination[pos] != 0) return;
  const auto b = static_cast<std::size_t>(plan.outgoing_begin[pos]);
  const auto e = static_cast<std::size_t>(plan.outgoing_end[pos]);
  outgoing_shares_.resize(e - b);
  if (!it.pooled) {
    split_at(plan, pos, inflow_[static_cast<std::size_t>(pos_global)], outgoing_shares_);
    for (std::size_t m = b; m < e; ++m) reach(it.first_occurrence + static_cast<std::int32_t>(m), outgoing_shares_[m - b]);
    return;
  }
  const auto [s0, s1] = slots(it, pos);
  const auto k = static_cast<std::size_t>(s1 - s0);
  remaining_.assign(slot_inflow_.begin() + s0, slot_inflow_.begin() + s1);
  Amount total = 0;
  for (Amount x : remaining_) total += x;
  split_at(plan, pos, total, outgoing_shares_);
  row_.resize(k);
  for (std::size_t m = b; m < e; ++m) {
    highest_averages(outgoing_shares_[m - b], remaining_, row_);
    const auto first = member_occurrence(it, m);
    for (std::size_t i = 0; i < k; ++i) {
      remaining_[i] -= row_[i];
      reach(first + static_cast<std::int32_t>(i), row_[i]);
    }
  }
}

void LossyForwardingEvaluator::evaluate_barrier(std::int32_t ev) {
  const auto l = static_cast<std::size_t>(barrier_link_[static_cast<std::size_t>(ev)]);
  const auto occs = std::span<const std::int32_t>(net_.occurrences_by_edge).subspan(
      static_cast<std::size_t>(net_.edge_occurrence_offsets[l]), static_cast<std::size_t>(net_.edge_occurrence_offsets[l + 1] - net_.edge_occurrence_offsets[l]));
  loss_rule(net_, occs, net_.available_capacity[l], arrival_, carried_, loss_);
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
  // The state x: every carried amount at the group's barriers, by EdgeId then demand, source and
  // entry (the barriers ascend by EdgeId; a link's sharing occurrences are in that order; a pooled
  // context's keep one coordinate per head).
  std::vector<std::int32_t> occs;
  for (auto b : barriers) {
    const auto l = static_cast<std::size_t>(barrier_link_[static_cast<std::size_t>(b)]);
    for (auto i = net_.edge_occurrence_offsets[l]; i < net_.edge_occurrence_offsets[l + 1]; ++i) occs.push_back(net_.occurrences_by_edge[static_cast<std::size_t>(i)]);
  }
  // What the group's routers hold before the rounds (from outside the group): a one-root split
  // route's inflow, a pooled router's slots.
  std::vector<std::pair<std::size_t, Amount>> inflow0;
  std::vector<std::pair<std::size_t, Amount>> slots0;
  for (auto r : routers) {
    const auto pos_global = static_cast<std::size_t>(r - num_barriers_);
    const auto& it = net_.entry_occurrences[static_cast<std::size_t>(net_.router_entry_positions[pos_global])];
    if (!it.pooled) {
      inflow0.emplace_back(pos_global, inflow_[pos_global]);
      continue;
    }
    const auto [s0, s1] = slots(it, pos_global - static_cast<std::size_t>(it.first_router_position));
    for (auto s = s0; s < s1; ++s) slots0.emplace_back(static_cast<std::size_t>(s), slot_inflow_[static_cast<std::size_t>(s)]);
  }
  auto restore = [&] {
    for (const auto& [i, x] : inflow0) inflow_[i] = x;
    for (const auto& [i, x] : slots0) slot_inflow_[i] = x;
  };
  std::vector<Amount> x(occs.size(), 0);
  std::vector<Amount> next(occs.size(), 0);
  std::vector<Amount> before(occs.size(), 0);  // x(r - 2)
  const auto limit = static_cast<std::int64_t>(64) * static_cast<std::int64_t>(barriers.size());
  auto replay = [&]() {
    restore();
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
      loss_rule(net_, span, net_.available_capacity[l], arrival_, carried_, loss_);
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
#if NETGRAPH_CORE_PRIVATE_HOOKS
  if (thread_lossy_counters().shift) {
    // The second planted fault of the check's test: one unit moved between two occurrences of one
    // barrier, from one whose tunnel meets no later barrier to one that dropped something, so the
    // link's total and every continuation stay consistent and only its loss rule is wrong.
    auto last_barrier = [&](std::int32_t o) {
      const auto& it = item_of(o);
      if (it.tunnel_position < 0) return false;
      for (auto k = o + 1; k < it.first_occurrence + it.num_occurrences; ++k) {
        if (event_of_occ(k) >= 0) return false;
      }
      return true;
    };
    std::size_t at = 0;
    bool moved = false;
    for (auto b : barriers) {
      const auto l = static_cast<std::size_t>(barrier_link_[static_cast<std::size_t>(b)]);
      const auto n = static_cast<std::size_t>(net_.edge_occurrence_offsets[l + 1] - net_.edge_occurrence_offsets[l]);
      for (std::size_t i = at; i < at + n && !moved; ++i) {
        for (std::size_t j = at; j < at + n && !moved; ++j) {
          if (i == j || x[i] == 0 || !last_barrier(occs[i]) || arrival_[static_cast<std::size_t>(occs[j])] <= x[j]) continue;
          --x[i];
          ++x[j];
          moved = true;
        }
      }
      at += n;
    }
  }
#endif
  // The solved state's forwarding, now also into the group's successors.
  iterating_ = false;
  current_component_ = -1;
  restore();
  for (std::size_t i = 0; i < occs.size(); ++i) carry(occs[i], x[i]);
  for (auto r : routers) split(r);
}

// The independent check of an iterated result (docs/semantics.md, "Losses on shared edges"):
// every entry's injection, every router's split of what reaches it (a pooled router's attribution
// to its heads too), every dark occurrence's drop, every link's loss rule and capacity, every
// continuation, and the conservation of each entry (each head of a pooled context: delivered +
// drops = its share).
void LossyForwardingEvaluator::verify() const {
  const auto fail = [] { throw std::runtime_error("lossy forwarding fixed-point verification failed"); };
  auto check_occ = [&](std::size_t o, Amount arrives) {
    if (arrival_[o] != arrives || carried_[o] > arrival_[o] || carried_[o] < 0) fail();
    if (net_.dark(o) && carried_[o] != 0) fail();
  };
  // Continuations and splits.
  std::vector<Amount> in(static_cast<std::size_t>(net_.num_router_positions), 0);
  std::vector<Amount> slot_in(net_.slot_heads.size(), 0);
  std::vector<Amount> part;
  std::vector<Amount> remaining;
  std::vector<Amount> row;
  for (std::size_t i = 0; i < net_.entry_occurrences.size(); ++i) {
    const auto& it = net_.entry_occurrences[i];
    const auto& share = shares_[static_cast<std::size_t>(it.demand)];
    if (it.tunnel_position >= 0) {
      const Amount s = share[static_cast<std::size_t>(it.entry_position)];
      Amount prev = s;
      Amount dropped = 0;
      for (std::int32_t j = 0; j < it.num_occurrences; ++j) {
        const auto o = static_cast<std::size_t>(it.first_occurrence + j);
        check_occ(o, prev);
        dropped += arrival_[o] - carried_[o];
        prev = carried_[o];
      }
      if (delivered_[i] != prev || delivered_[i] + dropped != s) fail();
      continue;
    }
    const auto& plan = *it.prepared_route;
    const auto base = static_cast<std::size_t>(it.first_router_position);
    if (!it.pooled) {
      const Amount s = share[static_cast<std::size_t>(it.entry_position)];
      in[base + static_cast<std::size_t>(it.source_position)] += s;
      Amount dropped = 0;
      Amount reached = 0;
      for (std::size_t pos = 0; pos < plan.topological_order.size(); ++pos) {
        const Amount amount = in[base + pos];
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
          check_occ(o, part[m - b]);
          dropped += arrival_[o] - carried_[o];
          in[base + static_cast<std::size_t>(plan.edge_destination_positions[m])] += carried_[o];
        }
      }
      if (delivered_[i] != reached || reached + dropped != s) fail();
      continue;
    }
    // A pooled context: per head its injection at its root's slot, its drops and what reaches the
    // destinations.
    const ForwardingView& view = *net_.views[static_cast<std::size_t>(it.demand)];
    const auto roots = pooled_root_positions(view);
    std::vector<Amount> dropped(roots.size(), 0);
    std::vector<Amount> reached(roots.size(), 0);
    auto slot_of = [&](std::size_t pos, std::int32_t head) {
      const auto [s0, s1] = slots(it, pos);
      const auto b = net_.slot_heads.begin();
      const auto at = std::lower_bound(b + s0, b + s1, head);
      if (at == b + s1 || *at != head) fail();
      return static_cast<std::size_t>(at - b);
    };
    for (std::size_t h = 0; h < roots.size(); ++h) {
      if (roots[h] < 0) continue;
      slot_in[slot_of(static_cast<std::size_t>(roots[h]), static_cast<std::int32_t>(h))] +=
          share[static_cast<std::size_t>(view.entries_for_source(h).front())];
    }
    for (std::size_t pos = 0; pos < plan.topological_order.size(); ++pos) {
      if (plan.is_destination[pos] != 0) continue;
      const auto [s0, s1] = slots(it, pos);
      remaining.assign(slot_in.begin() + s0, slot_in.begin() + s1);
      Amount total = 0;
      for (Amount x : remaining) total += x;
      const auto b = static_cast<std::size_t>(plan.outgoing_begin[pos]);
      const auto e = static_cast<std::size_t>(plan.outgoing_end[pos]);
      part.resize(e - b);
      split_at(plan, pos, total, part);
      row.resize(remaining.size());
      for (std::size_t m = b; m < e; ++m) {
        highest_averages(part[m - b], remaining, row);
        const auto first = static_cast<std::size_t>(member_occurrence(it, m));
        const auto q = static_cast<std::size_t>(plan.edge_destination_positions[m]);
        for (std::size_t k = 0; k < row.size(); ++k) {
          remaining[k] -= row[k];
          const auto o = first + k;
          const auto head = static_cast<std::size_t>(net_.slot_heads[static_cast<std::size_t>(s0) + k]);
          check_occ(o, row[k]);
          dropped[head] += arrival_[o] - carried_[o];
          if (plan.is_destination[q] != 0) {
            reached[head] += carried_[o];
          } else {
            slot_in[static_cast<std::size_t>(net_.occurrence_slots[o])] += carried_[o];
          }
        }
      }
    }
    const auto& delivered = head_delivered_[static_cast<std::size_t>(it.demand)];
    for (std::size_t h = 0; h < roots.size(); ++h) {
      if (roots[h] < 0) continue;
      const Amount s = share[static_cast<std::size_t>(view.entries_for_source(h).front())];
      if (delivered[h] != reached[h] || reached[h] + dropped[h] != s) fail();
    }
  }
  // Every link: its loss rule over the sharing arrivals (a transparent link passes everything on)
  // and its capacity.
  std::vector<Amount> expect(carried_.size(), 0);
  LossScratch scratch;
  for (std::size_t l = 0; l < net_.edges.size(); ++l) {
    const auto occs = std::span<const std::int32_t>(net_.occurrences_by_edge).subspan(
        static_cast<std::size_t>(net_.edge_occurrence_offsets[l]), static_cast<std::size_t>(net_.edge_occurrence_offsets[l + 1] - net_.edge_occurrence_offsets[l]));
    Amount sum = 0;
    for (auto o : occs) sum = saturating_add(sum, carried_[static_cast<std::size_t>(o)]);
    if (sum > net_.available_capacity[l]) fail();
    if (transparent_[l] != 0) {
      for (auto o : occs) {
        if (carried_[static_cast<std::size_t>(o)] != arrival_[static_cast<std::size_t>(o)]) fail();
      }
      continue;
    }
    loss_rule(net_, occs, net_.available_capacity[l], arrival_, expect, scratch);
    for (auto o : occs) {
      if (expect[static_cast<std::size_t>(o)] != carried_[static_cast<std::size_t>(o)]) fail();
    }
  }
}

std::vector<ProgramFlow> LossyForwardingEvaluator::run(std::span<const std::span<const Amount>> injections) {
  const std::size_t D = net_.views.size();
  std::vector<ProgramFlow> outs(D);
  shares_.resize(D);
  head_delivered_.resize(D);
  std::vector<Amount> source_drops_by_source;
  std::vector<Amount> w;
  std::vector<Amount> part;
  for (std::size_t d = 0; d < D; ++d) {
    const ForwardingView& view = *net_.views[d];
    const std::size_t k = view.num_sources();
    auto& out = outs[d];
    out = empty_flow(view, true);
    for (std::size_t h = 0; h < k; ++h) {
      out.offered_by_source[h] = injections[d][h];
      out.offered += injections[d][h];
    }
    out.admitted = out.offered;
    auto& share = shares_[d];
    entry_shares(view, injections[d], share, source_drops_by_source, w, part);
    for (std::size_t h = 0; h < k; ++h) {
      out.source_drops_by_source[h] = source_drops_by_source[h];
      out.source_drops += source_drops_by_source[h];
    }
    const auto& positions = net_.entry_occurrence_positions[d];
    for (std::size_t e = 0; e < share.size(); ++e) {
      out.entries[e].assigned = share[e];
      if (share[e] != 0 && positions[e] == -2) {
        out.entries[e].source_drops = share[e];
        out.source_drops += share[e];
      }
    }
    if (view.pooled) head_delivered_[d].assign(k, 0);
  }
  // A+: the no-drop loads first, then the barriers.
  std::vector<Amount> U0;
  std::vector<Amount> scratch;
  no_drop(net_, shares_, U0, &possible_, scratch);
  occurrence_source_positions_.assign(net_.occurrence_edge_positions.size(), -1);
  root_position_.assign(static_cast<std::size_t>(net_.num_router_positions), 0);
  for (const auto& it : net_.entry_occurrences) {
    if (it.tunnel_position >= 0) continue;
    const auto base = static_cast<std::size_t>(it.first_router_position);
    if (!it.pooled) {
      root_position_[base + static_cast<std::size_t>(it.source_position)] = 1;
      for (std::size_t pos = 0; pos < it.prepared_route->topological_order.size(); ++pos) {
        for (auto m = it.prepared_route->outgoing_begin[pos]; m < it.prepared_route->outgoing_end[pos]; ++m) occurrence_source_positions_[static_cast<std::size_t>(it.first_occurrence + m)] = static_cast<std::int32_t>(pos);
      }
      continue;
    }
    // A head whose share is 0 injects nothing.
    const ForwardingView& view = *net_.views[static_cast<std::size_t>(it.demand)];
    const auto roots = pooled_root_positions(view);
    for (std::size_t h = 0; h < roots.size(); ++h) {
      const Amount s = shares_[static_cast<std::size_t>(it.demand)][static_cast<std::size_t>(view.entries_for_source(h).front())];
      if (roots[h] >= 0 && s != 0) root_position_[base + static_cast<std::size_t>(roots[h])] = 1;
    }
    for (std::size_t m = 0; m < it.prepared_route->forwarding_edges.size(); ++m) {
      const auto tail = net_.member_tail_positions[it.first_member + m];
      for (auto o = member_occurrence(it, m); o < member_occurrence(it, m + 1); ++o) occurrence_source_positions_[static_cast<std::size_t>(o)] = tail;
    }
  }
  normalize_serial_links(U0);
  const auto n = net_.occurrence_edge_positions.size();
  arrival_.assign(n, 0);
  carried_.assign(n, 0);
  inflow_.assign(static_cast<std::size_t>(net_.num_router_positions), 0);
  slot_inflow_.assign(net_.slot_heads.size(), 0);
  delivered_.assign(net_.entry_occurrences.size(), 0);
  build_events();
  find_strong_components();
  // The injections: a tunnel's share enters its first hop, a split route's its root router, a
  // pooled context's heads' shares their roots' slots.
  for (std::size_t i = 0; i < net_.entry_occurrences.size(); ++i) {
    const auto& it = net_.entry_occurrences[i];
    const auto& share = shares_[static_cast<std::size_t>(it.demand)];
    if (it.tunnel_position >= 0) {
      reach(it.first_occurrence, share[static_cast<std::size_t>(it.entry_position)]);
    } else if (!it.pooled) {
      inflow_[static_cast<std::size_t>(it.first_router_position + it.source_position)] += share[static_cast<std::size_t>(it.entry_position)];
    } else {
      const ForwardingView& view = *net_.views[static_cast<std::size_t>(it.demand)];
      const auto roots = pooled_root_positions(view);
      for (std::size_t h = 0; h < roots.size(); ++h) {
        if (roots[h] < 0) continue;
        const auto [s0, s1] = slots(it, static_cast<std::size_t>(roots[h]));
        const auto b = net_.slot_heads.begin();
        const auto at = std::lower_bound(b + s0, b + s1, static_cast<std::int32_t>(h));
        slot_inflow_[static_cast<std::size_t>(at - b)] += share[static_cast<std::size_t>(view.entries_for_source(h).front())];
      }
    }
  }
  // The groups in dependency order, the smallest key first among the ready ones. Each group's
  // events, ascending, in one array (a group per offset range: a pooled plan has a group per router).
  const auto ns = static_cast<std::size_t>(num_scc_);
  std::vector<std::int32_t> member_offsets(ns + 1, 0);
  for (std::int32_t ev = 0; ev < num_events_; ++ev) ++member_offsets[static_cast<std::size_t>(scc_[static_cast<std::size_t>(ev)]) + 1];
  for (std::size_t g = 0; g < ns; ++g) member_offsets[g + 1] += member_offsets[g];
  std::vector<std::int32_t> member_events(static_cast<std::size_t>(num_events_));
  {
    std::vector<std::int32_t> fill(member_offsets.begin(), member_offsets.end() - 1);
    for (std::int32_t ev = 0; ev < num_events_; ++ev) member_events[static_cast<std::size_t>(fill[static_cast<std::size_t>(scc_[static_cast<std::size_t>(ev)])]++)] = ev;
  }
  auto members = [&](std::int32_t g) {
    const auto b = static_cast<std::size_t>(member_offsets[static_cast<std::size_t>(g)]);
    return std::span<const std::int32_t>(member_events).subspan(b, static_cast<std::size_t>(member_offsets[static_cast<std::size_t>(g) + 1]) - b);
  };
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
    const auto m = members(s);
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
    const auto m = members(s);
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
  // Split routes deliver what reaches their destinations (a pooled context per head).
  for (std::size_t i = 0; i < net_.entry_occurrences.size(); ++i) {
    const auto& it = net_.entry_occurrences[i];
    if (it.tunnel_position >= 0) continue;
    const auto& plan = *it.prepared_route;
    if (it.pooled) {
      auto& delivered = head_delivered_[static_cast<std::size_t>(it.demand)];
      for (std::int32_t j = 0; j < it.num_occurrences; ++j) {
        const auto o = static_cast<std::size_t>(it.first_occurrence + j);
        const auto m = static_cast<std::size_t>(net_.occurrence_route_positions[o]);
        if (plan.is_destination[static_cast<std::size_t>(plan.edge_destination_positions[m])] != 0) {
          delivered[static_cast<std::size_t>(net_.occurrence_heads[o])] += carried_[o];
        }
      }
      continue;
    }
    Amount d = 0;
    for (std::size_t pos = 0; pos < plan.topological_order.size(); ++pos) {
      if (plan.is_destination[pos] != 0) d += inflow_[static_cast<std::size_t>(it.first_router_position) + pos];
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
    verify();
  }
  for (std::size_t i = 0; i < net_.entry_occurrences.size(); ++i) {
    const auto& it = net_.entry_occurrences[i];
    auto& out = outs[static_cast<std::size_t>(it.demand)];
    if (it.pooled) {
      const ForwardingView& view = *net_.views[static_cast<std::size_t>(it.demand)];
      const auto& delivered = head_delivered_[static_cast<std::size_t>(it.demand)];
      for (std::size_t h = 0; h < view.num_sources(); ++h) {
        auto& f = out.entries[static_cast<std::size_t>(view.entries_for_source(h).front())];
        f.delivered = delivered[h];
        out.delivered += delivered[h];
      }
      // By member (ascending EdgeId), each occurrence to its head's entry.
      for (std::int32_t j = 0; j < it.num_occurrences; ++j) {
        const auto o = static_cast<std::size_t>(it.first_occurrence + j);
        auto& f = out.entries[static_cast<std::size_t>(view.entries_for_source(static_cast<std::size_t>(net_.occurrence_heads[o])).front())];
        const EdgeId e = net_.edges[static_cast<std::size_t>(net_.occurrence_edge_positions[o])];
        if (carried_[o] != 0) f.edge_flows.emplace_back(e, carried_[o]);
        if (arrival_[o] != carried_[o]) f.edge_drops.emplace_back(e, arrival_[o] - carried_[o]);
      }
      continue;
    }
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
  return outs;
}

// The joint evaluation of demands (one network, one evaluator), the hooks' counters reset.
std::vector<ProgramFlow> joint_lossy(std::span<const ForwardingView* const> views, std::span<const std::span<const Amount>> injections,
                                     const AvailableCapacityFn& free) {
  const ForwardingNetwork net = build_forwarding_network(views, free, true);
  LossyForwardingEvaluator run(net);
  return run.run(injections);
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
  thread_lossy_counters().event_structures = 0;
  thread_lossy_counters().barriers.clear();
  thread_lossy_counters().dependencies = 0;
#endif
  if (view.pooled) return pooled_lossy(view, free, injections);
  const ForwardingView* const views[] = {&view};
  const std::span<const Amount> inject[] = {injections};
  return std::move(joint_lossy(views, inject, free).front());
}

std::vector<ProgramFlow> lossy(std::span<const ForwardingDemand> demands, const AvailableCapacityFn& free) {
  std::vector<const ForwardingView*> views;
  std::vector<std::span<const Amount>> injections;
  views.reserve(demands.size());
  injections.reserve(demands.size());
  for (const auto& d : demands) {
    if (d.view == nullptr || d.injections.size() != d.view->num_sources()) {
      throw std::logic_error("forwarding: a demand's injections are not head-aligned");
    }
    views.push_back(d.view);
    injections.push_back(d.injections);
  }
#if NETGRAPH_CORE_PRIVATE_HOOKS
  thread_lossy_counters().groups = 0;
  thread_lossy_counters().rounds = 0;
  thread_lossy_counters().event_structures = 0;
  thread_lossy_counters().barriers.clear();
  thread_lossy_counters().dependencies = 0;
#endif
  return joint_lossy(views, injections, free);
}

} // namespace netgraph::core::detail
