/*
  The forwarding engine of programs (forwarding.hpp): views, the injection, the lossless bound
  and the loss rule.

  A pooled program (one split context) runs on its plan's kernels: allocate and the T_max search
  for the bound, forward() with its attribution for the losses. Every other program (tunnels and
  supplied split routes, each split route a one-root context) runs on the general network of one
  call: the union of its live entries' links, each link with its occurrences (a tunnel's hop, a
  split route's member), and per split route its routers. Each occurrence belongs to one source
  and one entry, so the loss rule at a link is two levels of highest averages.

  The evaluation (A+, 9.2): the no-drop loads U0 first; serial normalization marks the second link
  of every transparent pair (equal free capacity, every occurrence of the first continuing on the
  second and every occurrence of the second coming from the first, with nothing injected, absorbed,
  merged or split between them) as a pass-through; every other link with U0 > V is a barrier. The
  events are the barriers and the routers; an event waits for its predecessors (a router for its
  in-members' carried amounts, a barrier for all its arrivals, a tunnel's next barrier for its
  previous one). Tarjan's algorithm groups them into strongly connected components (SCCs), which
  run in dependency order, the smallest key first among the ready ones ((0, smallest barrier
  EdgeId) for a group with a barrier, else (1, its first event)). A cyclic group of b barriers is
  iterated from zero carried (9.4): one round replays the group's private forwarding from the
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
#include "netgraph/core/detail/pred_dag_builder.hpp"

namespace netgraph::core::detail {

#if NETGRAPH_CORE_PRIVATE_HOOKS
LossyCounters& thread_lossy_counters() noexcept {
  thread_local LossyCounters c;
  return c;
}
#endif

namespace {

// ================================================================ the view

bool has_black_holes(const Eligibility* dark) noexcept { return dark != nullptr && dark->uid() != 0; }

// The roster slots that reach an allowed destination over their basis under the mask: one
// reverse search from the allowed destinations over every allowed edge (topology), or per slot a
// forward search over its own allowed basis edges.
std::vector<std::uint8_t> replay_roster(const Graph& g, const ProgramData& p, const std::vector<std::uint8_t>& edge_ok,
                                        const std::vector<std::uint8_t>& node_ok,
                                        const std::vector<std::uint8_t>& is_dst) {
  const std::size_t k = p.roster.size();
  std::vector<std::uint8_t> keep(k, 0);
  const auto edge_src = g.edge_src();
  const auto edge_dst = g.edge_dst();
  if (p.topology) {
    const auto in_row = g.in_row_offsets();
    const auto in_adj = g.in_adj_edge_index();
    std::vector<std::uint8_t> seen(static_cast<std::size_t>(g.num_nodes()), 0);
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
        if (edge_ok[e] == 0) continue;
        const auto u = static_cast<std::size_t>(edge_src[e]);
        // A path stops at the first destination: a destination is not passed through.
        if (seen[u] != 0 || is_dst[u] != 0) continue;
        seen[u] = 1;
        stack.push_back(static_cast<NodeId>(u));
      }
    }
    for (std::size_t s = 0; s < k; ++s) {
      const auto h = static_cast<std::size_t>(p.roster[s]);
      keep[s] = node_ok[h] != 0 && seen[h] != 0 ? 1 : 0;
    }
    return keep;
  }
  std::vector<std::uint8_t> seen(static_cast<std::size_t>(g.num_nodes()), 0);
  std::vector<NodeId> touched;
  std::vector<NodeId> stack;
  std::vector<std::pair<std::span<const EdgeId>, bool>> parts;
  const auto row = g.row_offsets();
  for (std::size_t s = 0; s < k; ++s) {
    const NodeId h = p.roster[s];
    if (node_ok[static_cast<std::size_t>(h)] == 0) continue;
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
      if (edge_ok[e] == 0) return;
      const auto v = static_cast<std::size_t>(edge_dst[e]);
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

// An edge set under the mask, pruned to the walks from `src` to the allowed destinations, as a
// record bundle with its cost (an empty DAG and kUnreachable when nothing remains).
std::pair<std::shared_ptr<const PredDAG>, Cost> masked_route(const Graph& g, std::span<const EdgeId> edges, NodeId src,
                                                             std::span<const NodeId> dests, const Eligibility* elig) {
  const PredDAG dag = PredDAGBuilder::of_edges(g, edges);
  auto kept = prune(g, dag, src, dests, elig);
  if (kept.empty()) return {std::make_shared<const PredDAG>(PredDAGBuilder::empty(g)), kUnreachable};
  auto out = std::make_shared<const PredDAG>(PredDAGBuilder::of_edges(g, kept));
  return {out, route_cost(g, *out, src, dests)};
}

void finish_view(ForwardingView& v) {
  const ProgramData& p = *v.data;
  // The surviving entries per view head, ascending index.
  v.head_begin_.assign(v.slots_.size() + 1, 0);
  v.head_entries_.clear();
  for (std::size_t h = 0; h < v.slots_.size(); ++h) {
    for (std::int32_t e : p.entries_of(static_cast<std::size_t>(v.slots_[h]))) {
      if (v.status_[static_cast<std::size_t>(e)] != TunnelStatus::Failed) v.head_entries_.push_back(e);
    }
    v.head_begin_[h + 1] = static_cast<std::int32_t>(v.head_entries_.size());
  }
  // Black holes: the black-holed heads, and every up tunnel whose path reaches one.
  v.dark_head_.clear();
  if (v.dark != nullptr) {
    v.dark_head_.assign(v.slots_.size(), 0);
    for (std::size_t h = 0; h < v.slots_.size(); ++h) v.dark_head_[h] = v.dark->node_ok(v.head(h)) ? 0 : 1;
    for (std::size_t t = 0; t < p.tunnels.size(); ++t) {
      const auto e = static_cast<std::size_t>(p.tunnel_entry[t]);
      if (v.status_[e] != TunnelStatus::Up) continue;
      // An edge is black-holed when the Eligibility disallows it (itself or an endpoint).
      for (EdgeId x : path_of(p.tunnels[t])) {
        if (!v.dark->edge_ok(x)) {
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
  v.slots_.resize(p.roster.size());
  for (std::size_t s = 0; s < p.roster.size(); ++s) v.slots_[s] = static_cast<std::int32_t>(s);
  v.weights_ = p.head_weight;
  v.status_.resize(n);
  v.bundles_.resize(n);
  v.costs_.resize(n);
  v.root_of_.resize(n);
  for (std::size_t e = 0; e < n; ++e) {
    v.status_[e] = v.status(e);
    v.bundles_[e] = p.entries[e].bundle;
    v.costs_[e] = p.entries[e].cost;
    v.root_of_[e] = v.root_of(e);
  }
  v.plans_.clear();
  for (const auto& c : p.contexts) v.plans_.push_back(c.plan);
  v.borrowed = false;
}

} // namespace

ForwardingView make_view(const Graph& g, std::shared_ptr<const ForwardingProgram> program, const Eligibility* mask,
                         std::span<const NodeId> dests, const Eligibility* dark, bool records) {
  ForwardingView v;
  v.data = &ProgramAccess::data(*program);
  v.program = std::move(program);
  v.dark = has_black_holes(dark) ? dark : nullptr;
  const ProgramData& p = *v.data;
  const std::size_t n = p.entries.size();
  v.pooled = p.pooled;
  if (mask == nullptr || mask->uid() == 0) {
    // The program's own state: its roster, its routes and plans, borrowed from the program (a
    // tunnel's record bundle is made here; black holes change the statuses).
    v.borrowed = true;
    if (v.dark == nullptr && !(records && !p.tunnels.empty())) return v;
    materialize(v);
    for (std::size_t e = 0; records && e < n; ++e) {
      if (p.entries[e].tunnel < 0) continue;
      const auto path = path_of(p.tunnels[static_cast<std::size_t>(p.entries[e].tunnel)]);
      v.bundles_[e] = std::make_shared<const PredDAG>(path.empty() ? PredDAGBuilder::empty(g) : PredDAGBuilder::of_edges(g, path));
    }
    finish_view(v);
    return v;
  }
  v.status_.assign(n, TunnelStatus::Up);
  v.bundles_.resize(n);
  v.costs_.assign(n, kUnreachable);
  v.root_of_.assign(n, -1);
  // The replay (8.4): the mask and the program's own exclusions together.
  const auto num_nodes = static_cast<std::size_t>(g.num_nodes());
  const auto num_edges = static_cast<std::size_t>(g.num_edges());
  std::vector<std::uint8_t> node_ok(num_nodes, 0);
  std::vector<std::uint8_t> edge_ok(num_edges, 0);
  for (std::size_t x = 0; x < num_nodes; ++x) node_ok[x] = mask->node_ok(static_cast<NodeId>(x)) ? 1 : 0;
  for (NodeId x : p.excluded_nodes) node_ok[static_cast<std::size_t>(x)] = 0;
  for (std::size_t e = 0; e < num_edges; ++e) edge_ok[e] = mask->edge_ok(static_cast<EdgeId>(e)) ? 1 : 0;
  for (EdgeId e : p.excluded_edges) edge_ok[static_cast<std::size_t>(e)] = 0;
  std::vector<NodeId> live_dests;
  std::vector<std::uint8_t> is_dst(num_nodes, 0);
  for (NodeId d : dests) {
    if (node_ok[static_cast<std::size_t>(d)] == 0) continue;
    live_dests.push_back(d);
    is_dst[static_cast<std::size_t>(d)] = 1;
  }
  // The program's exclusions as an Eligibility of the replay (null when the program has none).
  std::optional<Eligibility> both;
  const Eligibility* replay = mask;
  if (!p.excluded_nodes.empty() || !p.excluded_edges.empty()) {
    both = Eligibility::intersect(*mask, g.eligibility(p.excluded_nodes, p.excluded_edges));
    replay = &*both;
  }
  const auto keep = replay_roster(g, p, edge_ok, node_ok, is_dst);
  for (std::size_t s = 0; s < p.roster.size(); ++s) {
    if (keep[s] == 0) continue;
    v.slots_.push_back(static_cast<std::int32_t>(s));
    if (!p.head_weight.empty()) v.weights_.push_back(p.head_weight[s]);
  }
  std::vector<std::int32_t> view_head(p.roster.size(), -1);
  for (std::size_t h = 0; h < v.slots_.size(); ++h) view_head[static_cast<std::size_t>(v.slots_[h])] = static_cast<std::int32_t>(h);
  // Entries: a tunnel fails when it touches an excluded node or edge (an edge is excluded with its
  // endpoints); a split route keeps its allowed edges on the walks from its head to an allowed
  // destination and fails when none remains. Every entry of a head outside the replay roster fails.
  // A tunnel that setup left down stays down with its status (replay never brings one up).
  for (std::size_t e = 0; e < n; ++e) {
    const auto& entry = p.entries[e];
    bool failed = view_head[static_cast<std::size_t>(entry.slot)] < 0;
    if (entry.tunnel >= 0) {
      const auto& t = p.tunnels[static_cast<std::size_t>(entry.tunnel)];
      if (t.status != TunnelStatus::Up) {
        v.status_[e] = t.status;
        if (records) v.bundles_[e] = std::make_shared<const PredDAG>(PredDAGBuilder::empty(g));
        continue;
      }
      for (EdgeId x : path_of(t)) failed = failed || edge_ok[static_cast<std::size_t>(x)] == 0;
      if (!failed) {
        if (records) v.bundles_[e] = std::make_shared<const PredDAG>(PredDAGBuilder::of_edges(g, path_of(t)));
        v.costs_[e] = entry.cost;
      }
    } else if (!failed) {
      // A split route's route: its bundle (its edges before weight pruning).
      const auto edges = PredDAGBuilder::ascending(*entry.bundle);
      auto [bundle, cost] = masked_route(g, edges, entry.head, live_dests, replay);
      failed = cost == kUnreachable;
      v.bundles_[e] = std::move(bundle);
      v.costs_[e] = cost;
    }
    if (failed) {
      v.status_[e] = TunnelStatus::Failed;
      v.bundles_[e] = records ? std::make_shared<const PredDAG>(PredDAGBuilder::empty(g)) : nullptr;
      v.costs_[e] = kUnreachable;
    }
  }
  // The split contexts restricted to the mask, rooted at their surviving entries' heads.
  v.plans_.reserve(p.contexts.size());
  std::vector<NodeId> roots;
  std::vector<std::int32_t> entries;
  for (const auto& c : p.contexts) {
    roots.clear();
    entries.clear();
    for (std::size_t r = 0; r < c.root_entry.size(); ++r) {
      const auto e = static_cast<std::size_t>(c.root_entry[r]);
      if (v.status_[e] == TunnelStatus::Failed) continue;
      roots.push_back(c.plan->roots[r]);
      entries.push_back(c.root_entry[r]);
    }
    v.plans_.push_back(std::make_shared<const KernelPlan>(restrict_plan(g, *c.plan, roots, live_dests, replay)));
    for (std::size_t r = 0; r < entries.size(); ++r) v.root_of_[static_cast<std::size_t>(entries[r])] = static_cast<std::int32_t>(r);
  }
  finish_view(v);
  return v;
}

ForwardingView with_black_holes(const ForwardingView& view, const Eligibility* dark) {
  ForwardingView v = view;
  v.dark = has_black_holes(dark) ? dark : nullptr;
  if (v.borrowed && v.dark == nullptr) return v;
  materialize(v);
  // The statuses again from the replay's (a black hole never fails a tunnel; a mask does).
  for (std::size_t t = 0; t < v.data->tunnels.size(); ++t) {
    const auto e = static_cast<std::size_t>(v.data->tunnel_entry[t]);
    if (v.status_[e] == TunnelStatus::Blackholed) v.status_[e] = TunnelStatus::Up;
  }
  finish_view(v);
  return v;
}

// ================================================================ the injection

namespace {

// x mod d for a 128-bit x and d > 0.
std::uint64_t mod_u128(const Counter128& x, std::uint64_t d) noexcept {
  const std::uint64_t q = div_u128(U128 {x.hi % d, x.lo}, d);
  return x.lo - q * d;
}

} // namespace

void head_shares(const ForwardingView& view, Amount x, std::span<Amount> out) {
  std::fill(out.begin(), out.end(), 0);
  if (out.empty()) return;
  if (view.weights().empty()) {
    (void)round_robin(x, out);
  } else {
    highest_averages(x, view.weights(), out);
  }
}

void head_increments(const ForwardingView& view, const Counter128& X, Amount v, std::span<Amount> out) {
  const std::size_t k = out.size();
  std::fill(out.begin(), out.end(), 0);
  if (k == 0) return;
  const auto weights = view.weights();
  if (weights.empty()) {
    // RR(X + v) - RR(X) depends only on r = X mod k: floor((r + v) / k) + [i < (r + v) mod k] - [i < r].
    const auto kk = static_cast<std::uint64_t>(k);
    const std::uint64_t r = mod_u128(X, kk);
    const auto total = static_cast<std::uint64_t>(v) + r;
    const auto base = static_cast<Amount>(total / kk);
    const std::uint64_t extra = total % kk;
    for (std::size_t i = 0; i < k; ++i) out[i] = base + (i < extra ? 1 : 0) - (i < r ? 1 : 0);
    return;
  }
  // H(X + c W) = H(X) + c w (the divisor rule's thresholds shift by c), so only r = X mod W and
  // r + v enter: H(X + v) - H(X) = floor((r + v) / W) w + H((r + v) mod W) - H(r). Captured
  // weights are a controller's solved supplies: W is its delivery, below the graph's capacity sum
  // (< 2^62), so W, r + v and c w fit 64 bits.
  std::uint64_t W = 0;
  for (Amount w : weights) W += static_cast<std::uint64_t>(w);
  if (W == 0) return;
  const std::uint64_t r = mod_u128(X, W);
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
// goes to `head_lost`.
void entry_shares(const ForwardingView& view, std::span<const Amount> heads, std::vector<Amount>& share,
                  std::vector<Amount>& head_lost, std::vector<Amount>& w, std::vector<Amount>& part) {
  const ProgramData& p = *view.data;
  share.assign(p.entries.size(), 0);
  head_lost.assign(heads.size(), 0);
  for (std::size_t h = 0; h < heads.size(); ++h) {
    const Amount x = heads[h];
    if (x == 0) continue;
    if (view.dark_head(h)) {
      head_lost[h] = x;
      continue;
    }
    const auto entries = view.head_entries(h);
    w.resize(entries.size());
    bool positive = false;  // a weight sum may exceed 2^63: only whether one weight is positive
    for (std::size_t i = 0; i < entries.size(); ++i) {
      w[i] = p.entries[static_cast<std::size_t>(entries[i])].weight;
      positive = positive || w[i] > 0;
    }
    if (!positive) {
      head_lost[h] = x;
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
  out.head_offered.assign(view.num_heads(), 0);
  if (lossy) out.head_lost.assign(view.num_heads(), 0);
  return out;
}

// ================================================================ the pooled form

// The free capacity of the pooled plan's members with the black holes at 0.
void pooled_view(const KernelPlan& plan, const ForwardingView& view, const FreeFn& free, std::vector<Amount>& out) {
  out.resize(plan.members.size());
  free(plan.members, out);
  if (view.dark == nullptr) return;
  for (std::size_t m = 0; m < plan.members.size(); ++m) {
    if (!view.dark->edge_ok(plan.members[m])) out[m] = 0;
  }
}

// The pooled form's per-root amounts as entry flows (root r carries entry root_entry).
void pooled_entries(const ForwardingView& view, const KernelPlan& plan, const Attribution* attr,
                    std::span<const Amount> alloc, std::span<const Amount> drops, ProgramFlow& out) {
  const ProgramData& p = *view.data;
  for (std::size_t h = 0; h < view.num_heads(); ++h) {
    const auto e = static_cast<std::size_t>(view.head_entries(h).front());
    auto& f = out.entries[e];
    (void)p;
    if (attr == nullptr) {
      // One root holds the whole forwarding (the attribution's degenerate form).
      f.carried = positive_pairs(plan.members, alloc);
      if (!drops.empty()) f.dropped = positive_pairs(plan.members, drops);
      continue;
    }
    for (const auto& x : attr->entries[h]) {
      const EdgeId edge = plan.members[static_cast<std::size_t>(x.member)];
      if (x.carried != 0) f.carried.emplace_back(edge, x.carried);
      if (x.dropped != 0) f.dropped.emplace_back(edge, x.dropped);
    }
  }
}

// The first black-holed head's index (the equal injection gives head i a unit from total i + 1),
// or the number of heads.
Amount first_dark(const ForwardingView& view) {
  for (std::size_t h = 0; h < view.num_heads(); ++h) {
    if (view.dark_head(h)) return static_cast<Amount>(h);
  }
  return std::numeric_limits<Amount>::max();
}

ProgramFlow pooled_lossless(const Graph& g, const ForwardingView& view, const FreeFn& free, Amount lo,
                            std::optional<Amount> hi) {
  const KernelPlan& plan = *view.plan(0);
  ProgramFlow out = empty_flow(view, false);
  const std::size_t k = view.num_heads();
  // A black-holed head may receive nothing: the total stays at most its index.
  const Amount cap = std::min(hi.value_or(std::numeric_limits<Amount>::max()), first_dark(view));
  std::vector<Amount> V;
  pooled_view(plan, view, free, V);
  std::vector<Amount> own(plan.members.size(), 0);
  if (lo != 0 && !allocate(plan, lo, own)) throw std::logic_error("forwarding: the admitted total does not fit");
  Amount total = lo;
  if (cap > lo && !plan.members.empty()) {
    const std::optional<Amount> request =
        cap == std::numeric_limits<Amount>::max() ? std::nullopt : std::optional<Amount>(cap - lo);
    auto step = lossless_step(g, plan, lo, V, own, request);
    total = step.total;
    own = std::move(step.alloc);
  }
  out.total = total;
  head_shares(view, total, out.head_offered);
  out.offered = total;
  out.delivered = total;
  for (std::size_t h = 0; h < k; ++h) {
    const auto e = static_cast<std::size_t>(view.head_entries(h).front());
    out.entries[e].assigned = out.head_offered[h];
    out.entries[e].delivered = out.head_offered[h];
  }
  if (k > 1 && total > 0) {
    // Several roots: forward the total again to attribute it.
    Attribution attr;
    const auto fw = forward(plan, Injection {out.head_offered, false}, {}, &attr);
    if (fw.carried != own) throw std::logic_error("forwarding: the attributed loads differ from the admitted ones");
    pooled_entries(view, plan, &attr, {}, {}, out);
  } else if (k == 1) {
    pooled_entries(view, plan, nullptr, own, {}, out);
  }
  return out;
}

ProgramFlow pooled_lossy(const ForwardingView& view, const FreeFn& free, std::span<const Amount> injections) {
  const KernelPlan& plan = *view.plan(0);
  ProgramFlow out = empty_flow(view, true);
  const std::size_t k = view.num_heads();
  std::vector<Amount> inject(injections.begin(), injections.end());
  for (std::size_t h = 0; h < k; ++h) {
    out.head_offered[h] = injections[h];
    out.offered += injections[h];
    if (view.dark_head(h)) {
      // A black-holed head loses its own share at the source.
      out.head_lost[h] = injections[h];
      out.source_loss += injections[h];
      inject[h] = 0;
    }
  }
  out.total = out.offered;
  std::vector<Amount> V;
  pooled_view(plan, view, free, V);
  Attribution attr;
  const auto fw = forward(plan, Injection {inject, false}, V, k > 1 ? &attr : nullptr);
  out.source_loss += fw.source_loss;
  out.delivered = fw.delivered;
  for (std::size_t h = 0; h < k; ++h) {
    const auto e = static_cast<std::size_t>(view.head_entries(h).front());
    auto& f = out.entries[e];
    f.assigned = inject[h];
    if (k > 1) {
      f.delivered = attr.delivered[h];
      f.lost = attr.source_loss[h];
    } else {
      f.delivered = fw.delivered;
      f.lost = fw.source_loss;
    }
  }
  pooled_entries(view, plan, k > 1 ? &attr : nullptr, fw.carried, fw.drops, out);
  return out;
}

// ================================================================ the general network

// The live entries of one call (surviving, positive weight, a usable route) and the union of their
// links. A tunnel's occurrences are its hops in path order, a split route's its plan's members.
struct Net {
  struct Item {
    std::int32_t entry {0};
    std::int32_t head {0};      // view head
    std::int32_t tunnel {-1};
    const KernelPlan* plan {nullptr};
    std::int32_t root_pos {-1};
    std::int32_t occ0 {0};      // its first occurrence
    std::int32_t nocc {0};
    std::int32_t pos0 {0};      // its first router (split routes: one per plan position)
  };
  std::vector<Item> items;
  std::vector<std::int32_t> item_of;  // per program entry: its item, or -1 (or -2: a routeless entry)
  std::vector<EdgeId> links;          // ascending
  std::vector<Amount> V;              // per link
  std::vector<std::int32_t> occ_link; // per occurrence
  std::vector<std::int32_t> occ_item;
  std::vector<std::int32_t> occ_pos;  // tunnel: hop; split route: member
  std::vector<std::int32_t> link_begin;  // per link + 1 into by_link
  std::vector<std::int32_t> by_link;     // occurrences, by (head NodeId, entry index)
  std::int32_t num_pos {0};
  std::vector<std::int32_t> pos_item;    // per router position: its item
  std::vector<std::int32_t> in_member;   // per router position: its one in-member, -1 none, -2 several
};

std::int32_t link_index(const Net& net, EdgeId e) {
  return static_cast<std::int32_t>(std::lower_bound(net.links.begin(), net.links.end(), e) - net.links.begin());
}

Net build_net(const ForwardingView& view, const FreeFn& free) {
  const ProgramData& p = *view.data;
  Net net;
  net.item_of.assign(p.entries.size(), -1);
  for (std::size_t h = 0; h < view.num_heads(); ++h) {
    for (const std::int32_t ei : view.head_entries(h)) {
      const auto e = static_cast<std::size_t>(ei);
      const auto& entry = p.entries[e];
      if (entry.weight == 0) continue;
      Net::Item it;
      it.entry = static_cast<std::int32_t>(e);
      it.head = static_cast<std::int32_t>(h);
      if (entry.tunnel >= 0) {
        it.tunnel = entry.tunnel;
        it.nocc = static_cast<std::int32_t>(path_of(p.tunnels[static_cast<std::size_t>(entry.tunnel)]).size());
      } else {
        it.plan = view.plan(static_cast<std::size_t>(entry.context)).get();
        const auto r = view.root_of(e);
        it.root_pos = r < 0 ? -1 : it.plan->root_pos[static_cast<std::size_t>(r)];
        if (it.root_pos < 0) {
          net.item_of[e] = -2;  // its share has no route: lost at the source
          continue;
        }
        it.nocc = static_cast<std::int32_t>(it.plan->members.size());
        it.pos0 = net.num_pos;
        net.num_pos += static_cast<std::int32_t>(it.plan->topo.size());
      }
      net.item_of[e] = static_cast<std::int32_t>(net.items.size());
      net.items.push_back(it);
    }
  }
  std::int32_t nocc = 0;
  for (auto& it : net.items) {
    it.occ0 = nocc;
    nocc += it.nocc;
    if (it.tunnel >= 0) {
      const auto edges = path_of(p.tunnels[static_cast<std::size_t>(it.tunnel)]);
      net.links.insert(net.links.end(), edges.begin(), edges.end());
    } else {
      net.links.insert(net.links.end(), it.plan->members.begin(), it.plan->members.end());
    }
  }
  std::sort(net.links.begin(), net.links.end());
  net.links.erase(std::unique(net.links.begin(), net.links.end()), net.links.end());
  net.V.resize(net.links.size());
  free(net.links, net.V);
  for (std::size_t l = 0; view.dark != nullptr && l < net.links.size(); ++l) {
    if (!view.dark->edge_ok(net.links[l])) net.V[l] = 0;
  }
  const auto n = static_cast<std::size_t>(nocc);
  net.occ_link.resize(n);
  net.occ_item.resize(n);
  net.occ_pos.resize(n);
  for (std::size_t i = 0; i < net.items.size(); ++i) {
    const auto& it = net.items[i];
    for (std::int32_t j = 0; j < it.nocc; ++j) {
      const auto o = static_cast<std::size_t>(it.occ0 + j);
      const EdgeId e = it.tunnel >= 0 ? path_of(p.tunnels[static_cast<std::size_t>(it.tunnel)])[static_cast<std::size_t>(j)]
                                      : it.plan->members[static_cast<std::size_t>(j)];
      net.occ_link[o] = link_index(net, e);
      net.occ_item[o] = static_cast<std::int32_t>(i);
      net.occ_pos[o] = j;
    }
  }
  net.pos_item.resize(static_cast<std::size_t>(net.num_pos));
  net.in_member.assign(static_cast<std::size_t>(net.num_pos), -1);
  for (std::size_t i = 0; i < net.items.size(); ++i) {
    const auto& it = net.items[i];
    if (it.tunnel >= 0) continue;
    for (std::size_t pos = 0; pos < it.plan->topo.size(); ++pos) net.pos_item[static_cast<std::size_t>(it.pos0) + pos] = static_cast<std::int32_t>(i);
    for (std::int32_t m = 0; m < it.nocc; ++m) {
      auto& slot = net.in_member[static_cast<std::size_t>(it.pos0 + it.plan->head_pos[static_cast<std::size_t>(m)])];
      slot = slot == -1 ? m : -2;
    }
  }
  net.link_begin.assign(net.links.size() + 1, 0);
  for (auto l : net.occ_link) ++net.link_begin[static_cast<std::size_t>(l) + 1];
  for (std::size_t l = 0; l < net.links.size(); ++l) net.link_begin[l + 1] += net.link_begin[l];
  net.by_link.resize(n);
  std::vector<std::int32_t> fill(net.link_begin.begin(), net.link_begin.end() - 1);
  // Occurrences are numbered by item, items by (head, entry index): each link's list comes out in
  // the loss rule's order (sources ascending, then entries ascending).
  for (std::size_t o = 0; o < n; ++o) net.by_link[static_cast<std::size_t>(fill[static_cast<std::size_t>(net.occ_link[o])]++)] = static_cast<std::int32_t>(o);
  return net;
}

// The no-drop loads of the entries' shares per link (and per occurrence when `occ` is given).
// A split route allocates its share over its plan from its root.
void no_drop(const Net& net, const ProgramData& p, std::span<const Amount> share, std::vector<Amount>& load,
             std::vector<Amount>* occ, std::vector<Amount>& scratch) {
  load.assign(net.links.size(), 0);
  if (occ != nullptr) occ->assign(net.occ_link.size(), 0);
  for (const auto& it : net.items) {
    const Amount s = share[static_cast<std::size_t>(it.entry)];
    if (s == 0) continue;
    if (it.tunnel >= 0) {
      for (std::int32_t j = 0; j < it.nocc; ++j) {
        const auto o = static_cast<std::size_t>(it.occ0 + j);
        load[static_cast<std::size_t>(net.occ_link[o])] += s;
        if (occ != nullptr) (*occ)[o] = s;
      }
      continue;
    }
    scratch.resize(it.plan->members.size());
    (void)allocate(*it.plan, s, scratch);
    for (std::int32_t j = 0; j < it.nocc; ++j) {
      const auto o = static_cast<std::size_t>(it.occ0 + j);
      load[static_cast<std::size_t>(net.occ_link[o])] += scratch[static_cast<std::size_t>(j)];
      if (occ != nullptr) (*occ)[o] = scratch[static_cast<std::size_t>(j)];
    }
  }
  (void)p;
}

// The router split of a split route's plan position p: round robin or highest averages over
// its out-members.
void split_at(const KernelPlan& plan, std::size_t pos, Amount amount, std::span<Amount> out) {
  const auto b = static_cast<std::size_t>(plan.out_begin[pos]);
  const auto e = static_cast<std::size_t>(plan.out_end[pos]);
  std::fill(out.begin(), out.end(), 0);
  if (amount == 0 || e == b) return;
  if (plan.placement == FlowPlacement::Weighted) {
    highest_averages(amount, std::span<const Amount>(plan.weights.data() + b, e - b), out);
  } else {
    (void)round_robin(amount, out);
  }
}

// The loss rule at one link: T = min(V, sum of the arrivals); highest averages over the sources'
// arrivals, then within a source over its entries'. `occs` lists the link's occurrences by
// (source, entry); writes their carried amounts.
void loss_rule(const Net& net, std::span<const std::int32_t> occs, Amount V, std::span<const Amount> arrival,
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
    const auto head = net.items[static_cast<std::size_t>(net.occ_item[static_cast<std::size_t>(occs[i])])].head;
    if (i == 0 || head != net.items[static_cast<std::size_t>(net.occ_item[static_cast<std::size_t>(occs[i - 1])])].head) {
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

ProgramFlow general_lossless(const ForwardingView& view, const FreeFn& free, Amount lo, std::optional<Amount> hi) {
  const ProgramData& p = *view.data;
  const Net net = build_net(view, free);
  const std::size_t k = view.num_heads();
  std::vector<Amount> heads(k);
  std::vector<Amount> share;
  std::vector<Amount> head_lost;
  std::vector<Amount> w;
  std::vector<Amount> part;
  std::vector<Amount> load;
  std::vector<Amount> scratch;
  auto fits = [&](Amount t) {
    head_shares(view, t, heads);
    // Units no head receives (every captured weight 0) have no route.
    Amount given = 0;
    for (Amount x : heads) given += x;
    if (given != t) return false;
    entry_shares(view, heads, share, head_lost, w, part);
    for (Amount x : head_lost) {
      if (x != 0) return false;
    }
    for (std::size_t e = 0; e < share.size(); ++e) {
      if (share[e] != 0 && net.item_of[e] == -2) return false;  // a share with no route
    }
    no_drop(net, p, share, load, nullptr, scratch);
    for (std::size_t l = 0; l < load.size(); ++l) {
      if (load[l] > net.V[l]) return false;
    }
    return true;
  };
  // An upper bound: every unit leaves on an entry's first link, which carries it without a drop;
  // with no head, nothing enters.
  Amount top = 0;
  {
    std::vector<std::uint8_t> first(net.links.size(), 0);
    for (const auto& it : net.items) {
      if (it.tunnel >= 0) {
        first[static_cast<std::size_t>(net.occ_link[static_cast<std::size_t>(it.occ0)])] = 1;
      } else {
        const auto rp = static_cast<std::size_t>(it.root_pos);
        for (auto m = it.plan->out_begin[rp]; m < it.plan->out_end[rp]; ++m) {
          first[static_cast<std::size_t>(net.occ_link[static_cast<std::size_t>(it.occ0 + m)])] = 1;
        }
      }
    }
    for (std::size_t l = 0; l < first.size(); ++l) {
      if (first[l] != 0) top += net.V[l];
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
  out.total = best;
  head_shares(view, best, out.head_offered);
  entry_shares(view, out.head_offered, share, head_lost, w, part);
  std::vector<Amount> occ;
  no_drop(net, p, share, load, &occ, scratch);
  for (std::size_t l = 0; l < load.size(); ++l) {
    if (load[l] > net.V[l]) throw std::logic_error("forwarding: the admitted total does not fit");
  }
  for (Amount x : out.head_offered) out.offered += x;
  out.delivered = out.offered;
  for (const auto& it : net.items) {
    auto& f = out.entries[static_cast<std::size_t>(it.entry)];
    f.assigned = share[static_cast<std::size_t>(it.entry)];
    f.delivered = f.assigned;
    for (std::int32_t j = 0; j < it.nocc; ++j) {
      const auto o = static_cast<std::size_t>(it.occ0 + j);
      if (occ[o] != 0) f.carried.emplace_back(net.links[static_cast<std::size_t>(net.occ_link[o])], occ[o]);
    }
    if (it.tunnel >= 0) sort_pairs(f.carried);
  }
  return out;
}

// ---------------------------------------------------------------- the general loss rule (A+)

class Lossy {
public:
  Lossy(const ForwardingView& view, const Net& net) : view_(view), p_(*view.data), net_(net) {}

  ProgramFlow run(std::span<const Amount> injections);

private:
  // Event ids: barriers [0, nb), routers [nb, nb + num_pos).
  [[nodiscard]] std::int32_t router(const Net::Item& it, std::int32_t pos) const noexcept {
    return nb_ + it.pos0 + pos;
  }
  [[nodiscard]] std::int32_t event_of_occ(std::int32_t o) const noexcept {
    return barrier_[static_cast<std::size_t>(net_.occ_link[static_cast<std::size_t>(o)])];
  }
  // A push into event `ev` is suppressed while a group iterates and ev lies outside it.
  [[nodiscard]] bool blocked(std::int32_t ev) const noexcept { return iterating_ && scc_[static_cast<std::size_t>(ev)] != current_; }

  void normalize(std::span<const Amount> U0);
  void build_events();
  void tarjan();
  void carry(std::int32_t o, Amount c);
  void split(std::int32_t ev);
  void evaluate_barrier(std::int32_t ev);
  void iterate(std::int32_t scc);
  void verify(std::span<const Amount> share) const;

  const ForwardingView& view_;
  const ProgramData& p_;
  const Net& net_;
  std::vector<std::int32_t> tail_;         // per occurrence of a split route: its member's tail position
  std::vector<std::uint8_t> transparent_;  // per link: the second link of a transparent pair
  std::vector<std::int32_t> barrier_;      // per link: its event, or -1 (a pass-through)
  std::vector<std::int32_t> barrier_link_;
  std::int32_t nb_ {0};
  std::int32_t num_events_ {0};
  std::vector<std::int32_t> adj_begin_;
  std::vector<std::int32_t> adj_;
  std::vector<std::int32_t> scc_;
  std::int32_t num_scc_ {0};
  bool iterating_ {false};
  std::int32_t current_ {-1};
  std::vector<Amount> arrival_;
  std::vector<Amount> carried_;
  std::vector<Amount> inflow_;    // per router position
  std::vector<Amount> delivered_; // per item
  std::vector<Amount> part_;
  std::vector<Amount> wa_;
  std::vector<Amount> wb_;
  bool iterated_ {false};
};

void Lossy::normalize(std::span<const Amount> U0) {
  const std::size_t L = net_.links.size();
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
  for (std::size_t o = 0; o < net_.occ_link.size(); ++o) {
    const auto l = static_cast<std::size_t>(net_.occ_link[o]);
    const auto& it = net_.items[static_cast<std::size_t>(net_.occ_item[o])];
    const std::int32_t j = net_.occ_pos[o];
    std::int32_t succ = -2;
    std::int32_t pred = -2;
    if (it.tunnel >= 0) {
      if (j + 1 < it.nocc) succ = net_.occ_link[static_cast<std::size_t>(it.occ0 + j + 1)];
      if (j > 0) pred = net_.occ_link[static_cast<std::size_t>(it.occ0 + j - 1)];
    } else {
      const auto& plan = *it.plan;
      // Continues: its head is no destination and no root, holds this member alone as its
      // in-member and has one out-member.
      const auto q = static_cast<std::size_t>(plan.head_pos[static_cast<std::size_t>(j)]);
      if (plan.terminal[q] == 0 && static_cast<std::int32_t>(q) != it.root_pos &&
          plan.out_end[q] - plan.out_begin[q] == 1 && net_.in_member[static_cast<std::size_t>(it.pos0) + q] == j) {
        succ = net_.occ_link[static_cast<std::size_t>(it.occ0 + plan.out_begin[q])];
      }
      // Comes from: its tail is no root, has one in-member and this one out-member.
      const auto t = static_cast<std::size_t>(tail_[static_cast<std::size_t>(it.occ0 + j)]);
      const auto in = net_.in_member[static_cast<std::size_t>(it.pos0) + t];
      if (static_cast<std::int32_t>(t) != it.root_pos && in >= 0 && plan.out_end[t] - plan.out_begin[t] == 1) {
        pred = net_.occ_link[static_cast<std::size_t>(it.occ0 + in)];
      }
    }
    merge(next[l], succ);
    merge(prev[l], pred);
  }
  for (std::size_t l = 0; l < L; ++l) {
    const auto n = next[l];
    if (n < 0) continue;
    const auto ni = static_cast<std::size_t>(n);
    if (prev[ni] != static_cast<std::int32_t>(l) || net_.V[l] != net_.V[ni]) continue;
    transparent_[ni] = 1;
  }
  barrier_.assign(L, -1);
  barrier_link_.clear();
  for (std::size_t l = 0; l < L; ++l) {
    if (transparent_[l] == 0 && U0[l] > net_.V[l]) {
      barrier_[l] = static_cast<std::int32_t>(barrier_link_.size());
      barrier_link_.push_back(static_cast<std::int32_t>(l));
    }
  }
  nb_ = static_cast<std::int32_t>(barrier_link_.size());
  num_events_ = nb_ + net_.num_pos;
}

void Lossy::build_events() {
  std::vector<std::pair<std::int32_t, std::int32_t>> edges;
  for (const auto& it : net_.items) {
    if (it.tunnel >= 0) {
      std::int32_t last = -1;
      for (std::int32_t j = 0; j < it.nocc; ++j) {
        const auto ev = event_of_occ(it.occ0 + j);
        if (ev < 0) continue;
        if (last >= 0) edges.emplace_back(last, ev);
        last = ev;
      }
      continue;
    }
    const auto& plan = *it.plan;
    for (std::size_t pos = 0; pos < plan.topo.size(); ++pos) {
      for (auto m = plan.out_begin[pos]; m < plan.out_end[pos]; ++m) {
        const auto o = it.occ0 + m;
        const auto from = router(it, static_cast<std::int32_t>(pos));
        const auto to = router(it, plan.head_pos[static_cast<std::size_t>(m)]);
        const auto ev = event_of_occ(o);
        if (ev >= 0) {
          edges.emplace_back(from, ev);
          edges.emplace_back(ev, to);
        } else {
          edges.emplace_back(from, to);
        }
      }
    }
  }
  adj_begin_.assign(static_cast<std::size_t>(num_events_) + 1, 0);
  for (const auto& [a, b] : edges) ++adj_begin_[static_cast<std::size_t>(a) + 1];
  for (std::size_t v = 0; v < static_cast<std::size_t>(num_events_); ++v) adj_begin_[v + 1] += adj_begin_[v];
  adj_.resize(edges.size());
  std::vector<std::int32_t> fill(adj_begin_.begin(), adj_begin_.end() - 1);
  for (const auto& [a, b] : edges) adj_[static_cast<std::size_t>(fill[static_cast<std::size_t>(a)]++)] = b;
}

// Tarjan's strongly connected components, iterative (no recursion depth).
void Lossy::tarjan() {
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
    work.emplace_back(static_cast<std::int32_t>(s), adj_begin_[s]);
    index[s] = low[s] = counter++;
    stack.push_back(static_cast<std::int32_t>(s));
    on[s] = 1;
    while (!work.empty()) {
      auto& [v, k] = work.back();
      const auto vi = static_cast<std::size_t>(v);
      if (k < adj_begin_[vi + 1]) {
        const auto w = static_cast<std::size_t>(adj_[static_cast<std::size_t>(k)]);
        ++k;
        if (index[w] < 0) {
          index[w] = low[w] = counter++;
          stack.push_back(static_cast<std::int32_t>(w));
          on[w] = 1;
          work.emplace_back(static_cast<std::int32_t>(w), adj_begin_[w]);
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
void Lossy::carry(std::int32_t o, Amount c) {
  carried_[static_cast<std::size_t>(o)] = c;
  const auto& it = net_.items[static_cast<std::size_t>(net_.occ_item[static_cast<std::size_t>(o)])];
  const std::int32_t j = net_.occ_pos[static_cast<std::size_t>(o)];
  if (it.tunnel < 0) {
    const auto q = it.plan->head_pos[static_cast<std::size_t>(j)];
    const auto ev = router(it, q);
    if (!blocked(ev)) inflow_[static_cast<std::size_t>(it.pos0 + q)] += c;
    return;
  }
  for (std::int32_t h = j + 1;; ++h) {
    if (h == it.nocc) {
      // The tunnel's end is no event of a group.
      if (!iterating_) delivered_[static_cast<std::size_t>(net_.occ_item[static_cast<std::size_t>(o)])] = c;
      return;
    }
    const auto o2 = it.occ0 + h;
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
void Lossy::split(std::int32_t ev) {
  const auto pos_global = ev - nb_;
  const auto& it = net_.items[static_cast<std::size_t>(net_.pos_item[static_cast<std::size_t>(pos_global)])];
  const auto& plan = *it.plan;
  const auto pos = static_cast<std::size_t>(pos_global - it.pos0);
  if (plan.terminal[pos] != 0) return;
  const auto b = static_cast<std::size_t>(plan.out_begin[pos]);
  const auto e = static_cast<std::size_t>(plan.out_end[pos]);
  part_.resize(e - b);
  split_at(plan, pos, inflow_[static_cast<std::size_t>(pos_global)], part_);
  for (std::size_t m = b; m < e; ++m) {
    const auto o = it.occ0 + static_cast<std::int32_t>(m);
    const Amount x = part_[m - b];
    const auto bev = event_of_occ(o);
    if (bev >= 0) {
      if (!blocked(bev)) arrival_[static_cast<std::size_t>(o)] = x;
      continue;
    }
    arrival_[static_cast<std::size_t>(o)] = x;
    carry(o, x);
  }
}

void Lossy::evaluate_barrier(std::int32_t ev) {
  const auto l = static_cast<std::size_t>(barrier_link_[static_cast<std::size_t>(ev)]);
  const auto occs = std::span<const std::int32_t>(net_.by_link).subspan(
      static_cast<std::size_t>(net_.link_begin[l]), static_cast<std::size_t>(net_.link_begin[l + 1] - net_.link_begin[l]));
  loss_rule(net_, occs, net_.V[l], arrival_, carried_, wa_, wb_);
  for (auto o : occs) carry(o, carried_[static_cast<std::size_t>(o)]);
}

// A cyclic group (9.4): rounds from zero carried at its barriers.
void Lossy::iterate(std::int32_t scc) {
  iterated_ = true;
  current_ = scc;
  std::vector<std::int32_t> barriers;
  std::vector<std::int32_t> routers;
  for (std::int32_t ev = 0; ev < num_events_; ++ev) {
    if (scc_[static_cast<std::size_t>(ev)] != scc) continue;
    (ev < nb_ ? barriers : routers).push_back(ev);
  }
  // The group's routers in an order that respects router -> router edges (barriers cut).
  {
    std::vector<std::int32_t> indeg(static_cast<std::size_t>(num_events_), 0);
    for (auto r : routers) {
      for (auto k = adj_begin_[static_cast<std::size_t>(r)]; k < adj_begin_[static_cast<std::size_t>(r) + 1]; ++k) {
        const auto w = adj_[static_cast<std::size_t>(k)];
        if (w >= nb_ && scc_[static_cast<std::size_t>(w)] == scc) ++indeg[static_cast<std::size_t>(w)];
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
      for (auto k = adj_begin_[static_cast<std::size_t>(r)]; k < adj_begin_[static_cast<std::size_t>(r) + 1]; ++k) {
        const auto w = adj_[static_cast<std::size_t>(k)];
        if (w >= nb_ && scc_[static_cast<std::size_t>(w)] == scc && --indeg[static_cast<std::size_t>(w)] == 0) ready.push_back(w);
      }
    }
    routers.swap(order);
  }
  // The state x: every carried amount at the group's barriers, by EdgeId then source then entry
  // (the barriers ascend by EdgeId; a link's occurrences are in that order).
  std::vector<std::int32_t> occs;
  for (auto b : barriers) {
    const auto l = static_cast<std::size_t>(barrier_link_[static_cast<std::size_t>(b)]);
    for (auto i = net_.link_begin[l]; i < net_.link_begin[l + 1]; ++i) occs.push_back(net_.by_link[static_cast<std::size_t>(i)]);
  }
  std::vector<Amount> inflow0(routers.size());
  for (std::size_t i = 0; i < routers.size(); ++i) inflow0[i] = inflow_[static_cast<std::size_t>(routers[i] - nb_)];
  std::vector<Amount> x(occs.size(), 0);
  std::vector<Amount> next(occs.size(), 0);
  std::vector<Amount> before(occs.size(), 0);  // x(r - 2)
  const auto limit = static_cast<std::int64_t>(64) * static_cast<std::int64_t>(barriers.size());
  auto replay = [&]() {
    for (std::size_t i = 0; i < routers.size(); ++i) inflow_[static_cast<std::size_t>(routers[i] - nb_)] = inflow0[i];
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
      const auto span = std::span<const std::int32_t>(net_.by_link).subspan(
          static_cast<std::size_t>(net_.link_begin[l]), static_cast<std::size_t>(net_.link_begin[l + 1] - net_.link_begin[l]));
      loss_rule(net_, span, net_.V[l], arrival_, carried_, wa_, wb_);
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
  current_ = -1;
  for (std::size_t i = 0; i < routers.size(); ++i) inflow_[static_cast<std::size_t>(routers[i] - nb_)] = inflow0[i];
  for (std::size_t i = 0; i < occs.size(); ++i) carry(occs[i], x[i]);
  for (auto r : routers) split(r);
}

// The independent check of an iterated result (9.4 step 5): every entry's injection, every
// router's split of what reaches it, every link's loss rule and capacity, every continuation, and
// each entry's conservation (delivered + drops = its share).
void Lossy::verify(std::span<const Amount> share) const {
  const auto fail = [] { throw std::runtime_error("lossy forwarding fixed-point verification failed"); };
  // Continuations and splits.
  std::vector<Amount> in(static_cast<std::size_t>(net_.num_pos), 0);
  std::vector<Amount> part;
  for (std::size_t i = 0; i < net_.items.size(); ++i) {
    const auto& it = net_.items[i];
    const Amount s = share[static_cast<std::size_t>(it.entry)];
    if (it.tunnel >= 0) {
      Amount prev = s;
      Amount dropped = 0;
      for (std::int32_t j = 0; j < it.nocc; ++j) {
        const auto o = static_cast<std::size_t>(it.occ0 + j);
        if (arrival_[o] != prev || carried_[o] > arrival_[o] || carried_[o] < 0) fail();
        dropped += arrival_[o] - carried_[o];
        prev = carried_[o];
      }
      if (delivered_[i] != prev || delivered_[i] + dropped != s) fail();
      continue;
    }
    const auto& plan = *it.plan;
    in[static_cast<std::size_t>(it.pos0 + it.root_pos)] += s;
    Amount dropped = 0;
    Amount reached = 0;
    for (std::size_t pos = 0; pos < plan.topo.size(); ++pos) {
      const Amount amount = in[static_cast<std::size_t>(it.pos0) + pos];
      if (plan.terminal[pos] != 0) {
        reached += amount;
        continue;
      }
      const auto b = static_cast<std::size_t>(plan.out_begin[pos]);
      const auto e = static_cast<std::size_t>(plan.out_end[pos]);
      part.resize(e - b);
      split_at(plan, pos, amount, part);
      for (std::size_t m = b; m < e; ++m) {
        const auto o = static_cast<std::size_t>(it.occ0) + m;
        if (arrival_[o] != part[m - b] || carried_[o] > arrival_[o] || carried_[o] < 0) fail();
        dropped += arrival_[o] - carried_[o];
        in[static_cast<std::size_t>(it.pos0 + plan.head_pos[m])] += carried_[o];
      }
    }
    if (delivered_[i] != reached || reached + dropped != s) fail();
  }
  // Every link: its loss rule over the arrivals (a transparent link passes everything on) and
  // its capacity.
  std::vector<Amount> expect(carried_.size(), 0);
  std::vector<Amount> wa;
  std::vector<Amount> wb;
  for (std::size_t l = 0; l < net_.links.size(); ++l) {
    const auto occs = std::span<const std::int32_t>(net_.by_link).subspan(
        static_cast<std::size_t>(net_.link_begin[l]), static_cast<std::size_t>(net_.link_begin[l + 1] - net_.link_begin[l]));
    Amount sum = 0;
    for (auto o : occs) sum += carried_[static_cast<std::size_t>(o)];
    if (sum > net_.V[l]) fail();
    if (transparent_[l] != 0) {
      for (auto o : occs) {
        if (carried_[static_cast<std::size_t>(o)] != arrival_[static_cast<std::size_t>(o)]) fail();
      }
      continue;
    }
    loss_rule(net_, occs, net_.V[l], arrival_, expect, wa, wb);
    for (auto o : occs) {
      if (expect[static_cast<std::size_t>(o)] != carried_[static_cast<std::size_t>(o)]) fail();
    }
  }
}

ProgramFlow Lossy::run(std::span<const Amount> injections) {
  const std::size_t k = view_.num_heads();
  ProgramFlow out = empty_flow(view_, true);
  for (std::size_t h = 0; h < k; ++h) {
    out.head_offered[h] = injections[h];
    out.offered += injections[h];
  }
  out.total = out.offered;
  std::vector<Amount> share;
  std::vector<Amount> head_lost;
  std::vector<Amount> w;
  std::vector<Amount> part;
  entry_shares(view_, injections, share, head_lost, w, part);
  for (std::size_t h = 0; h < k; ++h) {
    out.head_lost[h] = head_lost[h];
    out.source_loss += head_lost[h];
  }
  for (std::size_t e = 0; e < share.size(); ++e) {
    out.entries[e].assigned = share[e];
    if (share[e] != 0 && net_.item_of[e] == -2) {
      out.entries[e].lost = share[e];
      out.source_loss += share[e];
    }
  }
  // A+: the no-drop loads first, then the barriers.
  std::vector<Amount> U0;
  std::vector<Amount> scratch;
  no_drop(net_, p_, share, U0, nullptr, scratch);
  tail_.assign(net_.occ_link.size(), -1);
  for (const auto& it : net_.items) {
    if (it.tunnel >= 0) continue;
    for (std::size_t pos = 0; pos < it.plan->topo.size(); ++pos) {
      for (auto m = it.plan->out_begin[pos]; m < it.plan->out_end[pos]; ++m) tail_[static_cast<std::size_t>(it.occ0 + m)] = static_cast<std::int32_t>(pos);
    }
  }
  normalize(U0);
  const auto n = net_.occ_link.size();
  arrival_.assign(n, 0);
  carried_.assign(n, 0);
  inflow_.assign(static_cast<std::size_t>(net_.num_pos), 0);
  delivered_.assign(net_.items.size(), 0);
  build_events();
  tarjan();
  // The injections: a tunnel's share enters its first hop, a split route's its root router.
  for (std::size_t i = 0; i < net_.items.size(); ++i) {
    const auto& it = net_.items[i];
    const Amount s = share[static_cast<std::size_t>(it.entry)];
    if (it.tunnel >= 0) {
      const auto o = it.occ0;
      arrival_[static_cast<std::size_t>(o)] = s;
      if (event_of_occ(o) < 0) carry(o, s);
    } else {
      inflow_[static_cast<std::size_t>(it.pos0 + it.root_pos)] += s;
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
    for (auto kk = adj_begin_[static_cast<std::size_t>(ev)]; kk < adj_begin_[static_cast<std::size_t>(ev) + 1]; ++kk) {
      const auto b = scc_[static_cast<std::size_t>(adj_[static_cast<std::size_t>(kk)])];
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
      if (ev < nb_) edge = std::min<std::int64_t>(edge, net_.links[static_cast<std::size_t>(barrier_link_[static_cast<std::size_t>(ev)])]);
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
    } else if (m.front() < nb_) {
      evaluate_barrier(m.front());
    } else {
      split(m.front());
    }
    for (auto ev : m) {
      for (auto kk = adj_begin_[static_cast<std::size_t>(ev)]; kk < adj_begin_[static_cast<std::size_t>(ev) + 1]; ++kk) {
        const auto b = scc_[static_cast<std::size_t>(adj_[static_cast<std::size_t>(kk)])];
        if (b != s && --indeg[static_cast<std::size_t>(b)] == 0) ready.push(key(b));
      }
    }
  }
  // Split routes deliver what reaches their destinations.
  for (std::size_t i = 0; i < net_.items.size(); ++i) {
    const auto& it = net_.items[i];
    if (it.tunnel >= 0) continue;
    Amount d = 0;
    for (std::size_t pos = 0; pos < it.plan->topo.size(); ++pos) {
      if (it.plan->terminal[pos] != 0) d += inflow_[static_cast<std::size_t>(it.pos0) + pos];
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
  for (std::size_t i = 0; i < net_.items.size(); ++i) {
    const auto& it = net_.items[i];
    auto& f = out.entries[static_cast<std::size_t>(it.entry)];
    f.delivered = delivered_[i];
    out.delivered += delivered_[i];
    for (std::int32_t j = 0; j < it.nocc; ++j) {
      const auto o = static_cast<std::size_t>(it.occ0 + j);
      const EdgeId e = net_.links[static_cast<std::size_t>(net_.occ_link[o])];
      if (carried_[o] != 0) f.carried.emplace_back(e, carried_[o]);
      if (arrival_[o] != carried_[o]) f.dropped.emplace_back(e, arrival_[o] - carried_[o]);
    }
    if (it.tunnel >= 0) {
      sort_pairs(f.carried);
      sort_pairs(f.dropped);
    }
  }
  return out;
}

} // namespace

ProgramFlow lossless(const Graph& g, const ForwardingView& view, const FreeFn& free, Amount lo, std::optional<Amount> hi) {
  if (view.pooled) return pooled_lossless(g, view, free, lo, hi);
  return general_lossless(view, free, lo, hi);
}

ProgramFlow lossy(const Graph& g, const ForwardingView& view, const FreeFn& free, std::span<const Amount> injections) {
  (void)g;
  if (injections.size() != view.num_heads()) throw std::logic_error("forwarding: the injections are not head-aligned");
#if NETGRAPH_CORE_PRIVATE_HOOKS
  thread_lossy_counters().groups = 0;
  thread_lossy_counters().rounds = 0;
#endif
  if (view.pooled) return pooled_lossy(view, free, injections);
  const Net net = build_net(view, free);
  Lossy run(view, net);
  return run.run(injections);
}

} // namespace netgraph::core::detail
