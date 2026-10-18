/*
  ForwardingProgram: the value, from_routes (8.2), the programs of the computing models and the
  comparison of the change flags (10.4).

  Every program keeps its entries in ascending index, per roster slot its entries, its split
  contexts (one plan each; a CostOnly program's is pooled over its roster) and its replay basis.
  Nothing in a program changes after construction; policies share it by reference.
*/
#include "netgraph/core/program.hpp"

#include <algorithm>
#include <stdexcept>
#include <string>
#include <utility>

#include "netgraph/core/detail/edge_amounts.hpp"
#include "netgraph/core/detail/endpoints.hpp"
#include "netgraph/core/detail/enum_domain.hpp"
#include "netgraph/core/detail/program_data.hpp"
#include "netgraph/core/detail/pred_dag_builder.hpp"

namespace netgraph::core {

ForwardingProgram::ForwardingProgram(const detail::ProgramData* data) noexcept : data_(data) {}
ForwardingProgram::~ForwardingProgram() = default;

std::uint64_t ForwardingProgram::graph_uid() const noexcept { return data_->graph_uid; }
std::span<const NodeId> ForwardingProgram::sources() const noexcept { return data_->sources; }
std::span<const NodeId> ForwardingProgram::destinations() const noexcept { return data_->destinations; }
std::span<const Tunnel> ForwardingProgram::tunnels() const noexcept { return data_->tunnels; }
const std::optional<std::vector<std::pair<NodeId, Amount>>>& ForwardingProgram::source_weights() const noexcept {
  return data_->source_weights;
}

namespace detail {

namespace {

[[noreturn]] void invalid(const std::string& what) { throw std::invalid_argument(what); }

std::string route_at(std::size_t i) { return "routes[" + std::to_string(i) + "]: "; }

} // namespace

void ProgramData::lay_out(std::span<const NodeId> src, std::span<const NodeId> dst,
                          std::optional<std::span<const NodeId>> given, std::span<const NodeId> xn,
                          std::span<const EdgeId> xe) {
  const std::size_t n = entries.size();
  const std::size_t room = given.has_value() ? given->size() : n;  // the roster's room
  const std::size_t lists = pooled ? 0 : room + 2 + n;              // slot_first (counting sort) + slot_entries
  std::size_t roots = 0;
  for (const auto& c : contexts) roots += c.plan->roots.size();
  // The block, in this order: sources, destinations, roster, excluded nodes and edges, the per-slot
  // lists, the per-tunnel entries, the contexts' root entries.
  std::vector<std::int32_t> store(src.size() + dst.size() + room + xn.size() + xe.size() + lists + tunnels.size() +
                                  roots);
  std::int32_t* at = store.data();
  auto put = [&at](auto list) {
    std::int32_t* first = at;
    at = std::copy(list.begin(), list.end(), at);
    return first;
  };
  NodeId* s0 = put(src);
  NodeId* d0 = put(dst);
  NodeId* r0 = at;
  std::size_t k = 0;
  if (given.has_value()) {
    put(*given);
    k = given->size();
  } else {
    for (const auto& entry : entries) *at++ = entry.head;
    std::sort(r0, at);
    k = static_cast<std::size_t>(std::unique(r0, at) - r0);
    at = r0 + room;
  }
  NodeId* xn0 = put(xn);
  EdgeId* xe0 = put(xe);
  std::int32_t* first = at;
  std::int32_t* order = first + (pooled ? 0 : room + 2);
  std::int32_t* by_tunnel = first + lists;
  const std::span<const NodeId> ros(r0, k);
  for (auto& entry : entries) {
    entry.slot = static_cast<std::int32_t>(std::lower_bound(ros.begin(), ros.end(), entry.head) - ros.begin());
  }
  if (pooled) {
    // A pooled program's slots carry one entry each (entries_of): no lists, local index 0.
    for (auto& entry : entries) entry.local_index = 0;
  } else {
    // A stable counting sort by slot: counted at slot + 2, summed, placed through slot + 1 (which
    // then holds the slot's end, the next slot's start); ascending index within a slot.
    for (const auto& entry : entries) ++first[static_cast<std::size_t>(entry.slot) + 2];
    for (std::size_t x = 2; x < k + 2; ++x) first[x] += first[x - 1];
    for (std::size_t e = 0; e < n; ++e) order[first[static_cast<std::size_t>(entries[e].slot) + 1]++] = static_cast<std::int32_t>(e);
    for (std::size_t x = 0; x < k; ++x) {
      for (auto j = first[x]; j < first[x + 1]; ++j) {
        auto& entry = entries[static_cast<std::size_t>(order[j])];
        entry.local_index = j - first[x];
        if (entry.tunnel >= 0) tunnels[static_cast<std::size_t>(entry.tunnel)].local_index = entry.local_index;
      }
    }
  }
  std::fill(by_tunnel, by_tunnel + tunnels.size() + roots, -1);
  for (std::size_t e = 0; e < n; ++e) {
    if (entries[e].tunnel >= 0) by_tunnel[entries[e].tunnel] = static_cast<std::int32_t>(e);
  }
  // Each context's root entries: the entry of that context whose head is the root (a CostOnly
  // plan's roots are its roster, each with its entry; a supplied split route's is its head).
  std::int32_t* root0 = by_tunnel + tunnels.size();
  for (std::size_t c = 0, at_root = 0; c < contexts.size(); ++c) {
    contexts[c].root_entry = {root0 + at_root, contexts[c].plan->roots.size()};
    at_root += contexts[c].plan->roots.size();
  }
  for (std::size_t e = 0; e < n; ++e) {
    const auto c = entries[e].context;
    if (c < 0) continue;
    const auto& ctx = contexts[static_cast<std::size_t>(c)];
    const auto& plan_roots = ctx.plan->roots;
    const auto r = std::lower_bound(plan_roots.begin(), plan_roots.end(), entries[e].head) - plan_roots.begin();
    root0[(ctx.root_entry.data() - root0) + r] = static_cast<std::int32_t>(e);
  }
  // The views, into the block now owned (a vector's move keeps its elements where they are).
  store_ = std::move(store);
  sources = {s0, src.size()};
  destinations = {d0, dst.size()};
  roster = ros;
  excluded_nodes = {xn0, xn.size()};
  excluded_edges = {xe0, xe.size()};
  slot_first = pooled ? std::span<const std::int32_t> {} : std::span<const std::int32_t>(first, k + 1);
  slot_entries = pooled ? std::span<const std::int32_t> {} : std::span<const std::int32_t>(order, n);
  tunnel_entry = {by_tunnel, tunnels.size()};
}

std::shared_ptr<const ForwardingProgram> ProgramAccess::make(ProgramData&& data) {
  auto block = std::make_shared<const ProgramBlock>(std::move(data));
  return {block, &block->program};
}

void check_route_options(bool any_dag, std::optional<FlowPlacement> placement, std::optional<WeightBasis> weights,
                         bool edge_weights, bool residual) {
  if (!any_dag) {
    if (placement.has_value()) invalid("placement applies only to PredDAG routes");
    if (weights.has_value()) invalid("weights applies only to PredDAG routes");
    if (edge_weights) invalid("edge_weights applies only to PredDAG routes");
    if (residual) invalid("residual applies only to PredDAG routes");
    return;
  }
  if (placement.has_value()) check_member(*placement, "placement");
  if (weights.has_value()) check_member(*weights, "weights");
  if (placement == FlowPlacement::MaxFlow) invalid("from_routes placement must be EQUAL_BALANCED or WEIGHTED");
  const bool weighted = placement == FlowPlacement::Weighted;
  if (weights.has_value() && !weighted) invalid("weights require placement=WEIGHTED");
  if (weights == WeightBasis::Explicit && !edge_weights) invalid("weights=EXPLICIT requires edge_weights");
  if (edge_weights && weights != WeightBasis::Explicit) invalid("edge_weights require weights=EXPLICIT");
  if (weights == WeightBasis::ResidualSnapshot && !residual) invalid("weights=RESIDUAL_SNAPSHOT requires residual");
  if (residual && weights != WeightBasis::ResidualSnapshot) invalid("residual requires weights=RESIDUAL_SNAPSHOT");
}

// ---------------------------------------------------------------- the models' programs

std::shared_ptr<const ForwardingProgram> controller_program(
    const Graph& g, std::span<const NodeId> sources, std::span<const NodeId> destinations,
    std::span<const NodeId> roster, std::vector<ControllerPath> paths, std::span<const Amount> supplies,
    ControllerBasis basis) {
  ProgramData data;
  ProgramData* p = &data;
  p->graph_uid = g.uid();
  if (!supplies.empty()) {
    p->head_weight.assign(supplies.begin(), supplies.end());
    p->source_weights.emplace();
    p->source_weights->reserve(roster.size());
    for (std::size_t i = 0; i < roster.size(); ++i) p->source_weights->emplace_back(roster[i], supplies[i]);
  }
  const auto cost = g.cost();
  const auto edge_dst = g.edge_dst();
  p->tunnels.reserve(paths.size());
  p->entries.reserve(paths.size());
  for (std::size_t i = 0; i < paths.size(); ++i) {
    auto& path = paths[i];
    Tunnel t;
    t.index = static_cast<std::int32_t>(i);
    t.head = roster[static_cast<std::size_t>(path.slot)];
    t.destination = edge_dst[static_cast<std::size_t>(path.edges.back())];
    t.route_cost = 0;
    for (EdgeId e : path.edges) t.route_cost += cost[static_cast<std::size_t>(e)];
    t.weight = path.amount;
    ProgramEntry entry;
    entry.index = t.index;
    entry.head = t.head;
    entry.slot = path.slot;
    entry.weight = path.amount;
    entry.tunnel = static_cast<std::int32_t>(i);
    entry.cost = t.route_cost;
    t.edges = std::move(path.edges);
    p->tunnels.push_back(std::move(t));
    p->entries.push_back(std::move(entry));
  }
  p->topology = basis.topology;
  p->slot_routes = std::move(basis.per_source);
  p->common_basis = std::move(basis.common);
  p->lay_out(sources, destinations, roster, basis.excluded_nodes, basis.excluded_edges);
  return ProgramAccess::make(std::move(data));
}

std::shared_ptr<const ForwardingProgram> cost_only_program(
    const Graph& g, std::span<const NodeId> sources, std::span<const NodeId> destinations,
    std::span<const NodeId> roster, std::shared_ptr<const KernelPlan> plan, std::vector<ProgramEntry> entries) {
  ProgramData data;
  ProgramData* p = &data;
  p->graph_uid = g.uid();
  SplitContext ctx;
  ctx.plan = std::move(plan);
  p->entries = std::move(entries);
  for (std::size_t i = 0; i < roster.size(); ++i) {
    auto& entry = p->entries[i];
    entry.index = static_cast<std::int32_t>(i);
    entry.head = roster[i];
    entry.slot = static_cast<std::int32_t>(i);
    entry.weight = 1;
    entry.context = 0;
  }
  p->contexts.push_back(std::move(ctx));
  p->pooled = true;
  p->lay_out(sources, destinations, roster);
  return ProgramAccess::make(std::move(data));
}

std::shared_ptr<const ForwardingProgram> empty_program(std::uint64_t graph_uid, std::span<const NodeId> sources,
                                                       std::span<const NodeId> destinations, bool captured) {
  ProgramData data;
  ProgramData* p = &data;
  p->graph_uid = graph_uid;
  if (captured) p->source_weights.emplace();
  p->lay_out(sources, destinations, std::span<const NodeId> {});
  return ProgramAccess::make(std::move(data));
}

// ---------------------------------------------------------------- the change flags

namespace {

bool same_plan_structure(const KernelPlan& a, const KernelPlan& b) noexcept {
  return a.placement == b.placement && a.members == b.members && std::equal(a.roots.begin(), a.roots.end(), b.roots.begin(), b.roots.end());
}

} // namespace

void ProgramData::basis_parts(std::size_t slot, std::vector<std::pair<std::span<const EdgeId>, bool>>& out) const {
  if (common_basis) {
    out.emplace_back(*common_basis, true);
    return;
  }
  if (!slot_routes.empty()) {
    out.emplace_back(PredDAGBuilder::ascending(*slot_routes[slot]), true);
    return;
  }
  for (std::int32_t e : entries_of(slot)) {
    const auto& entry = entries[static_cast<std::size_t>(e)];
    if (entry.tunnel >= 0) {
      out.emplace_back(tunnels[static_cast<std::size_t>(entry.tunnel)].edges, false);
    } else {
      out.emplace_back(PredDAGBuilder::ascending(*entry.bundle), true);
    }
  }
}

ProgramChange compare_programs(const ForwardingProgram* before, const ForwardingProgram* after) noexcept {
  static const ProgramData kNone {};
  const ProgramData& a = before != nullptr ? ProgramAccess::data(*before) : kNone;
  const ProgramData& b = after != nullptr ? ProgramAccess::data(*after) : kNone;
  if (a.entries.empty() && b.entries.empty()) return {};
  ProgramChange out;
  // Routes: the ordered heads, each entry's kind and route (a tunnel's path, a split route's
  // context), duplicates included (entries are compared one by one).
  out.routes = !std::ranges::equal(a.roster, b.roster) || a.entries.size() != b.entries.size() ||
               a.contexts.size() != b.contexts.size();
  for (std::size_t i = 0; !out.routes && i < a.entries.size(); ++i) {
    const auto& x = a.entries[i];
    const auto& y = b.entries[i];
    if (x.index != y.index || x.head != y.head || (x.tunnel < 0) != (y.tunnel < 0)) {
      out.routes = true;
    } else if (x.tunnel >= 0) {
      out.routes = a.tunnels[static_cast<std::size_t>(x.tunnel)].edges != b.tunnels[static_cast<std::size_t>(y.tunnel)].edges;
    } else {
      out.routes = x.context != y.context;
    }
  }
  for (std::size_t c = 0; !out.routes && c < a.contexts.size(); ++c) {
    out.routes = !same_plan_structure(*a.contexts[c].plan, *b.contexts[c].plan);
  }
  // Weights: the injection (equal or the captured weights) and every entry and split-table weight.
  out.weights = a.source_weights.has_value() != b.source_weights.has_value() || a.head_weight != b.head_weight ||
                a.entries.size() != b.entries.size() || a.contexts.size() != b.contexts.size();
  for (std::size_t i = 0; !out.weights && i < a.entries.size(); ++i) out.weights = a.entries[i].weight != b.entries[i].weight;
  for (std::size_t c = 0; !out.weights && c < a.contexts.size(); ++c) {
    out.weights = a.contexts[c].plan->weights != b.contexts[c].plan->weights;
  }
  return out;
}

} // namespace detail

// ---------------------------------------------------------------- from_routes

std::shared_ptr<const ForwardingProgram> ForwardingProgram::from_routes(const Graph& g, std::span<const NodeId> sources,
                                                                        std::span<const NodeId> destinations,
                                                                        std::span<const Route> routes,
                                                                        const RouteOptions& options) {
  using detail::invalid;
  using detail::route_at;
  // 1. The endpoint sets.
  auto src = detail::node_set(sources, g.num_nodes(), "sources");
  auto dst = detail::node_set(destinations, g.num_nodes(), "destinations");
  detail::check_disjoint(src, dst, "sources", "destinations");
  // 2. The list.
  if (routes.empty()) invalid("routes must not be empty");
  // 3 to 5. The options.
  bool any_dag = false;
  for (const auto& r : routes) any_dag = any_dag || std::holds_alternative<std::shared_ptr<const PredDAG>>(r.bundle);
  detail::check_route_options(any_dag, options.placement, options.weights, options.edge_weights.has_value(),
                              options.residual.has_value());
  // 6. The vectors' lengths and entries.
  const auto E = static_cast<std::size_t>(g.num_edges());
  if (options.edge_weights.has_value()) {
    const auto& w = *options.edge_weights;
    if (w.size() != E) invalid("edge_weights: expected " + std::to_string(E) + " entries");
    for (std::size_t e = 0; e < E; ++e) {
      if (w[e] < 0 || w[e] >= kDomainLimit) invalid("edge_weights[" + std::to_string(e) + "] must be an integer in [0, 2^62)");
    }
  }
  if (options.residual.has_value()) detail::check_edge_amounts(g, *options.residual, "residual");
  // 7. Each route in list order, its entry made as soon as it is checked (a later route's failure
  // discards the program being made). The roster is the sources with a route.
  const auto edge_dst = g.edge_dst();
  const FlowPlacement placement = options.placement.value_or(FlowPlacement::EqualBalanced);
  const WeightBasis basis = options.weights.value_or(WeightBasis::Capacity);
  // Weighted reads its weights per edge: the explicit weights or the residual snapshot.
  std::span<const Amount> edge_values;
  if (placement == FlowPlacement::Weighted && basis == WeightBasis::Explicit) edge_values = *options.edge_weights;
  if (placement == FlowPlacement::Weighted && basis == WeightBasis::ResidualSnapshot) edge_values = *options.residual;
  detail::ProgramData data;
  detail::ProgramData* p = &data;
  p->graph_uid = g.uid();
  std::size_t num_dags = 0;
  for (const auto& r : routes) num_dags += std::holds_alternative<std::shared_ptr<const PredDAG>>(r.bundle) ? 1 : 0;
  p->entries.reserve(routes.size());
  p->contexts.reserve(num_dags);
  p->tunnels.reserve(routes.size() - num_dags);
  // The destination flags, made for the first Path (a tunnel's end and its earlier hops).
  std::vector<std::uint8_t> is_dst;
  for (std::size_t i = 0; i < routes.size(); ++i) {
    const auto& r = routes[i];
    const auto* path = std::get_if<std::shared_ptr<const Path>>(&r.bundle);
    const auto* dag = std::get_if<std::shared_ptr<const PredDAG>>(&r.bundle);
    if ((path != nullptr && !*path) || (dag != nullptr && !*dag)) invalid(route_at(i) + "bundle is null");
    const std::uint64_t uid = path != nullptr ? (*path)->graph_uid() : (*dag)->graph_uid();
    if (uid != g.uid()) invalid(route_at(i) + "bundle belongs to another graph");
    if (!std::binary_search(src.begin(), src.end(), r.head)) {
      invalid(route_at(i) + "head " + std::to_string(r.head) + " is not a source");
    }
    detail::ProgramEntry entry;
    entry.index = static_cast<std::int32_t>(i);
    entry.head = r.head;
    entry.weight = r.weight;
    if (path != nullptr) {
      if (is_dst.empty()) {
        is_dst.assign(static_cast<std::size_t>(g.num_nodes()), 0);
        for (NodeId d : dst) is_dst[static_cast<std::size_t>(d)] = 1;
      }
      const auto& hops = (*path)->hops();
      for (const auto& hop : hops) {
        if (hop.size() != 1) invalid(route_at(i) + "a tunnel has one edge per hop");
      }
      if ((*path)->src() != r.head) {
        invalid(route_at(i) + "the path starts at " + std::to_string((*path)->src()) + ", not at head " + std::to_string(r.head));
      }
      if (is_dst[static_cast<std::size_t>((*path)->dst())] == 0) {
        invalid(route_at(i) + "the path ends at " + std::to_string((*path)->dst()) + ", not at a destination");
      }
      for (std::size_t h = 0; h + 1 < hops.size(); ++h) {
        const NodeId v = edge_dst[static_cast<std::size_t>(hops[h].front())];
        if (is_dst[static_cast<std::size_t>(v)] != 0) {
          invalid(route_at(i) + "the path meets destination " + std::to_string(v) + " before its end");
        }
      }
      if (r.weight < 0 || r.weight >= kDomainLimit) invalid(route_at(i) + "weight must be an integer in [0, 2^62)");
      Tunnel t;
      t.index = entry.index;
      t.head = r.head;
      t.destination = (*path)->dst();
      t.route_cost = (*path)->cost();
      t.weight = r.weight;
      t.edges.reserve(hops.size());
      for (const auto& hop : hops) t.edges.push_back(hop.front());
      entry.tunnel = static_cast<std::int32_t>(p->tunnels.size());
      entry.cost = t.route_cost;
      p->tunnels.push_back(std::move(t));
    } else {
      // The plan in one pass: its route (the DAG's edges on the walks from the head to a
      // destination, before any weight pruning: the walk check) and the route's cost.
      Cost cost = kUnreachable;
      auto plan = detail::make_route_plan(g, placement, basis, edge_values, r.head, dst, **dag, &cost);
      if (plan.route.empty()) invalid(route_at(i) + "no walk from head " + std::to_string(r.head) + " to a destination");
      if (r.weight < 0 || r.weight >= kDomainLimit) invalid(route_at(i) + "weight must be an integer in [0, 2^62)");
      // Its own split context; the record bundle is the supplied DAG itself when nothing was
      // pruned (shared, never copied).
      const auto route_dag = plan.route.size() == detail::PredDAGBuilder::ascending(**dag).size()
                                 ? *dag
                                 : std::make_shared<const PredDAG>(detail::PredDAGBuilder::of_edges(g, plan.route));
      detail::SplitContext ctx;
      ctx.plan = std::make_shared<const detail::KernelPlan>(std::move(plan));
      entry.context = static_cast<std::int32_t>(p->contexts.size());
      entry.bundle = route_dag;
      entry.cost = cost;
      p->contexts.push_back(std::move(ctx));
    }
    p->entries.push_back(std::move(entry));
  }
  // One split route of positive weight forwards in its pooled form (its plan's kernels). The
  // roster: the heads with a route.
  p->pooled = p->entries.size() == 1 && p->tunnels.empty() && p->entries.front().weight > 0;
  p->lay_out(src, dst, std::nullopt);
  return detail::ProgramAccess::make(std::move(data));
}

} // namespace netgraph::core
