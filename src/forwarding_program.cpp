/*
  ForwardingProgram: the value, from_routes, the programs of the computing models and the
  comparison of the change flags.

  Every program keeps its entries in ascending index, per roster slot its entries, its split
  contexts (one plan each; a ShortestPaths program's is pooled over its roster) and its replay basis.
  Nothing in a program changes after construction; policies share it by reference.
*/
#include "netgraph/core/forwarding_program.hpp"

#include <algorithm>
#include <stdexcept>
#include <string>
#include <utility>

#include "netgraph/core/detail/edge_amounts.hpp"
#include "netgraph/core/detail/endpoints.hpp"
#include "netgraph/core/detail/enum_domain.hpp"
#include "netgraph/core/detail/program_data.hpp"
#include "netgraph/core/detail/route_builders.hpp"

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

void ProgramData::initialize_storage(std::span<const NodeId> source, std::span<const NodeId> destination,
                          std::optional<std::span<const NodeId>> given, std::span<const NodeId> xn,
                          std::span<const EdgeId> xe) {
  const std::size_t n = entries.size();
  const std::size_t room = given.has_value() ? given->size() : n;  // the roster's room
  const std::size_t lists = pools_source_inflows ? 0 : room + 2 + n;              // source_entry_offsets (counting sort) + source_entry_positions
  std::size_t roots = 0;
  for (const auto& c : split_contexts) roots += c.prepared_route->sources.size();
  // The block, in this order: sources, destinations, roster, excluded nodes and edges, the per-slot
  // lists, the per-tunnel entries, the contexts' root entries.
  std::vector<std::int32_t> store(source.size() + destination.size() + room + xn.size() + xe.size() + lists + tunnels.size() +
                                  roots);
  std::int32_t* at = store.data();
  auto put = [&at](auto list) {
    std::int32_t* first = at;
    at = std::copy(list.begin(), list.end(), at);
    return first;
  };
  NodeId* s0 = put(source);
  NodeId* d0 = put(destination);
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
  std::int32_t* order = first + (pools_source_inflows ? 0 : room + 2);
  std::int32_t* by_tunnel = first + lists;
  const std::span<const NodeId> ros(r0, k);
  for (auto& entry : entries) {
    entry.roster_index = static_cast<std::int32_t>(std::lower_bound(ros.begin(), ros.end(), entry.head) - ros.begin());
  }
  if (pools_source_inflows) {
    // A pooled program's slots carry one entry each (entries_for_source): no lists, local index 0.
    for (auto& entry : entries) entry.local_index = 0;
  } else {
    // A stable counting sort by slot: counted at slot + 2, summed, placed through slot + 1 (which
    // then holds the slot's end, the next slot's start); ascending index within a slot.
    for (const auto& entry : entries) ++first[static_cast<std::size_t>(entry.roster_index) + 2];
    for (std::size_t x = 2; x < k + 2; ++x) first[x] += first[x - 1];
    for (std::size_t e = 0; e < n; ++e) order[first[static_cast<std::size_t>(entries[e].roster_index) + 1]++] = static_cast<std::int32_t>(e);
    for (std::size_t x = 0; x < k; ++x) {
      for (auto j = first[x]; j < first[x + 1]; ++j) {
        auto& entry = entries[static_cast<std::size_t>(order[j])];
        entry.local_index = j - first[x];
        if (entry.tunnel_position >= 0) tunnels[static_cast<std::size_t>(entry.tunnel_position)].local_index = entry.local_index;
      }
    }
  }
  std::fill(by_tunnel, by_tunnel + tunnels.size() + roots, -1);
  for (std::size_t e = 0; e < n; ++e) {
    if (entries[e].tunnel_position >= 0) by_tunnel[entries[e].tunnel_position] = static_cast<std::int32_t>(e);
  }
  // Each context's root entries: the entry of that context whose head is the root (a ShortestPaths
  // plan's roots are its roster, each with its entry; a supplied split route's is its head).
  std::int32_t* root0 = by_tunnel + tunnels.size();
  for (std::size_t c = 0, at_root = 0; c < split_contexts.size(); ++c) {
    split_contexts[c].source_entry_positions = {root0 + at_root, split_contexts[c].prepared_route->sources.size()};
    at_root += split_contexts[c].prepared_route->sources.size();
  }
  for (std::size_t e = 0; e < n; ++e) {
    const auto c = entries[e].split_context_position;
    if (c < 0) continue;
    const auto& ctx = split_contexts[static_cast<std::size_t>(c)];
    const auto& plan_roots = ctx.prepared_route->sources;
    const auto r = std::lower_bound(plan_roots.begin(), plan_roots.end(), entries[e].head) - plan_roots.begin();
    root0[(ctx.source_entry_positions.data() - root0) + r] = static_cast<std::int32_t>(e);
  }
  // The views, into the block now owned (a vector's move keeps its elements where they are).
  store_ = std::move(store);
  sources = {s0, source.size()};
  destinations = {d0, destination.size()};
  reachable_sources = ros;
  excluded_nodes = {xn0, xn.size()};
  excluded_edges = {xe0, xe.size()};
  source_entry_offsets = pools_source_inflows ? std::span<const std::int32_t> {} : std::span<const std::int32_t>(first, k + 1);
  source_entry_positions = pools_source_inflows ? std::span<const std::int32_t> {} : std::span<const std::int32_t>(order, n);
  tunnel_entry_positions = {by_tunnel, tunnels.size()};
}

std::shared_ptr<const ForwardingProgram> ProgramAccess::make(ProgramData&& data) {
  auto block = std::make_shared<const ProgramBlock>(std::move(data));
  return {block, &block->program};
}

void check_route_options(const char* function, bool any_dag, std::optional<RouteAllocation> allocation,
                         std::optional<WeightBasis> weight_basis, bool edge_weights, bool residual) {
  if (!any_dag) {
    if (allocation.has_value()) invalid("allocation applies only to PredecessorDAG routes");
    if (weight_basis.has_value()) invalid("weight_basis applies only to PredecessorDAG routes");
    if (edge_weights) invalid("edge_weights applies only to PredecessorDAG routes");
    if (residual) invalid("residual applies only to PredecessorDAG routes");
    return;
  }
  if (allocation.has_value()) check_member(*allocation, "allocation");
  if (weight_basis.has_value()) check_member(*weight_basis, "weight_basis");
  if (allocation == RouteAllocation::Optimized) invalid(std::string(function) + " allocation must be EQUAL or WEIGHTED");
  const bool weighted = allocation == RouteAllocation::Weighted;
  if (weight_basis.has_value() && !weighted) invalid("weight_basis requires allocation=WEIGHTED");
  if (weight_basis == WeightBasis::Explicit && !edge_weights) invalid("weight_basis=EXPLICIT requires edge_weights");
  if (edge_weights && weight_basis != WeightBasis::Explicit) invalid("edge_weights require weight_basis=EXPLICIT");
  if (weight_basis == WeightBasis::ResidualSnapshot && !residual) invalid("weight_basis=RESIDUAL_SNAPSHOT requires residual");
  if (residual && weight_basis != WeightBasis::ResidualSnapshot) invalid("residual requires weight_basis=RESIDUAL_SNAPSHOT");
}

// ---------------------------------------------------------------- the models' programs

std::shared_ptr<const ForwardingProgram> optimized_flow_program(
    const Graph& graph, std::span<const NodeId> sources, std::span<const NodeId> destinations,
    std::span<const NodeId> reachable_sources, std::vector<OptimizedFlowPath> paths, std::span<const Amount> supplies,
    ReplayBasis basis) {
  ProgramData data;
  ProgramData* p = &data;
  p->graph_uid = graph.uid();
  if (!supplies.empty()) {
    p->source_weight_values.assign(supplies.begin(), supplies.end());
    p->source_weights.emplace();
    p->source_weights->reserve(reachable_sources.size());
    for (std::size_t i = 0; i < reachable_sources.size(); ++i) p->source_weights->emplace_back(reachable_sources[i], supplies[i]);
  }
  const auto cost = graph.edge_costs();
  const auto edge_destinations = graph.edge_destinations();
  p->tunnels.reserve(paths.size());
  p->entries.reserve(paths.size());
  for (std::size_t i = 0; i < paths.size(); ++i) {
    auto& path = paths[i];
    Tunnel t;
    t.entry_id = static_cast<std::int32_t>(i);
    t.head = reachable_sources[static_cast<std::size_t>(path.roster_index)];
    t.destination = edge_destinations[static_cast<std::size_t>(path.edges.back())];
    Cost route_cost = 0;
    for (EdgeId e : path.edges) route_cost += cost[static_cast<std::size_t>(e)];
    t.path_cost = route_cost;
    t.weight = path.weight;
    ProgramEntry entry;
    entry.entry_id = t.entry_id;
    entry.head = t.head;
    entry.roster_index = path.roster_index;
    entry.weight = path.weight;
    entry.tunnel_position = static_cast<std::int32_t>(i);
    entry.min_route_cost = route_cost;
    t.edges = std::move(path.edges);
    p->tunnels.push_back(std::move(t));
    p->entries.push_back(std::move(entry));
  }
  p->uses_topology_basis = basis.uses_topology_basis;
  p->source_supports = std::move(basis.source_supports);
  p->common_replay_support = std::move(basis.common_support);
  p->initialize_storage(sources, destinations, reachable_sources, basis.excluded_nodes, basis.excluded_edges);
  return ProgramAccess::make(std::move(data));
}

std::shared_ptr<const ForwardingProgram> shortest_path_program(
    const Graph& graph, std::span<const NodeId> sources, std::span<const NodeId> destinations,
    std::span<const NodeId> reachable_sources, std::shared_ptr<const PreparedRoute> plan, std::vector<ProgramEntry> entries) {
  ProgramData data;
  ProgramData* p = &data;
  p->graph_uid = graph.uid();
  SplitContext ctx;
  ctx.prepared_route = std::move(plan);
  p->entries = std::move(entries);
  for (std::size_t i = 0; i < reachable_sources.size(); ++i) {
    auto& entry = p->entries[i];
    entry.entry_id = static_cast<std::int32_t>(i);
    entry.head = reachable_sources[i];
    entry.roster_index = static_cast<std::int32_t>(i);
    entry.weight = 1;
    entry.split_context_position = 0;
  }
  p->split_contexts.push_back(std::move(ctx));
  p->pools_source_inflows = true;
  p->initialize_storage(sources, destinations, reachable_sources);
  return ProgramAccess::make(std::move(data));
}

std::shared_ptr<const ForwardingProgram> pooled_dag_program(
    const Graph& graph, std::span<const NodeId> sources, std::span<const NodeId> destinations,
    std::span<const NodeId> roster, const PredecessorDAG& dag, RouteAllocation allocation,
    const WeightSource& weights, std::vector<ProgramEntry> entries) {
  // Equal and Weighted take either admission: the plan is the same for both.
  auto plan = std::make_shared<const PreparedRoute>(
      prepare_route(graph, allocation, Admission::Lossless, weights, roster, destinations, dag));
  return shortest_path_program(graph, sources, destinations, roster, std::move(plan), std::move(entries));
}

std::shared_ptr<const ForwardingProgram> distributed_te_program(
    const Graph& g, std::span<const NodeId> sources, std::span<const NodeId> destinations,
    std::span<const NodeId> roster, std::vector<TunnelSetupResult> tunnels, std::span<const NodeId> excluded_nodes,
    std::span<const EdgeId> excluded_edges) {
  ProgramData data;
  ProgramData* p = &data;
  p->graph_uid = g.uid();
  const auto edge_destinations = g.edge_destinations();
  const auto H = static_cast<std::int32_t>(roster.size());
  p->tunnels.reserve(tunnels.size());
  p->entries.reserve(tunnels.size());
  for (std::size_t i = 0; i < tunnels.size(); ++i) {
    auto& s = tunnels[i];
    const bool up = s.status == TunnelStatus::Up;
    Tunnel t;
    t.entry_id = s.local_index * H + s.roster_index;
    t.head = roster[static_cast<std::size_t>(s.roster_index)];
    t.status = s.status;
    t.requested_bandwidth = s.requested_bandwidth;
    t.reserved_bandwidth = up ? s.requested_bandwidth : 0;
    t.planned_bandwidth = s.planned;
    t.weight = up ? 1 : 0;
    ProgramEntry entry;
    entry.entry_id = t.entry_id;
    entry.head = t.head;
    entry.roster_index = s.roster_index;
    entry.weight = t.weight;
    entry.tunnel_position = static_cast<std::int32_t>(i);
    entry.min_route_cost = up ? s.path_cost : kUnreachable;
    if (up) {
      t.destination = edge_destinations[static_cast<std::size_t>(s.edges.back())];
      t.path_cost = s.path_cost;
      t.edges = std::move(s.edges);
    }
    p->tunnels.push_back(std::move(t));
    p->entries.push_back(std::move(entry));
  }
  p->uses_topology_basis = true;
  p->initialize_storage(sources, destinations, roster, excluded_nodes, excluded_edges);
  return ProgramAccess::make(std::move(data));
}

std::shared_ptr<const ForwardingProgram> empty_program(std::uint64_t graph_uid, std::span<const NodeId> sources,
                                                       std::span<const NodeId> destinations, bool captured) {
  ProgramData data;
  ProgramData* p = &data;
  p->graph_uid = graph_uid;
  if (captured) p->source_weights.emplace();
  p->initialize_storage(sources, destinations, std::span<const NodeId> {});
  return ProgramAccess::make(std::move(data));
}

// ---------------------------------------------------------------- the change flags

namespace {

bool same_plan_structure(const PreparedRoute& a, const PreparedRoute& b) noexcept {
  return a.allocation == b.allocation && a.forwarding_edges == b.forwarding_edges && std::equal(a.sources.begin(), a.sources.end(), b.sources.begin(), b.sources.end());
}

} // namespace

void ProgramData::basis_parts(std::size_t slot, std::vector<std::pair<std::span<const EdgeId>, bool>>& out) const {
  if (common_replay_support) {
    out.emplace_back(*common_replay_support, true);
    return;
  }
  if (!source_supports.empty()) {
    out.emplace_back(PredecessorDAGBuilder::ascending(*source_supports[slot]), true);
    return;
  }
  for (std::int32_t e : entries_for_source(slot)) {
    const auto& entry = entries[static_cast<std::size_t>(e)];
    if (entry.tunnel_position >= 0) {
      out.emplace_back(path_of(tunnels[static_cast<std::size_t>(entry.tunnel_position)]), false);
    } else {
      out.emplace_back(PredecessorDAGBuilder::ascending(*entry.route), true);
    }
  }
}

ProgramChange compare_programs(const ForwardingProgram* before, const ForwardingProgram* after,
                               bool reservations) noexcept {
  static const ProgramData kNone {};
  const ProgramData& a = before != nullptr ? ProgramAccess::data(*before) : kNone;
  const ProgramData& b = after != nullptr ? ProgramAccess::data(*after) : kNone;
  if (a.entries.empty() && b.entries.empty()) return {};
  ProgramChange out;
  // Routes: the ordered heads, each entry's kind and route (a tunnel's path, a split route's
  // context), duplicates included (entries are compared one by one).
  out.routes = !std::ranges::equal(a.reachable_sources, b.reachable_sources) || a.entries.size() != b.entries.size() ||
               a.split_contexts.size() != b.split_contexts.size();
  for (std::size_t i = 0; !out.routes && i < a.entries.size(); ++i) {
    const auto& x = a.entries[i];
    const auto& y = b.entries[i];
    if (x.entry_id != y.entry_id || x.head != y.head || (x.tunnel_position < 0) != (y.tunnel_position < 0)) {
      out.routes = true;
    } else if (x.tunnel_position >= 0) {
      out.routes = a.tunnels[static_cast<std::size_t>(x.tunnel_position)].edges != b.tunnels[static_cast<std::size_t>(y.tunnel_position)].edges;
    } else {
      out.routes = x.split_context_position != y.split_context_position;
    }
  }
  for (std::size_t c = 0; !out.routes && c < a.split_contexts.size(); ++c) {
    out.routes = !same_plan_structure(*a.split_contexts[c].prepared_route, *b.split_contexts[c].prepared_route);
  }
  // Weights: the injection (equal or the captured weights) and every entry and split-table weight,
  // an entry absent from one program counting as weight 0 there.
  out.weights = a.source_weights.has_value() != b.source_weights.has_value() || a.source_weight_values != b.source_weight_values ||
                a.split_contexts.size() != b.split_contexts.size();
  const std::size_t n = std::max(a.entries.size(), b.entries.size());
  auto weight = [](const ProgramData& p, std::size_t i) { return i < p.entries.size() ? p.entries[i].weight : Amount {0}; };
  for (std::size_t i = 0; !out.weights && i < n; ++i) out.weights = weight(a, i) != weight(b, i);
  for (std::size_t c = 0; !out.weights && c < a.split_contexts.size(); ++c) {
    out.weights = a.split_contexts[c].prepared_route->weights != b.split_contexts[c].prepared_route->weights;
  }
  // Reservations: every tunnel's (requested, reserved) by index (both lists ascend), a tunnel
  // absent from one program counting as (0, 0) there (only T1 tunnels have them).
  if (!reservations) return out;
  auto held = [](const Tunnel& t) { return t.requested_bandwidth != 0 || t.reserved_bandwidth != 0; };
  std::size_t i = 0;
  std::size_t j = 0;
  while (!out.reservations && (i < a.tunnels.size() || j < b.tunnels.size())) {
    const Tunnel* x = i < a.tunnels.size() ? &a.tunnels[i] : nullptr;
    const Tunnel* y = j < b.tunnels.size() ? &b.tunnels[j] : nullptr;
    if (x != nullptr && (y == nullptr || x->entry_id < y->entry_id)) {
      out.reservations = held(*x);
      ++i;
    } else if (y != nullptr && (x == nullptr || y->entry_id < x->entry_id)) {
      out.reservations = held(*y);
      ++j;
    } else {
      out.reservations = x->requested_bandwidth != y->requested_bandwidth || x->reserved_bandwidth != y->reserved_bandwidth;
      ++i;
      ++j;
    }
  }
  return out;
}

} // namespace detail

// ---------------------------------------------------------------- from_routes

namespace {

// The last option step of from_routes and from_dag: edge_weights' length and entries, then residual's (the ledger's baseline messages).
void check_route_vectors(const Graph& graph, const RouteOptions& options) {
  const auto E = static_cast<std::size_t>(graph.num_edges());
  if (options.edge_weights.has_value()) {
    const auto& w = *options.edge_weights;
    if (w.size() != E) detail::invalid("edge_weights: expected " + std::to_string(E) + " entries");
    for (std::size_t e = 0; e < E; ++e) {
      if (w[e] < 0 || w[e] >= kDomainLimit) detail::invalid("edge_weights[" + std::to_string(e) + "] must be an integer in [0, 2^62)");
    }
  }
  if (options.residual.has_value()) detail::check_edge_amounts(graph, *options.residual, "residual");
}

} // namespace

std::shared_ptr<const ForwardingProgram> ForwardingProgram::from_routes(const Graph& graph, std::span<const NodeId> sources,
                                                                        std::span<const NodeId> destinations,
                                                                        std::span<const RouteDefinition> routes,
                                                                        const RouteOptions& options) {
  using detail::invalid;
  using detail::route_at;
  // 1. The endpoint sets.
  auto src = detail::node_set(sources, graph.num_nodes(), "sources");
  auto dst = detail::node_set(destinations, graph.num_nodes(), "destinations");
  detail::check_disjoint(src, dst, "sources", "destinations");
  // 2. The list.
  if (routes.empty()) invalid("routes must not be empty");
  // 3 to 5. The options.
  bool any_dag = false;
  for (const auto& r : routes) any_dag = any_dag || std::holds_alternative<std::shared_ptr<const PredecessorDAG>>(r.route);
  detail::check_route_options("from_routes", any_dag, options.allocation, options.weight_basis,
                              options.edge_weights.has_value(), options.residual.has_value());
  // 6. The vectors' lengths and entries.
  check_route_vectors(graph, options);
  // 7. Each route in list order, its entry made as soon as it is checked (a later route's failure
  // discards the program being made). The roster is the sources with a route.
  const auto edge_destinations = graph.edge_destinations();
  const RouteAllocation allocation = options.allocation.value_or(RouteAllocation::Equal);
  const WeightBasis basis = options.weight_basis.value_or(WeightBasis::Capacity);
  // Weighted reads its weights per edge: the explicit weights or the residual snapshot.
  std::span<const Amount> edge_values;
  if (allocation == RouteAllocation::Weighted && basis == WeightBasis::Explicit) edge_values = *options.edge_weights;
  if (allocation == RouteAllocation::Weighted && basis == WeightBasis::ResidualSnapshot) edge_values = *options.residual;
  detail::ProgramData data;
  detail::ProgramData* p = &data;
  p->graph_uid = graph.uid();
  std::size_t num_dags = 0;
  for (const auto& r : routes) num_dags += std::holds_alternative<std::shared_ptr<const PredecessorDAG>>(r.route) ? 1 : 0;
  p->entries.reserve(routes.size());
  p->split_contexts.reserve(num_dags);
  p->tunnels.reserve(routes.size() - num_dags);
  // The destination flags, made for the first Path (a tunnel's end and its earlier hops).
  std::vector<std::uint8_t> is_dst;
  for (std::size_t i = 0; i < routes.size(); ++i) {
    const auto& r = routes[i];
    const auto* path = std::get_if<std::shared_ptr<const Path>>(&r.route);
    const auto* dag = std::get_if<std::shared_ptr<const PredecessorDAG>>(&r.route);
    if ((path != nullptr && !*path) || (dag != nullptr && !*dag)) invalid(route_at(i) + "route is null");
    const std::uint64_t uid = path != nullptr ? (*path)->graph_uid() : (*dag)->graph_uid();
    if (uid != graph.uid()) invalid(route_at(i) + "route belongs to another graph");
    if (!std::binary_search(src.begin(), src.end(), r.head)) {
      invalid(route_at(i) + "head " + std::to_string(r.head) + " is not a source");
    }
    detail::ProgramEntry entry;
    entry.entry_id = static_cast<std::int32_t>(i);
    entry.head = r.head;
    entry.weight = r.weight;
    if (path != nullptr) {
      if (is_dst.empty()) {
        is_dst.assign(static_cast<std::size_t>(graph.num_nodes()), 0);
        for (NodeId d : dst) is_dst[static_cast<std::size_t>(d)] = 1;
      }
      const auto& hops = (*path)->hops();
      for (const auto& hop : hops) {
        if (hop.size() != 1) invalid(route_at(i) + "a tunnel has one edge per hop");
      }
      if ((*path)->source() != r.head) {
        invalid(route_at(i) + "the path starts at " + std::to_string((*path)->source()) + ", not at head " + std::to_string(r.head));
      }
      if (is_dst[static_cast<std::size_t>((*path)->destination())] == 0) {
        invalid(route_at(i) + "the path ends at " + std::to_string((*path)->destination()) + ", not at a destination");
      }
      for (std::size_t h = 0; h + 1 < hops.size(); ++h) {
        const NodeId v = edge_destinations[static_cast<std::size_t>(hops[h].front())];
        if (is_dst[static_cast<std::size_t>(v)] != 0) {
          invalid(route_at(i) + "the path meets destination " + std::to_string(v) + " before its end");
        }
      }
      if (r.weight < 0 || r.weight >= kDomainLimit) invalid(route_at(i) + "weight must be an integer in [0, 2^62)");
      Tunnel t;
      t.entry_id = entry.entry_id;
      t.head = r.head;
      t.destination = (*path)->destination();
      t.path_cost = (*path)->cost();
      t.weight = r.weight;
      auto& edges = t.edges.emplace();
      edges.reserve(hops.size());
      for (const auto& hop : hops) edges.push_back(hop.front());
      entry.tunnel_position = static_cast<std::int32_t>(p->tunnels.size());
      entry.min_route_cost = (*path)->cost();
      p->tunnels.push_back(std::move(t));
    } else {
      // The plan in one pass: its route (the DAG's edges on the walks from the head to a
      // destination, before any weight pruning: the walk check) and the route's cost.
      Cost cost = kUnreachable;
      auto plan = detail::prepare_supplied_route(graph, allocation, basis, edge_values, r.head, dst, **dag, &cost);
      if (plan.support_edges.empty()) invalid(route_at(i) + "no walk from head " + std::to_string(r.head) + " to a destination");
      if (r.weight < 0 || r.weight >= kDomainLimit) invalid(route_at(i) + "weight must be an integer in [0, 2^62)");
      // Its own split context; the record route is the supplied DAG itself when nothing was
      // pruned (shared, never copied).
      const auto route_dag = plan.support_edges.size() == detail::PredecessorDAGBuilder::ascending(**dag).size()
                                 ? *dag
                                 : std::make_shared<const PredecessorDAG>(detail::PredecessorDAGBuilder::of_edges(graph, plan.support_edges));
      detail::SplitContext ctx;
      ctx.prepared_route = std::make_shared<const detail::PreparedRoute>(std::move(plan));
      entry.split_context_position = static_cast<std::int32_t>(p->split_contexts.size());
      entry.route = route_dag;
      entry.min_route_cost = cost;
      p->split_contexts.push_back(std::move(ctx));
    }
    p->entries.push_back(std::move(entry));
  }
  // One split route of positive weight forwards in its pooled form (its plan's kernels). The
  // roster: the heads with a route.
  p->pools_source_inflows = p->entries.size() == 1 && p->tunnels.empty() && p->entries.front().weight > 0;
  p->initialize_storage(src, dst, std::nullopt);
  return detail::ProgramAccess::make(std::move(data));
}

std::shared_ptr<const ForwardingProgram> ForwardingProgram::from_dag(const Graph& graph, std::span<const NodeId> sources,
                                                                     std::span<const NodeId> destinations,
                                                                     const std::shared_ptr<const PredecessorDAG>& dag,
                                                                     const RouteOptions& options) {
  auto src = detail::node_set(sources, graph.num_nodes(), "sources");
  auto dst = detail::node_set(destinations, graph.num_nodes(), "destinations");
  detail::check_disjoint(src, dst, "sources", "destinations");
  if (!dag || dag->graph_uid() != graph.uid()) detail::invalid("dag belongs to another graph");
  detail::check_route_options("from_dag", true, options.allocation, options.weight_basis, options.edge_weights.has_value(),
                              options.residual.has_value());
  check_route_vectors(graph, options);
  // The roster: the sources with a walk to a destination in the DAG (trimmed at the first
  // destinations), each with its own walks as its record route and replay basis.
  std::vector<NodeId> roster;
  std::vector<detail::ProgramEntry> entries;
  std::vector<EdgeId> all;
  for (NodeId s : src) {
    auto part = detail::prune(graph, *dag, s, dst, nullptr);
    if (part.empty()) continue;
    roster.push_back(s);
    detail::ProgramEntry entry;
    entry.route = std::make_shared<const PredecessorDAG>(detail::PredecessorDAGBuilder::of_edges(graph, part));
    entry.min_route_cost = detail::min_route_cost(graph, *entry.route, s, dst);
    entries.push_back(std::move(entry));
    all.insert(all.end(), part.begin(), part.end());
  }
  if (roster.empty()) return detail::empty_program(graph.uid(), src, dst, false);
  std::sort(all.begin(), all.end());
  all.erase(std::unique(all.begin(), all.end()), all.end());
  const auto supplied = detail::PredecessorDAGBuilder::ascending(*dag);
  const auto pooled = std::equal(all.begin(), all.end(), supplied.begin(), supplied.end())
                          ? dag
                          : std::make_shared<const PredecessorDAG>(detail::PredecessorDAGBuilder::of_edges(graph, all));
  // Weighted's weights: Explicit edge-aligned; the snapshot aligned with the pooled DAG's edges.
  const RouteAllocation allocation = options.allocation.value_or(RouteAllocation::Equal);
  detail::WeightSource weights {options.weight_basis.value_or(WeightBasis::Capacity), {}};
  std::vector<Amount> snapshot;
  if (allocation == RouteAllocation::Weighted && weights.basis == WeightBasis::Explicit) weights.values = *options.edge_weights;
  if (allocation == RouteAllocation::Weighted && weights.basis == WeightBasis::ResidualSnapshot) {
    snapshot.reserve(all.size());
    for (EdgeId e : all) snapshot.push_back((*options.residual)[static_cast<std::size_t>(e)]);
    weights.values = snapshot;
  }
  return detail::pooled_dag_program(graph, src, dst, roster, *pooled, allocation, weights, std::move(entries));
}

} // namespace netgraph::core
