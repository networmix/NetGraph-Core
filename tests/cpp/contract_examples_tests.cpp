// The worked examples run on the library: each example's own inputs and expected values,
// one test per example (ContractExamples.<id> with '-' as '_'). The examples state edges as
// (source, destination, capacity, cost) in input order and amounts in units; every per-edge value
// below is in that order (translated through the graph's input positions). An example whose input
// is a Python value (a float factor, a NumPy scalar, a wrong type) is checked here through the C++
// value it stands for; the conversion itself is the Python suite's.

#include <gtest/gtest.h>

#include <algorithm>
#include <cstdint>
#include <map>
#include <memory>
#include <optional>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include "netgraph/core/detail/forwarding.hpp"
#include "netgraph/core/detail/route_kernels.hpp"
#include "netgraph/core/evaluators.hpp"
#include "netgraph/core/flow_ledger.hpp"
#include "netgraph/core/flow_policy.hpp"
#include "netgraph/core/forwarding_program.hpp"
#include "netgraph/core/search.hpp"
#include "netgraph/core/shortest_path_cache.hpp"
#include "test_utils.hpp"

using namespace netgraph::core;
using namespace netgraph::core::test;

namespace {

using Amounts = std::vector<Amount>;
using Ints = std::vector<std::int64_t>;
using Routes = std::vector<Ints>;
using Texts = std::vector<std::string>;
using DesignArg = std::optional<std::optional<Amount>>;

constexpr Admission LL = Admission::Lossless;
constexpr Admission LY = Admission::Lossy;
constexpr Amount kBig = Amount {1} << 40;  // a link that never limits anything

DesignArg finite(Amount x) { return DesignArg {std::in_place, x}; }
const DesignArg kUnlimited {std::in_place, std::nullopt};

Amount lo(const Counter128& c) {
  EXPECT_EQ(c.hi, 0U);
  return static_cast<Amount>(c.lo);
}
Amount lo(const std::optional<Counter128>& c) { return c.has_value() ? lo(*c) : -1; }

std::string name_of(TunnelStatus s) {
  switch (s) {
    case TunnelStatus::Up: return "UP";
    case TunnelStatus::InsufficientBandwidth: return "INSUFFICIENT_BANDWIDTH";
    case TunnelStatus::CostCeilingExceeded: return "COST_CEILING_EXCEEDED";
    case TunnelStatus::Failed: return "FAILED";
    case TunnelStatus::Blackholed: return "BLACKHOLED";
  }
  return "?";
}
std::string name_of(Termination t) {
  switch (t) {
    case Termination::Complete: return "COMPLETE";
    case Termination::NoRoute: return "NO_ROUTE";
    case Termination::Exhausted: return "EXHAUSTED";
    case Termination::NoActiveWeight: return "NO_ACTIVE_WEIGHT";
  }
  return "?";
}

// A graph of the example's edges (input order) with its input <-> EdgeId maps.
struct Net {
  std::shared_ptr<const Graph> g;
  std::vector<EdgeId> id;  // input position -> EdgeId

  explicit Net(const std::vector<InEdge>& edges, std::int32_t n = 0) {
    for (const auto& e : edges) n = std::max({n, e.u + 1, e.v + 1});
    g = shared(make_graph(n, edges));
    id = id_of_input(*g);
  }
  [[nodiscard]] std::int64_t pos(EdgeId e) const { return g->input_positions()[static_cast<std::size_t>(e)]; }
  [[nodiscard]] std::size_t m() const { return id.size(); }
  [[nodiscard]] Ints positions(const std::vector<EdgeId>& edges) const {
    Ints out;
    for (EdgeId e : edges) out.push_back(pos(e));
    return out;
  }
  [[nodiscard]] std::vector<EdgeId> ids(const Ints& positions) const {
    std::vector<EdgeId> out;
    for (auto p : positions) out.push_back(id[static_cast<std::size_t>(p)]);
    return out;
  }
  // A per-EdgeId array in input order.
  [[nodiscard]] Amounts in_order(const std::vector<Amount>& by_id) const { return to_input(*g, by_id); }
  // An (EdgeId, amount) list as a per-edge array in input order.
  [[nodiscard]] Amounts per_edge(const std::vector<std::pair<EdgeId, Amount>>& pairs) const {
    Amounts out(m(), 0);
    for (const auto& [e, x] : pairs) out[static_cast<std::size_t>(pos(e))] += x;
    return out;
  }
  [[nodiscard]] Amounts per_edge(const std::vector<std::pair<EdgeId, Counter128>>& pairs) const {
    Amounts out(m(), 0);
    for (const auto& [e, x] : pairs) out[static_cast<std::size_t>(pos(e))] += lo(x);
    return out;
  }
  // A per-input array in EdgeId order.
  [[nodiscard]] Amounts by_id(const Amounts& by_input) const { return to_ids(*g, by_input); }
  [[nodiscard]] std::shared_ptr<const Eligibility> mask(const std::vector<NodeId>& nodes, const Ints& edges) const {
    if (nodes.empty() && edges.empty()) return nullptr;
    return std::make_shared<const Eligibility>(g->eligibility(nodes, ids(edges)));
  }
  [[nodiscard]] std::shared_ptr<const PredecessorDAG> dag(const Ints& edges) const { return dag_of(*g, ids(edges)); }
  [[nodiscard]] std::shared_ptr<const Path> path(const Ints& edges) const {
    std::vector<std::vector<EdgeId>> hops;
    for (EdgeId e : ids(edges)) hops.push_back({e});
    return std::make_shared<const Path>(*g, hops);
  }
};

struct Pol {
  std::shared_ptr<FlowLedger> ledger;
  std::shared_ptr<const Eligibility> mask;
  std::unique_ptr<FlowPolicy> p;
  FlowPolicy* operator->() const { return p.get(); }
};
Pol policy(const Net& n, const std::vector<NodeId>& s, const std::vector<NodeId>& d, FlowPolicyConfig cfg,
           std::shared_ptr<const Eligibility> mask = nullptr, std::shared_ptr<FlowLedger> ledger = nullptr) {
  Pol out {ledger ? std::move(ledger) : std::make_shared<FlowLedger>(n.g), std::move(mask), nullptr};
  out.p = std::make_unique<FlowPolicy>(*out.ledger, s, d, std::move(cfg), out.mask.get());
  return out;
}

const OutputOptions kFlows {.edge_flows = true};

// ---------------------------------------------------------------- shared integer rules

Amounts ha(Amount n, const Amounts& w) {
  Amounts out(w.size(), 0);
  detail::highest_averages(n, w, out);
  return out;
}

// A reduced fraction num / den (den > 0) and its comparison by cross products.
struct Frac {
  std::int64_t num;
  std::int64_t den;
};
Frac reduce(std::int64_t num, std::int64_t den) {
  const auto g = std::gcd(num < 0 ? -num : num, den);
  return {num / g, den / g};
}
bool less(const Frac& a, const Frac& b) { return static_cast<__int128>(a.num) * b.den < static_cast<__int128>(b.num) * a.den; }
std::string text(const Frac& f) { return std::to_string(f.num) + "/" + std::to_string(f.den); }

// Each split's largest excess and deficit against its exact share n w / W.
struct SplitError {
  Amounts assigned;
  Frac max_excess;
  Frac max_deficit;
  bool below_k;
  bool below_one;
};
SplitError split_error(Amount n, const Amounts& w) {
  const Amounts got = ha(n, w);
  Amount total = 0;
  for (Amount x : w) total += x;
  SplitError out {got, {-1, 1}, {-1, 1}, true, true};
  bool first = true;
  const auto k = static_cast<std::int64_t>(w.size());
  for (std::size_t i = 0; i < w.size(); ++i) {
    const Frac excess = reduce(got[i] * total - n * w[i], total);
    const Frac deficit = reduce(n * w[i] - got[i] * total, total);
    if (first || less(out.max_excess, excess)) out.max_excess = excess;
    if (first || less(out.max_deficit, deficit)) out.max_deficit = deficit;
    first = false;
    const std::int64_t gap = excess.num < 0 ? -excess.num : excess.num;
    out.below_k = out.below_k && gap < k * excess.den;
    out.below_one = out.below_one && gap < excess.den;
  }
  return out;
}

// ---------------------------------------------------------------- searches

struct SearchOut {
  std::optional<Cost> distance;
  Ints dag;
  std::vector<std::pair<Cost, Ints>> paths;
  std::string error;
};
SearchOut search(const std::vector<InEdge>& edges, NodeId source, NodeId destination, std::optional<Amounts> residual = std::nullopt,
                 std::optional<std::int64_t> k = std::nullopt, std::optional<CostStretch> stretch = std::nullopt) {
  const Net n(edges);
  SearchOut out;
  Amounts res;
  if (residual.has_value()) res = residual->size() == n.m() ? n.by_id(*residual) : *residual;
  if (k.has_value()) {
    KspOptions o(*k);
    if (stretch.has_value()) o.cost_ceiling = *stretch;
    if (residual.has_value()) o.residual = std::span<const Amount>(res);
    out.error = error_of([&] {
      for (const auto& [cost, p] : ksp(*n.g, source, destination, o)) {
        std::vector<EdgeId> es;
        for (const auto& hop : p.hops()) es.push_back(hop.front());
        out.paths.emplace_back(cost, n.positions(es));
      }
    });
    return out;
  }
  SpfOptions o(Direction::Forward);
  o.targets = std::vector<NodeId> {destination};
  if (residual.has_value()) o.residual = std::span<const Amount>(res);
  const NodeId roots[1] = {source};
  const auto r = spf(*n.g, roots, o);
  out.distance = r.target_cost;
  out.dag = n.positions(ascending_edges(r.dag));
  std::sort(out.dag.begin(), out.dag.end());
  return out;
}

// The single path spf chooses under a tie rule (multipath false), on `available` as the residual.
Ints tie_choice(const std::vector<InEdge>& edges, NodeId source, NodeId destination, RouteTieBreak tie,
                std::optional<Amounts> available = std::nullopt) {
  const Net n(edges);
  SpfOptions o(Direction::Forward);
  o.targets = std::vector<NodeId> {destination};
  o.tie_break = tie;
  o.multipath = false;
  Amounts res;
  if (available.has_value()) {
    res = n.by_id(*available);
    o.residual = std::span<const Amount>(res);
  }
  const NodeId roots[1] = {source};
  const auto r = spf(*n.g, roots, o);
  const auto via = ascending_edges(r.dag);
  // One path: order its edges from the source.
  Ints out;
  NodeId at = source;
  while (at != destination) {
    const auto it = std::find_if(via.begin(), via.end(), [&](EdgeId e) { return n.g->edge_sources()[static_cast<std::size_t>(e)] == at; });
    if (it == via.end()) break;
    out.push_back(n.pos(*it));
    at = n.g->edge_destinations()[static_cast<std::size_t>(*it)];
  }
  return out;
}

// ---------------------------------------------------------------- exact flows (CentralizedTE)

struct FlowOut {
  Amount delivered;
  Amount cost;
  Amounts edge_flows;
  std::optional<Ints> min_cut;
};
FlowOut exact_flow(const std::vector<InEdge>& edges, const std::vector<NodeId>& heads, const std::vector<NodeId>& dests,
                   std::optional<Amount> limit, SourceAllocation s) {
  const Net n(edges);
  auto p = policy(n, heads, dests, {.routing = centralized(s)});
  const auto r = p->place(limit, OutputOptions {.edge_flows = true, .min_cut = true});
  FlowOut out {r.call_delivered, lo(r.total_cost), n.in_order(*r.edge_flows), std::nullopt};
  if (r.min_cut.has_value()) {
    out.min_cut = n.positions(*r.min_cut);
    std::sort(out.min_cut->begin(), out.min_cut->end());
  }
  return out;
}

// ---------------------------------------------------------------- T1

struct T1In {
  std::vector<InEdge> edges;
  std::vector<NodeId> heads;
  std::vector<NodeId> destinations;
  Amount design_demand;
  std::int64_t per_head;
  std::optional<Amount> offer {};
  Ints excluded {};
  std::optional<Cost> absolute {};
  std::optional<CostStretch> stretch {};
};
struct T1Out {
  std::string error;
  Amounts requested;
  Routes routes;
  Texts status;
  std::int64_t up {0};
  Amount reserved {0};
  Amounts edge_reservations;
  Amount bound {0};
  Amounts offered_by_entry;
  Amount source_drops {0};
  Amounts no_drop_load;
  Amount lossy_delivery {0};
  Amounts lossy_edge_drops;
  Amount lossless_admitted {0};
  Amount lossless_unadmitted {0};
};
T1Out t1(const T1In& in) {
  const Net n(in.edges);
  DistributedTE dte {TunnelCount {in.per_head}};
  if (in.absolute.has_value()) dte.cost_ceiling = AbsoluteCost {*in.absolute};
  if (in.stretch.has_value()) dte.cost_ceiling = *in.stretch;
  const auto mask = n.mask({}, in.excluded);
  T1Out out;
  out.error = error_of([&] { (void)policy(n, in.heads, in.destinations, {.routing = dte, .admission = LL}, mask); });
  if (!out.error.empty()) return out;
  auto bound = policy(n, in.heads, in.destinations, {.routing = dte, .admission = LL}, mask);
  const auto b = bound->replace(std::nullopt, finite(in.design_demand), OutputOptions {.edge_reservations = true});
  out.bound = b.call_delivered;
  out.reserved = b.plan->total_reserved_bandwidth;
  out.edge_reservations = n.in_order(*b.edge_reservations);
  Routes paths;
  for (const auto& t : b.program->tunnels()) {
    out.requested.push_back(t.requested_bandwidth);
    out.status.push_back(name_of(t.status));
    out.up += t.status == TunnelStatus::Up ? 1 : 0;
    paths.push_back(t.edges.has_value() ? n.positions(*t.edges) : Ints {});
  }
  out.routes = paths;
  if (!in.offer.has_value()) return out;
  auto lossy = policy(n, in.heads, in.destinations, {.routing = dte, .admission = LY}, mask);
  const auto r = lossy->replace(*in.offer, finite(in.design_demand));
  out.no_drop_load.assign(n.m(), 0);
  for (std::size_t i = 0; i < r.by_tunnel.size(); ++i) {
    const Amount a = lo(r.by_tunnel[i].assigned);
    out.offered_by_entry.push_back(a);
    for (auto e : paths[i]) out.no_drop_load[static_cast<std::size_t>(e)] += a;
  }
  out.source_drops = r.call_source_drops;
  out.lossy_delivery = r.call_delivered;
  out.lossy_edge_drops = n.per_edge(r.call_edge_drops);
  auto lossless = policy(n, in.heads, in.destinations, {.routing = dte, .admission = LL}, mask);
  out.lossless_admitted = lossless->replace(*in.offer, finite(in.design_demand)).call_delivered;
  out.lossless_unadmitted = *in.offer - out.lossless_admitted;
  return out;
}

// ---------------------------------------------------------------- T2

struct T2In {
  std::vector<InEdge> edges;
  std::vector<NodeId> heads;
  std::vector<NodeId> destinations;
  std::int64_t per_head;
  std::optional<Amount> design_demand {};
  std::optional<Amounts> cap {};
  bool prefer_headroom {false};
  std::optional<Amount> lossy_offer {};
  std::optional<Ints> support {};
  Ints excluded {};
  bool up_to {false};
};
struct T2Out {
  TunnelPlan plan;
  Ints counts;
  Routes paths;
  std::vector<std::int32_t> entry_ids;
  Amount delivered {0};  // Lossless: the design's delivery (the design, or the bound without one)
  Amount lossy_delivery {0};
  Amount lossy_source_drops {0};
};
std::shared_ptr<FlowLedger> ledger_of(const Net& n, const std::optional<Amounts>& service) {
  if (!service.has_value()) return std::make_shared<FlowLedger>(n.g);
  const Amounts v = n.by_id(*service);
  return std::make_shared<FlowLedger>(n.g, v);
}
T2Out t2(const T2In& in) {
  const Net n(in.edges);
  const auto support = in.support.has_value() ? n.dag(*in.support) : nullptr;
  const auto cte = centralized(SourceAllocation::Equal, support, TunnelCount {in.per_head, in.up_to}, in.prefer_headroom);
  const auto mask = n.mask({}, in.excluded);
  const DesignArg design = in.design_demand.has_value() ? finite(*in.design_demand) : kUnlimited;
  auto p = policy(n, in.heads, in.destinations, {.routing = cte, .admission = LL}, mask, ledger_of(n, in.cap));
  const auto r = p->replace(in.design_demand, design);
  T2Out out;
  out.plan = *r.plan;
  out.delivered = r.call_delivered;
  out.counts.assign(n.m(), 0);
  for (const auto& t : r.program->tunnels()) {
    const Ints ps = n.positions(*t.edges);
    for (auto e : ps) ++out.counts[static_cast<std::size_t>(e)];
    out.paths.push_back(ps);
    out.entry_ids.push_back(t.entry_id);
  }
  if (in.lossy_offer.has_value()) {
    auto ly = policy(n, in.heads, in.destinations, {.routing = cte, .admission = LY}, mask, ledger_of(n, in.cap));
    const auto l = ly->replace(*in.lossy_offer, design);
    out.lossy_delivery = l.call_delivered;
    out.lossy_source_drops = l.call_source_drops;
  }
  return out;
}
std::pair<Amount, Amount> fraction(Amount num, Amount den = 1) { return {num, den}; }
Counter128 c128(std::uint64_t x) { return Counter128 {0, x}; }

// ---------------------------------------------------------------- programs of explicit routes

// Routes over abstract links: each abstract link e of capacity cap[e] gets its own two nodes, and
// each route runs from its head through its links (joined by links that never limit anything) to
// the one destination. head_of[i] is route i's head (0, 1, ...); the destination follows the heads.
struct Routed {
  Net net;
  std::vector<EdgeId> link;  // abstract link -> EdgeId
  std::vector<std::vector<EdgeId>> paths;
  std::vector<NodeId> heads;
  NodeId destination;
};
Routed routed(const Amounts& cap, const Routes& routes, const std::vector<NodeId>& head_of) {
  const NodeId h = *std::max_element(head_of.begin(), head_of.end()) + 1;
  const NodeId dest = h;
  auto tail = [&](std::int64_t e) { return static_cast<NodeId>(h + 1 + 2 * e); };
  auto headn = [&](std::int64_t e) { return static_cast<NodeId>(h + 2 + 2 * e); };
  std::vector<InEdge> rows;
  for (std::size_t e = 0; e < cap.size(); ++e) rows.push_back({tail(static_cast<std::int64_t>(e)), headn(static_cast<std::int64_t>(e)), cap[e], 1});
  std::vector<std::vector<std::size_t>> rows_of(routes.size());
  for (std::size_t i = 0; i < routes.size(); ++i) {
    const auto& r = routes[i];
    rows_of[i].push_back(rows.size());
    rows.push_back({head_of[i], tail(r.front()), kBig, 1});
    for (std::size_t j = 0; j < r.size(); ++j) {
      rows_of[i].push_back(static_cast<std::size_t>(r[j]));
      if (j + 1 < r.size()) {
        rows_of[i].push_back(rows.size());
        rows.push_back({headn(r[j]), tail(r[j + 1]), kBig, 1});
      }
    }
    rows_of[i].push_back(rows.size());
    rows.push_back({headn(r.back()), dest, kBig, 1});
  }
  Routed out {Net(rows, h + 1 + 2 * static_cast<NodeId>(cap.size())), {}, {}, {}, dest};
  for (std::size_t e = 0; e < cap.size(); ++e) out.link.push_back(out.net.id[e]);
  for (const auto& rs : rows_of) {
    std::vector<EdgeId> p;
    for (auto r : rs) p.push_back(out.net.id[r]);
    out.paths.push_back(p);
  }
  for (NodeId v = 0; v < h; ++v) out.heads.push_back(v);
  return out;
}
std::shared_ptr<const ForwardingProgram> routed_program(const Routed& r, const std::vector<NodeId>& head_of, const Amounts& weights) {
  std::vector<RouteDefinition> defs;
  for (std::size_t i = 0; i < r.paths.size(); ++i) {
    std::vector<std::vector<EdgeId>> hops;
    for (EdgeId e : r.paths[i]) hops.push_back({e});
    defs.push_back(RouteDefinition {head_of[i], std::make_shared<const Path>(*r.net.g, hops), weights[i]});
  }
  const NodeId t[1] = {r.destination};
  return ForwardingProgram::from_routes(*r.net.g, r.heads, t, defs);
}

// The one Lossy forwarding of routes with offers: one tunnel per route, weight = its offer, under
// head label[i] (the hierarchy's middle level; heads share the offer equally, so a head's routes
// must offer equal totals). Per route: what it carries and drops at each of its links, and what it
// delivers; per link, what it carries.
struct LossyOut {
  std::string error;
  Routes carried;
  Routes edge_drops;
  Amounts delivered;
  Amounts edge_flows;
};
LossyOut lossy(const Amounts& offers, const Routes& routes, const Amounts& cap, std::vector<NodeId> head_of = {}) {
  if (head_of.empty()) head_of.assign(routes.size(), 0);
  const Routed r = routed(cap, routes, head_of);
  const auto prog = routed_program(r, head_of, offers);
  auto p = policy(r.net, r.heads, {r.destination}, {.routing = Pinned {prog, nullptr}, .admission = LY});
  Amount total = 0;
  for (Amount x : offers) total += x;
  LossyOut out;
  PlacementResult res;
  out.error = error_of([&] { res = p->place(total, kFlows); });
  if (!out.error.empty()) return out;
  std::map<std::int32_t, FlowRecord> rec;
  for (auto& x : p->records()) rec[*x.entry_id] = x;
  for (std::size_t i = 0; i < routes.size(); ++i) {
    Ints c;
    Ints d;
    for (auto e : routes[i]) {
      const EdgeId id = r.link[static_cast<std::size_t>(e)];
      Amount cf = 0;
      Amount df = 0;
      if (rec.count(static_cast<std::int32_t>(i)) != 0) {
        for (const auto& [x, a] : rec[static_cast<std::int32_t>(i)].edge_flows) cf += x == id ? a : 0;
        for (const auto& [x, a] : rec[static_cast<std::int32_t>(i)].edge_drops) df += x == id ? lo(a) : 0;
      }
      c.push_back(cf);
      d.push_back(df);
    }
    out.carried.push_back(c);
    out.edge_drops.push_back(d);
    out.delivered.push_back(res.by_tunnel[i].delivered);
  }
  for (EdgeId id : r.link) out.edge_flows.push_back((*res.edge_flows)[static_cast<std::size_t>(id)]);
  return out;
}

// The Lossless bound of one head's routes with weights.
Amount bound(const Amounts& cap, const Routes& routes, const Amounts& weights) {
  const std::vector<NodeId> head_of(routes.size(), 0);
  const Routed r = routed(cap, routes, head_of);
  auto p = policy(r.net, r.heads, {r.destination}, {.routing = Pinned {routed_program(r, head_of, weights), nullptr}});
  return p->place(std::nullopt).call_delivered;
}

// ---------------------------------------------------------------- stale forwarding

struct StaleIn {
  std::vector<InEdge> edges;
  std::vector<NodeId> sources;
  NodeId destination;
  Amount offer;
  std::vector<std::optional<Amount>> weights;
  Ints black_links {};
  std::vector<NodeId> black_nodes {};
  Ints mask_links {};
  std::vector<NodeId> mask_nodes {};
};
struct StaleOut {
  Amount delivered;
  Amounts edge_flows;
  Amounts edge_drops;
  Amount source_drops;
};
StaleOut stale_dag(const StaleIn& in) {
  const Net n(in.edges);
  Ints members;
  bool equal = true;
  Amounts w(n.m(), 0);
  for (std::size_t e = 0; e < in.weights.size(); ++e) {
    if (!in.weights[e].has_value()) continue;
    members.push_back(static_cast<std::int64_t>(e));
    equal = equal && *in.weights[e] == 1;
    w[e] = *in.weights[e];
  }
  RouteOptions o;
  if (!equal) {
    o.allocation = RouteAllocation::Weighted;
    o.weight_basis = WeightBasis::Explicit;
    o.edge_weights = n.by_id(w);
  }
  const NodeId t[1] = {in.destination};
  const auto prog = ForwardingProgram::from_dag(*n.g, in.sources, t, n.dag(members), o);
  auto p = policy(n, in.sources, {in.destination}, {.routing = Pinned {prog, n.mask(in.black_nodes, in.black_links)}, .admission = LY},
                  n.mask(in.mask_nodes, in.mask_links));
  const auto r = p->place(in.offer, kFlows);
  return {r.call_delivered, n.in_order(*r.edge_flows), n.per_edge(r.call_edge_drops), r.call_source_drops};
}

struct StalePathsIn {
  std::vector<InEdge> edges;
  std::vector<std::pair<NodeId, Ints>> tunnels;
  std::map<NodeId, Amount> offers;
  Ints black_links {};
  std::vector<NodeId> black_nodes {};
  Ints mask_links {};
  std::vector<NodeId> mask_nodes {};
};
struct StalePathsOut {
  Texts status;
  Amounts assigned;
  Amounts delivered;
  Amounts edge_drops;
  Amount source_drops;
  std::optional<Amount> lossless_admitted;
};
StalePathsOut stale_paths(const StalePathsIn& in) {
  const Net n(in.edges);
  std::vector<NodeId> heads;
  Amount total = 0;
  for (const auto& [h, x] : in.offers) {
    heads.push_back(h);
    total += x;
    EXPECT_EQ(x, in.offers.begin()->second);  // heads share the offer equally
  }
  std::vector<RouteDefinition> defs;
  NodeId dest = 0;
  for (const auto& [h, es] : in.tunnels) {
    const auto p = n.path(es);
    dest = p->destination();
    defs.push_back(RouteDefinition {h, p, 1});
  }
  const NodeId t[1] = {dest};
  const auto prog = ForwardingProgram::from_routes(*n.g, heads, t, defs);
  const FlowPolicyConfig ly {.routing = Pinned {prog, n.mask(in.black_nodes, in.black_links)}, .admission = LY};
  auto p = policy(n, heads, {dest}, ly, n.mask(in.mask_nodes, in.mask_links));
  const auto r = p->place(total);
  StalePathsOut out;
  for (const auto& row : r.by_tunnel) {
    out.status.push_back(name_of(row.status));
    out.assigned.push_back(lo(row.assigned));
    out.delivered.push_back(row.delivered);
  }
  out.edge_drops = n.per_edge(r.call_edge_drops);
  out.source_drops = r.call_source_drops;
  if (heads.size() == 1 && r.call_source_drops == 0) {
    FlowPolicyConfig ll = ly;
    ll.admission = LL;
    auto q = policy(n, heads, {dest}, ll, n.mask(in.mask_nodes, in.mask_links));
    out.lossless_admitted = q->place(total).call_delivered;
  }
  return out;
}

// ---------------------------------------------------------------- the configuration's examples

// The example graph of the configuration cases: 0 -> 1, 0 -> 2, 1 -> 2; one source [0] or two
// [0, 1]; destination 2; "dag" is the support {0 -> 2}.
const Net& config_net() {
  static const Net n({{0, 1, 5, 1}, {0, 2, 5, 1}, {1, 2, 5, 1}});
  return n;
}

std::string resolved(const Routing& r) {
  auto opt = [](const auto& x, auto f) { return x.has_value() ? f(*x) : std::string("None"); };
  auto alloc = [](RouteAllocation a) { return std::string(a == RouteAllocation::Equal ? "EQUAL" : a == RouteAllocation::Weighted ? "WEIGHTED" : "OPTIMIZED"); };
  auto src = [](SourceAllocation s) { return std::string(s == SourceAllocation::Equal ? "EQUAL" : "OPTIMIZED"); };
  auto basis = [](WeightBasis b) { return std::string(b == WeightBasis::Capacity ? "CAPACITY" : b == WeightBasis::Explicit ? "EXPLICIT" : "RESIDUAL_SNAPSHOT"); };
  auto tie = [](RouteTieBreak t) {
    switch (t) {
      case RouteTieBreak::FewerHops: return std::string("FEWER_HOPS");
      case RouteTieBreak::HigherCapacity: return std::string("HIGHER_CAPACITY");
      case RouteTieBreak::HigherAvailableCapacity: return std::string("HIGHER_AVAILABLE_CAPACITY");
      case RouteTieBreak::HigherHeadroom: return std::string("HIGHER_HEADROOM");
    }
    return std::string("?");
  };
  auto boolean = [](bool b) { return std::string(b ? "True" : "False"); };
  if (const auto* sp = std::get_if<ShortestPaths>(&r)) {
    return "allocation=" + alloc(sp->allocation) + " source_allocation=" + opt(sp->source_allocation, src) +
           " weight_basis=" + opt(sp->weight_basis, basis) + " tie_break=" + tie(sp->tie_break) + " multipath=" + boolean(sp->multipath) +
           " parallel_edges=" + boolean(sp->parallel_edges);
  }
  if (const auto* gl = std::get_if<CentralizedTE>(&r)) {
    auto count = [](const TunnelCount& t) { return std::to_string(t.per_head) + (t.up_to ? " up_to" : ""); };
    return "source_allocation=" + src(gl->source_allocation) + " support=" + (gl->support ? "dag" : "None") +
           " tunnel_count=" + opt(gl->tunnel_count, count) + " prefer_headroom=" + boolean(gl->prefer_headroom);
  }
  return "valid";
}

// One configuration case: the policy's first error, else the resolved fields of the model (a
// zero call placed first: a zero offer skips no validation).
std::string config_case(Routing routing, int sources = 1, Admission admission = LL, bool cache = false,
                        std::optional<std::pair<std::vector<NodeId>, std::vector<NodeId>>> endpoints = std::nullopt) {
  const Net& n = config_net();
  std::vector<NodeId> s = sources == 1 ? v({0}) : v({0, 1});
  std::vector<NodeId> d = v({2});
  if (endpoints.has_value()) std::tie(s, d) = *endpoints;
  ShortestPathCache c(*n.g);
  FlowLedger fg(n.g);
  std::string got;
  const auto e = error_of([&] {
    FlowPolicy p(fg, s, d, {.routing = std::move(routing), .admission = admission}, nullptr, cache ? &c : nullptr);
    const bool tunnels = std::holds_alternative<DistributedTE>(p.config().routing) ||
                         (std::holds_alternative<CentralizedTE>(p.config().routing) &&
                          std::get<CentralizedTE>(p.config().routing).tunnel_count.has_value());
    if (tunnels) {
      (void)p.replace(0, finite(0));
    } else {
      (void)p.place(0);
    }
    got = resolved(p.config().routing);
  });
  return e.empty() ? got : e;
}
ShortestPaths sp_with(RouteAllocation a = RouteAllocation::Equal, std::optional<SourceAllocation> s = std::nullopt,
                      std::optional<WeightBasis> b = std::nullopt, RouteTieBreak t = RouteTieBreak::FewerHops, bool multipath = true,
                      bool parallel = true) {
  ShortestPaths sp;
  sp.allocation = a;
  sp.source_allocation = s;
  sp.weight_basis = b;
  sp.tie_break = t;
  sp.multipath = multipath;
  sp.parallel_edges = parallel;
  return sp;
}
DistributedTE dte_with(std::int64_t count, bool up_to = false, std::optional<CostCeiling> ceiling = std::nullopt, Amount minimum = 0) {
  DistributedTE d {TunnelCount {count, up_to}};
  d.cost_ceiling = ceiling;
  d.min_unreserved_bandwidth = minimum;
  return d;
}

// One evaluation case: the evaluator's first error, else the offer and design it evaluated.
std::string evaluation_case(Routing routing, Admission admission, DesignArg offer, DesignArg design) {
  const Net& n = config_net();
  EvaluationOptions o {{.routing = std::move(routing), .admission = admission}};
  o.offer = offer;
  o.design_demand = design;
  std::string got;
  const auto e = error_of([&] {
    const auto r = evaluate_demand(*n.g, v({0}), v({2}), o);
    const std::string offered = !r.requested.has_value() ? (r.by_source[0].offered.has_value() ? "first-hop experiment" : "inf")
                                                         : std::to_string(lo(r.requested));
    const std::string d = !r.plan.has_value() ? "None"
                          : r.plan->design_demand.has_value() ? std::to_string(*r.plan->design_demand)
                                                              : "inf";
    got = "offer=" + offered + " design_demand=" + d;
  });
  return e.empty() ? got : e;
}
Pinned config_pinned(const Net& n, NodeId to) {
  const std::vector<RouteDefinition> defs = {RouteDefinition {0, n.path({to == 1 ? 0 : 1}), 1}};
  const NodeId t[1] = {to};
  return Pinned {ForwardingProgram::from_routes(*n.g, v({0}), t, defs), nullptr};
}

// ---------------------------------------------------------------- calls on one policy

struct EpochOut {
  Amount requested;
  Amount delivered;
  Amount edge_drops;
  Amounts edge_flows;
  Amounts residual;
  Amounts per_call;
};
EpochOut epoch(const Amounts& cap, const Amounts& calls, bool full, Admission admission) {
  std::vector<InEdge> rows;
  for (Amount c : cap) rows.push_back({0, 1, c, 1});
  const Net n(rows);
  auto p = policy(n, v({0}), v({1}), {.routing = ShortestPaths {}, .admission = admission});
  EpochOut out {};
  PlacementResult r;
  if (full) {
    Amount total = 0;
    for (Amount x : calls) total += x;
    r = p->place(total, kFlows);
    out.per_call.push_back(r.call_delivered);
  } else {
    for (Amount x : calls) {
      r = p->place(x, kFlows);
      out.per_call.push_back(r.call_delivered);
    }
  }
  out.requested = lo(r.requested);
  out.delivered = r.delivered;
  for (const auto& [e, x] : r.edge_drops) out.edge_drops += lo(x);
  out.edge_flows = n.in_order(*r.edge_flows);
  out.residual = n.in_order(p.ledger->residual());
  return out;
}

// Sources 0 .. k - 1, each with its own link to node k (equal sources, ShortestPaths ECMP); before
// call c the other owners' load on link i is outside[c][i] (one Pinned policy per link).
struct IncrementsOut {
  Amounts call_delivered;
  Amount requested;
  std::optional<Amount> admitted;
  std::optional<Amount> unadmitted;
  std::optional<Amount> offered;
  Amount delivered;
  Amounts edge_drops;
};
IncrementsOut increments(const Amounts& caps, const Amounts& calls, Admission admission, std::vector<Amounts> outside = {}) {
  const auto k = static_cast<NodeId>(caps.size());
  std::vector<InEdge> rows;
  for (NodeId i = 0; i < k; ++i) rows.push_back({i, k, caps[static_cast<std::size_t>(i)], 1});
  const Net n(rows);
  auto ledger = std::make_shared<FlowLedger>(n.g);
  std::vector<Pol> others;
  for (NodeId i = 0; i < k && !outside.empty(); ++i) {
    const std::vector<RouteDefinition> defs = {RouteDefinition {i, n.path({i}), 1}};
    const NodeId s[1] = {i};
    const NodeId t[1] = {k};
    others.push_back(policy(n, {i}, {k}, {.routing = Pinned {ForwardingProgram::from_routes(*n.g, s, t, defs), nullptr}}, nullptr, ledger));
  }
  std::vector<NodeId> sources;
  for (NodeId i = 0; i < k; ++i) sources.push_back(i);
  auto p = policy(n, sources, {k}, {.routing = ShortestPaths {}, .admission = admission}, nullptr, ledger);
  IncrementsOut out {};
  PlacementResult r;
  for (std::size_t c = 0; c < calls.size(); ++c) {
    for (std::size_t i = 0; i < others.size(); ++i) (void)others[i]->replace(outside[c][i]);
    r = p->place(calls[c]);
    out.call_delivered.push_back(r.call_delivered);
  }
  out.requested = lo(r.requested);
  out.admitted = admission == LL ? std::optional<Amount>(lo(r.admitted)) : std::nullopt;
  out.unadmitted = r.unadmitted.has_value() ? std::optional<Amount>(lo(r.unadmitted)) : std::nullopt;
  out.offered = r.offered.has_value() ? std::optional<Amount>(lo(r.offered)) : std::nullopt;
  out.delivered = r.delivered;
  out.edge_drops = n.per_edge(r.edge_drops);
  return out;
}

// The installed program of a setup forwarded at `offer`: the program's tunnel count, the Lossy and
// Lossless terminations and the roster's per-source drops; forward works on the installed program.
struct TerminationOut {
  std::size_t tunnels;
  std::string lossy;
  std::string lossless;
  Amounts source_drops;
  std::string forward_error;
};
TerminationOut termination(const Net& n, const std::vector<NodeId>& s, const std::vector<NodeId>& d, const Routing& routing,
                           Amount offer, DesignArg design) {
  auto ly = policy(n, s, d, {.routing = routing, .admission = LY});
  auto ll = policy(n, s, d, {.routing = routing, .admission = LL});
  const bool tunnel_model = !std::holds_alternative<Pinned>(routing);
  const auto a = tunnel_model ? ly->replace(offer, design) : ly->place(offer);
  const auto b = tunnel_model ? ll->replace(offer, design) : ll->place(offer);
  TerminationOut out {a.program ? a.program->tunnels().size() : 0, name_of(a.termination), name_of(b.termination), {}, ""};
  for (const auto& row : a.by_source) out.source_drops.push_back(lo(row.source_drops));
  if (tunnel_model) out.forward_error = error_of([&] { (void)ly->forward(offer); });
  return out;
}

// The edges of a T1 example with heads 0 .. h - 1 and one link each to node h.
std::vector<InEdge> fan(NodeId h, Amount cap) {
  std::vector<InEdge> rows;
  for (NodeId i = 0; i < h; ++i) rows.push_back({i, h, cap, 1});
  return rows;
}

const std::vector<InEdge> kX8 = {{0, 2, 0, 1}, {1, 2, 8, 3}};
const std::vector<InEdge> kG1 = {{0, 1, 8192, 1}, {0, 2, 40960, 1}, {1, 3, 8192, 1}, {2, 3, 40960, 1}};
const std::vector<InEdge> kX7 = {{0, 1, 6, 1}, {0, 2, 6, 3}, {1, 2, 4, 1}, {1, 3, 6, 3}, {2, 3, 6, 1}};
const std::vector<InEdge> kG10 = {
    {0, 1, 40960, 1}, {0, 5, 40960, 1}, {0, 3, 24576, 2}, {1, 0, 40960, 1}, {1, 2, 40960, 1}, {1, 4, 16384, 2},
    {2, 1, 40960, 1}, {2, 3, 40960, 1}, {2, 5, 20480, 2}, {3, 2, 40960, 1}, {3, 4, 40960, 1}, {3, 0, 24576, 2},
    {4, 3, 40960, 1}, {4, 5, 40960, 1}, {4, 1, 16384, 2}, {5, 0, 40960, 1}, {5, 4, 40960, 1}, {5, 2, 20480, 2}};
const std::vector<InEdge> kDP = {{0, 2, 10, 1}, {1, 2, 10, 1}, {2, 3, 1, 1}, {2, 3, 1, 1}};
const std::vector<InEdge> kSup = {{0, 1, 6, 1}, {0, 1, 6, 1}, {0, 1, 100, 2}};

}  // namespace

// ================================================================ shared integer rules

TEST(ContractExamples, E_configured) { EXPECT_EQ(ha(12, {1, 2}), (Amounts {4, 8})); }
TEST(ContractExamples, E_snapshot) { EXPECT_EQ(ha(12, {2, 10}), (Amounts {2, 10})); }
TEST(ContractExamples, E_source_survivor) { EXPECT_EQ(ha(16, {3, 1}), (Amounts {12, 4})); }
TEST(ContractExamples, H_not_one_unit) { EXPECT_EQ(ha(50, {100, 1, 1, 1}), (Amounts {50, 0, 0, 0})); }

// R3-H-large, H(2^126 - 1, [2^125, 2^125]) = [2^125, 2^125 - 1], is outside the domain (amounts
// and weights below 2^62): the library refuses its weight and its offer, and the same rounding at
// the largest inputs of the domain gives the rule's values (the Python test computes the same
// values with the rule written out and checks them on the library there too). Each analogue runs
// on the rule itself and as a LOSSY offer over tunnels of those weights on links of capacity 1.
TEST(ContractExamples, R3_H_large) {
  constexpr Amount kTop = (Amount {1} << 62) - 1;
  const Net n({{0, 1, 1, 1}, {0, 1, 1, 1}});
  const NodeId s[1] = {0};
  const NodeId t[1] = {1};
  const std::vector<RouteDefinition> beyond = {RouteDefinition {0, n.path({0}), kTop + 1}};
  EXPECT_THROW((void)ForwardingProgram::from_routes(*n.g, s, t, beyond), std::invalid_argument);
  const std::vector<RouteDefinition> unit = {RouteDefinition {0, n.path({0}), 1}};
  const auto one = ForwardingProgram::from_routes(*n.g, s, t, unit);
  auto p = policy(n, {0}, {1}, {.routing = Pinned {one, nullptr}, .admission = LY});
  EXPECT_THROW((void)p->place(kTop + 1), std::invalid_argument);

  constexpr Amount half = Amount {1} << 61;
  const std::vector<std::pair<std::pair<Amount, Amounts>, Amounts>> analogues = {
      {{kTop, {half, half}}, {half, half - 1}},
      {{kTop, {kTop, kTop}}, {half, half - 1}},
      {{kTop, {kTop, kTop - 1, 1}}, {half, half - 1, 0}},
      {{kTop - 1, {kTop, half, 3}}, {3074457345618258601, 1537228672809129300, 1}},
      {{kTop, {3, 5, 7}}, {922337203685477580, 1537228672809129301, 2152120141932781022}},
  };
  for (const auto& [in, want] : analogues) {
    const auto& [total, weights] = in;
    EXPECT_EQ(ha(total, weights), want);
    std::vector<InEdge> links;
    std::vector<RouteDefinition> defs;
    for (std::size_t i = 0; i < weights.size(); ++i) links.push_back({0, 1, 1, 1});
    const Net m(links);
    for (std::size_t i = 0; i < weights.size(); ++i) {
      defs.push_back(RouteDefinition {0, m.path({static_cast<std::int64_t>(i)}), weights[i]});
    }
    auto q = policy(m, {0}, {1}, {.routing = Pinned {ForwardingProgram::from_routes(*m.g, s, t, defs), nullptr}, .admission = LY});
    Amounts got;
    for (const auto& row : q->place(total).by_tunnel) got.push_back(lo(row.assigned));
    EXPECT_EQ(got, want);
  }
}

// The 128-bit cumulative offers near 2^126 (beyond any one call, reached by increments): from
// X = 2^126 - 2 (split equally) one more unit goes to source 0; from 2^126 - 4, two more units
// split (1, 1).
TEST(CumulativeOffers, IncrementsNear2To126) {
  const Net n({{0, 2, 10, 1}, {1, 2, 10, 1}});
  const std::vector<RouteDefinition> defs = {RouteDefinition {0, n.path({0}), 1}, RouteDefinition {1, n.path({1}), 1}};
  const NodeId t[1] = {2};
  const auto prog = ForwardingProgram::from_routes(*n.g, v({0, 1}), t, defs);
  const auto view = detail::make_forwarding_view(*n.g, prog, nullptr, v({2}), nullptr, false);
  const Counter128 just_below {(std::uint64_t {1} << 62) - 1, ~std::uint64_t {0} - 1};  // 2^126 - 2
  Amounts out(2, 0);
  detail::source_offer_increments(view, just_below, 1, out);
  EXPECT_EQ(out, (Amounts {1, 0}));
  const Counter128 four_below {(std::uint64_t {1} << 62) - 1, ~std::uint64_t {0} - 3};  // 2^126 - 4
  detail::source_offer_increments(view, four_below, 2, out);
  EXPECT_EQ(out, (Amounts {1, 1}));
}

TEST(ContractExamples, H_local_bound) {
  Amounts w(11, 1);
  w[0] = 1000000;
  const auto e = split_error(1000000, w);
  EXPECT_EQ(e.assigned, (Amounts {1000000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}));
  EXPECT_EQ(text(e.max_excess), "1000000/100001");
  EXPECT_EQ(text(e.max_deficit), "100000/100001");
  EXPECT_TRUE(e.below_k);
}

TEST(ContractExamples, RR_local_bound) {
  const auto e = split_error(10, {1, 1, 1});
  EXPECT_EQ(e.assigned, (Amounts {4, 3, 3}));
  EXPECT_EQ(text(e.max_excess), "2/3");
  EXPECT_EQ(text(e.max_deficit), "1/3");
  EXPECT_TRUE(e.below_one);
}

// Every split with total <= 24, 1 to 4 members and weights 0 .. 4 (not all zero): each member within
// k of its exact share, within 1 when the weights are equal.
TEST(ContractExamples, H_bound_sweep) {
  std::int64_t cases = 0;
  std::int64_t violations = 0;
  std::int64_t equal_violations = 0;
  for (std::size_t k = 1; k <= 4; ++k) {
    std::size_t combos = 1;
    for (std::size_t i = 0; i < k; ++i) combos *= 5;
    for (std::size_t c = 0; c < combos; ++c) {
      Amounts w(k);
      std::size_t x = c;
      for (std::size_t i = k; i-- > 0;) {
        w[i] = static_cast<Amount>(x % 5);
        x /= 5;
      }
      if (std::all_of(w.begin(), w.end(), [](Amount y) { return y == 0; })) continue;
      const bool equal = std::all_of(w.begin(), w.end(), [&](Amount y) { return y == w[0]; });
      for (Amount n = 0; n <= 24; ++n) {
        ++cases;
        const auto e = split_error(n, w);
        violations += e.below_k ? 0 : 1;
        if (equal) equal_violations += e.below_one ? 0 : 1;
      }
    }
  }
  EXPECT_EQ(cases, 19400);
  EXPECT_EQ(violations, 0);
  EXPECT_EQ(equal_violations, 0);
}

// Members in prepared order (ascending EdgeId): member 1 (input order) has the smaller EdgeId (its
// lower cost sorts it first), so one unit goes to it.
TEST(ContractExamples, R3_member_order) {
  const Net n({{0, 1, 5, 2}, {0, 1, 5, 1}});
  ASSERT_LT(n.id[1], n.id[0]);
  const std::vector<RouteDefinition> defs = {RouteDefinition {0, n.dag({0, 1}), 1}};
  const NodeId t[1] = {1};
  auto p = policy(n, v({0}), v({1}), {.routing = Pinned {ForwardingProgram::from_routes(*n.g, v({0}), t, defs), nullptr}, .admission = LY});
  EXPECT_EQ(n.in_order(*p->place(1, kFlows).edge_flows), (Amounts {0, 1}));
}

// Two sources of one unit each meet at router 2, which splits the pooled two units once over its
// two members (one each), attributed back by H: member 0's unit to source 0, member 1's to source 1.
// Splitting each source separately puts both units on the first member.
TEST(ContractExamples, B_pooled) {
  const Net n({{0, 2, 10, 1}, {1, 2, 10, 1}, {2, 3, 10, 1}, {2, 3, 10, 1}});
  const NodeId t[1] = {3};
  auto pooled = policy(n, v({0, 1}), v({3}), {.routing = Pinned {ForwardingProgram::from_dag(*n.g, v({0, 1}), t, n.dag({0, 1, 2, 3})), nullptr}, .admission = LY});
  const auto r = pooled->place(2, kFlows);
  EXPECT_EQ(n.in_order(*r.edge_flows), (Amounts {1, 1, 1, 1}));  // columns [1, 1]
  std::vector<std::vector<Amount>> by_column(2, Amounts(2, 0));
  for (const auto& rec : pooled->records()) {
    const auto src = rec.sources.at(0);
    for (const auto& [e, x] : rec.edge_flows) {
      const auto c = n.pos(e) - 2;
      if (c >= 0) by_column[static_cast<std::size_t>(c)][static_cast<std::size_t>(src)] += x;
    }
  }
  EXPECT_EQ(by_column, (std::vector<Amounts> {{1, 0}, {0, 1}}));
  const std::vector<RouteDefinition> defs = {RouteDefinition {0, n.dag({0, 2, 3}), 1}, RouteDefinition {1, n.dag({1, 2, 3}), 1}};
  auto separate = policy(n, v({0, 1}), v({3}), {.routing = Pinned {ForwardingProgram::from_routes(*n.g, v({0, 1}), t, defs), nullptr}, .admission = LY});
  const auto flows = n.in_order(*separate->place(2, kFlows).edge_flows);
  EXPECT_EQ((Amounts {flows[2], flows[3]}), (Amounts {2, 0}));  // independent columns
}

// ================================================================ searches

const std::vector<InEdge> kZeroEdge = {{0, 1, 0, 1}, {0, 2, 5, 1}, {2, 1, 5, 1}};

TEST(ContractExamples, SPF_zero_edge) {
  const auto r = search(kZeroEdge, 0, 1);
  EXPECT_EQ(r.distance, std::optional<Cost>(1));
  EXPECT_EQ(r.dag, (Ints {0}));
}
TEST(ContractExamples, SPF_zero_edge_residual) {
  const auto r = search(kZeroEdge, 0, 1, Amounts {0, 5, 5});
  EXPECT_EQ(r.distance, std::optional<Cost>(2));
  EXPECT_EQ(r.dag, (Ints {1, 2}));
}
TEST(ContractExamples, KSP_zero_edge) {
  const auto r = search(kZeroEdge, 0, 1, std::nullopt, 2);
  EXPECT_EQ(r.paths, (std::vector<std::pair<Cost, Ints>> {{1, {0}}, {2, {1, 2}}}));
}
TEST(ContractExamples, KSP_zero_edge_residual) {
  const auto r = search(kZeroEdge, 0, 1, Amounts {0, 5, 5}, 2);
  EXPECT_EQ(r.paths, (std::vector<std::pair<Cost, Ints>> {{2, {1, 2}}}));
}

// A ceiling of another type is a Python value; in C++ the ceiling is a variant of the two types.
// The valid part: with a residual of the wrong length, ksp refuses the residual (the type check
// that comes first has nothing to refuse here).
TEST(ContractExamples, KSP_ceiling_type) {
  static_assert(std::variant_size_v<CostCeiling> == 2);
  const auto r = search({{0, 1, 10, 5}}, 0, 1, Amounts {10, 10}, 2);
  EXPECT_EQ(r.error, "residual has length 2, expected 1 (one entry per edge)");
}
TEST(ContractExamples, KSP_stretch_below_one) {
  const auto r = search({{0, 1, 10, 5}}, 0, 1, std::nullopt, 2, CostStretch(1, 2));
  EXPECT_EQ(r.error, "cost_ceiling CostStretch factor must be >= 1 (got 1/2)");
}
TEST(ContractExamples, KSP_ceiling_before_residual) {
  const auto r = search({{0, 1, 10, 5}}, 0, 1, Amounts {10, 10}, 2, CostStretch(1, 2));
  EXPECT_EQ(r.error, "cost_ceiling CostStretch factor must be >= 1 (got 1/2)");
}
TEST(ContractExamples, KSP_residual_length) {
  const auto r = search({{0, 1, 10, 5}}, 0, 1, Amounts {10, 10}, 2, CostStretch(3, 2));
  EXPECT_EQ(r.error, "residual has length 2, expected 1 (one entry per edge)");
}
// numpy.float32(1.2) is the fraction 6/5 (its own shortest decimal).
TEST(ContractExamples, KSP_stretch_numpy) {
  const auto r = search({{0, 1, 0, 100000000}, {0, 1, 10, 120000001}}, 0, 1, Amounts {0, 10}, 2, CostStretch(6, 5));
  EXPECT_EQ(r.error, "");
  EXPECT_TRUE(r.paths.empty());
}
TEST(ContractExamples, KSP_stretch_numpy_equal) {
  const auto r = search({{0, 1, 0, 100000000}, {0, 1, 10, 120000000}}, 0, 1, Amounts {0, 10}, 2, CostStretch(6, 5));
  EXPECT_EQ(r.paths, (std::vector<std::pair<Cost, Ints>> {{120000000, {1}}}));
}
// The float 1.2 is the fraction 6/5.
TEST(ContractExamples, KSP_stretch_boundary) {
  const auto r = search({{0, 1, 0, 5}, {0, 1, 10, 6}}, 0, 1, Amounts {0, 10}, 2, CostStretch(6, 5));
  EXPECT_EQ(r.paths, (std::vector<std::pair<Cost, Ints>> {{6, {1}}}));
}
TEST(ContractExamples, KSP_stretch_below) {
  const auto r = search({{0, 1, 0, 5}, {0, 1, 10, 6}}, 0, 1, Amounts {0, 10}, 2, CostStretch(11, 10));
  EXPECT_EQ(r.error, "");
  EXPECT_TRUE(r.paths.empty());
}

TEST(ContractExamples, TIE_parallel) {
  const std::vector<InEdge> e = {{0, 1, 4, 1}, {0, 1, 10, 1}, {0, 1, 100, 1}, {0, 1, 1000, 2}};
  const Amounts avail = {4, 8, 5, 1000};
  EXPECT_EQ(tie_choice(e, 0, 1, RouteTieBreak::FewerHops, avail), (Ints {0}));
  EXPECT_EQ(tie_choice(e, 0, 1, RouteTieBreak::HigherCapacity, avail), (Ints {2}));
  EXPECT_EQ(tie_choice(e, 0, 1, RouteTieBreak::HigherAvailableCapacity, avail), (Ints {1}));
  EXPECT_EQ(tie_choice(e, 0, 1, RouteTieBreak::HigherHeadroom, avail), (Ints {0}));
}
TEST(ContractExamples, TIE_routes) {
  const std::vector<InEdge> e = {{0, 1, 10, 1}, {0, 2, 3, 2}, {1, 2, 10, 1}};
  EXPECT_EQ(tie_choice(e, 0, 2, RouteTieBreak::FewerHops), (Ints {1}));
  EXPECT_EQ(tie_choice(e, 0, 2, RouteTieBreak::HigherCapacity), (Ints {0, 2}));
  EXPECT_EQ(tie_choice(e, 0, 2, RouteTieBreak::HigherAvailableCapacity), (Ints {0, 2}));
  EXPECT_EQ(tie_choice(e, 0, 2, RouteTieBreak::HigherHeadroom), (Ints {1}));
}
TEST(ContractExamples, TIE_bottleneck) {
  const std::vector<InEdge> e = {{0, 1, 10, 1}, {0, 2, 3, 2}, {1, 2, 1, 1}};
  for (const auto t : {RouteTieBreak::FewerHops, RouteTieBreak::HigherCapacity, RouteTieBreak::HigherAvailableCapacity,
                       RouteTieBreak::HigherHeadroom}) {
    EXPECT_EQ(tie_choice(e, 0, 2, t), (Ints {1}));
  }
}

// ================================================================ exact flows

const std::vector<InEdge> kB1 = {{0, 3, 2, 1}, {1, 2, 10, 1}, {2, 3, 10, 2}};
const std::vector<InEdge> kB2 = {{0, 1, 4, 1}, {0, 2, 2, 1}, {1, 2, 10, 1}};

void expect_flow(const FlowOut& got, Amount delivered, Amount cost, const Amounts& flows) {
  EXPECT_EQ(got.delivered, delivered);
  EXPECT_EQ(got.cost, cost);
  EXPECT_EQ(got.edge_flows, flows);
}
TEST(ContractExamples, B1_equal) { expect_flow(exact_flow(kB1, v({0, 1}), v({3}), 12, SourceAllocation::Equal), 4, 8, {2, 2, 2}); }
TEST(ContractExamples, B1_optimized) {
  expect_flow(exact_flow(kB1, v({0, 1}), v({3}), 12, SourceAllocation::Optimized), 12, 32, {2, 10, 10});
}
TEST(ContractExamples, MINCUT_B1) {
  const auto r = exact_flow(kB1, v({0, 1}), v({3}), std::nullopt, SourceAllocation::Optimized);
  EXPECT_EQ(r.delivered, 12);
  EXPECT_EQ(r.min_cut, std::optional<Ints>(Ints {0, 1}));
}
TEST(ContractExamples, MINCUT_zero_edge) {
  const auto r = exact_flow({{0, 2, 0, 1}, {1, 2, 10, 1}}, v({0, 1}), v({2}), std::nullopt, SourceAllocation::Optimized);
  EXPECT_EQ(r.delivered, 10);
  EXPECT_EQ(r.min_cut, std::optional<Ints>(Ints {0, 1}));
}
TEST(ContractExamples, B2_equal_global) { expect_flow(exact_flow(kB2, v({0, 1}), v({2}), 12, SourceAllocation::Equal), 12, 16, {4, 2, 10}); }
TEST(ContractExamples, B2_shortest) {
  expect_flow(exact_flow({{0, 2, 2, 1}, {1, 2, 10, 1}}, v({0, 1}), v({2}), 12, SourceAllocation::Equal), 4, 4, {2, 2});
}
TEST(ContractExamples, B2_optimized) {
  expect_flow(exact_flow(kB2, v({0, 1}), v({2}), 12, SourceAllocation::Optimized), 12, 12, {0, 2, 10});
}
TEST(ContractExamples, B3_zero_first) {
  expect_flow(exact_flow({{0, 2, 0, 1}, {1, 2, 10, 1}}, v({0, 1}), v({2}), 10, SourceAllocation::Equal), 0, 0, {0, 0});
}
TEST(ContractExamples, B3_positive_first) {
  expect_flow(exact_flow({{0, 2, 10, 1}, {1, 2, 0, 1}}, v({0, 1}), v({2}), 10, SourceAllocation::Equal), 1, 1, {1, 0});
}
TEST(ContractExamples, B3_delete_zero) {
  expect_flow(exact_flow({{0, 2, 0, 1}, {1, 2, 10, 1}}, v({1}), v({2}), 10, SourceAllocation::Equal), 10, 10, {0, 10});
}
TEST(ContractExamples, E_support) {
  expect_flow(exact_flow({{0, 1, 10, 1}, {0, 2, 10, 1}, {1, 3, 2, 1}, {2, 3, 10, 2}}, v({0}), v({3}), 12, SourceAllocation::Equal), 12, 34,
              {2, 10, 2, 10});
}
TEST(ContractExamples, E_support_shortest) {
  expect_flow(exact_flow({{0, 1, 10, 1}, {1, 3, 2, 1}}, v({0}), v({3}), 12, SourceAllocation::Equal), 2, 4, {2, 2});
}
TEST(ContractExamples, T2_OPT_global) {
  expect_flow(exact_flow({{0, 2, 2, 1}, {1, 2, 10, 1}}, v({0, 1}), v({2}), 12, SourceAllocation::Optimized), 12, 12, {2, 10});
}

// Two sources with private links (capacities 2 and 10, costs 9 and 1), equal sources, Lossless
// increments of 3: the second call can add only 1 (both sources must carry the same).
TEST(ContractExamples, R4_lossless_equal) {
  const Net n({{0, 2, 2, 9}, {1, 2, 10, 1}});
  auto p = policy(n, v({0, 1}), v({2}), {.routing = centralized(SourceAllocation::Equal)});
  const auto a = p->place(3);
  const auto b = p->place(3, kFlows);
  EXPECT_EQ((Amounts {a.call_delivered, b.call_delivered}), (Amounts {3, 1}));
  EXPECT_EQ(n.in_order(*b.edge_flows), (Amounts {2, 2}));
  EXPECT_EQ(b.delivered, 4);
  EXPECT_EQ(lo(b.total_cost), 20);
}

// One unit from 0 to 2 over 0 -> 2 (cost 1) or 0 -> 1 -> 0 -> 2 (a cycle 0 -> 1 -> 0 of cost 2): the
// library's flow is the minimum-cost one, [0, 1, 0] at cost 1, with no positive-flow cycle; the
// cyclic flow [1, 1, 1] costs 3, above that minimum.
const std::vector<InEdge> kWitness = {{0, 1, 1, 1}, {0, 2, 1, 1}, {1, 0, 1, 1}};
TEST(ContractExamples, R3_residual_optimal) {
  const auto r = exact_flow(kWitness, v({0}), v({2}), 1, SourceAllocation::Optimized);
  EXPECT_EQ(r.edge_flows, (Amounts {0, 1, 0}));
  EXPECT_EQ(r.cost, 1);
}
TEST(ContractExamples, R3_residual_cyclic) {
  const auto r = exact_flow(kWitness, v({0}), v({2}), 1, SourceAllocation::Optimized);
  EXPECT_EQ(r.cost, 1);  // the minimum cost at value 1
  const Amounts cyclic = {1, 1, 1};
  Amount cost = 0;
  for (std::size_t e = 0; e < cyclic.size(); ++e) cost += cyclic[e] * kWitness[e].cost;
  EXPECT_EQ(cost, 3);
  EXPECT_NE(r.edge_flows, cyclic);  // not the minimum-cost flow
}

// ================================================================ T1

void expect_t1_base(const T1Out& r, const Amounts& requested, const Routes& routes, std::int64_t up, Amount reserved,
                    const Amounts& reservations, Amount bound_) {
  EXPECT_EQ(r.error, "");
  EXPECT_EQ(r.requested, requested);
  EXPECT_EQ(r.routes, routes);
  EXPECT_EQ(r.up, up);
  EXPECT_EQ(r.reserved, reserved);
  EXPECT_EQ(r.edge_reservations, reservations);
  EXPECT_EQ(r.bound, bound_);
}
void expect_t1_offer(const T1Out& r, const Amounts& offered, Amount source_drops, const Amounts& no_drop_load) {
  EXPECT_EQ(r.offered_by_entry, offered);
  EXPECT_EQ(r.source_drops, source_drops);
  EXPECT_EQ(r.no_drop_load, no_drop_load);
}

TEST(ContractExamples, A_X1) {
  const auto r = t1({{{0, 1, 6, 1}, {0, 1, 6, 2}}, v({0}), v({1}), 12, 4, 13});
  expect_t1_base(r, {3, 3, 3, 3}, {{0}, {0}, {1}, {1}}, 4, 12, {6, 6}, 12);
  expect_t1_offer(r, {4, 3, 3, 3}, 0, {7, 6});
  EXPECT_EQ(r.lossy_delivery, 12);
  EXPECT_EQ(r.lossy_edge_drops, (Amounts {1, 0}));
  EXPECT_EQ(r.lossless_admitted, 12);
  EXPECT_EQ(r.lossless_unadmitted, 1);
}
TEST(ContractExamples, A_X2) {
  const auto r = t1({{{0, 1, 3, 1}, {0, 1, 2, 2}}, v({0}), v({1}), 5, 2, 5});
  expect_t1_base(r, {3, 2}, {{0}, {1}}, 2, 5, {3, 2}, 5);
  expect_t1_offer(r, {3, 2}, 0, {3, 2});
}
TEST(ContractExamples, A_X2_small_share) {
  const auto r = t1({{{0, 1, 6, 1}, {0, 1, 6, 2}}, v({0}), v({1}), 1, 2, 1});
  expect_t1_base(r, {1, 0}, {{0}, {0}}, 2, 1, {1, 0}, 6);
  expect_t1_offer(r, {1, 0}, 0, {1, 0});
  EXPECT_EQ(r.status, (Texts {"UP", "UP"}));
}
TEST(ContractExamples, T1_over_reservation) {
  const auto r = t1({{{0, 1, 5, 1}, {0, 1, 5, 1}}, v({0}), v({1}), 8, 2, 10});
  EXPECT_EQ(r.requested, (Amounts {4, 4}));
  EXPECT_EQ(r.routes, (Routes {{0}, {1}}));
  EXPECT_EQ(r.edge_reservations, (Amounts {4, 4}));
  EXPECT_EQ(r.offered_by_entry, (Amounts {5, 5}));
  EXPECT_EQ(r.no_drop_load, (Amounts {5, 5}));
  EXPECT_EQ(r.lossy_delivery, 10);
  EXPECT_EQ(r.lossy_edge_drops, (Amounts {0, 0}));
  EXPECT_EQ(r.bound, 10);
}
TEST(ContractExamples, T1_zero_link_drop) {
  const auto r = t1({{{0, 1, 0, 1}}, v({0}), v({1}), 0, 1, 1});
  EXPECT_EQ(r.requested, (Amounts {0}));
  EXPECT_EQ(r.routes, (Routes {{0}}));
  EXPECT_EQ(r.status, (Texts {"UP"}));
  EXPECT_EQ(r.offered_by_entry, (Amounts {1}));
  EXPECT_EQ(r.source_drops, 0);
  EXPECT_EQ(r.lossy_delivery, 0);
  EXPECT_EQ(r.lossy_edge_drops, (Amounts {1}));
  EXPECT_EQ(r.bound, 0);
  EXPECT_EQ(r.lossless_admitted, 0);
  EXPECT_EQ(r.lossless_unadmitted, 1);
}
TEST(ContractExamples, A_X3) {
  const auto r = t1({{{0, 1, 2, 1}, {0, 1, 10, 2}}, v({0}), v({1}), 12, 4, 12});
  expect_t1_base(r, {3, 3, 3, 3}, {{1}, {1}, {1}, {}}, 3, 9, {0, 9}, 10);
  expect_t1_offer(r, {4, 4, 4, 0}, 0, {0, 12});
  EXPECT_EQ(r.status, (Texts {"UP", "UP", "UP", "INSUFFICIENT_BANDWIDTH"}));
}
TEST(ContractExamples, A_X4) {
  const auto r = t1({{{0, 2, 6, 1}, {1, 2, 6, 1}, {2, 3, 6, 1}}, v({0, 1}), v({3}), 12, 2, 12});
  expect_t1_base(r, {3, 3, 3, 3}, {{0, 2}, {1, 2}, {}, {}}, 2, 6, {3, 3, 6}, 6);
  expect_t1_offer(r, {6, 6, 0, 0}, 0, {6, 6, 12});
}
TEST(ContractExamples, A_X5) {
  const auto r = t1({fan(3, 20), v({0, 1, 2}), v({3}), 2, 1, 3});
  expect_t1_base(r, {1, 1, 0}, {{0}, {1}, {2}}, 3, 2, {1, 1, 0}, 60);
  expect_t1_offer(r, {1, 1, 1}, 0, {1, 1, 1});
}
TEST(ContractExamples, A_X6) {
  const auto r = t1({{{0, 1, 5, 1}, {0, 1, 5, 1}}, v({0}), v({1}), 8, 2, 8});
  expect_t1_base(r, {4, 4}, {{0}, {1}}, 2, 8, {4, 4}, 10);
  expect_t1_offer(r, {4, 4}, 0, {4, 4});
}
TEST(ContractExamples, A_X8) {
  const auto r = t1({kX8, v({0, 1}), v({2}), 4, 1, 4});
  expect_t1_base(r, {2, 2}, {{}, {1}}, 1, 2, {0, 2}, 0);
  expect_t1_offer(r, {0, 2}, 2, {0, 2});
}
TEST(ContractExamples, A_X9) {
  const auto r = t1({kG1, v({0}), v({3}), 49152, 16, 49152});
  Routes routes = {{0, 2}, {0, 2}};
  for (int i = 0; i < 13; ++i) routes.push_back({1, 3});
  routes.push_back({});
  expect_t1_base(r, Amounts(16, 3072), routes, 15, 46080, {6144, 39936, 6144, 39936}, 47262);
  Amounts offered(12, 3277);
  for (int i = 0; i < 3; ++i) offered.push_back(3276);
  offered.push_back(0);
  expect_t1_offer(r, offered, 0, {6554, 42598, 6554, 42598});
}
TEST(ContractExamples, A_X9_large_design) {
  const auto r = t1({kG1, v({0}), v({3}), 409600, 16});
  Routes routes = {{1, 3}};
  for (int i = 0; i < 15; ++i) routes.push_back({});
  expect_t1_base(r, Amounts(16, 25600), routes, 1, 25600, {0, 25600, 0, 25600}, 40960);
}
TEST(ContractExamples, A_X8_delete_zero_edge) {
  T1In in {kX8, v({0, 1}), v({2}), 4, 1};
  in.excluded = {0};
  expect_t1_base(t1(in), {4}, {{1}}, 1, 4, {0, 4}, 8);
}
TEST(ContractExamples, A_X8_Stretch2) {
  T1In in {kX8, v({0, 1}), v({2}), 4, 1};
  in.stretch = CostStretch(2);
  const auto r = t1(in);
  expect_t1_base(r, {2, 2}, {{}, {}}, 0, 0, {0, 0}, 0);
  EXPECT_EQ(r.status, (Texts {"INSUFFICIENT_BANDWIDTH", "COST_CEILING_EXCEEDED"}));
}
TEST(ContractExamples, A_X8_Absolute3) {
  T1In in {kX8, v({0, 1}), v({2}), 4, 1};
  in.absolute = 3;
  const auto r = t1(in);
  expect_t1_base(r, {2, 2}, {{}, {1}}, 1, 2, {0, 2}, 0);
  EXPECT_EQ(r.status, (Texts {"INSUFFICIENT_BANDWIDTH", "UP"}));
}
// The float factor 1.2 is the fraction 6/5: cost 6 against 5 * 6/5 = 6 is admitted.
TEST(ContractExamples, A_X8_Stretch_exact) {
  T1In in {{{0, 1, 0, 5}, {0, 1, 1, 6}}, v({0}), v({1}), 1, 1, 1};
  in.stretch = CostStretch(6, 5);
  const auto r = t1(in);
  EXPECT_EQ(r.status, (Texts {"UP"}));
  EXPECT_EQ(r.routes, (Routes {{1}}));
  EXPECT_EQ(r.up, 1);
  EXPECT_EQ(r.lossy_delivery, 1);
  EXPECT_EQ(r.source_drops, 0);
  EXPECT_EQ(r.lossless_admitted, 1);
}
TEST(ContractExamples, A_X8_Stretch_cost7) {
  T1In in {{{0, 1, 0, 5}, {0, 1, 1, 7}}, v({0}), v({1}), 1, 1, 1};
  in.stretch = CostStretch(6, 5);
  const auto r = t1(in);
  EXPECT_EQ(r.status, (Texts {"COST_CEILING_EXCEEDED"}));
  EXPECT_EQ(r.routes, (Routes {{}}));
  EXPECT_EQ(r.up, 0);
  EXPECT_EQ(r.lossy_delivery, 0);
  EXPECT_EQ(r.source_drops, 1);
  EXPECT_EQ(r.lossless_admitted, 0);
}
TEST(ContractExamples, A_X7_8) { expect_t1_base(t1({kX7, v({0}), v({3}), 8, 2}), {4, 4}, {{0, 2, 4}, {}}, 1, 4, {4, 0, 4, 0, 4}, 4); }
TEST(ContractExamples, A_X7_10) { expect_t1_base(t1({kX7, v({0}), v({3}), 10, 2}), {5, 5}, {{0, 3}, {1, 4}}, 2, 10, {5, 5, 0, 5, 5}, 12); }
TEST(ContractExamples, A_X7_14) { expect_t1_base(t1({kX7, v({0}), v({3}), 14, 2}), {7, 7}, {{}, {}}, 0, 0, {0, 0, 0, 0, 0}, 0); }
TEST(ContractExamples, T1_two_heads) {
  const auto r = t1({{{0, 2, 2, 1}, {1, 2, 10, 1}}, v({0, 1}), v({2}), 12, 2, 12});
  expect_t1_base(r, {3, 3, 3, 3}, {{}, {1}, {}, {1}}, 2, 6, {0, 6}, 0);
  EXPECT_EQ(r.status, (Texts {"INSUFFICIENT_BANDWIDTH", "UP", "INSUFFICIENT_BANDWIDTH", "UP"}));
  expect_t1_offer(r, {0, 3, 0, 3}, 6, {0, 6});
  EXPECT_EQ(r.lossy_delivery, 6);
}
TEST(ContractExamples, T1_idle) {
  const auto r = t1({{{0, 1, 10, 1}}, v({0}), v({1}), 3, 4, 3});
  expect_t1_base(r, {1, 1, 1, 0}, {{0}, {0}, {0}, {0}}, 4, 3, {3}, 10);
  expect_t1_offer(r, {1, 1, 1, 0}, 0, {3});
  EXPECT_EQ(r.lossy_delivery, 3);
}
TEST(ContractExamples, T1_idle_zero_link) {
  const auto r = t1({{{0, 1, 0, 1}, {0, 1, 10, 2}}, v({0}), v({1}), 1, 2, 2});
  expect_t1_base(r, {1, 0}, {{1}, {0}}, 2, 1, {0, 1}, 1);
  EXPECT_EQ(r.status, (Texts {"UP", "UP"}));
  expect_t1_offer(r, {1, 1}, 0, {1, 1});
  EXPECT_EQ(r.lossy_delivery, 1);
}
TEST(ContractExamples, T1_zero_design) {
  const auto r = t1({{{0, 1, 10, 1}}, v({0}), v({1}), 0, 2, 4});
  expect_t1_base(r, {0, 0}, {{0}, {0}}, 2, 0, {0}, 10);
  expect_t1_offer(r, {2, 2}, 0, {4});
  EXPECT_EQ(r.lossy_delivery, 4);
}
TEST(ContractExamples, T1_stretch_below_one) {
  T1In in {{{0, 1, 10, 5}}, v({0}), v({1}), 1, 1, 1};
  in.stretch = CostStretch(1, 2);
  EXPECT_EQ(t1(in).error, "cost_ceiling CostStretch factor must be >= 1 (got 1/2)");
}
// Two heads with two tunnels each, design 5: RR(5, 2) = (3, 2) per head, then per tunnel in round
// order: (2, 1) for head 0's rounds and (1, 1) for head 1's.
TEST(ContractExamples, T1_split) {
  const auto r = t1({fan(2, 100), v({0, 1}), v({2}), 5, 2});
  EXPECT_EQ(r.requested, (Amounts {2, 1, 1, 1}));
}

// ================================================================ T2

void expect_t2_core(const T2Out& r, Amount U, const Ints& counts, std::pair<Amount, Amount> fluid, Amount max_lossless) {
  EXPECT_EQ(r.plan.max_bandwidth_per_tunnel, U);
  EXPECT_EQ(r.counts, counts);
  EXPECT_EQ(r.plan.max_fractional_bandwidth_per_tunnel, fluid);
  EXPECT_EQ(r.plan.max_lossless_offer, max_lossless);
}
Routes repeat(const Routes& parts, const std::vector<int>& times) {
  Routes out;
  for (std::size_t i = 0; i < parts.size(); ++i) {
    for (int t = 0; t < times[i]; ++t) out.push_back(parts[i]);
  }
  return out;
}

TEST(ContractExamples, C_G1) {
  T2In in {kG1, v({0}), v({3}), 16, 49152};
  in.prefer_headroom = true;
  in.lossy_offer = 49152;
  const auto r = t2(in);
  expect_t2_core(r, 2925, {2, 14, 2, 14}, fraction(20480, 7), 46812);
  EXPECT_EQ(r.lossy_delivery, 47104);
  EXPECT_EQ(r.paths, repeat({{0, 2}, {1, 3}}, {2, 14}));
  EXPECT_EQ(r.plan.headroom_penalty, c128(73610035200));
  EXPECT_EQ(r.plan.total_planned_bandwidth, 46800);
  EXPECT_EQ(r.plan.total_tunnel_cost, c128(32));
  EXPECT_EQ(r.plan.planned_bandwidth_per_tunnel, 2925);
}
TEST(ContractExamples, C_G1_finite) {
  T2In in {kG1, v({0}), v({3}), 16, 8192};
  in.prefer_headroom = true;
  const auto r = t2(in);
  expect_t2_core(r, 2925, {2, 14, 2, 14}, fraction(20480, 7), 46812);
  EXPECT_EQ(r.paths, repeat({{0, 2}, {1, 3}}, {2, 14}));
  EXPECT_EQ(r.plan.headroom_penalty, c128(12884901876));
  EXPECT_EQ(r.plan.total_planned_bandwidth, 8192);
  EXPECT_EQ(r.plan.total_tunnel_cost, c128(32));
  EXPECT_EQ(r.plan.planned_bandwidth_per_tunnel, 512);
}
TEST(ContractExamples, C_G4) {
  std::vector<InEdge> e;
  for (NodeId i = 1; i <= 15; ++i) e.push_back({0, i, 4096, 1});
  for (NodeId i = 1; i <= 15; ++i) e.push_back({i, 16, 4096, 1});
  T2In in {e, v({0}), v({16}), 16, 61440};
  in.prefer_headroom = true;
  in.lossy_offer = 61440;
  const auto r = t2(in);
  Ints counts(30, 1);
  counts[0] = 2;
  counts[15] = 2;
  expect_t2_core(r, 2048, counts, fraction(2048), 32768);
  EXPECT_EQ(r.plan.planned_bandwidth_per_tunnel, 2048);
  Routes paths = {{0, 15}, {0, 15}};
  for (std::int64_t i = 1; i < 15; ++i) paths.push_back({i, 15 + i});
  EXPECT_EQ(r.paths, paths);
  EXPECT_EQ(r.lossy_delivery, 57856);
  EXPECT_EQ(r.plan.headroom_penalty, c128(73014444032));
  EXPECT_EQ(r.plan.total_planned_bandwidth, 32768);
  EXPECT_EQ(r.plan.total_tunnel_cost, c128(32));
}
TEST(ContractExamples, C_equal) {
  T2In in {{{0, 2, 3, 1}, {1, 2, 8, 1}}, v({0, 1}), v({2}), 2, 6};
  in.prefer_headroom = true;
  const auto r = t2(in);
  expect_t2_core(r, 1, {2, 2}, fraction(3, 2), 6);
  EXPECT_EQ(r.paths, (Routes {{0}, {1}, {0}, {1}}));
  EXPECT_EQ(r.plan.headroom_penalty, c128(5905580031));
  EXPECT_EQ(r.plan.total_planned_bandwidth, 4);
  EXPECT_EQ(r.plan.total_tunnel_cost, c128(4));
  EXPECT_EQ(r.plan.planned_bandwidth_per_tunnel, 1);
}
TEST(ContractExamples, C_finite_cost) {
  T2In in {{{0, 1, 2, 1}, {0, 1, 10, 2}}, v({0}), v({1}), 2, 2};
  in.prefer_headroom = true;
  const auto r = t2(in);
  expect_t2_core(r, 5, {2, 0}, fraction(5), 2);
  EXPECT_EQ(r.paths, (Routes {{0}, {0}}));
  EXPECT_EQ(r.plan.headroom_penalty, c128(6442450944));
  EXPECT_EQ(r.plan.total_planned_bandwidth, 2);
  EXPECT_EQ(r.plan.total_tunnel_cost, c128(2));
  EXPECT_EQ(r.plan.planned_bandwidth_per_tunnel, 1);
}
TEST(ContractExamples, C_G10_0) {
  T2In in {kG10, v({0}), v({3}), 16, 49152};
  in.prefer_headroom = true;
  in.cap = Amounts {40960, 40960, 24576, 40960, 40960, 16384, 40960, 40960, 20480, 40960, 40960, 24576, 40960, 40960, 16384, 40960, 40960, 20480};
  const auto r = t2(in);
  expect_t2_core(r, 6144, {4, 4, 8, 0, 4, 0, 0, 4, 0, 0, 0, 0, 4, 0, 0, 0, 4, 0}, fraction(6144), 49160);
  EXPECT_EQ(r.paths, repeat({{0, 4, 7}, {1, 16, 12}, {2}}, {4, 4, 8}));
  EXPECT_EQ(r.plan.headroom_penalty, c128(38654705652));
  EXPECT_EQ(r.plan.total_planned_bandwidth, 49152);
  EXPECT_EQ(r.plan.total_tunnel_cost, c128(40));
  EXPECT_EQ(r.plan.planned_bandwidth_per_tunnel, 3072);
}
TEST(ContractExamples, C_G10_1) {
  T2In in {kG10, v({1}), v({4}), 16, 49152};
  in.prefer_headroom = true;
  in.cap = Amounts {28672, 28672, 0, 40960, 28672, 16384, 40960, 28672, 20480, 40960, 40960, 24576, 28672, 40960, 16384, 40960, 28672, 20480};
  const auto r = t2(in);
  expect_t2_core(r, 4096, {0, 6, 0, 6, 5, 5, 0, 5, 0, 0, 5, 0, 0, 0, 0, 0, 6, 0}, fraction(4096), 52431);
  EXPECT_EQ(r.paths, repeat({{3, 1, 16}, {4, 7, 10}, {5}}, {6, 5, 5}));
  EXPECT_EQ(r.plan.total_planned_bandwidth, 49152);
  EXPECT_EQ(r.plan.total_tunnel_cost, c128(43));
  EXPECT_EQ(r.plan.planned_bandwidth_per_tunnel, 3072);
}
TEST(ContractExamples, C_G10_2) {
  T2In in {kG10, v({2}), v({5}), 16, 49152};
  in.prefer_headroom = true;
  in.cap = Amounts {28672, 10240, 0, 22528, 13312, 1024, 40960, 13312, 20480, 40960, 25600, 24576, 28672, 40960, 16384, 40960, 10240, 20480};
  const auto r = t2(in);
  expect_t2_core(r, 2560, {0, 4, 0, 4, 0, 0, 4, 4, 8, 0, 4, 0, 0, 4, 0, 0, 0, 0}, fraction(2560), 40960);
  EXPECT_EQ(r.paths, repeat({{6, 3, 1}, {7, 10, 13}, {8}}, {4, 4, 8}));
  EXPECT_EQ(r.plan.total_planned_bandwidth, 40960);
  EXPECT_EQ(r.plan.total_tunnel_cost, c128(40));
  EXPECT_EQ(r.plan.planned_bandwidth_per_tunnel, 2560);
}
// The plan's alg-internal supply budget (21) is no public field; every public one is checked.
TEST(ContractExamples, T2_dead_head) {
  T2In in {{{0, 2, 0, 1}, {1, 2, 10, 1}}, v({0, 1}), v({2}), 2, 10};
  in.prefer_headroom = true;
  in.lossy_offer = 10;
  const auto r = t2(in);
  EXPECT_EQ(r.entry_ids, (std::vector<std::int32_t> {1, 3}));
  EXPECT_EQ(r.plan.design_by_head, (std::vector<std::pair<NodeId, Amount>> {{0, 5}, {1, 5}}));
  EXPECT_EQ(r.plan.heads_with_capacity, std::optional<std::vector<NodeId>>(v({1})));
  expect_t2_core(r, 5, {0, 2}, fraction(5), 0);
  EXPECT_EQ(r.plan.planned_bandwidth_per_tunnel, 3);
  EXPECT_EQ(r.paths, (Routes {{1}, {1}}));
  EXPECT_EQ(r.plan.total_planned_bandwidth, 6);
  EXPECT_EQ(r.plan.total_tunnel_cost, c128(2));
  EXPECT_EQ(r.plan.headroom_penalty, c128(3865470565));
  EXPECT_EQ(r.lossy_delivery, 5);
  EXPECT_EQ(r.lossy_source_drops, 5);
}
TEST(ContractExamples, T2_two_heads_equal) {
  T2In in {{{0, 2, 2, 1}, {1, 2, 10, 1}}, v({0, 1}), v({2}), 2};
  in.prefer_headroom = true;
  const auto r = t2(in);
  EXPECT_EQ(r.entry_ids, (std::vector<std::int32_t> {0, 1, 2, 3}));
  EXPECT_TRUE(r.plan.design_by_head.empty());
  expect_t2_core(r, 1, {2, 2}, fraction(1), 4);
  EXPECT_EQ(r.plan.planned_bandwidth_per_tunnel, 1);
  EXPECT_EQ(r.paths, (Routes {{0}, {1}, {0}, {1}}));
  EXPECT_EQ(r.plan.total_planned_bandwidth, 4);
  EXPECT_EQ(r.plan.total_tunnel_cost, c128(4));
  EXPECT_EQ(r.plan.headroom_penalty, c128(7730941132));
}
TEST(ContractExamples, T2_idle) {
  T2In in {{{0, 1, 10, 1}}, v({0}), v({1}), 4, 3};
  in.prefer_headroom = true;
  in.lossy_offer = 3;
  const auto r = t2(in);
  expect_t2_core(r, 2, {4}, fraction(5, 2), 10);
  EXPECT_EQ(r.plan.planned_bandwidth_per_tunnel, 1);
  EXPECT_EQ(r.plan.total_planned_bandwidth, 4);
  EXPECT_EQ(r.plan.headroom_penalty, c128(4294967294));
  EXPECT_EQ(r.lossy_delivery, 3);
}
TEST(ContractExamples, T2_zero_design) {
  T2In in {{{0, 1, 10, 1}}, v({0}), v({1}), 2, 0};
  in.prefer_headroom = true;
  in.lossy_offer = 4;
  const auto r = t2(in);
  expect_t2_core(r, 5, {2}, fraction(5), 10);
  EXPECT_EQ(r.plan.planned_bandwidth_per_tunnel, 1);
  EXPECT_EQ(r.plan.total_planned_bandwidth, 2);
  EXPECT_EQ(r.lossy_delivery, 4);
}
TEST(ContractExamples, T2_remainder_first) {
  T2In in {{{0, 1, 3, 1}, {0, 1, 2, 1}}, v({0}), v({1}), 2, 5};
  in.prefer_headroom = true;
  const auto r = t2(in);
  expect_t2_core(r, 2, {1, 1}, fraction(2), 5);
  EXPECT_EQ(r.plan.planned_bandwidth_per_tunnel, 2);
  EXPECT_EQ(r.paths, (Routes {{0}, {1}}));
  EXPECT_EQ(r.plan.total_planned_bandwidth, 4);
}
TEST(ContractExamples, T2_remainder_second) {
  T2In in {{{0, 1, 2, 1}, {0, 1, 3, 1}}, v({0}), v({1}), 2, 5};
  in.prefer_headroom = true;
  const auto r = t2(in);
  expect_t2_core(r, 2, {1, 1}, fraction(2), 4);
  EXPECT_EQ(r.plan.planned_bandwidth_per_tunnel, 2);
  EXPECT_EQ(r.paths, (Routes {{0}, {1}}));
  EXPECT_EQ(r.plan.total_planned_bandwidth, 4);
}
TEST(ContractExamples, T2_tiny_head) {
  T2In in {{{0, 2, 1, 1}, {1, 2, 100, 1}}, v({0, 1}), v({2}), 16};
  in.prefer_headroom = true;
  const auto r = t2(in);
  EXPECT_EQ(r.plan.heads_with_capacity, std::optional<std::vector<NodeId>>(v({0, 1})));
  expect_t2_core(r, 0, {0, 0}, fraction(1, 16), 0);
  EXPECT_TRUE(r.paths.empty());
  EXPECT_EQ(r.plan.insufficient_resolution, std::optional<bool>(true));
}
TEST(ContractExamples, T2_prefer_headroom_off) {
  T2In in {{{0, 1, 10, 1}, {0, 1, 10, 1}}, v({0}), v({1}), 2, 2};
  in.lossy_offer = 12;
  const auto r = t2(in);
  EXPECT_EQ(r.plan.max_bandwidth_per_tunnel, 10);
  EXPECT_EQ(r.plan.planned_bandwidth_per_tunnel, 1);
  EXPECT_EQ(r.counts, (Ints {2, 0}));
  EXPECT_EQ(r.plan.max_lossless_offer, 10);
  EXPECT_EQ(r.lossy_delivery, 10);
}
TEST(ContractExamples, T2_prefer_headroom_on) {
  T2In in {{{0, 1, 10, 1}, {0, 1, 10, 1}}, v({0}), v({1}), 2, 2};
  in.prefer_headroom = true;
  in.lossy_offer = 12;
  const auto r = t2(in);
  EXPECT_EQ(r.plan.max_bandwidth_per_tunnel, 10);
  EXPECT_EQ(r.plan.planned_bandwidth_per_tunnel, 1);
  EXPECT_EQ(r.counts, (Ints {1, 1}));
  EXPECT_EQ(r.plan.headroom_penalty, c128(858993458));
  EXPECT_EQ(r.plan.max_lossless_offer, 20);
  EXPECT_EQ(r.lossy_delivery, 12);
}
TEST(ContractExamples, T2_support) {
  T2In in {kSup, v({0}), v({1}), 2};
  in.support = Ints {0, 1};
  const auto r = t2(in);
  EXPECT_EQ(r.plan.max_bandwidth_per_tunnel, 6);
  EXPECT_EQ(r.counts, (Ints {1, 1, 0}));
  EXPECT_EQ(r.plan.max_lossless_offer, 12);
}
TEST(ContractExamples, T2_support_masked) {
  T2In in {kSup, v({0}), v({1}), 2};
  in.support = Ints {0, 1};
  in.excluded = {1};
  const auto r = t2(in);
  EXPECT_EQ(r.plan.max_bandwidth_per_tunnel, 3);
  EXPECT_EQ(r.counts, (Ints {2, 0, 0}));
  EXPECT_EQ(r.paths, (Routes {{0}, {0}}));
  EXPECT_EQ(r.entry_ids, (std::vector<std::int32_t> {0, 1}));
  EXPECT_EQ(r.plan.max_lossless_offer, 6);
}
TEST(ContractExamples, T2_support_mask_outside) {
  T2In in {kSup, v({0}), v({1}), 2};
  in.support = Ints {0, 1};
  in.excluded = {2};
  const auto r = t2(in);
  EXPECT_EQ(r.plan.max_bandwidth_per_tunnel, 6);
  EXPECT_EQ(r.counts, (Ints {1, 1, 0}));
  EXPECT_EQ(r.paths, (Routes {{0}, {1}}));
  EXPECT_EQ(r.entry_ids, (std::vector<std::int32_t> {0, 1}));
  EXPECT_EQ(r.plan.max_lossless_offer, 12);
}

// "Up to L": the chosen count, the delivery of the design (the bound without one) and the paths.
// Each count's bound is the lossless maximum of that count's program (every head live here).
std::vector<Amount> upto_bounds(const T2In& base, std::int64_t limit) {
  std::vector<Amount> out;
  for (std::int64_t k = 1; k <= limit; ++k) {
    T2In in = base;
    in.per_head = k;
    in.up_to = false;
    T2In unlimited = in;
    unlimited.design_demand = std::nullopt;
    out.push_back(*t2(unlimited).plan.max_lossless_offer);
  }
  return out;
}
TEST(ContractExamples, C_G1_UP_TO) {
  T2In in {kG1, v({0}), v({3}), 6, 49152};
  in.prefer_headroom = true;
  in.up_to = true;
  const auto r = t2(in);
  EXPECT_EQ(r.plan.tunnels_per_head, 6);
  EXPECT_EQ(r.delivered, 49152);
  EXPECT_EQ(upto_bounds(in, 6), (Amounts {40960, 40960, 40960, 40960, 40960, 49152}));  // = the scores (design above)
}
TEST(ContractExamples, T2_UP_TO_common) {
  T2In in {{{0, 3, 6, 1}, {1, 3, 3, 1}, {1, 3, 3, 1}}, v({0, 1}), v({3}), 2, 12};
  in.prefer_headroom = true;
  in.up_to = true;
  const auto r = t2(in);
  EXPECT_EQ(r.plan.tunnels_per_head, 2);
  EXPECT_EQ(r.delivered, 12);
  EXPECT_EQ(r.paths, (Routes {{0}, {1}, {0}, {2}}));
  EXPECT_EQ(upto_bounds(in, 2), (Amounts {7, 12}));  // = the scores (design 12)
}
// With a dead head the count-selection scores (its supply budget 21 and live share 5) are internal;
// the public outcome: one tunnel per head, delivery 0 (the dead head's share blocks), path [1].
TEST(ContractExamples, T2_UP_TO_dead_head) {
  T2In in {{{0, 2, 0, 1}, {1, 2, 10, 1}}, v({0, 1}), v({2}), 2, 10};
  in.prefer_headroom = true;
  in.up_to = true;
  const auto r = t2(in);
  EXPECT_EQ(r.plan.tunnels_per_head, 1);
  EXPECT_EQ(r.delivered, 0);
  EXPECT_EQ(r.paths, (Routes {{1}}));
}
TEST(ContractExamples, T2_tiny_head_UP_TO) {
  T2In in {{{0, 2, 1, 1}, {1, 2, 100, 1}}, v({0, 1}), v({2}), 3};
  in.prefer_headroom = true;
  in.up_to = true;
  const auto r = t2(in);
  EXPECT_EQ(r.plan.tunnels_per_head, 1);
  EXPECT_EQ(r.delivered, 2);
  EXPECT_EQ(r.paths, (Routes {{0}, {1}}));
  EXPECT_EQ(upto_bounds(in, 3), (Amounts {2, 0, 0}));  // = the scores (no design)
}

// ================================================================ programs

TEST(ContractExamples, E_duplicates) { EXPECT_EQ(bound({4, 8}, {{0}, {1}, {1}}, {1, 1, 1}), 12); }
TEST(ContractExamples, E_deduplicated) { EXPECT_EQ(bound({4, 8}, {{0}, {1}}, {1, 1}), 8); }
TEST(ContractExamples, E_static_path_fill) { EXPECT_EQ(bound({6, 4}, {{0}, {1}}, {1, 1}), 9); }
TEST(ContractExamples, E_zero_first) { EXPECT_EQ(bound({0, 4}, {{0}, {1}}, {1, 1}), 0); }
TEST(ContractExamples, E_zero_last) { EXPECT_EQ(bound({4, 0}, {{0}, {1}}, {1, 1}), 1); }
TEST(ContractExamples, E_support_pinned_bound) { EXPECT_EQ(bound({10, 2, 10, 10}, {{0, 1}, {2, 3}}, {1, 1}), 4); }

// One head, two links of capacity 1 to the destination, two equal split routes over both: separate
// contexts each put their share on the first link; pooled (from_dag) they split once.
std::pair<Amount, Amounts> split_contexts(bool pool) {
  const Net n({{0, 1, 1, 1}, {0, 1, 1, 1}});
  const NodeId t[1] = {1};
  std::shared_ptr<const ForwardingProgram> prog;
  if (pool) {
    prog = ForwardingProgram::from_dag(*n.g, v({0}), t, n.dag({0, 1}));
  } else {
    const std::vector<RouteDefinition> defs = {RouteDefinition {0, n.dag({0, 1}), 1}, RouteDefinition {0, n.dag({0, 1}), 1}};
    prog = ForwardingProgram::from_routes(*n.g, v({0}), t, defs);
  }
  auto ll = policy(n, v({0}), v({1}), {.routing = Pinned {prog, nullptr}});
  auto ly = policy(n, v({0}), v({1}), {.routing = Pinned {prog, nullptr}, .admission = LY});
  const auto r = ly->place(2, kFlows);
  Amounts loads = n.in_order(*r.edge_flows);
  const Amounts drops = n.per_edge(r.call_edge_drops);
  for (std::size_t e = 0; e < loads.size(); ++e) loads[e] += drops[e];
  return {ll->place(std::nullopt).call_delivered, loads};
}
TEST(ContractExamples, E_split_contexts) { EXPECT_EQ(split_contexts(false), (std::pair<Amount, Amounts> {1, {2, 0}})); }
TEST(ContractExamples, E_split_contexts_pooled) { EXPECT_EQ(split_contexts(true), (std::pair<Amount, Amounts> {2, {1, 1}})); }

std::pair<Amount, Amount> dag_contexts(bool pooled, Amount offer) {
  const Net n(kDP);
  const NodeId t[1] = {3};
  std::shared_ptr<const ForwardingProgram> prog;
  if (pooled) {
    prog = ForwardingProgram::from_dag(*n.g, v({0, 1}), t, n.dag({0, 1, 2, 3}));
  } else {
    const std::vector<RouteDefinition> defs = {RouteDefinition {0, n.dag({0, 1, 2, 3}), 1}, RouteDefinition {1, n.dag({0, 1, 2, 3}), 1}};
    prog = ForwardingProgram::from_routes(*n.g, v({0, 1}), t, defs);
  }
  auto ll = policy(n, v({0, 1}), v({3}), {.routing = Pinned {prog, nullptr}});
  auto ly = policy(n, v({0, 1}), v({3}), {.routing = Pinned {prog, nullptr}, .admission = LY});
  return {ll->place(std::nullopt).call_delivered, ly->place(offer).call_delivered};
}
TEST(ContractExamples, DAG_pooled) { EXPECT_EQ(dag_contexts(true, 2), (std::pair<Amount, Amount> {2, 2})); }
TEST(ContractExamples, DAG_separate) { EXPECT_EQ(dag_contexts(false, 2), (std::pair<Amount, Amount> {1, 1})); }

// ================================================================ stale forwarding and repair

const std::vector<InEdge> kDiamond100 = {{0, 1, 100, 1}, {0, 2, 100, 1}, {1, 3, 100, 1}, {2, 3, 100, 1}};
std::vector<std::optional<Amount>> ones(std::size_t k) { return std::vector<std::optional<Amount>>(k, Amount {1}); }

void expect_stale(const StaleOut& r, Amount delivered, const Amounts& flows, const Amounts& drops, std::optional<Amount> source = std::nullopt) {
  EXPECT_EQ(r.delivered, delivered);
  EXPECT_EQ(r.edge_flows, flows);
  EXPECT_EQ(r.edge_drops, drops);
  if (source.has_value()) EXPECT_EQ(r.source_drops, *source);
}
TEST(ContractExamples, STALE_remote_tail) {
  StaleIn in {kDiamond100, v({0}), 3, 100, ones(4)};
  in.black_links = {3};
  expect_stale(stale_dag(in), 50, {50, 50, 50, 0}, {0, 0, 0, 50});
}
TEST(ContractExamples, STALE_remote_tail_repair) {
  StaleIn in {kDiamond100, v({0}), 3, 100, ones(4)};
  in.mask_links = {3};
  expect_stale(stale_dag(in), 100, {100, 0, 100, 0}, {0, 0, 0, 0});
}
TEST(ContractExamples, STALE_wcmp_member) {
  StaleIn in {{{0, 1, 10, 1}, {0, 1, 40, 1}}, v({0}), 1, 40, {Amount {10}, Amount {40}}};
  in.black_links = {0};
  expect_stale(stale_dag(in), 32, {0, 32}, {8, 0});
}
TEST(ContractExamples, STALE_wcmp_zero_capacity) {
  const StaleIn in {{{0, 1, 0, 1}, {0, 1, 40, 1}}, v({0}), 1, 40, {Amount {0}, Amount {40}}};
  expect_stale(stale_dag(in), 40, {0, 40}, {0, 0});
}
TEST(ContractExamples, STALE_pooled_original) {
  const StaleIn in {{{0, 2, 10, 1}, {1, 2, 10, 1}, {2, 3, 10, 1}, {2, 3, 10, 1}}, v({0, 1}), 3, 2, ones(4)};
  expect_stale(stale_dag(in), 2, {1, 1, 1, 1}, {0, 0, 0, 0}, 0);
}
TEST(ContractExamples, STALE_pooled_remainder) {
  StaleIn in {{{0, 2, 10, 1}, {1, 2, 10, 1}, {2, 3, 10, 1}, {2, 3, 10, 1}}, v({0, 1}), 3, 2, ones(4)};
  in.black_links = {0, 3};
  expect_stale(stale_dag(in), 1, {0, 1, 1, 0}, {1, 0, 0, 0}, 0);
}
TEST(ContractExamples, STALE_all_zero_weight) {
  const StaleIn in {{{0, 1, 0, 1}}, v({0}), 1, 2, {Amount {0}}};
  expect_stale(stale_dag(in), 0, {0}, {0}, 2);
}
TEST(ContractExamples, STALE_black_source) {
  StaleIn in {{{0, 1, 10, 1}}, v({0}), 1, 10, ones(1)};
  in.black_nodes = {0};
  expect_stale(stale_dag(in), 0, {0}, {0}, 10);
}

TEST(ContractExamples, STALE_relief) {
  StalePathsIn in {{{0, 2, 10, 1}, {1, 2, 10, 1}, {2, 3, 10, 1}}, {{0, {0, 2}}, {1, {1, 2}}}, {{0, 10}, {1, 10}}};
  in.black_links = {0};
  const auto r = stale_paths(in);
  EXPECT_EQ(r.status, (Texts {"BLACKHOLED", "UP"}));
  EXPECT_EQ(r.delivered, (Amounts {0, 10}));
  EXPECT_EQ(r.edge_drops, (Amounts {10, 0, 0}));
}
TEST(ContractExamples, STALE_relief_original) {
  const StalePathsIn in {{{0, 2, 10, 1}, {1, 2, 10, 1}, {2, 3, 10, 1}}, {{0, {0, 2}}, {1, {1, 2}}}, {{0, 10}, {1, 10}}};
  const auto r = stale_paths(in);
  EXPECT_EQ(r.status, (Texts {"UP", "UP"}));
  EXPECT_EQ(r.delivered, (Amounts {5, 5}));
  EXPECT_EQ(r.edge_drops, (Amounts {0, 0, 10}));
}
const std::vector<InEdge> kDiamond10 = {{0, 1, 10, 1}, {0, 2, 10, 1}, {1, 3, 10, 1}, {2, 3, 10, 1}};
TEST(ContractExamples, STALE_tunnel_node) {
  StalePathsIn in {kDiamond10, {{0, {0, 2}}, {0, {1, 3}}}, {{0, 10}}};
  in.black_nodes = {2};
  const auto r = stale_paths(in);
  EXPECT_EQ(r.status, (Texts {"UP", "BLACKHOLED"}));
  EXPECT_EQ(r.assigned, (Amounts {5, 5}));
  EXPECT_EQ(r.delivered, (Amounts {5, 0}));
  EXPECT_EQ(r.edge_drops, (Amounts {0, 5, 0, 0}));
  EXPECT_EQ(r.lossless_admitted, std::optional<Amount>(1));
}
TEST(ContractExamples, STALE_tunnel_node_repair) {
  StalePathsIn in {kDiamond10, {{0, {0, 2}}, {0, {1, 3}}}, {{0, 10}}};
  in.mask_nodes = {2};
  const auto r = stale_paths(in);
  EXPECT_EQ(r.status, (Texts {"UP", "FAILED"}));
  EXPECT_EQ(r.assigned, (Amounts {10, 0}));
  EXPECT_EQ(r.delivered, (Amounts {10, 0}));
  EXPECT_EQ(r.edge_drops, (Amounts {0, 0, 0, 0}));
  EXPECT_EQ(r.lossless_admitted, std::optional<Amount>(10));
}
TEST(ContractExamples, STALE_black_source_tunnel) {
  StalePathsIn in {{{0, 1, 10, 1}}, {{0, {0}}}, {{0, 10}}};
  in.black_nodes = {0};
  const auto r = stale_paths(in);
  EXPECT_EQ(r.status, (Texts {"BLACKHOLED"}));
  EXPECT_EQ(r.assigned, (Amounts {0}));
  EXPECT_EQ(r.delivered, (Amounts {0}));
  EXPECT_EQ(r.source_drops, 10);
  EXPECT_EQ(r.edge_drops, (Amounts {0}));
}
TEST(ContractExamples, STALE_mask_and_black_hole) {
  StalePathsIn in {{{0, 1, 10, 1}, {0, 2, 10, 1}, {0, 3, 10, 1}, {1, 4, 10, 1}, {2, 4, 10, 1}, {3, 4, 10, 1}},
                   {{0, {0, 3}}, {0, {1, 4}}, {0, {2, 5}}},
                   {{0, 10}}};
  in.mask_links = {3};
  in.black_links = {3, 4};
  const auto r = stale_paths(in);
  EXPECT_EQ(r.status, (Texts {"FAILED", "BLACKHOLED", "UP"}));
  EXPECT_EQ(r.assigned, (Amounts {0, 5, 5}));
  EXPECT_EQ(r.delivered, (Amounts {0, 0, 5}));
  EXPECT_EQ(r.edge_drops, (Amounts {0, 0, 0, 0, 5, 0}));
  EXPECT_EQ(r.source_drops, 0);
}

// ================================================================ the Lossy loss rule

void expect_lossy(const LossyOut& r, const Routes& carried, const std::optional<Routes>& drops, const Amounts& delivered,
                  const Amounts& flows) {
  EXPECT_EQ(r.error, "");
  EXPECT_EQ(r.carried, carried);
  if (drops.has_value()) EXPECT_EQ(r.edge_drops, *drops);
  EXPECT_EQ(r.delivered, delivered);
  EXPECT_EQ(r.edge_flows, flows);
}
TEST(ContractExamples, D_merge) { expect_lossy(lossy({8, 4}, {{0}, {0}}, {9}), {{6}, {3}}, Routes {{2}, {1}}, {6, 3}, {9}); }
TEST(ContractExamples, D_key_tie) { expect_lossy(lossy({1, 1}, {{0}, {0}}, {1}), {{1}, {0}}, Routes {{0}, {1}}, {1, 0}, {1}); }
TEST(ContractExamples, D_rounding) {
  expect_lossy(lossy({5, 1, 1}, {{0}, {0}, {0}}, {3}), {{3}, {0}, {0}}, Routes {{2}, {1}, {1}}, {3, 0, 0}, {3});
}
// Demand 0 (routes 0 and 1) and demand 1 (route 2): the demand level is the head level here (two
// heads offered 2 each), the route level below it.
TEST(ContractExamples, D_hierarchy) {
  expect_lossy(lossy({1, 1, 2}, {{0}, {0}, {0}}, {1}, {0, 0, 1}), {{1}, {0}, {0}}, Routes {{0}, {1}, {2}}, {1, 0, 0}, {1});
}
TEST(ContractExamples, D_cycle_compressed) {
  expect_lossy(lossy({1, 1}, {{0, 1}, {1, 0}}, {2, 2}), {{1, 1}, {1, 1}}, Routes {{0, 0}, {0, 0}}, {1, 1}, {2, 2});
}
TEST(ContractExamples, D_zero_service) { expect_lossy(lossy({2}, {{0}}, {0}), {{0}}, Routes {{2}}, {0}, {0}); }
TEST(ContractExamples, D_subdivision_before) {
  expect_lossy(lossy({1, 1, 1}, {{0, 1}, {0, 2}, {0, 3}}, {1, 1, 1, 1}), {{1, 1}, {0, 0}, {0, 0}}, std::nullopt, {1, 0, 0}, {1, 1, 0, 0});
}
TEST(ContractExamples, D_subdivision_after) {
  expect_lossy(lossy({1, 1, 1}, {{0, 1, 2}, {0, 3}, {0, 4}}, {1, 1, 1, 1, 1}), {{1, 1, 1}, {0, 0}, {0, 0}}, std::nullopt, {1, 0, 0},
               {1, 1, 1, 0, 0});
}
TEST(ContractExamples, E_support_pinned_loss) {
  expect_lossy(lossy({6, 6}, {{0, 1}, {2, 3}}, {10, 2, 10, 10}), {{6, 2}, {6, 6}}, Routes {{0, 4}, {0, 0}}, {2, 6}, {6, 2, 6, 6});
}
// The forwarding's alternating state is an error (the round count it took is internal).
TEST(ContractExamples, D_C2) {
  EXPECT_EQ(lossy({1, 2}, {{0, 1}, {1, 0}}, {1, 2}).error, "lossy forwarding iteration entered a period-2 cycle");
}

// ================================================================ ledger, epochs, increments

// A T1 tunnel reserving 6 of a link of 10 and carrying 4: the ledger's four views; a replace from 4
// to 2 has delivery delta -2.
TEST(ContractExamples, E_ledger) {
  const Net n({{0, 1, 10, 1}});
  auto p = policy(n, v({0}), v({1}), {.routing = dte_with(1), .admission = LY});
  (void)p->replace(4, finite(6));
  EXPECT_EQ(p.ledger->carried(), (Amounts {4}));
  EXPECT_EQ(p.ledger->residual(), (Amounts {6}));
  EXPECT_EQ(p.ledger->reserved(), (Amounts {6}));
  EXPECT_EQ(p.ledger->unreserved(), (Amounts {4}));
  EXPECT_EQ(p->replace(2, finite(6)).delivery_delta, -2);
}

void expect_epoch(const EpochOut& r, Amount requested, Amount delivered, Amount drops, const Amounts& flows, const Amounts& residual,
                  const Amounts& per_call) {
  EXPECT_EQ(r.requested, requested);
  EXPECT_EQ(r.delivered, delivered);
  EXPECT_EQ(r.edge_drops, drops);
  EXPECT_EQ(r.edge_flows, flows);
  EXPECT_EQ(r.residual, residual);
  EXPECT_EQ(r.per_call, per_call);
}
TEST(ContractExamples, RC02_increment) { expect_epoch(epoch({1, 1}, {1, 1}, false, LY), 2, 1, 1, {1, 0}, {0, 1}, {1, 0}); }
TEST(ContractExamples, RC02_increment_lossless) { expect_epoch(epoch({1, 1}, {1, 1}, false, LL), 2, 2, 0, {1, 1}, {0, 0}, {1, 1}); }
TEST(ContractExamples, RC02_full) { expect_epoch(epoch({1, 1}, {1, 1}, true, LY), 2, 2, 0, {1, 1}, {0, 0}, {2}); }

TEST(ContractExamples, INC_lossless) {
  const auto r = increments({1, 3}, {3, 1}, LL);
  EXPECT_EQ(r.call_delivered, (Amounts {2, 0}));
  EXPECT_EQ(r.requested, 4);
  EXPECT_EQ(r.admitted, std::optional<Amount>(2));
  EXPECT_EQ(r.unadmitted, std::optional<Amount>(2));
  EXPECT_EQ(r.offered, std::nullopt);
  EXPECT_EQ(r.delivered, 2);
}
TEST(ContractExamples, INC_lossy) {
  const auto r = increments({1, 3}, {3, 1}, LY);
  EXPECT_EQ(r.call_delivered, (Amounts {2, 1}));
  EXPECT_EQ(r.requested, 4);
  EXPECT_EQ(r.offered, std::optional<Amount>(4));
  EXPECT_EQ(r.edge_drops, (Amounts {1, 0}));
  EXPECT_EQ(r.unadmitted, std::optional<Amount>(0));
  EXPECT_EQ(r.delivered, 3);
}
TEST(ContractExamples, INC_lossless_refusal) {
  const auto r = increments({10}, {5, 1}, LL, {{8}, {5}});
  EXPECT_EQ(r.call_delivered, (Amounts {2, 1}));
  EXPECT_EQ(r.requested, 6);
  EXPECT_EQ(r.admitted, std::optional<Amount>(3));
  EXPECT_EQ(r.unadmitted, std::optional<Amount>(3));
  EXPECT_EQ(r.delivered, 3);
}
// The delivery after each call: a refused unit is not offered again (2 then 3), not 2 then 5.
TEST(ContractExamples, R4_refusal) {
  const auto r = increments({10}, {5, 1}, LL, {{8}, {5}});
  const Amounts correct = {r.call_delivered[0], r.call_delivered[0] + r.call_delivered[1]};
  EXPECT_EQ(correct, (Amounts {2, 3}));
  EXPECT_NE(correct, (Amounts {2, 5}));
}

// place(0) on a new policy: nothing prepared, nothing changed; the residual is the service baseline
// less the other owners' loads (1 and 2).
TEST(ContractExamples, RC02_zero_unprepared) {
  const Net n({{0, 1, 5, 1}, {0, 1, 8, 1}});
  auto ledger = std::make_shared<FlowLedger>(n.g, n.by_id({5, 8}));
  const NodeId s[1] = {0};
  const NodeId t[1] = {1};
  const std::vector<RouteDefinition> defs = {RouteDefinition {0, n.path({0}), 1}, RouteDefinition {0, n.path({1}), 2}};
  auto other = policy(n, v({0}), v({1}), {.routing = Pinned {ForwardingProgram::from_routes(*n.g, s, t, defs), nullptr}, .admission = LY},
                      nullptr, ledger);
  (void)other->place(3);
  ASSERT_EQ(n.in_order(ledger->carried()), (Amounts {1, 2}));
  auto p = policy(n, v({0}), v({1}), {.routing = ShortestPaths {}}, nullptr, ledger);
  const auto r = p->place(0, OutputOptions {.edge_flows = true, .residual = true, .min_cut = true, .edge_reservations = true});
  EXPECT_EQ(r.call_delivered, 0);
  EXPECT_EQ(r.call_undelivered, std::optional<Amount>(0));
  EXPECT_EQ(r.delivery_delta, 0);
  EXPECT_TRUE(r.call_edge_drops.empty());
  EXPECT_EQ(r.call_source_drops, 0);
  EXPECT_EQ(r.termination, Termination::Complete);
  EXPECT_FALSE(r.routes_changed || r.weights_changed || r.reservations_changed);
  EXPECT_EQ(r.program, nullptr);
  EXPECT_EQ(n.in_order(*r.edge_flows), (Amounts {0, 0}));
  EXPECT_EQ(n.in_order(*r.edge_reservations), (Amounts {0, 0}));
  EXPECT_EQ(n.in_order(*r.residual), (Amounts {4, 6}));
  EXPECT_EQ(r.min_cut, std::nullopt);
  EXPECT_EQ(r.marginal_cost_curve, std::nullopt);
}

// The change flags of the second and third calls on a one-path program (path capacity 20): a
// controller re-solving changes its path's weight; a fixed program (Pinned) changes nothing.
TEST(ContractExamples, FLAGS_controller_increment) {
  const Net n({{0, 1, 20, 1}});
  auto p = policy(n, v({0}), v({1}), {.routing = centralized()});
  (void)p->place(5);
  const auto b = p->place(5);
  const auto c = p->replace(10);
  EXPECT_EQ((std::vector<bool> {b.routes_changed, c.routes_changed}), (std::vector<bool> {false, false}));
  EXPECT_EQ((std::vector<bool> {b.weights_changed, c.weights_changed}), (std::vector<bool> {true, false}));
  EXPECT_EQ((std::vector<bool> {b.reservations_changed, c.reservations_changed}), (std::vector<bool> {false, false}));
}
TEST(ContractExamples, FLAGS_fixed_offer) {
  const Net n({{0, 1, 20, 1}});
  const NodeId s[1] = {0};
  const NodeId t[1] = {1};
  const std::vector<RouteDefinition> defs = {RouteDefinition {0, n.path({0}), 1}};
  auto p = policy(n, v({0}), v({1}), {.routing = Pinned {ForwardingProgram::from_routes(*n.g, s, t, defs), nullptr}});
  (void)p->replace(5);
  for (int i = 0; i < 2; ++i) {
    const auto r = p->replace(10);
    EXPECT_FALSE(r.routes_changed || r.weights_changed || r.reservations_changed);
  }
}

// ================================================================ terminations

TEST(ContractExamples, TERM_T1_all_down) {
  const Net n(kX7);
  const auto r = termination(n, v({0}), v({3}), dte_with(2), 14, finite(14));
  EXPECT_EQ(r.tunnels, 2U);
  EXPECT_EQ(r.lossy, "COMPLETE");
  EXPECT_EQ(r.lossless, "EXHAUSTED");
  EXPECT_EQ(r.source_drops, (Amounts {14}));
  EXPECT_EQ(r.forward_error, "");
}
TEST(ContractExamples, TERM_T2_all_dead) {
  const Net n({{0, 2, 0, 1}, {1, 2, 0, 1}});
  const auto r = termination(n, v({0, 1}), v({2}), centralized(SourceAllocation::Equal, {}, TunnelCount {2}), 10, finite(10));
  EXPECT_EQ(r.tunnels, 0U);
  EXPECT_EQ(r.lossy, "COMPLETE");
  EXPECT_EQ(r.lossless, "EXHAUSTED");
  EXPECT_EQ(r.source_drops, (Amounts {5, 5}));
  EXPECT_EQ(r.forward_error, "");
}
// One head, one tunnel of weight 0: no head can forward.
TEST(ContractExamples, TERM_zero_weight_route) {
  const Net n({{0, 1, 10, 1}});
  const NodeId s[1] = {0};
  const NodeId t[1] = {1};
  const std::vector<RouteDefinition> defs = {RouteDefinition {0, n.path({0}), 0}};
  const auto r = termination(n, v({0}), v({1}), Pinned {ForwardingProgram::from_routes(*n.g, s, t, defs), nullptr}, 4, std::nullopt);
  EXPECT_EQ(r.tunnels, 1U);
  EXPECT_EQ(r.lossy, "NO_ACTIVE_WEIGHT");
  EXPECT_EQ(r.lossless, "NO_ACTIVE_WEIGHT");
  EXPECT_EQ(r.source_drops, (Amounts {4}));
}
TEST(ContractExamples, T2_no_tunnels_forward) {
  const Net n({{0, 2, 1, 1}, {1, 2, 100, 1}});
  const auto r = termination(n, v({0, 1}), v({2}), centralized(SourceAllocation::Equal, {}, TunnelCount {16}), 10, kUnlimited);
  EXPECT_EQ(r.tunnels, 0U);
  EXPECT_EQ(r.lossy, "COMPLETE");
  EXPECT_EQ(r.lossless, "EXHAUSTED");
  EXPECT_EQ(r.source_drops, (Amounts {5, 5}));
  EXPECT_EQ(r.forward_error, "");
}

// ================================================================ texts

TEST(ContractExamples, ERR_ownership_texts) {
  const Net n(kDP);
  const Net other(kDP);
  const NodeId t[1] = {3};
  const auto foreign_program = ForwardingProgram::from_dag(*other.g, v({0, 1}), t, other.dag({0, 1, 2, 3}));
  const auto foreign_mask = std::make_shared<const Eligibility>(other.g->eligibility({}, std::vector<EdgeId> {0}));
  ShortestPathCache foreign_cache(*other.g);
  FlowLedger fg(n.g);
  auto make = [&](FlowPolicyConfig c, const Eligibility* m = nullptr, ShortestPathCache* cache = nullptr) {
    return error_of([&] { FlowPolicy p(fg, v({0, 1}), v({3}), std::move(c), m, cache); });
  };
  const Texts got = {
      make({.routing = Pinned {foreign_program, nullptr}}),
      make({.routing = centralized(SourceAllocation::Optimized, other.dag({0, 2}))}),
      make({.routing = Pinned {ForwardingProgram::from_dag(*n.g, v({0, 1}), t, n.dag({0, 1, 2, 3})), foreign_mask}}),
      make({.routing = ShortestPaths {}}, foreign_mask.get()),
      make({.routing = ShortestPaths {}}, nullptr, &foreign_cache),
      error_of([&] { (void)ForwardingProgram::from_dag(*n.g, v({0, 1}), t, other.dag({0, 1, 2, 3})); }),
  };
  EXPECT_EQ(got, (Texts {"program belongs to another graph", "support belongs to another graph", "black_holes belongs to another graph",
                         "foreign Eligibility: built for another graph", "foreign ShortestPathCache: bound to another graph",
                         "dag belongs to another graph"}));
}

// The T1 search minimum and the ledger's reservation baseline (its length, then each entry).
TEST(ContractExamples, ERR_unreserved_texts) {
  const Net n({{0, 1, 10, 1}, {0, 1, 10, 1}});
  auto minimum = [&](Amount x) { return error_of([&] { (void)policy(n, v({0}), v({1}), {.routing = dte_with(1, false, std::nullopt, x)}); }); };
  auto baseline = [&](Amounts r) {
    const Amounts service = {10, 10};
    return error_of([&] { FlowLedger fg(n.g, service, r); });
  };
  EXPECT_EQ(minimum(Amount {1} << 62), "min_unreserved_bandwidth must be an integer in [0, 2^62)");
  EXPECT_EQ(baseline({1, 2, 3}), "reservation_baseline has length 3, expected 2 (one entry per edge)");
  EXPECT_EQ(minimum(-1), "min_unreserved_bandwidth must be an integer in [0, 2^62)");
  EXPECT_EQ(baseline(n.by_id({5, 11})), "reservation_baseline[1] = 11 is outside [0, edge_capacity_units[1] = 10]");
  EXPECT_EQ(minimum(0), "");
  EXPECT_EQ(baseline({10, 0}), "");
}

// ================================================================ the configuration

TEST(ContractExamples, CFG_defaults) {
  EXPECT_EQ(config_case(ShortestPaths {}),
            "allocation=EQUAL source_allocation=EQUAL weight_basis=None tie_break=FEWER_HOPS multipath=True parallel_edges=True");
  EXPECT_EQ(config_case(sp_with(RouteAllocation::Equal)),
            "allocation=EQUAL source_allocation=EQUAL weight_basis=None tie_break=FEWER_HOPS multipath=True parallel_edges=True");
  EXPECT_EQ(config_case(sp_with(RouteAllocation::Weighted)),
            "allocation=WEIGHTED source_allocation=EQUAL weight_basis=CAPACITY tie_break=FEWER_HOPS multipath=True parallel_edges=True");
  EXPECT_EQ(config_case(sp_with(RouteAllocation::Optimized)),
            "allocation=OPTIMIZED source_allocation=OPTIMIZED weight_basis=None tie_break=FEWER_HOPS multipath=True parallel_edges=True");
  EXPECT_EQ(config_case(centralized(SourceAllocation::Equal, {}, TunnelCount {16})),
            "source_allocation=EQUAL support=None tunnel_count=16 prefer_headroom=False");
}

// Every case is a construction TypeError of the Python surface (a missing or None source
// allocation, a keyword of another model). In C++ the same shapes cannot be written: CentralizedTE
// has no default constructor, and the models have only their own fields.
TEST(ContractExamples, CFG_construction) {
  static_assert(!std::is_default_constructible_v<CentralizedTE>);
  static_assert(std::is_constructible_v<CentralizedTE, SourceAllocation>);
  static_assert(std::is_aggregate_v<DistributedTE> && std::is_aggregate_v<ShortestPaths> && std::is_aggregate_v<FlowPolicyConfig>);
}

TEST(ContractExamples, CFG_matrix) {
  const Net& n = config_net();
  const Texts got = {
      config_case(sp_with(RouteAllocation::Equal, SourceAllocation::Optimized)),
      config_case(sp_with(RouteAllocation::Optimized), 1, LY),
      config_case(sp_with(RouteAllocation::Equal, std::nullopt, std::nullopt, RouteTieBreak::HigherCapacity)),
      config_case(sp_with(RouteAllocation::Equal, std::nullopt, std::nullopt, RouteTieBreak::HigherCapacity, false), 2),
      config_case(sp_with(RouteAllocation::Equal, std::nullopt, std::nullopt, RouteTieBreak::HigherCapacity, false)),
      config_case(sp_with(RouteAllocation::Equal, std::nullopt, std::nullopt, RouteTieBreak::HigherHeadroom, true, false), 2),
      config_case(sp_with(RouteAllocation::Equal, std::nullopt, WeightBasis::Capacity)),
      config_case(sp_with(RouteAllocation::Weighted, std::nullopt, WeightBasis::Explicit)),
      config_case(centralized(SourceAllocation::Optimized, {}, TunnelCount {2}), 2),
      config_case(centralized(SourceAllocation::Optimized, {}, TunnelCount {2})),
      config_case(centralized(SourceAllocation::Equal), 1, LY),
      config_case(centralized(SourceAllocation::Equal, {}, std::nullopt, true)),
      config_case(centralized(SourceAllocation::Equal, n.dag({1}), TunnelCount {2}), 1, LY),
      config_case(dte_with(4, true)),
  };
  EXPECT_EQ(got, (Texts {
                     "ShortestPaths with EQUAL or WEIGHTED requires source_allocation=EQUAL",
                     "ShortestPaths OPTIMIZED requires admission=LOSSLESS",
                     "ShortestPaths tie_break requires multipath=False or parallel_edges=False",
                     "ShortestPaths with multiple input sources requires tie_break=FEWER_HOPS when multipath=False",
                     "allocation=EQUAL source_allocation=EQUAL weight_basis=None tie_break=HIGHER_CAPACITY multipath=False parallel_edges=True",
                     "allocation=EQUAL source_allocation=EQUAL weight_basis=None tie_break=HIGHER_HEADROOM multipath=True parallel_edges=False",
                     "weight_basis requires allocation=WEIGHTED",
                     "weight_basis=EXPLICIT requires edge_weights",
                     "CentralizedTE with tunnel_count and several sources requires source_allocation=EQUAL",
                     "source_allocation=OPTIMIZED support=None tunnel_count=2 prefer_headroom=False",
                     "CentralizedTE without tunnel_count requires admission=LOSSLESS",
                     "prefer_headroom requires tunnel_count",
                     "source_allocation=EQUAL support=dag tunnel_count=2 prefer_headroom=False",
                     "DistributedTE requires an exact tunnel count",
                 }));
}

TEST(ContractExamples, CFG_first_error) {
  const Texts got = {
      config_case(centralized(SourceAllocation::Equal, {}, std::nullopt, true), 1, LY),
      config_case(sp_with(RouteAllocation::Equal, SourceAllocation::Optimized, WeightBasis::Capacity)),
      config_case(sp_with(RouteAllocation::Optimized, std::nullopt, std::nullopt, RouteTieBreak::HigherCapacity), 1, LY),
      config_case(sp_with(RouteAllocation::Equal, std::nullopt, WeightBasis::Capacity, RouteTieBreak::HigherCapacity)),
      config_case(sp_with(RouteAllocation::Equal, std::nullopt, WeightBasis::Capacity, RouteTieBreak::HigherCapacity, false), 2),
      config_case(centralized(SourceAllocation::Optimized, {}, TunnelCount {0}, true), 2),
      config_case(dte_with(0, true, CostStretch(1, 2))),
      config_case(dte_with(4, true, CostStretch(1, 2))),
      config_case(dte_with(4, false, CostStretch(1, 2), -1)),
      config_case(dte_with(4, false, std::nullopt, -1)),
  };
  EXPECT_EQ(got, (Texts {
                     "CentralizedTE without tunnel_count requires admission=LOSSLESS",
                     "ShortestPaths with EQUAL or WEIGHTED requires source_allocation=EQUAL",
                     "ShortestPaths OPTIMIZED requires admission=LOSSLESS",
                     "ShortestPaths tie_break requires multipath=False or parallel_edges=False",
                     "ShortestPaths with multiple input sources requires tie_break=FEWER_HOPS when multipath=False",
                     "CentralizedTE with tunnel_count and several sources requires source_allocation=EQUAL",
                     "tunnel count must be in [1, 10000]",
                     "DistributedTE requires an exact tunnel count",
                     "cost_ceiling CostStretch factor must be >= 1 (got 1/2)",
                     "min_unreserved_bandwidth must be an integer in [0, 2^62)",
                 }));
}

TEST(ContractExamples, CFG_cache_order) {
  const Net& n = config_net();
  const auto to2 = config_pinned(n, 2);
  const auto to1 = config_pinned(n, 1);
  const std::pair<std::vector<NodeId>, std::vector<NodeId>> ends {v({0}), v({1})};
  const Texts got = {
      config_case(dte_with(0), 1, LL, true),
      config_case(to2, 1, LL, true, ends),
      config_case(centralized(SourceAllocation::Equal, {}, TunnelCount {0}), 1, LL, true),
      config_case(centralized(SourceAllocation::Equal, {}, std::nullopt, true), 1, LL, true),
      config_case(sp_with(RouteAllocation::Weighted), 1, LL, true),
      config_case(dte_with(0)),
      config_case(to2, 1, LL, false, ends),
      config_case(to1, 1, LL, false, ends),
  };
  EXPECT_EQ(got, (Texts {
                     "DistributedTE does not use ShortestPathCache",
                     "Pinned does not use ShortestPathCache",
                     "CentralizedTE does not use ShortestPathCache",
                     "prefer_headroom requires tunnel_count",
                     "allocation=WEIGHTED source_allocation=EQUAL weight_basis=CAPACITY tie_break=FEWER_HOPS multipath=True parallel_edges=True",
                     "tunnel count must be in [1, 10000]",
                     "program endpoints do not match policy endpoints",
                     "valid",
                 }));
}

// A ceiling of another type (a float, an int) is a Python TypeError: the C++ ceiling is a variant of
// AbsoluteCost and CostStretch. The valid case: the fraction 3/2.
TEST(ContractExamples, CFG_ceiling_type) {
  static_assert(!std::is_constructible_v<CostCeiling, double>);
  EXPECT_EQ(config_case(dte_with(4, false, CostStretch(3, 2))), "valid");
}

TEST(ContractExamples, EVAL_first_error) {
  const Net& n = config_net();
  const DesignArg inf {std::in_place, std::nullopt};
  const Texts got = {
      evaluation_case(centralized(SourceAllocation::Equal, {}, std::nullopt, true), LL, finite(-1), std::nullopt),
      evaluation_case(dte_with(4), LY, finite(-1), inf),
      evaluation_case(dte_with(4), LY, finite(5), inf),
      evaluation_case(dte_with(4), LY, std::nullopt, std::nullopt),
      evaluation_case(dte_with(4), LY, std::nullopt, finite(10)),
      evaluation_case(ShortestPaths {}, LL, inf, finite(5)),
      evaluation_case(ShortestPaths {}, LL, std::nullopt, finite(5)),
      evaluation_case(centralized(SourceAllocation::Equal, {}, TunnelCount {16}), LY, std::nullopt, std::nullopt),
      evaluation_case(centralized(SourceAllocation::Equal, {}, TunnelCount {16}), LY, std::nullopt, finite(10)),
      evaluation_case(config_pinned(n, 2), LY, std::nullopt, std::nullopt),
      evaluation_case(ShortestPaths {}, LY, std::nullopt, std::nullopt),
      evaluation_case(ShortestPaths {}, LY, finite(4), std::nullopt),
  };
  EXPECT_EQ(got, (Texts {
                     "prefer_headroom requires tunnel_count",
                     "offer must be an integer in [0, 2^62)",
                     "DistributedTE requires a finite design demand",
                     "DistributedTE evaluation requires design_demand",
                     "DistributedTE LOSSY evaluation requires a finite offer",
                     "offer must be an integer in [0, 2^62)",
                     "ShortestPaths does not accept design_demand",
                     "CentralizedTE LOSSY evaluation requires a finite offer",
                     "offer=10 design_demand=10",
                     "Pinned LOSSY evaluation requires a finite offer",
                     "offer=first-hop experiment design_demand=None",
                     "offer=4 design_demand=None",
                 }));
}

// An explicit +inf offer is refused by the evaluator after the configuration. A string offer is a
// Python TypeError with no C++ form (the C++ offer is an amount or +inf).
TEST(ContractExamples, EVAL_offer_infinite) {
  const DesignArg inf {std::in_place, std::nullopt};
  EXPECT_EQ(evaluation_case(centralized(SourceAllocation::Equal, {}, std::nullopt, true), LL, inf, std::nullopt),
            "prefer_headroom requires tunnel_count");
  EXPECT_EQ(evaluation_case(ShortestPaths {}, LL, inf, std::nullopt), "offer must be an integer in [0, 2^62)");
}

// ECMP Lossy over parallel links from 0 to 1: the baseline and each single-link deletion at a fixed
// offer, or without one (the first-hop experiment: each first-hop link injects its own capacity).
// The baseline and each deletion's reduction, and the baseline's accounting: requested (-1 for an
// infinite request), offered and the drop on each link.
struct OfferOut {
  Amount baseline;
  Amounts reductions;
  Amount requested;
  Amount offered;
  Amounts edge_drops;
  bool operator==(const OfferOut&) const = default;
};
OfferOut evaluate_offer(const Amounts& caps, DesignArg offer) {
  std::vector<InEdge> rows;
  for (Amount c : caps) rows.push_back({0, 1, c, 1});
  const Net n(rows);
  EvaluationOptions o {{.routing = ShortestPaths {}, .admission = LY}};
  o.offer = offer;
  const auto s = evaluate_edge_failures(*n.g, v({0}), v({1}), o, 1);
  Amounts reductions(n.m(), 0);
  for (const auto& e : s.entries) reductions[static_cast<std::size_t>(n.pos(e.edge))] = e.delivery_reduction;
  return {s.baseline.call_delivered, reductions, lo(s.baseline.requested), lo(s.baseline.offered),
          n.per_edge(s.baseline.edge_drops)};
}
TEST(ContractExamples, EVAL_fixed_offer) {
  const auto r = evaluate_offer({10, 10}, finite(4));
  EXPECT_EQ(r.baseline, 4);
  EXPECT_EQ(r.reductions, (Amounts {0, 0}));
}
TEST(ContractExamples, EVAL_first_hop) {
  const auto r = evaluate_offer({10, 10}, std::nullopt);
  EXPECT_EQ(r.baseline, 20);
  EXPECT_EQ(r.reductions, (Amounts {10, 10}));
}
// Each link injects its own free capacity (here its capacity): 14 delivered, nothing dropped.
TEST(ContractExamples, EVAL_first_hop_unequal) {
  EXPECT_EQ(evaluate_offer({4, 10}, std::nullopt), (OfferOut {14, {4, 10}, -1, 14, {0, 0}}));
}
// A finite offer splits by round robin: 6 per link, 2 dropped on the link of 4.
TEST(ContractExamples, EVAL_fixed_offer_unequal) {
  EXPECT_EQ(evaluate_offer({4, 10}, finite(12)), (OfferOut {10, {0, 6}, 12, 12, {2, 0}}));
}

// ================================================================ CostStretch's exact fraction

std::string reduced(const CostStretch& f) { return std::to_string(f.numerator()) + "/" + std::to_string(f.denominator()); }

// The C++ value of each factor: a fraction reduced; an integer n/1; the float 1.2 is 6/5, 1.1 is
// 11/10, 1.5 is 3/2. The Python conversions (a float's shortest decimal, a non-finite float, an
// integer or a component beyond 64 bits, a bool or a string) have no C++ form; a denominator below 1
// is the C++ refusal.
TEST(ContractExamples, STRETCH_values) {
  EXPECT_EQ(reduced(CostStretch(6, 5)), "6/5");
  EXPECT_EQ(reduced(CostStretch(11, 10)), "11/10");
  EXPECT_EQ(reduced(CostStretch(3, 2)), "3/2");
  EXPECT_EQ(reduced(CostStretch(7500000000000001, 25000000000000000)), "7500000000000001/25000000000000000");
  EXPECT_EQ(reduced(CostStretch(12, 10)), "6/5");
  EXPECT_EQ(reduced(CostStretch(Amount {1} << 62, Amount {1} << 60)), "4/1");
  EXPECT_EQ(reduced(CostStretch(2)), "2/1");
  EXPECT_EQ(reduced(CostStretch(3)), "3/1");
  EXPECT_EQ(error_of([] { (void)CostStretch(1, 0); }), "CostStretch fraction denominator must be in [1, 2^63)");
  EXPECT_TRUE(CostStretch(15, 10) == CostStretch(3, 2));
  EXPECT_TRUE(CostStretch(6, 5) == CostStretch(12, 10));
  EXPECT_FALSE(CostStretch(6, 5) == CostStretch(5404319552844595, 4503599627370496));
}
// NumPy factors are Python values; their C++ values: float32/float16/float64 1.2 are 6/5, float32
// 0.30000001192092896 is 3/10, int64 2 is 2/1.
TEST(ContractExamples, STRETCH_numpy) {
  EXPECT_EQ(reduced(CostStretch(6, 5)), "6/5");
  EXPECT_EQ(reduced(CostStretch(3, 10)), "3/10");
  EXPECT_EQ(reduced(CostStretch(2)), "2/1");
  EXPECT_TRUE(CostStretch(6, 5) == CostStretch(6, 5));
}
