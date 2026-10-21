// The worked examples run on the library: each example's own inputs and expected values,
// one test per example (ContractExamples.<id> with '-' as '_'). The examples state edges as
// (source, destination, capacity, cost) in input order and amounts in units; every per-edge value
// below is in that order (translated through the graph's input positions). An example whose input
// is a Python value (a float factor, a NumPy scalar, a wrong type) is checked here through the C++
// value it stands for; the conversion itself is the Python suite's.

#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <iomanip>
#include <limits>
#include <random>
#include <cstdint>
#include <functional>
#include <map>
#include <memory>
#include <optional>
#include <set>
#include <string>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

#include "netgraph/core/detail/forwarding.hpp"
#include "netgraph/core/detail/concurrency_models.hpp"
#include "netgraph/core/detail/ledger_access.hpp"
#include "netgraph/core/detail/policy_engine.hpp"
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
bool less(const Frac& a, const Frac& b) { return static_cast<Int128>(a.num) * b.den < static_cast<Int128>(b.num) * a.den; }
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
// Each route a tunnel, or (where dags[i] is set) a supplied DAG route of the same links.
std::shared_ptr<const ForwardingProgram> routed_program(const Routed& r, const std::vector<NodeId>& head_of, const Amounts& weights,
                                                        const std::vector<bool>& dags = {}) {
  std::vector<RouteDefinition> defs;
  for (std::size_t i = 0; i < r.paths.size(); ++i) {
    if (i < dags.size() && dags[i]) {
      defs.push_back(RouteDefinition {head_of[i], dag_of(*r.net.g, r.paths[i]), weights[i]});
      continue;
    }
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
// delivers; per link, what it carries. `dark`: the abstract links the call black-holes; `dags[i]`:
// route i a supplied DAG route instead of a tunnel; `total`: the call's offer (default: the sum of
// the offers). The structure (hooks): the abstract links that wait, the cyclic groups, their rounds
// and the dependency arcs.
struct LossyOut {
  std::string error;
  std::int64_t rounds {-1};  // hooks only, as the three below
  std::int64_t groups {-1};
  Ints barrier_links;
  std::int64_t dependencies {-1};
  Routes carried;
  Routes edge_drops;
  Amounts delivered;
  Amounts edge_flows;
};
LossyOut lossy(const Amounts& offers, const Routes& routes, const Amounts& cap, std::vector<NodeId> head_of = {}, const Ints& dark = {},
               const std::vector<bool>& dags = {}, std::optional<Amount> total = std::nullopt) {
  if (head_of.empty()) head_of.assign(routes.size(), 0);
  const Routed r = routed(cap, routes, head_of);
  const auto prog = routed_program(r, head_of, offers, dags);
  auto p = policy(r.net, r.heads, {r.destination}, {.routing = Pinned {prog, dark.empty() ? nullptr : r.net.mask({}, dark)}, .admission = LY});
  if (!total.has_value()) {
    total = 0;
    for (Amount x : offers) *total += x;
  }
  LossyOut out;
  PlacementResult res;
  out.error = error_of([&] { res = p->place(*total, kFlows); });
#if NETGRAPH_CORE_PRIVATE_HOOKS
  const auto& counters = detail::thread_lossy_counters();
  out.rounds = counters.rounds;
  out.groups = counters.groups;
  out.dependencies = counters.dependencies;
  for (std::size_t e = 0; e < cap.size(); ++e) {
    if (std::find(counters.barriers.begin(), counters.barriers.end(), r.link[e]) != counters.barriers.end()) out.barrier_links.push_back(static_cast<std::int64_t>(e));
  }
#endif
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
    // by_tunnel has the tunnels only (a DAG route's delivery is not reported: -1).
    const auto tunnel = static_cast<std::size_t>(std::count(dags.begin(), dags.begin() + static_cast<std::ptrdiff_t>(std::min(i, dags.size())), true));
    out.delivered.push_back(i < dags.size() && dags[i] ? -1 : res.by_tunnel[i - tunnel].delivered);
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

// ================================================================ several demands at once (apply_actions)

namespace {

// The list form of the loss examples: each route is its own demand, a Pinned LOSSY policy (created
// in route order), all in one SIMULTANEOUS list. Each abstract link e gets its own two nodes; each
// route runs from its own head through its links to one destination, consecutive links joined by a
// link that never limits anything, except where e is fused to f (fused_links): there e's head is f's
// tail (no join), so the examples' serial pairs stay pairs. `dark[i]`: the abstract links route i's
// demand black-holes. `forms[i]` (default kTunnel) realizes route i as:
// - kTunnel: one tunnel;
// - kDagRoute: the same links as a supplied DAG route alone (a pooled context);
// - kDagBesideExit: that DAG route beside a tunnel on its head's own exit to the destination (a split
//   route of a program that is not pooled), offered twice the route's offer (each entry gets it);
// - kZeroBranch: a DAG route of the route's links plus that exit, beside the exit tunnel, offered 2:
//   each entry gets 1 and round robin sends the DAG route's unit to the exit (a lower EdgeId than
//   the route's first link), so every member on the route has a whole-unit share of 0 (the route's
//   offer is not used);
// - kEntryTunnel, kEntryDag: a tunnel or a DAG route in route 0's program, from route 0's head, of
//   weight its offer (0: an empty entry); no participant of its own, its black holes not used; route
//   0 (a tunnel) then has weight its offer and is offered `offer0` (default: the sum of its
//   program's weights): with every other route an entry, the list is one scalar call.
// Per participant (a route of its own) its rows on its route's links and its delivery (a split route
// beside an exit: without the exit's units); the structure (hooks): the abstract links that wait
// (after serial normalization), the cyclic groups, their rounds and the dependency arcs.
enum Form { kTunnel, kDagRoute, kDagBesideExit, kZeroBranch, kEntryTunnel, kEntryDag };
struct ListOut {
  std::string error;
  std::int64_t rounds {-1};  // hooks only, as the three below
  std::int64_t groups {-1};
  Ints barrier_links;
  std::int64_t dependencies {-1};
  Routes arrivals;
  Routes carried;
  Routes edge_drops;
  Amounts delivered;
  Amounts edge_flows;
  Texts status;
};
bool black_holed(const std::vector<Ints>& dark, std::size_t i, std::int64_t e) {
  return i < dark.size() && std::find(dark[i].begin(), dark[i].end(), e) != dark[i].end();
}
// Per link the link it is fused to (-1: none): every occurrence of e that can carry something
// continues on f and every such occurrence of f comes from e. A zero offer carries nothing, and a
// route's traffic stops at its own black hole (which still arrives from the link before it). A
// pair that would make a route visit a node twice (a dead occurrence on e going elsewhere and one
// on f) is joined instead: a tunnel is a simple path.
Ints fused_links(const Amounts& offers, const Routes& routes, std::size_t L, const std::vector<Ints>& dark) {
  std::vector<std::set<std::int64_t>> next(L);
  std::vector<std::set<std::int64_t>> prev(L);
  for (std::size_t i = 0; i < routes.size(); ++i) {
    if (offers[i] == 0) continue;
    const auto& r = routes[i];
    for (std::size_t j = 0; j < r.size(); ++j) {
      prev[static_cast<std::size_t>(r[j])].insert(j > 0 ? r[j - 1] : -1);
      if (black_holed(dark, i, r[j])) break;
      next[static_cast<std::size_t>(r[j])].insert(j + 1 < r.size() ? r[j + 1] : -1);
    }
  }
  Ints out(L, -1);
  for (std::size_t e = 0; e < L; ++e) {
    const auto& s = next[e];
    if (s.size() == 1 && *s.begin() >= 0 && prev[static_cast<std::size_t>(*s.begin())] == std::set<std::int64_t> {static_cast<std::int64_t>(e)}) out[e] = *s.begin();
  }
  for (bool changed = true; changed;) {
    changed = false;
    for (const auto& r : routes) {
      // Its nodes: a tail is e, a head L + e (a fused head is the next link's tail).
      std::vector<std::int64_t> nodes {r.front()};
      for (std::size_t j = 0; j < r.size(); ++j) {
        const auto f = out[static_cast<std::size_t>(r[j])];
        nodes.push_back(f >= 0 ? f : static_cast<std::int64_t>(L) + r[j]);
        if (j + 1 < r.size() && f != r[j + 1]) nodes.push_back(r[j + 1]);
      }
      if (std::set<std::int64_t>(nodes.begin(), nodes.end()).size() == nodes.size()) continue;
      for (std::size_t j = 0; j < r.size(); ++j) {
        auto& f = out[static_cast<std::size_t>(r[j])];
        if (f >= 0 && (j + 1 == r.size() || f != r[j + 1])) f = -1;  // e without its f
        for (auto& g : out) {
          if (g == r[j] && (j == 0 || out[static_cast<std::size_t>(r[j - 1])] != r[j])) g = -1;  // f without its e
        }
      }
      changed = true;
    }
  }
  return out;
}
ListOut list_lossy(const Amounts& offers, const Routes& routes, const Amounts& cap, const std::vector<Ints>& dark = {},
                   const std::optional<Ints>& fusing = std::nullopt, const std::vector<Form>& forms = {},
                   std::optional<Amount> offer0 = std::nullopt) {
  const auto L = static_cast<std::int64_t>(cap.size());
  const Ints fuse = fusing.value_or(fused_links(offers, routes, cap.size(), dark));
  auto fused = [&](std::int64_t e, std::int64_t f) { return fuse[static_cast<std::size_t>(e)] == f; };
  auto form = [&](std::size_t i) { return i < forms.size() ? forms[i] : kTunnel; };
  auto own = [&](std::size_t i) { return form(i) != kEntryTunnel && form(i) != kEntryDag; };
  // Nodes: heads 0..n-1, the destination n, then per link a tail and a head (a fused head is the
  // next link's tail).
  const auto n = static_cast<NodeId>(routes.size());
  const NodeId dest = n;
  std::vector<NodeId> tail(cap.size());
  std::vector<NodeId> head(cap.size());
  NodeId next_node = n + 1;
  for (std::int64_t e = 0; e < L; ++e) tail[static_cast<std::size_t>(e)] = next_node++;
  for (std::int64_t e = 0; e < L; ++e) {
    const auto f = fuse[static_cast<std::size_t>(e)];
    head[static_cast<std::size_t>(e)] = f >= 0 ? tail[static_cast<std::size_t>(f)] : next_node++;
  }
  std::vector<InEdge> rows;
  for (std::int64_t e = 0; e < L; ++e) rows.push_back({tail[static_cast<std::size_t>(e)], head[static_cast<std::size_t>(e)], cap[static_cast<std::size_t>(e)], 1});
  std::vector<std::vector<std::size_t>> rows_of(routes.size());
  std::vector<std::size_t> exit_row(routes.size(), 0);
  for (std::size_t i = 0; i < routes.size(); ++i) {
    const auto& r = routes[i];
    rows_of[i].push_back(rows.size());
    rows.push_back({own(i) ? static_cast<NodeId>(i) : NodeId {0}, tail[static_cast<std::size_t>(r.front())], kBig, 1});
    for (std::size_t j = 0; j < r.size(); ++j) {
      rows_of[i].push_back(static_cast<std::size_t>(r[j]));
      if (j + 1 < r.size() && !fused(r[j], r[j + 1])) {
        rows_of[i].push_back(rows.size());
        rows.push_back({head[static_cast<std::size_t>(r[j])], tail[static_cast<std::size_t>(r[j + 1])], kBig, 1});
      }
    }
    rows_of[i].push_back(rows.size());
    rows.push_back({head[static_cast<std::size_t>(r.back())], dest, kBig, 1});
    if (form(i) == kDagBesideExit || form(i) == kZeroBranch) {
      exit_row[i] = rows.size();
      rows.push_back({static_cast<NodeId>(i), dest, kBig, 1});
    }
  }
  const Net net(rows, next_node);
  auto edges_of = [&](const std::vector<std::size_t>& rs) {
    std::vector<EdgeId> out;
    for (auto row : rs) out.push_back(net.id[row]);
    return out;
  };
  auto tunnel_of = [&](const std::vector<std::size_t>& rs) {
    std::vector<std::vector<EdgeId>> hops;
    for (EdgeId e : edges_of(rs)) hops.push_back({e});
    return std::make_shared<const Path>(*net.g, hops);
  };
  auto ledger = std::make_shared<FlowLedger>(net.g);
  std::vector<std::shared_ptr<const Eligibility>> black;
  std::vector<std::unique_ptr<FlowPolicy>> policies;
  std::vector<std::size_t> route_of;  // per participant
  std::vector<Action> list;
  for (std::size_t i = 0; i < routes.size(); ++i) {
    if (!own(i)) continue;
    const auto h = static_cast<NodeId>(i);
    std::vector<RouteDefinition> defs;
    Amount offer = offers[i];
    switch (form(i)) {
      case kTunnel:
        defs.push_back(RouteDefinition {h, tunnel_of(rows_of[i]), i == 0 ? offers[0] : 1});
        break;
      case kDagRoute:
        defs.push_back(RouteDefinition {h, dag_of(*net.g, edges_of(rows_of[i])), 1});
        break;
      case kDagBesideExit:
      case kZeroBranch: {
        auto rs = rows_of[i];
        if (form(i) == kZeroBranch) rs.push_back(exit_row[i]);
        defs.push_back(RouteDefinition {h, dag_of(*net.g, edges_of(rs)), 1});
        defs.push_back(RouteDefinition {h, tunnel_of({exit_row[i]}), 1});
        offer = form(i) == kZeroBranch ? 2 : 2 * offers[i];
        break;
      }
      default:
        break;
    }
    if (i == 0) {  // route 0's program holds the entries
      for (std::size_t k = 1; k < routes.size(); ++k) {
        if (form(k) == kEntryTunnel) defs.push_back(RouteDefinition {0, tunnel_of(rows_of[k]), offers[k]});
        if (form(k) == kEntryDag) defs.push_back(RouteDefinition {0, dag_of(*net.g, edges_of(rows_of[k])), offers[k]});
        if (!own(k)) offer += offers[k];
      }
      offer = offer0.value_or(offer);
    }
    const NodeId s[1] = {h};
    const NodeId t[1] = {dest};
    const auto prog = ForwardingProgram::from_routes(*net.g, s, t, defs);
    black.push_back(i < dark.size() ? net.mask({}, dark[i]) : nullptr);
    policies.push_back(std::make_unique<FlowPolicy>(*ledger, std::vector<NodeId> {h}, std::vector<NodeId> {dest},
                                                    FlowPolicyConfig {.routing = Pinned {prog, black.back()}, .admission = LY}));
    route_of.push_back(i);
    list.push_back(Replace {*policies.back(), offer});
  }
  ListOut out;
  std::vector<PlacementResult> res;
  out.error = error_of([&] { res = apply_actions(*ledger, list, BatchMode::Simultaneous); });
#if NETGRAPH_CORE_PRIVATE_HOOKS
  const auto& counters = detail::thread_lossy_counters();
  out.rounds = counters.rounds;
  out.groups = counters.groups;
  out.dependencies = counters.dependencies;
  for (std::int64_t e = 0; e < L; ++e) {
    const EdgeId id = net.id[static_cast<std::size_t>(e)];
    if (std::find(counters.barriers.begin(), counters.barriers.end(), id) != counters.barriers.end()) out.barrier_links.push_back(e);
  }
#endif
  if (!out.error.empty()) return out;
  for (std::size_t k = 0; k < policies.size(); ++k) {
    const std::size_t i = route_of[k];
    const auto recs = policies[k]->records();
    Ints a;
    Ints c;
    Ints d;
    for (auto e : routes[i]) {
      const EdgeId id = net.id[static_cast<std::size_t>(e)];
      Amount cf = 0;
      Amount df = 0;
      for (const auto& rec : recs) {
        for (const auto& [x, v] : rec.edge_flows) cf += x == id ? v : 0;
        for (const auto& [x, v] : rec.edge_drops) df += x == id ? lo(v) : 0;
      }
      a.push_back(cf + df);
      c.push_back(cf);
      d.push_back(df);
    }
    out.arrivals.push_back(a);
    out.carried.push_back(c);
    out.edge_drops.push_back(d);
    const Amount exit = form(i) == kDagBesideExit ? offers[i] : form(i) == kZeroBranch ? 2 : 0;
    out.delivered.push_back(res[k].call_delivered - exit);
    out.status.push_back(form(i) == kTunnel ? name_of(res[k].by_tunnel.front().status) : "SPLIT");
  }
  for (std::int64_t e = 0; e < L; ++e) out.edge_flows.push_back(ledger->carried()[static_cast<std::size_t>(net.id[static_cast<std::size_t>(e)])]);
  return out;
}

// The list form of the stale-forwarding cycle examples: each head is its own demand (the reference's
// labels), a Pinned LOSSY policy over its tunnels with the example's black holes, created in head
// order, all in one SIMULTANEOUS list; a head's destination is its tunnels' end.
ListOut stale_list(const StalePathsIn& in) {
  const Net n(in.edges);
  auto ledger = std::make_shared<FlowLedger>(n.g);
  const auto black = n.mask(in.black_nodes, in.black_links);
  std::vector<std::unique_ptr<FlowPolicy>> policies;
  std::vector<Action> list;
  for (const auto& [h, offer] : in.offers) {
    std::vector<RouteDefinition> defs;
    NodeId dest = 0;
    for (const auto& [th, es] : in.tunnels) {
      if (th != h) continue;
      const auto p = n.path(es);
      dest = p->destination();
      defs.push_back(RouteDefinition {h, p, 1});
    }
    const NodeId s[1] = {h};
    const NodeId t[1] = {dest};
    const auto prog = ForwardingProgram::from_routes(*n.g, s, t, defs);
    policies.push_back(std::make_unique<FlowPolicy>(*ledger, std::vector<NodeId> {h}, std::vector<NodeId> {dest},
                                                    FlowPolicyConfig {.routing = Pinned {prog, black}, .admission = LY}));
    list.push_back(Replace {*policies.back(), offer});
  }
  ListOut out;
  std::vector<PlacementResult> res;
  out.error = error_of([&] { res = apply_actions(*ledger, list, BatchMode::Simultaneous); });
#if NETGRAPH_CORE_PRIVATE_HOOKS
  out.rounds = detail::thread_lossy_counters().rounds;
#endif
  // The statuses, from each policy's own view of its program (the same whatever the forwarding).
  for (const auto& p : policies) {
    const auto& pe = detail::PolicyAccess::engine(*p);
    for (std::size_t e = 0; e < pe.view.data->entries.size(); ++e) out.status.push_back(name_of(pe.view.status(e)));
  }
  if (!out.error.empty()) return out;
  for (const auto& r : res) out.delivered.push_back(r.call_delivered);
  return out;
}

// The examples' list helper for one-tunnel LOSSY participants on a fresh ledger, created in list
// order: per entry its delivery, drops, tunnel status, residual and saturated_boundary_edges (from
// the entry's own source).
struct ListAction {
  NodeId source;
  Ints route;
  Amount offer;
  Ints black_holes {};
};
struct EntryOut {
  Amount delivered;
  std::vector<std::pair<std::int64_t, Amount>> call_edge_drops;
  std::string status;
  Amounts residual;
  Ints saturated_boundary_edges;
};
std::vector<EntryOut> action_list(const std::vector<InEdge>& edges, const std::vector<ListAction>& actions, BatchMode mode) {
  const Net n(edges);
  auto ledger = std::make_shared<FlowLedger>(n.g);
  std::vector<std::shared_ptr<const Eligibility>> black;
  std::vector<std::unique_ptr<FlowPolicy>> policies;
  std::vector<Action> list;
  for (const auto& a : actions) {
    const auto p = n.path(a.route);
    const NodeId s[1] = {a.source};
    const NodeId t[1] = {p->destination()};
    const std::vector<RouteDefinition> defs {RouteDefinition {a.source, p, 1}};
    black.push_back(n.mask({}, a.black_holes));
    policies.push_back(std::make_unique<FlowPolicy>(*ledger, std::vector<NodeId> {a.source}, std::vector<NodeId> {p->destination()},
                                                    FlowPolicyConfig {.routing = Pinned {ForwardingProgram::from_routes(*n.g, s, t, defs), black.back()},
                                                                      .admission = LY}));
    list.push_back(Replace {*policies.back(), a.offer});
  }
  const auto res = apply_actions(*ledger, list, mode, OutputOptions {.residual = true, .saturated_boundary_edges = true});
  std::vector<EntryOut> out;
  for (const auto& r : res) {
    EntryOut e {r.call_delivered, {}, name_of(r.by_tunnel.front().status), n.in_order(*r.residual), {}};
    for (const auto& [x, d] : r.call_edge_drops) e.call_edge_drops.emplace_back(n.pos(x), d);
    std::sort(e.call_edge_drops.begin(), e.call_edge_drops.end());
    for (EdgeId x : *r.saturated_boundary_edges) e.saturated_boundary_edges.push_back(n.pos(x));
    std::sort(e.saturated_boundary_edges.begin(), e.saturated_boundary_edges.end());
    out.push_back(e);
  }
  return out;
}

void expect_entry(const EntryOut& e, Amount delivered, const std::vector<std::pair<std::int64_t, Amount>>& drops, const std::string& status,
                  const Amounts& residual, const Ints& boundary) {
  EXPECT_EQ(e.delivered, delivered);
  EXPECT_EQ(e.call_edge_drops, drops);
  EXPECT_EQ(e.status, status);
  EXPECT_EQ(e.residual, residual);
  EXPECT_EQ(e.saturated_boundary_edges, boundary);
}

const std::vector<InEdge> kCycleEdges = {{0, 1, 12, 1},     {1, 6, 1000, 1},   {2, 3, 9, 1},      {3, 0, 1000, 1},
                                         {3, 10, 1000, 1},  {4, 5, 7, 1},      {5, 6, 1000, 1},   {6, 7, 13, 1},
                                         {7, 2, 1000, 1},   {7, 10, 1000, 1},  {8, 9, 1, 1},      {10, 11, 10, 1},
                                         {11, 2, 1000, 1},  {12, 4, 1000, 1},  {13, 6, 1000, 1},  {14, 0, 1000, 1}};
const std::vector<std::pair<NodeId, Ints>> kCycleTunnels = {
    {12, {13, 5, 6, 7, 8, 2, 4, 11}}, {13, {14, 7, 9, 11, 12, 2, 3, 0}}, {14, {15, 0, 1, 7}}};
const std::vector<InEdge> kCycleBlackEdges = {
    {0, 1, 12, 1},    {1, 6, 1000, 1},  {1, 16, 1000, 1}, {2, 3, 9, 1},     {3, 0, 1000, 1},  {3, 10, 1000, 1}, {4, 5, 7, 1},
    {5, 6, 1000, 1},  {6, 7, 13, 1},    {7, 2, 1000, 1},  {7, 10, 1000, 1}, {8, 9, 1, 1},     {10, 11, 10, 1},  {11, 2, 1000, 1},
    {12, 4, 1000, 1}, {13, 6, 1000, 1}, {14, 0, 1000, 1}, {15, 0, 1000, 1}, {16, 2, 1000, 1}};
const std::vector<std::pair<NodeId, Ints>> kCycleBlackTunnels = {
    {12, {14, 6, 7, 8, 9, 3, 5, 12}}, {13, {15, 8, 10, 12, 13, 3, 4, 0}}, {14, {16, 0, 1, 8}}, {15, {17, 0, 2, 18, 3, 5, 12}}};

const std::string kPeriodTwo = "lossy forwarding iteration entered a period-2 cycle";
const std::string kRoundLimit = "lossy forwarding iteration did not reach a fixed point within the round limit";

void expect_list(const ListOut& r, const Routes& arrivals, const Routes& carried, const std::optional<Routes>& drops, const Amounts& delivered,
                 const Amounts& flows) {
  EXPECT_EQ(r.error, "");
  EXPECT_EQ(r.arrivals, arrivals);
  EXPECT_EQ(r.carried, carried);
  if (drops.has_value()) EXPECT_EQ(r.edge_drops, *drops);
  EXPECT_EQ(r.delivered, delivered);
  EXPECT_EQ(r.edge_flows, flows);
}

// D-no-H-fixed-point's oracle, the test's own: unit-at-a-time highest averages by exact
// cross-multiplication (ties to the earlier label), sharing nothing with the library.
Amounts own_h(Amount x, const Amounts& w) {
  Amounts a(w.size(), 0);
  for (Amount u = 0; u < x; ++u) {
    std::size_t best = w.size();
    for (std::size_t i = 0; i < w.size(); ++i) {
      if (w[i] <= 0) continue;
      if (best == w.size() || w[i] * (a[best] + 1) > w[best] * (a[i] + 1)) best = i;
    }
    if (best == w.size()) return Amounts(w.size(), 0);
    ++a[best];
  }
  return a;
}

} // namespace

TEST(ContractExamples, D_crossing) {
  expect_list(list_lossy({8, 4, 4}, {{0, 1}, {0}, {1}}, {9, 7}), {{8, 6}, {4}, {4}}, {{6, 4}, {3}, {3}}, Routes {{2, 2}, {1}, {1}}, {4, 3, 3},
              {9, 7});
}
TEST(ContractExamples, D_C3) {
  const auto r = list_lossy({1, 1, 2}, {{0, 1}, {1, 2}, {2, 0}}, {1, 1, 2});
  EXPECT_EQ(r.error, kPeriodTwo);
#if NETGRAPH_CORE_PRIVATE_HOOKS
  EXPECT_EQ(r.rounds, 3);
#endif
}
TEST(ContractExamples, D_cycle_prefix_suffix) {
  expect_list(list_lossy({8, 4}, {{0, 1, 2}, {2, 1, 3}}, {2, 5, 5, 10}), {{8, 2, 2}, {4, 3, 3}}, {{2, 2, 2}, {3, 3, 3}},
              Routes {{6, 0, 0}, {1, 0, 0}}, {2, 3}, {2, 5, 5, 3});
}
TEST(ContractExamples, D_cycle_subdivision_before) {
  expect_list(list_lossy({2, 1, 4}, {{1, 0, 2}, {2, 0}, {2, 0, 1}}, {2, 4, 2}), {{2, 2, 1}, {1, 0}, {4, 2, 1}}, {{2, 1, 0}, {0, 0}, {2, 1, 1}},
              std::nullopt, {0, 0, 1}, {2, 3, 2});
}
TEST(ContractExamples, D_cycle_subdivision_after) {
  expect_list(list_lossy({2, 1, 4}, {{1, 0, 2, 3}, {2, 3, 0}, {2, 3, 0, 1}}, {2, 4, 2, 2}), {{2, 2, 1, 0}, {1, 0, 0}, {4, 2, 2, 1}},
              {{2, 1, 0, 0}, {0, 0, 0}, {2, 2, 1, 1}}, std::nullopt, {0, 0, 1}, {2, 3, 2, 2});
}
// A dead occurrence (a zero offer or share, traffic behind its own black hole) takes no part in the
// no-drop loads, the waits or serial normalization, so adding it changes no other output, structure
// included: with and without it, the other participants' rows and deliveries, every link flow, the
// error and (hooks) the links that wait, the groups, their rounds and the dependency arcs are equal.
// "Without" is the list without the dead participant, or (behind a black hole) with its route cut
// after its black hole, as the contract's property states it. `extra_arcs`: the dependency arcs of
// the added element's own live part that the comparison has none of.
void expect_same_without(const ListOut& with, const ListOut& without, std::size_t n, const std::string& where,
                         std::int64_t extra_arcs = 0) {
  EXPECT_EQ(with.error, without.error) << where;
  EXPECT_EQ(with.rounds, without.rounds) << where;
  EXPECT_EQ(with.groups, without.groups) << where;
  EXPECT_EQ(with.barrier_links, without.barrier_links) << where;
#if NETGRAPH_CORE_PRIVATE_HOOKS
  EXPECT_EQ(with.dependencies, without.dependencies + extra_arcs) << where;
#else
  (void)extra_arcs;
#endif
  auto first = [n](const auto& v) {
    using V = std::decay_t<decltype(v)>;
    return V(v.begin(), v.begin() + static_cast<std::ptrdiff_t>(std::min(n, v.size())));
  };
  EXPECT_EQ(first(with.arrivals), first(without.arrivals)) << where;
  EXPECT_EQ(first(with.carried), first(without.carried)) << where;
  EXPECT_EQ(first(with.edge_drops), first(without.edge_drops)) << where;
  EXPECT_EQ(first(with.delivered), first(without.delivered)) << where;
  EXPECT_EQ(with.edge_flows, without.edge_flows) << where;
}
std::string form_name(Form f) {
  static const char* const kNames[] = {"tunnel", "DAG route", "DAG route beside an exit", "zero branch", "entry tunnel", "entry DAG route"};
  return kNames[static_cast<std::size_t>(f)];
}
const Form kDeadForms[] = {kTunnel, kDagRoute, kDagBesideExit};

// D-dead-tail: Z's traffic dies on link 3, its own black hole, so its occurrences on links 0 and 1
// are dead: the period-2 error at round 3, as with Z cut after link 3 (and as without Z), with Z a
// tunnel, a DAG route alone or a DAG route beside an exit. The scalar call (one head, the routes
// its entries with weights = offers, the call's black hole on link 3, which only Z uses) gives the
// same, Z a tunnel or a DAG route, against Z cut after link 3.
void dead_tail(const Amounts& offers) {
  const Routes routes {{2, 0}, {1, 0, 2}, {3, 0, 1}};
  const Routes cut {{2, 0}, {1, 0, 2}, {3}};
  const Amounts cap {2, 2, 1, 9};
  const std::vector<Ints> dark {{}, {}, {3}};
  const auto alone = list_lossy({offers[0], offers[1]}, {routes[0], routes[1]}, cap);
  EXPECT_EQ(alone.error, kPeriodTwo);
  for (Form f : kDeadForms) {
    const std::vector<Form> forms {kTunnel, kTunnel, f};
    const auto with = list_lossy(offers, routes, cap, dark, std::nullopt, forms);
    EXPECT_EQ(with.error, kPeriodTwo) << form_name(f);
#if NETGRAPH_CORE_PRIVATE_HOOKS
    EXPECT_EQ(with.rounds, 3) << form_name(f);
#endif
    expect_same_without(with, list_lossy(offers, cut, cap, dark, std::nullopt, forms), 2, form_name(f) + " against Z cut");
    EXPECT_EQ(with.barrier_links, alone.barrier_links) << form_name(f);
    EXPECT_EQ(with.groups, alone.groups) << form_name(f);
  }
  for (bool dag : {false, true}) {
    const auto with = lossy(offers, routes, cap, {}, {3}, {false, false, dag});
    const auto without = lossy(offers, cut, cap, {}, {3}, {false, false, dag});
    EXPECT_EQ(with.error, kPeriodTwo) << "scalar, dag " << dag;
    EXPECT_EQ(with.error, without.error) << "scalar, dag " << dag;
    EXPECT_EQ(with.rounds, without.rounds) << "scalar, dag " << dag;
    EXPECT_EQ(with.groups, without.groups) << "scalar, dag " << dag;
    EXPECT_EQ(with.barrier_links, without.barrier_links) << "scalar, dag " << dag;
    EXPECT_EQ(with.dependencies, without.dependencies) << "scalar, dag " << dag;
#if NETGRAPH_CORE_PRIVATE_HOOKS
    EXPECT_EQ(with.rounds, 3) << "scalar, dag " << dag;
#endif
  }
}
TEST(ContractExamples, D_dead_tail) { dead_tail({1, 2, 1}); }
// D-dead-waits: link 1 is a barrier with or without Z; Z's dead occurrences add no wait, so links 0
// and 2 alone form the group.
TEST(ContractExamples, D_dead_waits) { dead_tail({1, 3, 3}); }
// D-zero-share-pair: D-cycle-subdivision-after plus a zero offer over link 2 alone, which ends there:
// links 2 and 3 stay a transparent pair, with Z a tunnel, a DAG route alone (a pooled head whose
// share is 0), a DAG route beside an exit (a split route whose share is 0) or a zero branch (a
// split route member whose whole-unit share is 0; its exit member's unit makes one arc of its own).
TEST(ContractExamples, D_zero_share_pair) {
  const Routes routes {{1, 0, 2, 3}, {2, 3, 0}, {2, 3, 0, 1}, {2}};
  const Amounts cap {2, 4, 2, 2};
  const auto without = list_lossy({2, 1, 4}, {routes[0], routes[1], routes[2]}, cap);
  for (Form f : {kTunnel, kDagRoute, kDagBesideExit, kZeroBranch}) {
    const auto r = list_lossy({2, 1, 4, 0}, routes, cap, {}, std::nullopt, {kTunnel, kTunnel, kTunnel, f});
    expect_list(r, {{2, 2, 1, 0}, {1, 0, 0}, {4, 2, 2, 1}, {0}}, {{2, 1, 0, 0}, {0, 0, 0}, {2, 2, 1, 1}, {0}}, std::nullopt, {0, 0, 1, 0},
                {2, 3, 2, 2});
#if NETGRAPH_CORE_PRIVATE_HOOKS
    EXPECT_EQ(r.rounds, 4) << form_name(f);
#endif
    expect_same_without(r, without, 3, form_name(f), f == kZeroBranch ? 1 : 0);
  }
}
// D-dead-subdivision: Z over links 2 and 3 with its own black hole on link 2 alone: its unit is
// dropped there and its occurrence on link 3 is dead, as before the split; with Z a tunnel, a DAG
// route alone or a DAG route beside an exit, against Z cut after link 2.
TEST(ContractExamples, D_dead_subdivision) {
  const Routes routes {{1, 0, 2, 3}, {2, 3, 0}, {2, 3, 0, 1}, {2, 3}};
  const Routes cut {{1, 0, 2, 3}, {2, 3, 0}, {2, 3, 0, 1}, {2}};
  const Amounts cap {2, 4, 2, 2};
  const std::vector<Ints> dark {{}, {}, {}, {2}};
  for (Form f : kDeadForms) {
    const std::vector<Form> forms {kTunnel, kTunnel, kTunnel, f};
    const auto r = list_lossy({2, 1, 4, 1}, routes, cap, dark, std::nullopt, forms);
    expect_list(r, {{2, 2, 1, 0}, {1, 0, 0}, {4, 2, 2, 1}, {1, 0}}, {{2, 1, 0, 0}, {0, 0, 0}, {2, 2, 1, 1}, {0, 0}},
                Routes {{0, 1, 1, 0}, {1, 0, 0}, {2, 0, 1, 0}, {1, 0}}, {0, 0, 1, 0}, {2, 3, 2, 2});
    if (f == kTunnel) EXPECT_EQ(r.status.back(), "BLACKHOLED");
    expect_same_without(r, list_lossy({2, 1, 4, 1}, cut, cap, dark, std::nullopt, forms), 3, form_name(f) + " against Z cut");
  }
}

// A dead occurrence loads no link, so it makes no barrier. A: links 1, 0 (offer 1); B: 0, 2, 1
// (offer 2); V = 2, 1, 2. With or without Z (links 3, 2, offer 2, its own black hole on link 3; a
// tunnel, a DAG route alone or a DAG route beside an exit) the period-2 error at round 3; Z's unit
// on link 2 counted in the no-drop loads would make link 2 a third barrier of the group and the
// iteration fail at the round limit (192 rounds).
TEST(DeadOccurrences, LoadNoLink) {
  const Amounts cap {2, 1, 2, 9};
  const auto without = list_lossy({1, 2}, {{1, 0}, {0, 2, 1}}, cap);
  EXPECT_EQ(without.error, kPeriodTwo);
#if NETGRAPH_CORE_PRIVATE_HOOKS
  EXPECT_EQ(without.rounds, 3);
#endif
  for (Form f : kDeadForms) {
    const std::vector<Form> forms {kTunnel, kTunnel, f};
    const auto with = list_lossy({1, 2, 2}, {{1, 0}, {0, 2, 1}, {3, 2}}, cap, {{}, {}, {3}}, std::nullopt, forms);
    expect_same_without(with, list_lossy({1, 2, 2}, {{1, 0}, {0, 2, 1}, {3}}, cap, {{}, {}, {3}}, std::nullopt, forms), 2,
                        form_name(f) + " against Z cut");
    EXPECT_EQ(with.error, kPeriodTwo) << form_name(f);
    EXPECT_EQ(with.barrier_links, without.barrier_links) << form_name(f);
  }
}
// An empty entry (weight 0: it gets no share and is not forwarded) changes nothing: D-cycle-
// subdivision-after with a tunnel or a DAG route of weight 0 in A's program, over link 2 alone (it
// would break the pair 2, 3) or over links 2, 3, 0.
TEST(DeadOccurrences, EmptyEntriesChangeNothing) {
  const Routes routes {{1, 0, 2, 3}, {2, 3, 0}, {2, 3, 0, 1}};
  const Amounts cap {2, 4, 2, 2};
  const auto without = list_lossy({2, 1, 4}, routes, cap);
  for (Form f : {kEntryTunnel, kEntryDag}) {
    for (const Ints& extra : {Ints {2}, Ints {2, 3, 0}}) {
      Routes with_routes = routes;
      with_routes.push_back(extra);
      const auto with = list_lossy({2, 1, 4, 0}, with_routes, cap, {}, std::nullopt, {kTunnel, kTunnel, kTunnel, f});
      expect_same_without(with, without, 3, form_name(f) + " over " + std::to_string(extra.size()) + " links");
      EXPECT_EQ(with.delivered, (Amounts {0, 0, 1}));
    }
  }
}
// The scalar call with a real serial pair: D-zero-share-pair's routes as the entries of one head's
// program on the list layout (links 2 and 3 joined directly), weights 2, 1, 4, 1 at offer 7, so
// highest averages give Z (link 2 alone, then its own exit) a share of 0. Z a tunnel or a DAG
// route: delivery 1 and the structure of the call without Z (weights 2, 1, 4 at offer 7); counted in
// normalization, Z would keep links 2 and 3 from passing traffic through as one.
TEST(DeadOccurrences, ScalarSerialPairWithAZeroShareEntry) {
  const Routes routes {{1, 0, 2, 3}, {2, 3, 0}, {2, 3, 0, 1}, {2}};
  const Amounts cap {2, 4, 2, 2};
  const Ints fuse = fused_links({2, 1, 4, 0}, routes, cap.size(), {});
  const auto without = list_lossy({2, 1, 4}, {routes[0], routes[1], routes[2]}, cap, {}, fuse, {kTunnel, kEntryTunnel, kEntryTunnel}, 7);
  EXPECT_EQ(without.error, "");
  EXPECT_EQ(without.delivered, (Amounts {1}));
  for (Form f : {kEntryTunnel, kEntryDag}) {
    const auto with = list_lossy({2, 1, 4, 1}, routes, cap, {}, fuse, {kTunnel, kEntryTunnel, kEntryTunnel, f}, 7);
    expect_same_without(with, without, 1, "scalar, Z an " + form_name(f));
#if NETGRAPH_CORE_PRIVATE_HOOKS
    EXPECT_EQ(with.barrier_links, (Ints {0, 1, 2})) << form_name(f);  // link 3 passes link 2's traffic through
#endif
  }
}
// The scalar call: a dead entry changes no other output, structure included. D-crossing (A over
// links 0, 1 with weight 8; B over 0, 4; C over 1, 4; V 9, 7) plus Z over its own link 2 (the
// call's black hole) then link 0, weight 4, offer 20: Z's share dies on link 2. And D-dead-tail's
// routes with weights 1, 2, 1 and offer 3 (no black hole): highest averages give Z a share of 0.
// Z a tunnel or a DAG route; the outputs equal those without Z (D-crossing: offer 16; D-dead-tail:
// weights 1, 2 and offer 3, the same shares), and D-crossing's dependency arcs those with Z cut after
// its black hole (a DAG route's live part keeps the arcs of its own routers).
TEST(DeadOccurrences, ScalarCallsWithoutThem) {
  for (bool dag : {false, true}) {
    const std::string where = dag ? "DAG route" : "tunnel";
    const auto with = lossy({8, 4, 4, 4}, {{0, 1}, {0}, {1}, {2, 0}}, {9, 7, 9}, {}, {2}, {false, false, false, dag});
    const auto without = lossy({8, 4, 4}, {{0, 1}, {0}, {1}}, {9, 7, 9});
    EXPECT_EQ(with.error, "") << where;
    EXPECT_EQ(Routes(with.carried.begin(), with.carried.begin() + 3), without.carried) << where;
    EXPECT_EQ(Routes(with.edge_drops.begin(), with.edge_drops.begin() + 3), without.edge_drops) << where;
    EXPECT_EQ(Amounts(with.delivered.begin(), with.delivered.begin() + 3), without.delivered) << where;
    EXPECT_EQ(with.delivered[0], 4) << where;
    EXPECT_EQ(with.edge_flows, without.edge_flows) << where;
    EXPECT_EQ(with.barrier_links, without.barrier_links) << where;
    EXPECT_EQ(with.groups, without.groups) << where;
    // Z's live part up to its black hole (a DAG route's routers) keeps its own arcs: against Z cut.
    const auto cut = lossy({8, 4, 4, 4}, {{0, 1}, {0}, {1}, {2}}, {9, 7, 9}, {}, {2}, {false, false, false, dag});
    EXPECT_EQ(with.dependencies, cut.dependencies) << where;
    EXPECT_EQ(with.barrier_links, cut.barrier_links) << where;
    const auto zero = lossy({1, 2, 1}, {{2, 0}, {1, 0, 2}, {3, 0, 1}}, {2, 2, 1, 9}, {}, {}, {false, false, dag}, 3);
    const auto none = lossy({1, 2}, {{2, 0}, {1, 0, 2}}, {2, 2, 1, 9}, {}, {}, {}, 3);
    EXPECT_EQ(zero.error, kPeriodTwo) << where;
    EXPECT_EQ(zero.error, none.error) << where;
    EXPECT_EQ(zero.rounds, none.rounds) << where;
    EXPECT_EQ(zero.groups, none.groups) << where;
    EXPECT_EQ(zero.barrier_links, none.barrier_links) << where;
    EXPECT_EQ(zero.dependencies, none.dependencies) << where;
  }
}
// The forwarding sizes follow the same loads: D-dead-tail has two barriers (links 0 and 2) and
// two dependency arcs (B: 0 before 2; A: 2 before 0; Z: none) with or without Z, so its event and
// dependency counts are those of the network without Z (counted in Z's no-drop load, link 1 would
// be a third barrier, with arcs to and from link 0).
TEST(DeadOccurrences, ForwardingSizesWithoutThem) {
#if NETGRAPH_CORE_PRIVATE_HOOKS
  using detail::ForwardingCounter;
  struct Restore {
    ~Restore() { detail::thread_forwarding_limits() = {}; }
  } restore;
  // The smallest limit of a counter that the call meets.
  auto smallest = [](ForwardingCounter counter, const std::function<ListOut()>& call) {
    for (std::int64_t n = 1; n < 20; ++n) {
      detail::thread_forwarding_limits() = {};
      detail::thread_forwarding_limits()[static_cast<std::size_t>(counter)] = n;
      if (call().error.find("limit exceeded") == std::string::npos) return n;
    }
    return std::int64_t {-1};
  };
  const std::function<ListOut()> with_z = [] { return list_lossy({1, 2, 1}, {{2, 0}, {1, 0, 2}, {3, 0, 1}}, {2, 2, 1, 9}, {{}, {}, {3}}); };
  const std::function<ListOut()> without_z = [] { return list_lossy({1, 2}, {{2, 0}, {1, 0, 2}}, {2, 2, 1, 9}); };
  EXPECT_EQ(smallest(ForwardingCounter::Event, with_z), 2);
  EXPECT_EQ(smallest(ForwardingCounter::Event, without_z), 2);
  EXPECT_EQ(smallest(ForwardingCounter::Dependency, with_z), 2);
  EXPECT_EQ(smallest(ForwardingCounter::Dependency, without_z), 2);
#else
  GTEST_SKIP() << "needs the private hooks";
#endif
}

// Python's random.Random(seed) for a small non-negative seed (MT19937 seeded from the seed's word, and
// the random(), randint, randrange and sample of CPython 3.11 and later on populations of at most
// 21), so the property below draws the contract's own inputs.
class PyRandom {
public:
  explicit PyRandom(std::uint32_t seed) {
    init(19650218U);
    std::size_t i = 1;
    for (std::size_t k = kN; k != 0; --k) {  // the key array's mixing, with `seed` its one word
      mt_[i] = (mt_[i] ^ ((mt_[i - 1] ^ (mt_[i - 1] >> 30)) * 1664525U)) + seed;
      if (++i >= kN) {
        mt_[0] = mt_[kN - 1];
        i = 1;
      }
    }
    for (std::size_t k = kN - 1; k != 0; --k) {
      mt_[i] = (mt_[i] ^ ((mt_[i - 1] ^ (mt_[i - 1] >> 30)) * 1566083941U)) - static_cast<std::uint32_t>(i);
      if (++i >= kN) {
        mt_[0] = mt_[kN - 1];
        i = 1;
      }
    }
    mt_[0] = 0x80000000U;
  }
  double random() {
    const std::uint32_t a = next() >> 5;
    const std::uint32_t b = next() >> 6;
    return (a * 67108864.0 + b) * (1.0 / 9007199254740992.0);
  }
  std::int64_t randbelow(std::int64_t n) {
    int k = 0;
    while ((std::int64_t {1} << k) <= n) ++k;
    std::int64_t r = next() >> (32 - k);
    while (r >= n) r = next() >> (32 - k);
    return r;
  }
  std::int64_t randint(std::int64_t a, std::int64_t b) { return a + randbelow(b - a + 1); }
  Ints sample(std::int64_t n, std::int64_t k) {  // sample(range(n), k)
    Ints pool(static_cast<std::size_t>(n));
    for (std::int64_t i = 0; i < n; ++i) pool[static_cast<std::size_t>(i)] = i;
    Ints out;
    for (std::int64_t i = 0; i < k; ++i) {
      const auto j = static_cast<std::size_t>(randbelow(n - i));
      out.push_back(pool[j]);
      pool[j] = pool[static_cast<std::size_t>(n - i - 1)];
    }
    return out;
  }

private:
  static constexpr std::size_t kN = 624;
  void init(std::uint32_t s) {
    mt_[0] = s;
    for (std::size_t i = 1; i < kN; ++i) mt_[i] = 1812433253U * (mt_[i - 1] ^ (mt_[i - 1] >> 30)) + static_cast<std::uint32_t>(i);
    index_ = kN;
  }
  std::uint32_t next() {
    if (index_ >= kN) {
      for (std::size_t i = 0; i < kN; ++i) {
        const std::uint32_t y = (mt_[i] & 0x80000000U) | (mt_[(i + 1) % kN] & 0x7fffffffU);
        mt_[i] = mt_[(i + 397) % kN] ^ (y >> 1) ^ ((y & 1U) != 0 ? 0x9908b0dfU : 0U);
      }
      index_ = 0;
    }
    std::uint32_t y = mt_[index_++];
    y ^= y >> 11;
    y ^= (y << 7) & 0x9d2c5680U;
    y ^= (y << 15) & 0xefc60000U;
    y ^= y >> 18;
    return y;
  }
  std::array<std::uint32_t, kN> mt_ {};
  std::size_t index_ {kN};
};

// D-dead-demand-property (seed 4, 2000 trials): one added dead demand Z changes nothing. Small random
// inputs of simple-path demands; Z has a zero offer over 1 to 3 of their links, or an offer of 1 to 3
// whose route reaches its own black hole (link k, a link of Z's own, inserted into its route, or one of
// its shared links) and continues behind it. The outcome with Z equals the outcome with Z's dead
// occurrences removed (without Z for a zero offer, Z's route cut after its black hole otherwise): every
// arrival, carried amount, drop and delivery of the others and of Z's prefix, every link flow, and an
// error with its round count (hooks); Z's dead occurrences show 0. Structure included (hooks): the
// links that wait, the groups and the dependency arcs are equal too. Z is a tunnel, a DAG route alone
// (a pooled context) and a DAG route beside an exit (a split route of a program that is not pooled),
// each on the same inputs. Both runs of a trial share one layout (fused_links of the input with Z).
TEST(ContractExamples, D_dead_demand_property) {
  PyRandom rng(4);
  std::array<std::int64_t, 3> agree {};
  std::array<std::int64_t, 3> errors {};
  std::int64_t zero = 0;
  std::int64_t tail = 0;
  for (int trial = 0; trial < 2000; ++trial) {
    const auto k = rng.randint(3, 5);
    Amounts cap;
    for (std::int64_t e = 0; e < k; ++e) cap.push_back(rng.randint(1, 3));
    cap.push_back(9);
    const auto n = static_cast<std::size_t>(rng.randint(2, 4));
    Routes routes;
    for (std::size_t i = 0; i < n; ++i) {
      const auto len = rng.randint(2, std::min<std::int64_t>(4, k));
      routes.push_back(rng.sample(k, len));
    }
    Amounts offers;
    for (std::size_t i = 0; i < n; ++i) offers.push_back(rng.randint(1, 4));
    const Ints links = rng.sample(k, rng.randint(1, std::min<std::int64_t>(3, k)));
    Amount z = 0;
    Ints zr = links;
    std::int64_t hole = -1;
    if (rng.random() < 0.5) {
      ++zero;
    } else {
      ++tail;
      z = rng.randint(1, 3);
      const auto p = rng.randbelow(static_cast<std::int64_t>(links.size()));
      if (rng.random() < 0.5) {
        zr.insert(zr.begin() + p, k);
        hole = k;
      } else {
        hole = links[static_cast<std::size_t>(p)];
      }
    }
    std::vector<Ints> dark(n);
    dark.push_back(hole < 0 ? Ints {} : Ints {hole});
    Amounts all_offers = offers;
    all_offers.push_back(z);
    Routes all_routes = routes;
    all_routes.push_back(zr);
    const Ints fuse = fused_links(all_offers, all_routes, cap.size(), dark);
    const std::size_t cut = hole < 0 ? 0 : static_cast<std::size_t>(std::find(zr.begin(), zr.end(), hole) - zr.begin()) + 1;
    Routes cut_routes = all_routes;
    if (cut != 0) cut_routes.back().resize(cut);
    // The comparison: the others' rows and Z's prefix, deliveries, link flows, Z's dead tail all 0,
    // the structure.
    auto view = [&](const ListOut& r) {
      std::vector<Routes> rows;
      bool quiet = true;
      const auto structure = std::make_tuple(r.rounds, r.groups, r.barrier_links, r.dependencies);
      if (!r.error.empty()) return std::make_tuple(r.error, structure, rows, Amounts {}, Amounts {}, quiet);
      for (const Routes* key : {&r.arrivals, &r.carried, &r.edge_drops}) {
        Routes x(key->begin(), key->begin() + static_cast<std::ptrdiff_t>(std::min(n, key->size())));
        if (cut != 0) x.push_back(Ints((*key)[n].begin(), (*key)[n].begin() + static_cast<std::ptrdiff_t>(cut)));
        if (key != &r.edge_drops && key->size() > n) {
          for (std::size_t j = cut; j < (*key)[n].size(); ++j) quiet = quiet && (*key)[n][j] == 0;
        }
        rows.push_back(x);
      }
      const Amounts delivered(r.delivered.begin(), r.delivered.begin() + static_cast<std::ptrdiff_t>(std::min(n, r.delivered.size())));
      return std::make_tuple(r.error, structure, rows, delivered, r.edge_flows, quiet);
    };
    for (std::size_t f = 0; f < 3; ++f) {
      std::vector<Form> forms(n, kTunnel);
      forms.push_back(kDeadForms[f]);
      const auto full = list_lossy(all_offers, all_routes, cap, dark, fuse, forms);
      const auto shorter = hole < 0 ? list_lossy(offers, routes, cap, {}, fuse) : list_lossy(all_offers, cut_routes, cap, dark, fuse, forms);
      const bool same = view(full) == view(shorter);
      EXPECT_TRUE(same) << "trial " << trial << ", Z a " << form_name(kDeadForms[f]);
      agree[f] += same ? 1 : 0;
      errors[f] += full.error.empty() ? 0 : 1;
    }
  }
  EXPECT_EQ(zero, 993);
  EXPECT_EQ(tail, 1007);
  for (std::size_t f = 0; f < 3; ++f) {
    EXPECT_EQ(agree[f], 2000) << form_name(kDeadForms[f]);
    EXPECT_EQ(errors[f], 55) << form_name(kDeadForms[f]);
  }
}
// The documented call: two ECMP LOSSY policies on one link of 10.
TEST(ContractExamples, D_simultaneous) {
  const Net n({{0, 1, 10, 1}});
  FlowLedger ledger(n.g);
  const FlowPolicyConfig ecmp {.routing = ShortestPaths {}, .admission = LY};
  FlowPolicy a(ledger, v({0}), v({1}), ecmp);
  FlowPolicy b(ledger, v({0}), v({1}), ecmp);
  const auto res = apply_actions(ledger, std::vector<Action> {Replace {a, 10}, Replace {b, 10}}, BatchMode::Simultaneous);
  EXPECT_EQ(res[0].call_delivered, 5);
  EXPECT_EQ(res[1].call_delivered, 5);
  EXPECT_EQ(res[0].call_edge_drops, (std::vector<std::pair<EdgeId, Amount>> {{0, 5}}));
  EXPECT_EQ(res[1].call_edge_drops, (std::vector<std::pair<EdgeId, Amount>> {{0, 5}}));
  EXPECT_EQ(ledger.carried(), (Amounts {10}));
}
// The C3 instance (A: 0, 1; B: 1, 2; C: 2, 0; offers 1, 1, 2; V 1, 1, 2): the test's own enumeration
// of every integer carried record within the continuation bounds and work conservation finds 5, none
// meeting the loss rule's link equations (each record's violating links named); the library's highest averages
// agrees with the test's own on every equation the enumeration evaluates; the list raises the
// period-2 error.
TEST(ContractExamples, D_no_H_fixed_point) {
  const Amounts cap {1, 1, 2};
  int records = 0;
  int solutions = 0;
  std::vector<std::vector<int>> violations;
  for (Amount a0 = 0; a0 <= 1; ++a0) {
    for (Amount a1 = 0; a1 <= a0; ++a1) {
      for (Amount b1 = 0; b1 <= 1; ++b1) {
        for (Amount b2 = 0; b2 <= b1; ++b2) {
          for (Amount c2 = 0; c2 <= 2; ++c2) {
            for (Amount c0 = 0; c0 <= c2; ++c0) {
              // Arrivals: A 1 then a0; B 1 then b1; C 2 then c2. Links: 0 (A0, C0), 1 (A1, B1), 2 (B2, C2).
              const Amounts arr0 {1, c2};
              const Amounts car0 {a0, c0};
              const Amounts arr1 {a0, 1};
              const Amounts car1 {a1, b1};
              const Amounts arr2 {b1, 2};
              const Amounts car2 {b2, c2};
              auto conserving = [&](const Amounts& arr, const Amounts& car, Amount v) {
                return car[0] + car[1] == std::min(v, arr[0] + arr[1]);
              };
              if (!conserving(arr0, car0, cap[0]) || !conserving(arr1, car1, cap[1]) || !conserving(arr2, car2, cap[2])) continue;
              ++records;
              std::vector<int> bad;
              const std::vector<std::pair<Amounts, Amounts>> links {{arr0, car0}, {arr1, car1}, {arr2, car2}};
              for (std::size_t l = 0; l < 3; ++l) {
                const Amount t = std::min(cap[l], links[l].first[0] + links[l].first[1]);
                const Amounts mine = own_h(t, links[l].first);
                Amounts lib(2, 0);
                detail::highest_averages(t, links[l].first, lib);
                EXPECT_EQ(lib, mine) << "link " << l;
                if (mine != links[l].second) bad.push_back(static_cast<int>(l));
              }
              solutions += bad.empty() ? 1 : 0;
              violations.push_back(bad);
            }
          }
        }
      }
    }
  }
  EXPECT_EQ(records, 5);
  EXPECT_EQ(solutions, 0);
  EXPECT_EQ(violations, (std::vector<std::vector<int>> {{2}, {0}, {0, 1, 2}, {1}, {0}}));
  EXPECT_EQ(list_lossy({1, 1, 2}, {{0, 1}, {1, 2}, {2, 0}}, {1, 1, 2}).error, kPeriodTwo);
}
TEST(ContractExamples, STALE_cycle_period_2) {
  StalePathsIn in {{{0, 1, 768, 1}, {1, 2, 25600, 1}, {2, 3, 768, 1}, {3, 0, 25600, 1}, {4, 2, 25600, 1}, {5, 0, 25600, 1},
                    {6, 2, 25600, 1}, {7, 8, 25600, 1}},
                   {{4, {4, 2, 3, 0}}, {5, {5, 0, 1, 2}}, {6, {6, 2}}, {7, {7}}},
                   {{4, 1024}, {5, 256}, {6, 256}, {7, 512}}};
  in.black_links = {7};
  const auto r = stale_list(in);
  EXPECT_EQ(r.error, kPeriodTwo);
#if NETGRAPH_CORE_PRIVATE_HOOKS
  EXPECT_EQ(r.rounds, 8);
#endif
}
TEST(ContractExamples, STALE_cycle_round_limit_baseline) {
  const StalePathsIn in {kCycleEdges, kCycleTunnels, {{12, 13}, {13, 11}, {14, 10}}};
  const auto r = stale_list(in);
  EXPECT_EQ(r.error, "");
  EXPECT_EQ(r.delivered, (Amounts {3, 4, 4}));
}
TEST(ContractExamples, STALE_cycle_round_limit) {
  StalePathsIn in {kCycleEdges, kCycleTunnels, {{12, 13}, {13, 11}, {14, 10}}};
  in.black_links = {5};
  const auto r = stale_list(in);
  EXPECT_EQ(r.error, kRoundLimit);
#if NETGRAPH_CORE_PRIVATE_HOOKS
  EXPECT_EQ(r.rounds, 256);
#endif
}
TEST(ContractExamples, STALE_black_hole_in_cycle) {
  StalePathsIn in {kCycleBlackEdges, kCycleBlackTunnels, {{12, 13}, {13, 11}, {14, 10}, {15, 5}}};
  in.black_links = {2, 6, 17};
  const auto r = stale_list(in);
  EXPECT_EQ(r.status, (Texts {"BLACKHOLED", "UP", "UP", "BLACKHOLED"}));
  EXPECT_EQ(r.error, kRoundLimit);
#if NETGRAPH_CORE_PRIVATE_HOOKS
  EXPECT_EQ(r.rounds, 256);
#endif
}
TEST(ContractExamples, LIST_black_hole_in_cycle) {
  std::vector<ListAction> actions;
  const Amounts offers {13, 11, 10, 5};
  for (std::size_t i = 0; i < kCycleBlackTunnels.size(); ++i) actions.push_back({kCycleBlackTunnels[i].first, kCycleBlackTunnels[i].second, offers[i], {2, 6, 17}});
  EXPECT_EQ(error_of([&] { (void)action_list(kCycleBlackEdges, actions, BatchMode::Simultaneous); }), kRoundLimit);
#if NETGRAPH_CORE_PRIVATE_HOOKS
  EXPECT_EQ(detail::thread_lossy_counters().rounds, 256);
#endif
}
TEST(ContractExamples, LIST_black_hole_own) {
  const auto r = action_list({{0, 1, 10, 1}}, {{0, {0}, 10, {0}}, {0, {0}, 10}}, BatchMode::Simultaneous);
  expect_entry(r[0], 0, {{0, 10}}, "BLACKHOLED", {0}, {0});
  expect_entry(r[1], 10, {}, "UP", {0}, {0});
}
TEST(ContractExamples, LIST_black_hole_upstream) {
  const auto r = action_list({{0, 1, 10, 1}, {1, 2, 10, 1}}, {{0, {0, 1}, 10, {1}}, {0, {0, 1}, 10}}, BatchMode::Simultaneous);
  expect_entry(r[0], 0, {{0, 5}, {1, 5}}, "BLACKHOLED", {0, 5}, {0});
  expect_entry(r[1], 5, {{0, 5}}, "UP", {0, 5}, {0});
}
TEST(ContractExamples, LIST_outputs_ordered) {
  const auto r = action_list({{0, 1, 10, 1}}, {{0, {0}, 7}, {0, {0}, 7}}, BatchMode::Ordered);
  expect_entry(r[0], 7, {}, "UP", {3}, {});
  expect_entry(r[1], 3, {{0, 4}}, "UP", {0}, {0});
}
TEST(ContractExamples, LIST_outputs_simultaneous) {
  const auto r = action_list({{0, 2, 5, 1}, {1, 2, 10, 1}, {2, 3, 10, 1}}, {{0, {0, 2}, 5}, {1, {1, 2}, 5}}, BatchMode::Simultaneous);
  expect_entry(r[0], 5, {}, "UP", {0, 5, 0}, {0});
  expect_entry(r[1], 5, {}, "UP", {0, 5, 0}, {2});
}
// The fourteen lists with real policies (a, b of the list's ledger; x of another; the cycling action
// is the cycle's Pinned replace(3)).
TEST(ContractExamples, LIST_first_error) {
  const Net n({{0, 1, 10, 1}});
  const Net other({{0, 1, 10, 1}});
  const Net cyc({{0, 1, 1, 1}, {1, 2, 1000, 1}, {1, 5, 1000, 1}, {2, 3, 2, 1}, {3, 0, 1000, 1}, {3, 5, 1000, 1}, {4, 0, 1000, 1}, {4, 2, 1000, 1},
                 {6, 7, 10, 1}});
  FlowLedger ledger(n.g);
  FlowLedger foreign(other.g);
  FlowLedger cl(cyc.g);
  const FlowPolicyConfig ly {.routing = ShortestPaths {}, .admission = LY};
  const FlowPolicyConfig ll {.routing = ShortestPaths {}, .admission = LL};
  FlowPolicy a(ledger, v({0}), v({1}), ly);
  FlowPolicy x(foreign, v({0}), v({1}), ly);
  FlowPolicy all(ledger, v({0}), v({1}), ll);
  FlowPolicy t(ledger, v({0}), v({1}), {.routing = dte_with(1), .admission = LL});
  FlowPolicy c(ledger, v({0}), v({1}), {.routing = centralized()});
  (void)c.place(1);
  const NodeId s4[1] = {4};
  const NodeId t5[1] = {5};
  const std::vector<RouteDefinition> two {RouteDefinition {4, cyc.path({6, 0, 1, 3, 5}), 1}, RouteDefinition {4, cyc.path({7, 3, 4, 0, 2}), 2}};
  const auto pc_prog = ForwardingProgram::from_routes(*cyc.g, s4, t5, two);
  const NodeId s6[1] = {6};
  const NodeId t7[1] = {7};
  const std::vector<RouteDefinition> one {RouteDefinition {6, cyc.path({8}), 1}};
  const auto pb_prog = ForwardingProgram::from_routes(*cyc.g, s6, t7, one);
  FlowPolicy pc(cl, v({4}), v({5}), {.routing = Pinned {pc_prog, nullptr}, .admission = LY});
  FlowPolicy pb(cl, v({6}), v({7}), {.routing = Pinned {pb_prog, nullptr}, .admission = LY});
  FlowPolicy p(ledger, v({0}), v({1}), {.routing = Pinned {ForwardingProgram::from_routes(*n.g, std::vector<NodeId> {0}, std::vector<NodeId> {1},
                                                                                          std::vector<RouteDefinition> {RouteDefinition {0, n.path({0}), 1}}),
                                                           nullptr},
                                        .admission = LY});
  FlowPolicy inf(ledger, v({0}), v({1}), ly);
  (void)inf.place(std::nullopt);
  FlowPolicy fw(ledger, v({0}), v({1}), ly);
  (void)fw.replace(1);
  const auto foreign_black = other.mask({}, {0});
  const BatchMode O = BatchMode::Ordered;
  const BatchMode S = BatchMode::Simultaneous;
  auto first = [&](FlowLedger& l, std::vector<Action> list, BatchMode mode) {
    const std::string e = error_of([&] { (void)apply_actions(l, list, mode); });
    return e.empty() ? std::string("valid") : e;
  };
  const Texts got {
      first(ledger, {Replace {a, -1}, Replace {x, 1}}, O),
      first(ledger, {Replace {a, 1}, Replace {a, -1}}, O),
      first(ledger, {Forward {fw, std::nullopt}}, S),
      first(ledger, {Place {all, 1}}, S),
      first(ledger, {Replace {t, 1}}, S),
      first(ledger, {Forward {c, std::nullopt}}, S),
      first(cl, {Replace {pc, 3}, Replace {pb, Amount {1} << 62}}, S),
      first(cl, {Replace {pc, 3}, Replace {pb, 1}}, S),
      first(ledger, {Place {t, -1}}, O),
      first(ledger, {Forward {p, std::nullopt}}, O),
      first(ledger, {Place {inf, 5}, Replace {x, 1}}, O),
      first(ledger, {Replace {a, 1, finite(5)}, Replace {x, 1}}, O),
      first(ledger, {Forward {fw, 1, foreign_black.get()}, Replace {x, 1}}, O),
      first(ledger, {}, S),
  };
  EXPECT_EQ(got, (Texts {"offer must be an integer in [0, 2^62) or permitted +inf", "duplicate policy in list", "SIMULTANEOUS requires finite offers",
                         "SIMULTANEOUS does not accept Place", "SIMULTANEOUS Replace requires ShortestPaths EQUAL or WEIGHTED, or Pinned",
                         "SIMULTANEOUS requires admission=LOSSY", "offer must be an integer in [0, 2^62) or permitted +inf", kPeriodTwo,
                         "DistributedTE does not support incremental placement", "Pinned forwards with replace",
                         "list policy belongs to another ledger", "ShortestPaths does not accept design_demand",
                         "black_holes belongs to another graph", "valid"}));
}
// The forwarding sizes' counting step on the example's sizes (per participant, groups of tunnels
// of one path length): the entry count, then the occurrence count, each in wide arithmetic against
// 2^31 - 1 (the network itself cannot be built at these sizes; the end-to-end refusals run under
// the hooks' small limits, ForwardingSize and ApplyActions tests).
TEST(ContractExamples, ERR_forwarding_limit) {
  using Groups = std::vector<std::pair<std::int64_t, std::int64_t>>;
  auto check = [](const std::vector<Groups>& participants) {
    std::int64_t entries = 0;
    std::int64_t occurrences = 0;
    for (const auto& g : participants) {
      for (const auto& [count, length] : g) {
        entries += count;
        occurrences += count * length;
      }
    }
    try {
      detail::check_forwarding_count(detail::ForwardingCounter::Entry, entries);
      detail::check_forwarding_count(detail::ForwardingCounter::Occurrence, occurrences);
    } catch (const std::invalid_argument& e) {
      return std::string(e.what());
    }
    return std::string("within the entry and occurrence limits");
  };
#if NETGRAPH_CORE_PRIVATE_HOOKS
  detail::thread_forwarding_limits() = {};
#endif
  EXPECT_EQ(check({{{10000, 214749}}}), "forwarding occurrence limit exceeded: 2147490000 > 2147483647");
  EXPECT_EQ(check({{{10000, 214748}, {1, 3647}}}), "within the entry and occurrence limits");
  EXPECT_EQ(check({{{10000, 214748}, {1, 3648}}}), "forwarding occurrence limit exceeded: 2147483648 > 2147483647");
  EXPECT_EQ(check({{{1, 1073741824}}}), "within the entry and occurrence limits");
  EXPECT_EQ(check({{{1, 1073741824}}, {{1, 1073741824}}}), "forwarding occurrence limit exceeded: 2147483648 > 2147483647");
  EXPECT_EQ(check({{{2147483648, 1}}}), "forwarding entry limit exceeded: 2147483648 > 2147483647");
}
// The list's offers as C++ values (the Python conversions, a bool, a string, a fraction, a float or
// an int beyond 64 bits, have no C++ form): -1 fails the call's domain, 5 and +inf are valid.
TEST(ContractExamples, ERR_class_list) {
  const Net n({{0, 1, 10, 1}});
  auto run = [&](std::optional<Amount> offer) {
    FlowLedger ledger(n.g);
    FlowPolicy p(ledger, v({0}), v({1}), {.routing = ShortestPaths {}, .admission = LY});
    const std::string e = error_of([&] { (void)apply_actions(ledger, std::vector<Action> {Place {p, offer}}, BatchMode::Ordered); });
    return e.empty() ? std::string("valid") : e;
  };
  EXPECT_EQ(run(-1), "offer must be an integer in [0, 2^62) or permitted +inf");
  EXPECT_EQ(run(5), "valid");
  EXPECT_EQ(run(std::nullopt), "valid");
}

TEST(ContractExamples, FPublicationLocked) {
  const auto out = detail::publication_model(true);
  EXPECT_EQ(out.mixed, 0U);
  EXPECT_EQ(out.states, 15U);
}
TEST(ContractExamples, FPublicationUnlocked) {
  const auto out = detail::publication_model(false);
  EXPECT_EQ(out.mixed, 2U);
  EXPECT_EQ(out.states, 23U);
}
namespace {
void expect_phases(std::string_view name, std::size_t states, std::size_t partial,
                   std::size_t python_under_lock, std::size_t late_hits) {
  const auto out = detail::phase_model(name);
  EXPECT_EQ(out.states, states);
  EXPECT_EQ(out.partial, partial);
  EXPECT_EQ(out.python_under_lock, python_under_lock);
  EXPECT_EQ(out.late_hits, late_hits);
}
} // namespace
TEST(ContractExamples, FPublishResultFirst) { expect_phases("publish-result-first", 22, 0, 0, 0); }
TEST(ContractExamples, FPublishResultAfter) { expect_phases("publish-result-after", 22, 1, 0, 0); }
TEST(ContractExamples, FReadResultFailure) { expect_phases("read-result-failure", 20, 0, 0, 0); }
TEST(ContractExamples, FReadResultUnderGate) { expect_phases("read-result-under-gate", 18, 0, 2, 0); }
TEST(ContractExamples, FFrozenHits) { expect_phases("frozen-hits", 25, 0, 0, 0); }
TEST(ContractExamples, FFrozenHitsAfterLease) { expect_phases("frozen-hits-after-lease", 28, 0, 0, 2); }

TEST(ContractExamples, LIST_state_second) {
  const Texts expected {
    "offer must be an integer in [0, 2^62) or permitted +inf",
    "black_holes belongs to another graph",
    "list policy belongs to another ledger",
    "offer must be an integer in [0, 2^62) or permitted +inf",
    "offer must be an integer in [0, 2^62) or permitted +inf",
    "requested limit exceeded",
    "forward requires an installed program",
    "finite placement after infinite placement requires a new epoch",
    "offer must be an integer in [0, 2^62) or permitted +inf"};
  Texts observed;
  for (std::size_t which = 0; which < expected.size(); ++which) {
    const Net net({{0, 1, 10, 1}}), foreign({{0, 1, 10, 1}});
    FlowLedger ledger(net.g), other(net.g);
    FlowPolicy a(ledger, v({0}), v({1}), {.routing = ShortestPaths{}, .admission = LY});
    FlowPolicy b(ledger, v({0}), v({1}), {.routing = ShortestPaths{}, .admission = LY});
    FlowPolicy x(other, v({0}), v({1}), {.routing = ShortestPaths{}, .admission = LY});
    const auto program = ForwardingProgram::from_routes(*net.g, v({0}), v({1}),
        std::vector<RouteDefinition>{{0, net.path({0}), 1}});
    FlowPolicy pin(ledger, v({0}), v({1}), {.routing = Pinned{program}, .admission = LY});
    const auto holes = foreign.mask({}, {0});
    if (which == 3 || which == 7) (void)a.place(std::nullopt);
    if (which == 6) (void)b.place(std::nullopt);
    if (which == 4 || which == 5) {
      const auto gate = detail::LedgerAccess::lock(ledger);
      detail::PolicyAccess::engine(a).epoch.requested = Counter128{(std::uint64_t{1} << 63) - 1, ~std::uint64_t{0}};
    }
    std::vector<Action> actions;
    switch (which) {
      case 0: actions = {Forward{a, -1}}; break;
      case 1: actions = {Forward{a, 1, holes.get()}}; break;
      case 2: actions = {Forward{a, 1}, Replace{x, 1}}; break;
      case 3: actions = {Place{a, 5}, Replace{b, -1}}; break;
      case 4: actions = {Place{a, 1}, Replace{b, -1}}; break;
      case 5: actions = {Place{a, 1}}; break;
      case 6: actions = {Forward{a, 1}, Place{b, 5}}; break;
      case 7: actions = {Place{a, 5}, Forward{b, 1}}; break;
      case 8: actions = {Forward{a, 1}, Replace{pin, kDomainLimit}}; break;
    }
    const auto carried = ledger.carried();
    const auto holds = ledger.reserved();
    const auto a_program = a.program(), b_program = b.program(), pin_program = pin.program();
    const auto record_id = detail::LedgerAccess::next_record_counter(ledger);
    observed.push_back(error_of([&] { (void)apply_actions(ledger, actions,
        which == 8 ? BatchMode::Simultaneous : BatchMode::Ordered); }));
    EXPECT_EQ(ledger.carried(), carried);
    EXPECT_EQ(ledger.reserved(), holds);
    EXPECT_EQ(a.program(), a_program);
    EXPECT_EQ(b.program(), b_program);
    EXPECT_EQ(pin.program(), pin_program);
    EXPECT_EQ(detail::LedgerAccess::next_record_counter(ledger), record_id);
  }
  EXPECT_EQ(observed, expected);
}

// The numeric/ownership forms of CALL-state-second. Wrong Python types and
// malformed pair shapes have no C++ value and are checked by its Python port.
TEST(ContractExamples, CALL_state_second) {
  const Net net({{0, 1, 10, 1}}), foreign({{0, 1, 10, 1}});
  FlowLedger ledger(net.g);
  FlowPolicy p(ledger, v({0}), v({1}), {.routing = ShortestPaths{}, .admission = LY});
  const auto holes = foreign.mask({}, {0});
  const auto program = ForwardingProgram::from_routes(*net.g, v({0}), v({1}),
      std::vector<RouteDefinition>{{0, net.path({0}), 1}});
  FlowPolicy pin(ledger, v({0}), v({1}), {.routing = Pinned{program}, .admission = LY});
  ShortestPathCache cache(*net.g), wrong(*foreign.g);
  const std::vector<std::pair<std::vector<NodeId>, std::vector<NodeId>>> pairs {{{0}, {1}}};
  EvaluationOptions options;
  options.cache = &cache;
  const std::string mode = "a supplied building ShortestPathCache has one writer: it is rejected unless workers == 1 (got workers = 2; a frozen cache may be shared by every worker)";
  Texts observed;
  observed.push_back(error_of([&] { (void)p.forward(-1); }));
  observed.push_back(error_of([&] { (void)p.forward(std::nullopt); }));
  observed.push_back(error_of([&] { (void)p.forward(1, holes.get()); }));
  observed.push_back(error_of([&] { (void)pin.forward(-1); }));
  observed.push_back(error_of([&] { (void)p.forward(1); }));
  options.config.routing = Pinned{program};
  observed.push_back(error_of([&] { (void)evaluate_pairs(*net.g, pairs, options, {}, 2); }));
  options.config.routing = ShortestPaths{};
  options.offer = kUnlimited;
  observed.push_back(error_of([&] { (void)evaluate_pairs(*net.g, pairs, options, {}, 2); }));
  options.offer.reset();
  options.cache = &wrong;
  observed.push_back(error_of([&] { (void)evaluate_pairs(*net.g, pairs, options, {}, 2); }));
  options.cache = &cache;
  observed.push_back(error_of([&] { (void)evaluate_pairs(*net.g, pairs, options, {}, 2); }));
  options.offer = kUnlimited;
  observed.push_back(error_of([&] { (void)evaluate_edge_failures(*net.g, v({0}), v({1}), options, 2); }));
  options.offer.reset();
  observed.push_back(error_of([&] { (void)evaluate_edge_failures(*net.g, v({0}), v({1}), options, 2); }));
  EXPECT_EQ(observed, (Texts{
      "offer must be an integer in [0, 2^62) or permitted +inf", "LOSSY forwarding requires a finite offer",
      "black_holes belongs to another graph", "Pinned forwards with replace", "forward requires an installed program",
      "evaluate_pairs does not accept Pinned", "offer must be an integer in [0, 2^62)",
      "foreign ShortestPathCache: bound to another graph", mode, "offer must be an integer in [0, 2^62)", mode}));
}

namespace {

// Independently decode b = significand * 2^exponent, then divide integers with a half-up
// remainder test. No production kernel, floating scaling, floor or floating +1/2 is used.
Amount binary_units_reference(double b, int k) {
  const auto bits = std::bit_cast<std::uint64_t>(b);
  const auto biased = static_cast<int>((bits >> 52) & 0x7ffU);
  const auto significand = (bits & ((std::uint64_t {1} << 52) - 1)) |
                          (biased == 0 ? 0 : (std::uint64_t {1} << 52));
  const int exponent = (biased == 0 ? -1074 : biased - 1023 - 52) - k;
  if (exponent >= 0) return static_cast<Amount>(significand << exponent);
  const int down = -exponent;
  if (down >= 64) return 0;
  const auto remainder = significand & ((std::uint64_t {1} << down) - 1);
  return static_cast<Amount>((significand >> down) + (remainder >= (std::uint64_t {1} << (down - 1))));
}

void unit_error(double b, double q, Amount expected, double error_in_units) {
  const auto u = to_units(b, q);
  EXPECT_EQ(u, expected);
  const auto product = as_float(u, q);
  EXPECT_EQ(product, static_cast<double>(u) * q);
  EXPECT_EQ((product - b) / q, error_in_units);
  EXPECT_GT((product - b) / q, -0.5);
  EXPECT_LE((product - b) / q, 0.5);
}

} // namespace

TEST(ContractExamples, UNITS_literal_tie) {
  const double q = 0x1p-21;
  unit_error(5501195.1796, q, 11536842473289LL, 0.5);
  EXPECT_EQ(to_units(std::nextafter(5501195.1796, 0.0), q), 11536842473288LL);
}

TEST(ContractExamples, UNITS_computed_tie) {
  const double q = 0x1p-18;
  unit_error(400.0 + 0x1p-19, q, 104857601, 0.5);
  EXPECT_EQ(to_units(400.0, q), 104857600);
}

TEST(ContractExamples, UNITS_float32) {
  const double q = 0x1p-32;
  // float32(0.1) = 13421773/134217728, not the decimal its short printing suggests.
  const double b = static_cast<double>(0.1F);
  EXPECT_EQ(b, 13421773.0 / 134217728.0);
  unit_error(b, q, 429496736, 0.0);
  EXPECT_EQ(to_units(0.1, q), 429496730);
}

TEST(ContractExamples, UNITS_ties) {
  EXPECT_EQ(to_units(0x1p-33, 0x1p-32), 1);
  EXPECT_EQ(to_units(std::nextafter(0x1p-33, 0.0), 0x1p-32), 0);
  EXPECT_EQ(to_units(0x1p31, 0x1p32), 1);
  for (Amount n = 0; n != 4; ++n) EXPECT_EQ(to_units(static_cast<double>(n) + 0.5, 1.0), n + 1);
  EXPECT_EQ(to_units(0x1p52 - 0.5, 1.0), Amount {1} << 52);
  EXPECT_EQ(to_units(0x1p52 + 1, 1.0), (Amount {1} << 52) + 1);
  EXPECT_EQ(to_units(0x1p53 + 2, 1.0), (Amount {1} << 53) + 2);
  EXPECT_EQ(to_units(-0.0, 1.0), 0);
  EXPECT_EQ(to_units(std::numeric_limits<double>::denorm_min(), 0x1p32), 0);
}

TEST(ContractExamples, UNITS_round_trip) {
  for (double q : {1.0, 0x1p-32, 0x1p32}) {
    for (Amount u : {Amount {1} << 53, (Amount {1} << 53) + 2, (Amount {1} << 62) - 512}) {
      EXPECT_EQ(to_units(as_float(u, q), q), u);
    }
    EXPECT_EQ(to_units(as_float((Amount {1} << 53) + 1, q), q), Amount {1} << 53);
    EXPECT_EQ(as_float((Amount {1} << 53) + 1, q), 0x1p53 * q);
    EXPECT_THROW((void)to_units(as_float((Amount {1} << 62) - 1, q), q), std::invalid_argument);
  }
  EXPECT_EQ(to_units(as_float(120401795907070192LL, 0x1p-27), 0x1p-27), 120401795907070192LL);
}

TEST(ContractExamples, UNITS_promoted_integer) {
  const Amount first = (Amount {1} << 53) + 1;
  const Amount second = (Amount {1} << 53) + 3;
  unit_error(static_cast<double>(first), 1.0, first - 1, 0.0);
  unit_error(static_cast<double>(second), 1.0, second + 1, 0.0);
  // Compare as integers to preserve the promotion error which a double subtraction would hide.
  EXPECT_EQ(to_units(static_cast<double>(first), 1.0) - first, -1);
  EXPECT_EQ(to_units(static_cast<double>(second), 1.0) - second, 1);
}

TEST(ContractExamples, UNITS_as_float_range) {
  EXPECT_EQ(as_float(-1, 0.5), -0.5);
  EXPECT_EQ(as_float((Amount {1} << 53) + 3, 1.0), 9007199254740996.0);
  // The native signature is int64, so the unbounded-integer finite-double boundary and
  // 2^70+1 case belong to the Python API tests. Native int64 extrema are always finite.
  for (int k = -32; k <= 32; ++k) {
    const double q = std::ldexp(1.0, k);
    EXPECT_EQ(as_float(std::numeric_limits<Amount>::max(), q), std::ldexp(1.0, 63 + k));
    EXPECT_EQ(as_float(std::numeric_limits<Amount>::min(), q), -std::ldexp(1.0, 63 + k));
    EXPECT_EQ(as_float(-(Amount {1} << 53), q), -0x1p53 * q);
    EXPECT_EQ(as_float(-((Amount {1} << 53) + 3), q), -(0x1p53 + 4) * q);
  }
}

TEST(ContractExamples, UNITS_sums) {
  EXPECT_EQ(to_units(0.5, 1.0), 1);
  EXPECT_EQ(to_units(0.5, 1.0) + to_units(0.5, 1.0), 2);
  EXPECT_EQ(to_units(0.5 + 0.5, 1.0), 1);
  EXPECT_EQ(to_units(0.25, 1.0), 0);
  EXPECT_EQ(to_units(0.25, 1.0) + to_units(0.25, 1.0), 0);
  EXPECT_EQ(to_units(0.25 + 0.25, 1.0), 1);
  EXPECT_EQ(as_float(2, 1.0) - (0.5 + 0.5), 1.0);
  EXPECT_EQ(as_float(0, 1.0) - (0.25 + 0.25), -0.5);
}

TEST(ContractExamples, UNITS_double_steps) {
  // Two hand-derived differentiators for each oracle branch: below half/at half, odd/even
  // ties, the last tie/first integral double, and exact/inexact integer promotion.
  EXPECT_EQ(binary_units_reference(0.25, 0), 0);
  EXPECT_EQ(binary_units_reference(0.5, 0), 1);
  EXPECT_EQ(binary_units_reference(1.5, 0), 2);
  EXPECT_EQ(binary_units_reference(2.5, 0), 3);
  EXPECT_EQ(binary_units_reference(0x1p52 - 0.5, 0), Amount {1} << 52);
  EXPECT_EQ(binary_units_reference(0x1p52 + 1, 0), (Amount {1} << 52) + 1);
  EXPECT_EQ(binary_units_reference(0x1p53, 0), Amount {1} << 53);
  EXPECT_EQ(binary_units_reference(0x1p53 + 2, 0), (Amount {1} << 53) + 2);
  EXPECT_EQ(binary_units_reference(0x1p-1074, 32), 0);
  EXPECT_EQ(binary_units_reference(0x1p-1074, -32), 0);
  std::mt19937_64 rng(10);
  int checked = 0;
  int returned = 0;
  for (int k = -32; k <= 32; ++k) {
    const double q = std::ldexp(1.0, k);
    const auto check = [&](double b) {
      const auto expected = binary_units_reference(b, k);
      const auto u = to_units(b, q);
      EXPECT_EQ(u, expected) << "b=" << std::hexfloat << b << " k=" << k;
      const auto value = as_float(u, q);
      EXPECT_EQ(value, static_cast<double>(u) * q);
      EXPECT_GT((value - b) / q, -0.5);
      EXPECT_LE((value - b) / q, 0.5);
      ++checked;
    };
    for (int i = 0; i < 68; ++i) {
      const Amount n = i == 0 ? 0 : i == 1 ? 1 : i == 2 ? Amount {1} << 26 :
                       i == 3 ? (Amount {1} << 52) - 1 : static_cast<Amount>(rng() % (std::uint64_t {1} << 52));
      const double tie = static_cast<double>(2 * n + 1) * q / 2;
      EXPECT_EQ(to_units(tie, q), n + 1);
      check(tie);
      check(std::nextafter(tie, 0.0));
      check(std::nextafter(tie, std::numeric_limits<double>::infinity()));
    }
    for (double edge : {0x1p52, 0x1p53}) {
      check(edge * q);
      check(std::nextafter(edge * q, 0.0));
    }
    check(std::numeric_limits<double>::denorm_min());
    for (int i = 0; i < 64; ++i) {
      // A 53-bit significand, with exponents covering tiny values through <2^62 scaled units.
      const double mantissa = static_cast<double>(rng() >> 11) * 0x1p-53;
      const int exponent = static_cast<int>(rng() % 103) - 40;
      check(std::ldexp(mantissa, exponent) * q);
    }
    for (int i = 0; i < 68; ++i) {
      const Amount u = i == 0 ? 0 : i == 1 ? 1 : i == 2 ? (Amount {1} << 53) - 1 :
                       i == 3 ? Amount {1} << 53 : static_cast<Amount>(rng() % ((std::uint64_t {1} << 53) + 1));
      EXPECT_EQ(to_units(as_float(u, q), q), u);
      ++returned;
    }
  }
  EXPECT_EQ(checked, 17745);
  EXPECT_EQ(returned, 4420);
}
