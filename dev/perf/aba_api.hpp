// The A/B/A row adapter for the current Core API. The rows (aba_rows.hpp) call only `abaapi::`;
// this file implements each call with the public Core API, so the same rows run unchanged on any
// two builds that share this API (README.md: a build with another public API needs this file
// adapted first).
//
// Row mapping of the search rows:
// The selection of a search is spelled in target terms (Sel: multipath, multi_edge,
// prefer_priority as PreferHigherPriority on the Residual priority); the default Sel is
// multipath and multi_edge with the Deterministic tie-break. A supplied residual always gates its
// edges at min_residual = 1 unit (through to_units), and static capacity >= 1 unit is
// unconditional.
// * spf (spf.full.*, except the storage row spf.full.random1m, which runs search_run below):
//   spf(g, {src}, Forward, no targets), the default selection, with or without the masks (as one
//   Eligibility) and the residual. Digest: 16 sampled costs, the DAG's entry count and its first
//   and last via.
// * spf_dag (spf.dag.*): spf(g, {src}, Forward, targets {dst}), the default selection: the DAG of
//   the shortest walks that end in the target. Digest: target_cost, entries, first and last via.
// * spf_to (spf.to.random16k): spf(g, {dst}, Reverse, no targets), the default selection.
// * targetset_* (spf.to.targetset_clos): the Clos graph and the destination set; the timed call is
//   a full reverse spf rooted at the destination set, targets absent.
// * spf_path (spf.single.canonical16k): spf(g, {src}, Forward, targets {dst}) with the row's
//   selection (multipath = false) and the residual: the canonical single path. Digest:
//   target_cost and the walk.
// * ksp (ksp.fabric3_k8): ksp(g, src, dst, KspOptions(k, PathSelection(Deterministic, Capacity,
//   multipath, multi_edge))), no residual, no ceiling, no masks. Digest: the path count, then per
//   path its cost, hop count and B' edge ids (the graph's edge ids, ordered by source, cost, target
//   and input position).
// The kernel rows place on fixed DAGs read with the public spf with targets (base_dag).
//
// The storage rows' interface (README "Rows"): construct_prepare and construct_run (the
// constructor on prepared unit capacities, then the sampled digest of construct_sample),
// queue_choice (the library's queue rule, read through the private detail::bucket_width and
// detail::uses_buckets: the one use of a non-public function here) and search_run (one full
// search plus the sampled digest of search_sample). The digests read only API-independent
// values: node-indexed counts and distances, costs, physical capacities and edges named by
// their input position, so they are the same for any API that builds the same graph.
#pragma once

#include "netgraph/core/dag_cache.hpp"
#include "netgraph/core/evaluators.hpp"
#include "netgraph/core/flow_ledger.hpp"
#include "netgraph/core/flow_policy.hpp"
#include "netgraph/core/graph.hpp"
#include "netgraph/core/search.hpp"
#include "netgraph/core/detail/traversal.hpp"

#include <algorithm>
#include <cstdint>
#include <cstring>
#include <memory>
#include <optional>
#include <span>
#include <string>
#include <utility>
#include <vector>

namespace abaapi {
namespace nc = netgraph::core;

// Digest: order-sensitive 64-bit mix (splitmix64 finaliser over a running state).
struct Digest {
  std::uint64_t h = 0x9e3779b97f4a7c15ULL;
  void add(std::uint64_t v) {
    std::uint64_t z = (h ^ v) + 0x9e3779b97f4a7c15ULL;
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    h = z ^ (z >> 31);
  }
  void add_i(std::int64_t v) { add(static_cast<std::uint64_t>(v)); }
  void add_d(double v) {
    std::uint64_t bits;
    std::memcpy(&bits, &v, sizeof bits);
    add(bits);
  }
};

// API-independent graph description: edges in input order, physical capacities.
struct Arrays {
  int n = 0;
  std::vector<std::int32_t> src, dst;
  std::vector<double> cap;
  std::vector<std::int64_t> cost;
  int num_edges() const { return static_cast<int>(src.size()); }
};

// The resolution a NetGraph scenario supplies by default: one integer unit is 1/4096 of a
// capacity unit.
inline constexpr double kResolution = 1.0 / 4096.0;

inline nc::Graph graph_of(const Arrays& a) {
  std::vector<nc::Amount> units(a.cap.size());
  for (std::size_t i = 0; i < a.cap.size(); ++i) units[i] = nc::to_units(a.cap[i], kResolution);
  return nc::Graph(a.n, a.src, a.dst, units, a.cost, kResolution);
}

// The graph plus the caller's inverse of input_index() (input edge index -> EdgeId), and the
// warm DagCache of the cache.warm row (one per graph, prepopulated by the untimed iteration).
struct Graph {
  std::shared_ptr<const nc::Graph> g;
  std::vector<std::int32_t> edge_of_input;
  std::shared_ptr<nc::DagCache> warm_cache;
  const nc::Graph& ref() const { return *g; }
};

// Construction from physical arrays: to_units on every capacity, then the one constructor.
inline std::uint64_t construct(const Arrays& a) {
  auto g = graph_of(a);
  Digest d;
  d.add_i(g.num_edges());
  d.add_i(g.max_cost());
  d.add_i(g.cost_gcd());
  d.add_i(g.row_offsets()[static_cast<std::size_t>(a.n / 2)]);
  return d.h;
}

inline Graph build(const Arrays& a) {
  Graph out;
  out.g = std::make_shared<const nc::Graph>(graph_of(a));
  out.edge_of_input.assign(a.src.size(), -1);
  auto ii = out.g->input_index();
  for (std::size_t e = 0; e < ii.size(); ++e) out.edge_of_input[static_cast<std::size_t>(ii[e])] = static_cast<std::int32_t>(e);
  return out;
}

// Per-edge values given in input-edge order, in id order.
inline std::vector<double> edge_values(const Graph& g, const std::vector<double>& by_input) {
  std::vector<double> out(by_input.size());
  for (std::size_t i = 0; i < by_input.size(); ++i) out[static_cast<std::size_t>(g.edge_of_input[i])] = by_input[i];
  return out;
}

// Static eligibility: the Eligibility value built from excluded ids (input edges translated).
struct Elig {
  std::optional<nc::Eligibility> el;
  const nc::Eligibility* get() const { return el ? &*el : nullptr; }
};

inline Elig eligibility(const Graph& g, const std::vector<int>& excluded_nodes,
                        const std::vector<int>& excluded_input_edges) {
  std::vector<nc::NodeId> xn(excluded_nodes.begin(), excluded_nodes.end());
  std::vector<nc::EdgeId> xe;
  xe.reserve(excluded_input_edges.size());
  for (int i : excluded_input_edges) xe.push_back(g.edge_of_input[static_cast<std::size_t>(i)]);
  Elig out;
  out.el = g.ref().eligibility(xn, xe);
  return out;
}

inline std::uint64_t eligibility_build_intersect(const Graph& g, const std::vector<int>& common,
                                                 const std::vector<std::vector<int>>& states) {
  const Elig base = eligibility(g, {}, common);
  Digest d;
  d.add(base.el->uid() != 0 ? 1U : 0U);
  for (const auto& st : states) {
    const Elig el = eligibility(g, {}, st);
    const auto both = nc::Eligibility::intersect(*base.el, *el.el);
    std::uint64_t allowed = 0;
    for (auto x : both.edge_allowed()) allowed += x;
    d.add(allowed);
  }
  return d.h;
}

// Search selection, spelled in target terms.
struct Sel {
  bool multipath = true;
  bool multi_edge = true;
  bool prefer_priority = false;
};

inline nc::PathSelection to_sel(const Sel& s) {
  return nc::PathSelection(s.prefer_priority ? nc::EdgeTieBreak::PreferHigherPriority : nc::EdgeTieBreak::Deterministic,
                           nc::EdgePriority::Residual, s.multipath, s.multi_edge);
}

// A residual given in physical units (id order) as units.
inline std::vector<nc::Amount> residual_units(const std::vector<double>& residual) {
  std::vector<nc::Amount> out(residual.size());
  for (std::size_t e = 0; e < residual.size(); ++e) out[e] = nc::to_units(residual[e], kResolution);
  return out;
}

inline nc::SpfResult run_spf(const Graph& g, std::span<const nc::NodeId> roots, nc::Direction dir,
                             std::optional<std::vector<nc::NodeId>> targets, const Sel& s, const Elig* el,
                             const std::vector<double>* residual) {
  nc::SpfOptions o(dir, to_sel(s));
  o.targets = std::move(targets);
  std::vector<nc::Amount> units;
  if (residual != nullptr) {
    units = residual_units(*residual);
    o.residual = std::span<const nc::Amount>(units);
  }
  o.elig = el ? el->get() : nullptr;
  return nc::spf(g.ref(), roots, o);
}

inline void fold_dag(Digest& d, const nc::PredDAG& dag) {
  const auto via = dag.via_edges();
  d.add_i(static_cast<std::int64_t>(via.size()));
  if (!via.empty()) {
    d.add_i(via.front());
    d.add_i(via.back());
  }
}

inline void fold_costs(Digest& d, const std::vector<nc::Cost>& dist) {
  const std::size_t n = dist.size();
  for (std::size_t k = 0; k < 16 && n; ++k) d.add_i(dist[(k * 7919) % n]);
}

inline std::uint64_t spf(const Graph& g, int src, const Elig* el, const std::vector<double>* residual) {
  const nc::NodeId roots[1] = {src};
  auto r = run_spf(g, roots, nc::Direction::Forward, std::nullopt, Sel {}, el, residual);
  Digest d;
  fold_costs(d, *r.costs);
  fold_dag(d, r.dag);
  return d.h;
}

inline std::uint64_t spf_dag(const Graph& g, int src, int dst) {
  const nc::NodeId roots[1] = {src};
  auto r = run_spf(g, roots, nc::Direction::Forward, std::vector<nc::NodeId> {dst}, Sel {}, nullptr, nullptr);
  Digest d;
  d.add_i(*r.target_cost);
  fold_dag(d, r.dag);
  return d.h;
}

inline std::uint64_t spf_to(const Graph& g, int dst) {
  const nc::NodeId roots[1] = {dst};
  auto r = run_spf(g, roots, nc::Direction::Reverse, std::nullopt, Sel {}, nullptr, nullptr);
  Digest d;
  fold_costs(d, *r.costs);
  fold_dag(d, r.dag);
  return d.h;
}

// ---- Destination, endpoint-set and canonical single-path rows (README "Rows")

// spf.to.targetset_clos: the Clos graph and the destination set; nothing else is built.
struct TargetSet {
  std::shared_ptr<Graph> g;
  std::vector<nc::NodeId> destinations;
};

inline TargetSet targetset_prepare(const Arrays& clos, const std::vector<int>& destinations) {
  TargetSet ts;
  ts.g = std::make_shared<Graph>(build(clos));
  ts.destinations.assign(destinations.begin(), destinations.end());
  return ts;
}

inline std::uint64_t targetset_run(const TargetSet& ts) {
  auto r = run_spf(*ts.g, ts.destinations, nc::Direction::Reverse, std::nullopt, Sel {}, nullptr, nullptr);
  Digest d;
  fold_costs(d, *r.costs);
  fold_dag(d, r.dag);
  return d.h;
}

// spf.single.canonical16k: the canonical single path (multipath = false) to a target.
inline std::uint64_t spf_path(const Graph& g, int src, int dst, const Sel& s, const std::vector<double>* residual) {
  const nc::NodeId roots[1] = {src};
  auto r = run_spf(g, roots, nc::Direction::Forward, std::vector<nc::NodeId> {dst}, s, nullptr, residual);
  Digest d;
  d.add_i(*r.target_cost);
  for (auto e : r.dag.via_edges()) d.add_i(e);
  return d.h;
}

// A fixed DAG for placement kernels: the cost-only destination search (setup, not timed).
struct Dag {
  std::shared_ptr<const nc::PredDAG> dag;
};
inline Dag base_dag(const Graph& g, int src, int dst) {
  const nc::NodeId roots[1] = {src};
  auto r = run_spf(g, roots, nc::Direction::Forward, std::vector<nc::NodeId> {dst}, Sel {}, nullptr, nullptr);
  return Dag {std::make_shared<const nc::PredDAG>(std::move(r.dag))};
}

// The placement kernels the rows time; Lossy is EqualBalanced placement with Lossy admission.
enum class Kernel { MaxFlow, EqualBalanced, Lossy, Weighted };

struct KernelAxes {
  nc::FlowPlacement placement = nc::FlowPlacement::MaxFlow;
  nc::Admission admission = nc::Admission::Lossless;
  nc::EdgePriority priority = nc::EdgePriority::Residual;
};

inline KernelAxes axes(Kernel k) {
  KernelAxes a;
  switch (k) {
    case Kernel::MaxFlow: break;
    case Kernel::EqualBalanced: a.placement = nc::FlowPlacement::EqualBalanced; break;
    case Kernel::Lossy: a.placement = nc::FlowPlacement::EqualBalanced; a.admission = nc::Admission::Lossy; break;
    case Kernel::Weighted: a.placement = nc::FlowPlacement::Weighted; a.priority = nc::EdgePriority::Capacity; break;
  }
  return a;
}

inline nc::FlowPolicyConfig config(nc::Routing routing, Kernel k, nc::PathSelection sel) {
  const auto a = axes(k);
  return nc::FlowPolicyConfig {std::move(routing), a.placement, a.admission, sel};
}

// Manual placements on fixed DAGs into one ledger: one Pinned policy with that one bundle per
// placement, destroyed after the call (its record stays in the ledger as an orphan).
struct Ledger {
  nc::FlowLedger fg;
  explicit Ledger(const Graph& g) : fg(g.g) {}
  double place(int src, int dst, const Dag& dag, double amount, Kernel k) {
    nc::Pinned pinned;
    pinned.bundles.push_back(dag.dag);
    const auto a = axes(k);
    const nc::NodeId s[1] = {src};
    const nc::NodeId t[1] = {dst};
    nc::FlowPolicy policy(fg, s, t,
                          nc::FlowPolicyConfig {std::move(pinned), a.placement, a.admission,
                                                nc::PathSelection(nc::EdgeTieBreak::Deterministic, a.priority)});
    const auto r = policy.place(nc::to_units(amount, kResolution));
    return nc::as_float(r.placed, kResolution);
  }
};

inline nc::PathSelection costonly_selection() {
  return nc::PathSelection(nc::EdgeTieBreak::Deterministic, nc::EdgePriority::Capacity, true, true);
}

// Scalar CostOnly policy capacity: place(+inf) of a CostOnly policy on a fresh ledger.
inline double costonly_capacity(const Graph& g, int src, int dst, Kernel k, const Elig& el) {
  const nc::NodeId s[1] = {src};
  const nc::NodeId t[1] = {dst};
  nc::CapacityOptions opts {config(nc::CostOnly {}, k, costonly_selection()), el.get(), nullptr, {}};
  return nc::as_float(nc::capacity(g.ref(), s, t, opts).placed, kResolution);
}

inline nc::FlowPolicyConfig global_config() {
  return nc::FlowPolicyConfig {nc::Global {}, nc::FlowPlacement::MaxFlow, nc::Admission::Lossless,
                               nc::PathSelection(nc::EdgeTieBreak::Deterministic, nc::EdgePriority::Capacity)};
}

// Exact maximum flow: capacity() with Global + MaxFlow, totals only or with
// every optional output (edge flows, residual, bottleneck with the Global min_cut).
inline std::uint64_t maxflow_exact(const Graph& g, int src, int dst, const Elig& el, bool detailed) {
  const nc::NodeId s[1] = {src};
  const nc::NodeId t[1] = {dst};
  nc::OutputOptions out;
  if (detailed) out = nc::OutputOptions {true, true, true};
  nc::CapacityOptions opts {global_config(), el.get(), nullptr, out};
  const auto r = nc::capacity(g.ref(), s, t, opts);
  Digest d;
  d.add_d(nc::as_float(r.placed, kResolution));
  d.add_i(r.min_cut ? static_cast<std::int64_t>(r.min_cut->size()) : -1);
  d.add_i(r.cost_distribution ? static_cast<std::int64_t>(r.cost_distribution->size()) : -1);
  if (detailed) {
    std::int64_t residual = 0;
    for (auto x : *r.residual) residual += x;
    d.add_i(residual);
    d.add_i(static_cast<std::int64_t>(r.bottleneck->size()));
  }
  return d.h;
}

inline std::uint64_t fold_sensitivity(const nc::SensitivityResult& r) {
  Digest d;
  d.add_d(nc::as_float(r.baseline.placed, kResolution));
  for (const auto& e : r.entries) {
    d.add_i(e.edge);
    d.add_d(nc::as_float(e.loss, kResolution));
  }
  return d.h;
}

// Exact-model sensitivity: Global + MaxFlow (candidates: the baseline's positive-flow edges).
inline std::uint64_t sensitivity_exact(const Graph& g, int src, int dst, const Elig& el) {
  const nc::NodeId s[1] = {src};
  const nc::NodeId t[1] = {dst};
  nc::SensitivityOptions opts {global_config(), el.get(), nullptr, {}};
  return fold_sensitivity(nc::sensitivity(g.ref(), s, t, opts, 1));
}

// Policy-model sensitivity over every allowed edge of at least 1 unit: CostOnly EqualBalanced.
inline std::uint64_t sensitivity_policy_all_edges(const Graph& g, int src, int dst, const Elig& el) {
  const nc::NodeId s[1] = {src};
  const nc::NodeId t[1] = {dst};
  const auto cfg = config(nc::CostOnly {}, Kernel::EqualBalanced, costonly_selection());
  nc::SensitivityOptions opts {cfg, el.get(), nullptr, {}};
  return fold_sensitivity(nc::sensitivity(g.ref(), s, t, opts, 1));
}

// Cohort placement: one TE_ECMP CapacityAware policy per pair on one shared ledger (NetGraph
// presets TE_ECMP_16_LSP: 16 single-path EB tunnels; TE_ECMP_UP_TO_256_LSP: 1..256),
// reoptimised, place(to_units(volume)).
inline std::uint64_t cohort_place(const Graph& g, const std::vector<std::pair<int, int>>& pairs, double volume,
                                  int min_flows, int max_flows) {
  nc::CapacityAware ca;
  ca.flow_count = min_flows;
  ca.max_flow_count = max_flows;
  ca.reoptimize = true;
  const nc::FlowPolicyConfig cfg {
      ca, nc::FlowPlacement::EqualBalanced, nc::Admission::Lossless,
      nc::PathSelection(nc::EdgeTieBreak::PreferHigherPriority, nc::EdgePriority::Residual, false, false)};
  nc::FlowLedger fg(g.g);
  std::vector<std::unique_ptr<nc::FlowPolicy>> policies;
  Digest d;
  const auto units = nc::to_units(volume, kResolution);
  for (auto [s, t] : pairs) {
    const nc::NodeId sv[1] = {s};
    const nc::NodeId tv[1] = {t};
    policies.push_back(std::make_unique<nc::FlowPolicy>(fg, sv, tv, cfg));
    const auto r = policies.back()->place(units);
    d.add_d(nc::as_float(r.placed, kResolution));
    d.add_i(static_cast<std::int64_t>(policies.back()->records().size()));
  }
  return d.h;
}

// One series of CostOnly EqualBalanced placements with one DagCache per series (the support of
// each destination set under the series Eligibility; NetGraph's hop-by-hop engine).
struct CostOnlySeries {
  const Graph& g;
  const Elig& el;
  nc::DagCache cache;
  CostOnlySeries(const Graph& graph, const Elig& elig) : g(graph), el(elig), cache(graph.ref()) {}
  // Places every demand on a fresh ledger; returns the placed amount per demand.
  std::vector<double> place_all(const std::vector<std::pair<int, int>>& demands, const std::vector<double>& volumes) {
    nc::FlowLedger fg(g.g);
    const auto cfg = config(nc::CostOnly {}, Kernel::EqualBalanced, costonly_selection());
    std::vector<double> placed(demands.size());
    for (std::size_t i = 0; i < demands.size(); ++i) {
      const nc::NodeId s[1] = {demands[i].first};
      const nc::NodeId t[1] = {demands[i].second};
      nc::FlowPolicy policy(fg, s, t, cfg, el.get(), &cache);
      placed[i] = nc::as_float(policy.place(nc::to_units(volumes[i], kResolution)).placed, kResolution);
    }
    return placed;
  }
};

inline double grid_resolution() { return kResolution; }

// DagCache series: lookups (roots {r}, Reverse, no eligibility) for every root; cold: a fresh
// cache per series (misses: one support search each); warm: the graph's cache, prepopulated by the
// first (untimed) iteration (hits).
inline std::uint64_t dag_cache_series(const Graph& g, const std::vector<int>& roots, bool warm) {
  std::shared_ptr<nc::DagCache> cache;
  if (warm) {
    auto& w = const_cast<Graph&>(g).warm_cache;
    if (!w) w = std::make_shared<nc::DagCache>(g.ref());
    cache = w;
  } else {
    cache = std::make_shared<nc::DagCache>(g.ref());
  }
  Digest d;
  for (int r : roots) {
    const nc::NodeId root[1] = {r};
    const auto v = cache->lookup(g.ref(), root, nc::Direction::Reverse, nullptr);
    d.add_i(static_cast<std::int64_t>(v->via_edges().size()));
  }
  d.add_i(static_cast<std::int64_t>(cache->stats().entries));
  return d.h;
}

// K shortest paths: Deterministic tie-break, Capacity priority, multipath and multi_edge.
inline std::uint64_t ksp(const Graph& g, int src, int dst, int k) {
  const nc::KspOptions opts(k, nc::PathSelection(nc::EdgeTieBreak::Deterministic, nc::EdgePriority::Capacity, true, true));
  const auto paths = nc::ksp(g.ref(), src, dst, opts);
  Digest d;
  d.add_i(static_cast<std::int64_t>(paths.size()));
  for (const auto& [cost, path] : paths) {
    d.add_i(cost);
    const auto& hops = path.hops();
    d.add_i(static_cast<std::int64_t>(hops.size()));
    for (const auto& hop : hops) {
      for (const auto e : hop) d.add_i(e);
    }
  }
  return d.h;
}

// ---- The storage rows' interface (README "Rows")

// Construction rows: the API's own input (capacities in units of kResolution), prepared untimed.
struct ConstructInput {
  const Arrays* a = nullptr;
  std::vector<nc::Amount> units;
};

inline ConstructInput construct_prepare(const Arrays& a) {
  ConstructInput in;
  in.a = &a;
  in.units.resize(a.cap.size());
  for (std::size_t i = 0; i < a.cap.size(); ++i) in.units[i] = nc::to_units(a.cap[i], kResolution);
  return in;
}

// The sampled construction digest (API-independent, see the file comment):
// N, E, max_cost, cost_gcd; then for the 16 nodes u_k = (k * 7919) mod N: the out-degree and the
// in-degree; for a non-empty out-row the input positions of its first and last edge in
// traversal order and the cost and physical capacity of the first; for a non-empty in-row the
// input positions of its first and last slot in reverse order.
inline std::uint64_t construct_sample(const nc::Graph& g) {
  Digest d;
  const std::int64_t n = g.num_nodes();
  d.add_i(n);
  d.add_i(g.num_edges());
  d.add_i(g.max_cost());
  d.add_i(g.cost_gcd());
  const auto row = g.row_offsets();
  const auto in_row = g.in_row_offsets();
  const auto in_adj = g.in_adj_edge_index();
  const auto input = g.input_index();
  for (std::int64_t k = 0; k < 16 && n > 0; ++k) {
    const auto u = static_cast<std::size_t>((k * 7919) % n);
    const std::int32_t b = row[u], e = row[u + 1], ib = in_row[u], ie = in_row[u + 1];
    d.add_i(e - b);
    d.add_i(ie - ib);
    if (e > b) {
      d.add_i(input[static_cast<std::size_t>(b)]);
      d.add_i(input[static_cast<std::size_t>(e - 1)]);
      d.add_i(g.cost()[static_cast<std::size_t>(b)]);
      d.add_d(nc::as_float(g.capacity_units()[static_cast<std::size_t>(b)], kResolution));
    }
    if (ie > ib) {
      d.add_i(input[static_cast<std::size_t>(in_adj[static_cast<std::size_t>(ib)])]);
      d.add_i(input[static_cast<std::size_t>(in_adj[static_cast<std::size_t>(ie - 1)])]);
    }
  }
  return d.h;
}

// Timed: the constructor on the prepared input, then the sampled digest.
inline std::uint64_t construct_run(const ConstructInput& in) {
  const auto& [n, tails, heads, caps, weights] = *in.a;  // the Arrays fields, in declaration order
  const nc::Graph g(n, tails, heads, in.units, weights, kResolution);
  return construct_sample(g);
}

// The queue a search with the default queue (Auto) takes on g: the production rule.
struct QueueChoice {
  std::int64_t width = 0;
  bool buckets = false;
};

inline QueueChoice queue_choice(const Graph& g) {
  return {nc::detail::bucket_width(g.ref()), nc::detail::uses_buckets(g.ref(), nc::SpfQueue::Auto)};
}

// The sampled search digest (API-independent, see the file comment): the
// number of DAG entries; the distances of the 16 nodes v_k = (k * 7919) mod N (-1 unreachable);
// for each of those nodes the number of its entries and the smallest and largest input position
// of their edges (an entry's node is its edge's head; the entries of a node are a set).
inline std::uint64_t search_sample(const Graph& g, const nc::SpfResult& r) {
  Digest d;
  const auto& [distances, dag, nearest] = r;  // costs, dag, target_cost
  const auto via = dag.via_edges();
  const auto dst = g.ref().edge_dst();
  const auto input = g.ref().input_index();
  const auto& costs = *distances;
  const std::int64_t n = static_cast<std::int64_t>(costs.size());
  d.add_i(static_cast<std::int64_t>(via.size()));
  for (std::int64_t k = 0; k < 16 && n > 0; ++k) {
    const auto v = static_cast<std::int32_t>((k * 7919) % n);
    const nc::Cost c = costs[static_cast<std::size_t>(v)];
    d.add_i(c == nc::kUnreachable ? -1 : c);
    // The entries are grouped by node in ascending order: the run of node v.
    const auto lo = std::partition_point(via.begin(), via.end(), [&](nc::EdgeId e) { return dst[static_cast<std::size_t>(e)] < v; });
    const auto hi = std::partition_point(lo, via.end(), [&](nc::EdgeId e) { return dst[static_cast<std::size_t>(e)] <= v; });
    std::int64_t lo_in = -1, hi_in = -1;
    for (auto it = lo; it != hi; ++it) {
      const std::int64_t x = input[static_cast<std::size_t>(*it)];
      lo_in = lo_in < 0 ? x : std::min(lo_in, x);
      hi_in = std::max(hi_in, x);
    }
    d.add_i(hi - lo);
    d.add_i(lo_in);
    d.add_i(hi_in);
  }
  return d.h;
}

// Timed: one full search (default selection: multipath, every parallel member) and its digest.
inline std::uint64_t search_run(const Graph& g, int root, bool reverse) {
  const nc::NodeId roots[1] = {root};
  return search_sample(g, run_spf(g, roots, reverse ? nc::Direction::Reverse : nc::Direction::Forward, std::nullopt,
                                  Sel {}, nullptr, nullptr));
}

}  // namespace abaapi
