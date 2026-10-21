// The A/B/A row adapter for the current Core API. The rows (aba_rows.hpp) call only `abaapi::`;
// this file implements each call with the public Core API, so the same rows run unchanged on any
// two builds that share this API (README.md: a build with another public API needs this file
// adapted first).
//
// Row mapping of the search rows:
// The route choice of a search is spelled in target terms (Sel: multipath, parallel_edges,
// higher_available_capacity as the HIGHER_AVAILABLE_CAPACITY tie_break); the default Sel is
// multipath and parallel_edges with the FEWER_HOPS tie_break. A supplied residual always gates its
// edges at min_residual = 1 unit (through to_units); capacity is no condition of a search (the
// row graphs have no edge below 1 unit, so this does not change their digests).
// * spf (spf.full.*, except the storage row spf.full.random1m, which runs search_run below):
//   spf(g, {src}, Forward, no targets), the default route choice, with or without the masks (as one
//   Eligibility) and the residual. Digest: 16 sampled costs, the DAG's entry count and its first
//   and last via.
// * spf_dag (spf.dag.*): spf(g, {src}, Forward, targets {dst}), the default route choice: the DAG of
//   the shortest walks that end in the target. Digest: target_cost, entries, first and last via.
// * spf_to (spf.to.random16k): spf(g, {dst}, Reverse, no targets), the default route choice.
// * targetset_* (spf.to.targetset_clos): the Clos graph and the destination set; the timed call is
//   a full reverse spf rooted at the destination set, targets absent.
// * spf_path (spf.single.canonical16k): spf(g, {src}, Forward, targets {dst}) with the row's
//   route choice (multipath = false) and the residual: the canonical single path. Digest:
//   target_cost and the walk.
// * ksp (ksp.fabric3_k8): ksp(g, src, dst, KspOptions(k)) (FEWER_HOPS, every parallel edge), no
//   residual, no ceiling, no masks. Digest: the path count, then per
//   path its cost, hop count and B' edge ids (the graph's edge ids, ordered by source, cost, target
//   and input position).
// The kernel rows place on fixed DAGs read with the public spf with targets (base_dag): Optimized
// as CentralizedTE restricted to the DAG (CentralizedTE(support = dag)), the other kernels as Pinned over a
// program of the one split route (from_routes per placement; Weighted with the capacity weights).
//
// The storage rows' interface (README "Rows"): construct_prepare and construct_run (the
// constructor on prepared unit capacities, then the sampled digest of construct_sample),
// queue_choice (the library's queue rule, read through the private detail::bucket_width and
// detail::uses_buckets: the one use of a non-public function here) and search_run (one full
// search plus the sampled digest of search_sample). The digests read only API-independent
// values: node-indexed counts and distances, costs, physical capacities and edges named by
// their input position, so they are the same for any API that builds the same graph.
#pragma once

#include "netgraph/core/shortest_path_cache.hpp"
#include "netgraph/core/evaluators.hpp"
#include "netgraph/core/flow_ledger.hpp"
#include "netgraph/core/flow_policy.hpp"
#include "netgraph/core/graph.hpp"
#include "netgraph/core/forwarding_program.hpp"
#include "netgraph/core/search.hpp"
#include "netgraph/core/detail/search_detail.hpp"

#include <algorithm>
#include <array>
#include <exception>
#include <stdexcept>
#include <thread>
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

// The graph plus the caller's inverse of input_positions() (input edge index -> EdgeId), and the
// warm ShortestPathCache of the cache.warm row (one per graph, prepopulated by the untimed iteration).
struct Graph {
  std::shared_ptr<const nc::Graph> g;
  std::vector<std::int32_t> edge_of_input;
  std::shared_ptr<nc::ShortestPathCache> warm_cache;
  const nc::Graph& ref() const { return *g; }
};

// Construction from physical arrays: to_units on every capacity, then the one constructor.
inline std::uint64_t construct(const Arrays& a) {
  auto g = graph_of(a);
  Digest d;
  d.add_i(g.num_edges());
  d.add_i(g.max_edge_cost());
  d.add_i(g.edge_cost_gcd());
  d.add_i(g.outgoing_offsets()[static_cast<std::size_t>(a.n / 2)]);
  return d.h;
}

inline Graph build(const Arrays& a) {
  Graph out;
  out.g = std::make_shared<const nc::Graph>(graph_of(a));
  out.edge_of_input.assign(a.src.size(), -1);
  auto ii = out.g->input_positions();
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

// A search's route choice, spelled in target terms (multipath: the DAG of every equal-cost walk).
struct Sel {
  bool multipath = true;
  bool parallel_edges = true;
  bool higher_available_capacity = false;
};

inline void choose(nc::SpfOptions& o, const Sel& s) {
  o.tie_break = s.higher_available_capacity ? nc::RouteTieBreak::HigherAvailableCapacity : nc::RouteTieBreak::FewerHops;
  o.multipath = s.multipath;
  o.parallel_edges = s.parallel_edges;
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
  nc::SpfOptions o(dir);
  choose(o, s);
  o.targets = std::move(targets);
  std::vector<nc::Amount> units;
  if (residual != nullptr) {
    units = residual_units(*residual);
    o.residual = std::span<const nc::Amount>(units);
  }
  o.eligibility = el ? el->get() : nullptr;
  return nc::spf(g.ref(), roots, o);
}

inline void fold_dag(Digest& d, const nc::PredecessorDAG& dag) {
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
  std::shared_ptr<const nc::PredecessorDAG> dag;
};
inline Dag base_dag(const Graph& g, int src, int dst) {
  const nc::NodeId roots[1] = {src};
  auto r = run_spf(g, roots, nc::Direction::Forward, std::vector<nc::NodeId> {dst}, Sel {}, nullptr, nullptr);
  return Dag {std::make_shared<const nc::PredecessorDAG>(std::move(r.dag))};
}

// The placement kernels the rows time; Lossy is Equal placement with Lossy admission.
enum class Kernel { Optimized, Equal, Lossy, Weighted };

struct KernelAxes {
  nc::RouteAllocation allocation = nc::RouteAllocation::Optimized;
  nc::Admission admission = nc::Admission::Lossless;
};

inline KernelAxes axes(Kernel k) {
  KernelAxes a;
  switch (k) {
    case Kernel::Optimized: break;
    case Kernel::Equal: a.allocation = nc::RouteAllocation::Equal; break;
    case Kernel::Lossy: a.allocation = nc::RouteAllocation::Equal; a.admission = nc::Admission::Lossy; break;
    case Kernel::Weighted: a.allocation = nc::RouteAllocation::Weighted; break;
  }
  return a;
}

// ShortestPaths `sp` with the kernel's allocation, and the kernel's admission.
inline nc::FlowPolicyConfig config(nc::ShortestPaths sp, Kernel k) {
  const auto a = axes(k);
  sp.allocation = a.allocation;
  return nc::FlowPolicyConfig {std::move(sp), a.admission};
}

// CentralizedTE with its source allocation, on `support` (nullptr: every eligible edge).
inline nc::CentralizedTE centralized(nc::SourceAllocation source_allocation,
                                     std::shared_ptr<const nc::PredecessorDAG> support = {}) {
  nc::CentralizedTE gl(source_allocation);
  gl.support = std::move(support);
  return gl;
}

// Manual placements on fixed DAGs into one ledger: one policy per placement on that one DAG
// (Optimized: CentralizedTE with the DAG as its support; the other kernels: Pinned over the program of the
// one split route), destroyed after the call (its records stay in the ledger as orphans).
struct Ledger {
  nc::FlowLedger fg;
  explicit Ledger(const Graph& g) : fg(g.g) {}
  double place(int src, int dst, const Dag& dag, double amount, Kernel k) {
    const auto a = axes(k);
    const nc::NodeId s[1] = {src};
    const nc::NodeId t[1] = {dst};
    nc::FlowPolicyConfig cfg;
    cfg.admission = a.admission;
    if (k == Kernel::Optimized) {
      cfg.routing = centralized(nc::SourceAllocation::Optimized, dag.dag);
    } else {
      nc::RouteOptions options;
      options.allocation = a.allocation;
      const nc::RouteDefinition routes[1] = {nc::RouteDefinition {src, dag.dag, 1}};
      cfg.routing = nc::Pinned {nc::ForwardingProgram::from_routes(*fg.graph(), s, t, routes, options), nullptr};
    }
    nc::FlowPolicy policy(fg, s, t, std::move(cfg));
    const auto r = policy.place(nc::to_units(amount, kResolution));
    return nc::as_float(r.call_delivered, kResolution);
  }
};

// Scalar ShortestPaths policy capacity: place(+inf) of a ShortestPaths policy on a fresh ledger.
inline double costonly_capacity(const Graph& g, int src, int dst, Kernel k, const Elig& el) {
  const nc::NodeId s[1] = {src};
  const nc::NodeId t[1] = {dst};
  nc::EvaluationOptions opts {config(nc::ShortestPaths {}, k), el.get(), nullptr, {}};
  return nc::as_float(nc::evaluate_demand(g.ref(), s, t, opts).call_delivered, kResolution);
}

inline nc::FlowPolicyConfig global_config() {
  return nc::FlowPolicyConfig {.routing = centralized(nc::SourceAllocation::Optimized)};
}

// Exact maximum flow: evaluate_demand() with CentralizedTE (optimized source allocation), totals only
// or with every optional output (edge flows, residual, saturated boundary edges, min_cut).
inline std::uint64_t maxflow_exact(const Graph& g, int src, int dst, const Elig& el, bool detailed) {
  const nc::NodeId s[1] = {src};
  const nc::NodeId t[1] = {dst};
  nc::OutputOptions out;
  if (detailed) out = nc::OutputOptions {true, true, true, true};
  nc::EvaluationOptions opts {global_config(), el.get(), nullptr, out};
  const auto r = nc::evaluate_demand(g.ref(), s, t, opts);
  Digest d;
  d.add_d(nc::as_float(r.call_delivered, kResolution));
  d.add_i(r.min_cut ? static_cast<std::int64_t>(r.min_cut->size()) : -1);
  d.add_i(r.marginal_cost_curve ? static_cast<std::int64_t>(r.marginal_cost_curve->size()) : -1);
  if (detailed) {
    std::int64_t residual = 0;
    for (auto x : *r.residual) residual += x;
    d.add_i(residual);
    d.add_i(static_cast<std::int64_t>(r.saturated_boundary_edges->size()));
  }
  return d.h;
}

inline std::uint64_t fold_edge_failures(const nc::SensitivityResult& r) {
  Digest d;
  d.add_d(nc::as_float(r.baseline.call_delivered, kResolution));
  for (const auto& e : r.entries) {
    d.add_i(e.edge);
    d.add_d(nc::as_float(e.delivery_reduction, kResolution));
  }
  return d.h;
}

// Exact-model failure study: CentralizedTE, optimized source allocation (candidates: the baseline's positive-flow edges).
inline std::uint64_t edge_failures_exact(const Graph& g, int src, int dst, const Elig& el) {
  const nc::NodeId s[1] = {src};
  const nc::NodeId t[1] = {dst};
  nc::EvaluationOptions opts {global_config(), el.get(), nullptr, {}};
  return fold_edge_failures(nc::evaluate_edge_failures(g.ref(), s, t, opts, 1));
}

// Policy-model failure study over every allowed edge: ShortestPaths Equal.
inline std::uint64_t edge_failures_policy_all_edges(const Graph& g, int src, int dst, const Elig& el) {
  const nc::NodeId s[1] = {src};
  const nc::NodeId t[1] = {dst};
  const auto cfg = config(nc::ShortestPaths {}, Kernel::Equal);
  nc::EvaluationOptions opts {cfg, el.get(), nullptr, {}};
  return fold_edge_failures(nc::evaluate_edge_failures(g.ref(), s, t, opts, 1));
}

// Pooled ShortestPaths placement: one Equal policy (Lossless, or Lossy) per demand of several
// sources to one destination on one shared ledger, place(to_units(offer)) each. The digest folds
// every call's delivered amount, its drop total and count, and every source row (source, admitted,
// delivered, source loss, dropped).
inline std::uint64_t pooled_place(const Graph& g, const std::vector<std::pair<std::vector<int>, int>>& demands,
                                  double offer, bool lossy) {
  const auto cfg = config(nc::ShortestPaths {}, lossy ? Kernel::Lossy : Kernel::Equal);
  nc::FlowLedger fg(g.g);
  std::vector<std::unique_ptr<nc::FlowPolicy>> policies;
  Digest d;
  const auto units = nc::to_units(offer, kResolution);
  for (const auto& [srcs, dst] : demands) {
    const std::vector<nc::NodeId> sv(srcs.begin(), srcs.end());
    const nc::NodeId tv[1] = {dst};
    policies.push_back(std::make_unique<nc::FlowPolicy>(fg, sv, tv, cfg));
    const auto r = policies.back()->place(units);
    d.add_i(r.call_delivered);
    nc::Amount drops = 0;
    for (const auto& [e, x] : r.call_edge_drops) drops += x;
    d.add_i(drops);
    d.add_i(static_cast<std::int64_t>(r.call_edge_drops.size()));
    for (const auto& row : r.by_source) {
      d.add_i(row.source);
      d.add(row.admitted.lo);
      d.add_i(row.delivered);
      d.add(row.source_drops.lo);
      d.add(row.total_edge_drops.lo);
    }
  }
  return d.h;
}

// Several demands at once: the pooled demands of pooled_place (ShortestPaths Equal Lossy, one
// ledger per iteration, the policies created in demand order) given a full offer each: as scalar
// replace(to_units(offer)) calls in demand order, or as one apply_actions list of the same Replace
// actions (ORDERED: the same work as the scalar calls, the same results; SIMULTANEOUS: forwarded at
// once). The digest folds every entry as pooled_place folds a call.
enum class Batch { Scalar, Ordered, Simultaneous };
inline std::uint64_t batch_place(const Graph& g, const std::vector<std::pair<std::vector<int>, int>>& demands, double offer, Batch mode) {
  const auto cfg = config(nc::ShortestPaths {}, Kernel::Lossy);
  nc::FlowLedger fg(g.g);
  std::vector<std::unique_ptr<nc::FlowPolicy>> policies;
  const auto units = nc::to_units(offer, kResolution);
  for (const auto& [srcs, dst] : demands) {
    const std::vector<nc::NodeId> sv(srcs.begin(), srcs.end());
    const nc::NodeId tv[1] = {dst};
    policies.push_back(std::make_unique<nc::FlowPolicy>(fg, sv, tv, cfg));
  }
  std::vector<nc::PlacementResult> results;
  if (mode == Batch::Scalar) {
    for (auto& p : policies) results.push_back(p->replace(units));
  } else {
    std::vector<nc::Action> actions;
    for (auto& p : policies) actions.push_back(nc::Replace {*p, units});
    results = nc::apply_actions(fg, actions, mode == Batch::Ordered ? nc::BatchMode::Ordered : nc::BatchMode::Simultaneous);
  }
  Digest d;
  for (const auto& r : results) {
    d.add_i(r.call_delivered);
    nc::Amount drops = 0;
    for (const auto& [e, x] : r.call_edge_drops) drops += x;
    d.add_i(drops);
    d.add_i(static_cast<std::int64_t>(r.call_edge_drops.size()));
    for (const auto& row : r.by_source) {
      d.add_i(row.source);
      d.add(row.admitted.lo);
      d.add_i(row.delivered);
      d.add(row.source_drops.lo);
      d.add(row.total_edge_drops.lo);
    }
  }
  return d.h;
}

// The ring: routers 0..63 with ring links i -> (i + 1) mod 64 of 3 units (EdgeIds 0..63), heads
// 64..127 with wide feeders head i -> router i, wide exits router i -> 128 (2^40 units each; built
// in units, not through a resolution). Tunnel i: its feeder, ring links i and i + 1, the exit of
// router i + 2. Every ring link carries two tunnels of 3 (load 6 > 3): 64 barriers in one cyclic
// group that settles at round 6, every head delivering 1. As one 64-head Pinned LOSSY program (one
// policy, replace(192): 3 per head) or as 64 one-tunnel Pinned LOSSY policies created in head
// order (one SIMULTANEOUS list of Replace(3) each), each on its own ledger, built once.
struct Ring {
  std::shared_ptr<const nc::Graph> g;
  std::unique_ptr<nc::FlowLedger> program_ledger;
  std::unique_ptr<nc::FlowLedger> list_ledger;
  std::unique_ptr<nc::FlowPolicy> program;
  std::vector<std::unique_ptr<nc::FlowPolicy>> policies;
};
inline std::shared_ptr<Ring> ring_fixture() {
  constexpr int n = 64;
  constexpr nc::Amount wide = nc::Amount {1} << 40;
  std::vector<std::int32_t> src;
  std::vector<std::int32_t> dst;
  std::vector<nc::Amount> units;
  for (int i = 0; i < n; ++i) src.push_back(i), dst.push_back((i + 1) % n), units.push_back(3);
  for (int i = 0; i < n; ++i) src.push_back(i), dst.push_back(2 * n), units.push_back(wide);
  for (int i = 0; i < n; ++i) src.push_back(n + i), dst.push_back(i), units.push_back(wide);
  auto ring = std::make_shared<Ring>();
  ring->g = std::make_shared<const nc::Graph>(2 * n + 1, src, dst, units, std::vector<std::int64_t>(src.size(), 1), kResolution);
  const auto inv = ring->g->input_positions();
  std::vector<std::int32_t> id(inv.size());
  for (std::size_t e = 0; e < inv.size(); ++e) id[static_cast<std::size_t>(inv[e])] = static_cast<std::int32_t>(e);
  auto tunnel = [&](int i) {
    const std::vector<std::vector<nc::EdgeId>> hops = {{id[static_cast<std::size_t>(2 * n + i)]}, {id[static_cast<std::size_t>(i)]},
                                                       {id[static_cast<std::size_t>((i + 1) % n)]}, {id[static_cast<std::size_t>(n + (i + 2) % n)]}};
    return nc::RouteDefinition {n + i, std::make_shared<const nc::Path>(*ring->g, hops), 1};
  };
  const nc::NodeId sink[1] = {2 * n};
  std::vector<nc::NodeId> heads;
  std::vector<nc::RouteDefinition> all;
  for (int i = 0; i < n; ++i) heads.push_back(n + i), all.push_back(tunnel(i));
  ring->program_ledger = std::make_unique<nc::FlowLedger>(ring->g);
  ring->list_ledger = std::make_unique<nc::FlowLedger>(ring->g);
  auto pinned = [](std::shared_ptr<const nc::ForwardingProgram> prog) {
    nc::FlowPolicyConfig cfg;
    cfg.routing = nc::Pinned {std::move(prog), nullptr};
    cfg.admission = nc::Admission::Lossy;
    return cfg;
  };
  ring->program = std::make_unique<nc::FlowPolicy>(*ring->program_ledger, heads, sink, pinned(nc::ForwardingProgram::from_routes(*ring->g, heads, sink, all)));
  for (int i = 0; i < n; ++i) {
    const nc::NodeId h[1] = {n + i};
    const std::vector<nc::RouteDefinition> one = {tunnel(i)};
    ring->policies.push_back(std::make_unique<nc::FlowPolicy>(*ring->list_ledger, h, sink, pinned(nc::ForwardingProgram::from_routes(*ring->g, h, sink, one))));
  }
  return ring;
}
// The digest of a ring forwarding: per head its delivery, then per ring link its carried load.
inline std::uint64_t ring_digest(const std::vector<nc::Amount>& delivered, std::span<const nc::Amount> carried) {
  Digest d;
  for (auto x : delivered) d.add_i(x);
  for (int e = 0; e < 64; ++e) d.add_i(carried[static_cast<std::size_t>(e)]);
  return d.h;
}
// One replace(192) of the ring's one program (a replace releases its own load: every iteration
// starts from the same state).
inline std::uint64_t ring_program(Ring& r) {
  const auto res = r.program->replace(192);
  std::vector<nc::Amount> delivered;
  for (const auto& t : res.by_tunnel) delivered.push_back(t.delivered);
  return ring_digest(delivered, r.program_ledger->carried());
}
// One SIMULTANEOUS list of the 64 one-tunnel policies, Replace(3) each.
inline std::uint64_t ring_list(Ring& r) {
  std::vector<nc::Action> actions;
  for (auto& p : r.policies) actions.push_back(nc::Replace {*p, 3});
  const auto res = nc::apply_actions(*r.list_ledger, actions, nc::BatchMode::Simultaneous);
  std::vector<nc::Amount> delivered;
  for (const auto& x : res) delivered.push_back(x.call_delivered);
  return ring_digest(delivered, r.list_ledger->carried());
}

// Pooled programs from supplied DAGs: per demand, ForwardingProgram::from_dag of its sources over
// the reverse shortest-path DAG to its destination (trimmed to the sources' walks: their selected
// supports, the DAG ShortestPaths Equal pools).
using PooledPrograms = std::vector<std::shared_ptr<const nc::ForwardingProgram>>;
inline PooledPrograms pooled_dag_programs(const Graph& g, const std::vector<std::pair<std::vector<int>, int>>& demands) {
  PooledPrograms out;
  for (const auto& [srcs, dst] : demands) {
    const std::vector<nc::NodeId> sv(srcs.begin(), srcs.end());
    const nc::NodeId tv[1] = {dst};
    auto dag = std::make_shared<const nc::PredecessorDAG>(run_spf(g, tv, nc::Direction::Reverse, std::nullopt, Sel {}, nullptr, nullptr).dag);
    out.push_back(nc::ForwardingProgram::from_dag(g.ref(), sv, tv, dag));
  }
  return out;
}

// Building the pooled programs: the digest folds each program's roster size and its Lossless bound
// is not computed (building only).
inline std::uint64_t pooled_dag_build(const Graph& g, const std::vector<std::pair<std::vector<int>, int>>& demands) {
  Digest d;
  for (const auto& p : pooled_dag_programs(g, demands)) {
    d.add_i(static_cast<std::int64_t>(p->sources().size()));
    d.add_i(static_cast<std::int64_t>(p->tunnels().size()));
  }
  return d.h;
}

// Pinned over prebuilt pooled programs: one policy per demand on one shared ledger,
// place(to_units(offer)) each; the digest as pooled_place's.
inline std::uint64_t pooled_pinned_place(const Graph& g, const std::vector<std::pair<std::vector<int>, int>>& demands,
                                         const PooledPrograms& programs, double offer, bool lossy) {
  nc::FlowLedger fg(g.g);
  std::vector<std::unique_ptr<nc::FlowPolicy>> policies;
  Digest d;
  const auto units = nc::to_units(offer, kResolution);
  for (std::size_t i = 0; i < demands.size(); ++i) {
    const auto& [srcs, dst] = demands[i];
    const std::vector<nc::NodeId> sv(srcs.begin(), srcs.end());
    const nc::NodeId tv[1] = {dst};
    const nc::FlowPolicyConfig cfg {.routing = nc::Pinned {programs[i], nullptr},
                                    .admission = lossy ? nc::Admission::Lossy : nc::Admission::Lossless};
    policies.push_back(std::make_unique<nc::FlowPolicy>(fg, sv, tv, cfg));
    const auto r = policies.back()->place(units);
    d.add_i(r.call_delivered);
    nc::Amount drops = 0;
    for (const auto& [e, x] : r.call_edge_drops) drops += x;
    d.add_i(drops);
    d.add_i(static_cast<std::int64_t>(r.call_edge_drops.size()));
    for (const auto& row : r.by_source) {
      d.add_i(row.source);
      d.add(row.admitted.lo);
      d.add_i(row.delivered);
      d.add(row.source_drops.lo);
      d.add(row.total_edge_drops.lo);
    }
  }
  return d.h;
}

// ---- forwarding programs

// One demand's program: the optimized CentralizedTE solve of `offer` from `srcs` to `dst` on a fresh
// ledger (its paths with their solved amounts, the captured source weights).
struct Program {
  std::shared_ptr<const nc::ForwardingProgram> p;
  std::vector<nc::NodeId> srcs;
  nc::NodeId dst {0};
};
inline Program global_program(const Graph& g, const std::vector<int>& srcs, int dst, double offer) {
  Program out;
  out.srcs.assign(srcs.begin(), srcs.end());
  out.dst = dst;
  nc::FlowLedger fg(g.g);
  const nc::NodeId tv[1] = {dst};
  nc::FlowPolicy policy(fg, out.srcs, tv, global_config());
  out.p = policy.replace(nc::to_units(offer, kResolution)).program;
  return out;
}

// The failure states as shared Eligibility values (a Pinned black-hole set is held by reference).
using Masks = std::vector<std::shared_ptr<const nc::Eligibility>>;
inline Masks masks(const Graph& g, const std::vector<std::vector<int>>& states) {
  Masks out;
  for (const auto& st : states) out.push_back(std::make_shared<const nc::Eligibility>(*eligibility(g, {}, st).el));
  return out;
}

// A status as a fixed code, so a digest does not move when the enum gains or reorders values:
// Up 0, Failed 1, Blackholed 2 (the codes of every older Core), InsufficientBandwidth 3,
// CostCeilingExceeded 4.
inline std::int64_t status_code(nc::TunnelStatus s) {
  switch (s) {
    case nc::TunnelStatus::Up: return 0;
    case nc::TunnelStatus::Failed: return 1;
    case nc::TunnelStatus::Blackholed: return 2;
    case nc::TunnelStatus::InsufficientBandwidth: return 3;
    case nc::TunnelStatus::CostCeilingExceeded: return 4;
  }
  return -1;
}

inline void add_rows(Digest& d, const nc::PlacementResult& r) {
  d.add_i(r.call_delivered);
  d.add(r.source_drops.lo);
  nc::Amount drops = 0;
  for (const auto& [e, x] : r.call_edge_drops) drops += x;
  d.add_i(drops);
  for (const auto& t : r.by_tunnel) {
    d.add_i(status_code(t.status));
    d.add(t.assigned.lo);
    d.add_i(t.delivered);
  }
}

// An installed program's policy: forward(offer) for each offer (no route computation). The policy
// and its ledger live for the whole row; each call retires and creates its records.
struct Installed {
  nc::FlowLedger fg;
  std::unique_ptr<nc::FlowPolicy> policy;
  Installed(const Graph& g, const std::vector<int>& srcs, int dst, double offer) : fg(g.g) {
    const std::vector<nc::NodeId> sv(srcs.begin(), srcs.end());
    const nc::NodeId tv[1] = {dst};
    policy = std::make_unique<nc::FlowPolicy>(fg, sv, tv, global_config());
    (void)policy->replace(nc::to_units(offer, kResolution));
  }
  std::uint64_t forward_all(const std::vector<double>& offers) {
    Digest d;
    for (double o : offers) add_rows(d, policy->forward(nc::to_units(o, kResolution)));
    return d.h;
  }
};

// Pinned Lossy over a program under each failure state, each on a fresh ledger: the state as the
// policy's mask (replay: failed tunnels, the heads split over the survivors) or as its black holes
// (stale forwarding: the program unchanged, traffic dropped at the failures); replace(offer).
inline std::uint64_t pinned_failures(const Graph& g, const Program& prog, const Masks& states, double offer,
                                     bool stale) {
  Digest d;
  const nc::NodeId tv[1] = {prog.dst};
  const auto units = nc::to_units(offer, kResolution);
  for (const auto& st : states) {
    nc::FlowLedger fg(g.g);
    nc::FlowPolicyConfig cfg;
    cfg.routing = nc::Pinned {prog.p, stale ? st : nullptr};
    cfg.admission = nc::Admission::Lossy;
    nc::FlowPolicy policy(fg, prog.srcs, tv, std::move(cfg), stale ? nullptr : st.get());
    add_rows(d, policy.replace(units));
  }
  return d.h;
}

// A ShortestPaths Equal demand replaced while an outside ShortestPaths demand changes its load: for
// each load, the outside policy replaces its offer, then the demand replace(offer) re-prepares its
// program on the new released view and retires its records. One ledger for the row's iteration.
inline std::uint64_t replace_under_load(const Graph& g, const std::vector<int>& srcs, int dst,
                                        const std::vector<int>& outside_srcs, int outside_dst,
                                        const std::vector<double>& loads, double offer) {
  const auto cfg = config(nc::ShortestPaths {}, Kernel::Equal);
  nc::FlowLedger fg(g.g);
  const std::vector<nc::NodeId> sv(srcs.begin(), srcs.end());
  const std::vector<nc::NodeId> ov(outside_srcs.begin(), outside_srcs.end());
  const nc::NodeId tv[1] = {dst};
  const nc::NodeId otv[1] = {outside_dst};
  nc::FlowPolicy demand(fg, sv, tv, cfg);
  nc::FlowPolicy outside(fg, ov, otv, cfg);
  Digest d;
  const auto units = nc::to_units(offer, kResolution);
  for (double load : loads) {
    (void)outside.replace(nc::to_units(load, kResolution));
    const auto r = demand.replace(units);
    d.add_i(r.call_delivered);
    d.add_i(r.delivery_delta);
    d.add_i(r.routes_changed ? 1 : 0);
    d.add_i(r.weights_changed ? 1 : 0);
  }
  return d.h;
}

// The equal-source-allocation CentralizedTE configuration (every roster source sends its round-robin share).
inline nc::FlowPolicyConfig global_equal_config() {
  return nc::FlowPolicyConfig {.routing = centralized(nc::SourceAllocation::Equal)};
}

// evaluate_demand() of equal-source-allocation CentralizedTE from several sources to one destination. The digest folds the
// delivered amount, the exact total cost (both halves) and every source row's admitted amount.
inline std::uint64_t global_equal_capacity(const Graph& g, const std::vector<int>& srcs, int dst) {
  const std::vector<nc::NodeId> sv(srcs.begin(), srcs.end());
  const nc::NodeId tv[1] = {dst};
  const auto r = nc::evaluate_demand(g.ref(), sv, tv, nc::EvaluationOptions {global_equal_config(), nullptr, nullptr, {}});
  Digest d;
  d.add_i(r.call_delivered);
  d.add(r.total_cost->hi);
  d.add(r.total_cost->lo);
  for (const auto& row : r.by_source) {
    d.add_i(row.source);
    d.add(row.admitted.lo);
  }
  return d.h;
}

// evaluate_edge_failures() of equal-source-allocation CentralizedTE (every allowed edge is a candidate, each solved afresh).
inline std::uint64_t edge_failures_equal_global(const Graph& g, const std::vector<int>& srcs, int dst) {
  const std::vector<nc::NodeId> sv(srcs.begin(), srcs.end());
  const nc::NodeId tv[1] = {dst};
  return fold_edge_failures(nc::evaluate_edge_failures(g.ref(), sv, tv, nc::EvaluationOptions {global_equal_config()}, 1));
}

// T1: one DistributedTE(TunnelCount(tunnels)) policy per pair on one shared ledger, Lossless, each head
// setting up its tunnels by CSPF on the unreserved view (the earlier pairs' holds count).
struct T1Fabric {
  nc::FlowLedger fg;
  std::vector<std::unique_ptr<nc::FlowPolicy>> policies;
  T1Fabric(const Graph& g, const std::vector<std::pair<int, int>>& pairs, int tunnels, nc::Admission admission)
      : fg(g.g) {
    const nc::FlowPolicyConfig cfg {.routing = nc::DistributedTE {nc::TunnelCount {tunnels}}, .admission = admission};
    for (auto [s, t] : pairs) {
      const nc::NodeId sv[1] = {s};
      const nc::NodeId tv[1] = {t};
      policies.push_back(std::make_unique<nc::FlowPolicy>(fg, sv, tv, cfg));
    }
  }
};

// The digest of a T1 call: its delivery, its up tunnels and its holds in all.
inline void add_t1(Digest& d, const nc::PlacementResult& r) {
  d.add_i(r.call_delivered);
  std::int64_t up = 0;
  for (const auto& t : r.program->tunnels()) up += t.status == nc::TunnelStatus::Up ? 1 : 0;
  d.add_i(up);
  d.add_i(r.plan->total_reserved_bandwidth);
}

// t1.setup.fabric: every pair's replace(to_units(offer), design = the same) on a fresh ledger, in
// pair order: CSPF setup, the holds' commit and the Lossless forwarding.
inline std::uint64_t t1_setup(const Graph& g, const std::vector<std::pair<int, int>>& pairs, double offer, int tunnels) {
  T1Fabric f(g, pairs, tunnels, nc::Admission::Lossless);
  Digest d;
  const auto units = nc::to_units(offer, kResolution);
  for (auto& p : f.policies) add_t1(d, p->replace(units, units));
  return d.h;
}

// The fabric of t1.forward.fabric: Lossy, every pair set up once by replace(to_units(offer),
// design = the same), outside the timed iterations.
inline std::shared_ptr<T1Fabric> t1_installed(const Graph& g, const std::vector<std::pair<int, int>>& pairs, double offer,
                                              int tunnels) {
  auto f = std::make_shared<T1Fabric>(g, pairs, tunnels, nc::Admission::Lossy);
  const auto units = nc::to_units(offer, kResolution);
  for (auto& p : f->policies) (void)p->replace(units, units);
  return f;
}

// t1.forward.fabric: the installed Lossy programs of every pair forwarded under each offer (no
// setup: forward keeps the tunnels and the holds).
inline std::uint64_t t1_forward(T1Fabric& f, const std::vector<double>& offers) {
  Digest d;
  for (double offer : offers) {
    const auto units = nc::to_units(offer, kResolution);
    for (auto& p : f.policies) add_t1(d, p->forward(units));
  }
  return d.h;
}

// T2: a CentralizedTE(source_allocation=EQUAL, tunnel_count=TunnelCount(count, up_to),
// prefer_headroom=headroom) configuration, Lossless.
inline nc::FlowPolicyConfig t2_config(int count, bool up_to, bool headroom) {
  auto t2 = centralized(nc::SourceAllocation::Equal);
  t2.tunnel_count = nc::TunnelCount {count, up_to};
  t2.prefer_headroom = headroom;
  return {.routing = std::move(t2)};
}

// The digest of a T2 call: its delivery, the chosen count, U, u, F, the route cost and the
// tunnels.
inline void add_t2(Digest& d, const nc::PlacementResult& r) {
  d.add_i(r.call_delivered);
  d.add_i(r.plan->tunnels_per_head);
  d.add_i(*r.plan->max_bandwidth_per_tunnel);
  d.add_i(*r.plan->planned_bandwidth_per_tunnel);
  d.add_i(*r.plan->max_lossless_offer);
  d.add(r.plan->total_tunnel_cost->lo);
  d.add_i(static_cast<std::int64_t>(r.program->tunnels().size()));
}

// t2.exact16.*, t2.upto256.headroom: one T2 policy per pair on one shared ledger, in pair order,
// each replace(to_units(offer), design = the same) on the residual the earlier pairs left: live
// heads, sizing (the quota flows), s*, the canonical augmentation (with Phi under Headroom), the
// cut into tunnels, the lossless bound and the forwarding.
inline std::uint64_t t2_fabric(const Graph& g, const std::vector<std::pair<int, int>>& pairs, double offer, int count,
                               bool up_to, bool headroom) {
  nc::FlowLedger fg(g.g);
  const auto cfg = t2_config(count, up_to, headroom);
  std::vector<std::unique_ptr<nc::FlowPolicy>> policies;
  Digest d;
  const auto units = nc::to_units(offer, kResolution);
  for (auto [s, t] : pairs) {
    const nc::NodeId sv[1] = {s};
    const nc::NodeId tv[1] = {t};
    policies.push_back(std::make_unique<nc::FlowPolicy>(fg, sv, tv, cfg));
    add_t2(d, policies.back()->replace(units, units));
  }
  return d.h;
}

// t2.size.mesh: one T2 policy of the heads to the destination on a ledger whose residual is
// `free` (input order, units), replace(+inf): an unlimited design, u = U.
inline std::uint64_t t2_mesh(const Graph& g, const std::vector<std::int64_t>& free, const std::vector<int>& heads,
                             int destination, int count, bool up_to, bool headroom) {
  std::vector<nc::Amount> residual(free.size());
  for (std::size_t i = 0; i < free.size(); ++i) residual[static_cast<std::size_t>(g.edge_of_input[i])] = free[i];
  nc::FlowLedger fg(g.g, residual);
  std::vector<nc::NodeId> sv(heads.begin(), heads.end());
  const nc::NodeId tv[1] = {destination};
  nc::FlowPolicy p(fg, sv, tv, t2_config(count, up_to, headroom));
  Digest d;
  add_t2(d, p.replace(std::nullopt));
  return d.h;
}

// One series of ShortestPaths Equal placements with one ShortestPathCache per series (the support of
// each destination set under the series Eligibility; NetGraph's hop-by-hop engine).
struct CostOnlySeries {
  const Graph& g;
  const Elig& el;
  nc::ShortestPathCache cache;
  CostOnlySeries(const Graph& graph, const Elig& elig) : g(graph), el(elig), cache(graph.ref()) {}
  // Places every demand on a fresh ledger; returns the delivered amount per demand.
  std::vector<double> place_all(const std::vector<std::pair<int, int>>& demands, const std::vector<double>& offers) {
    nc::FlowLedger fg(g.g);
    const auto cfg = config(nc::ShortestPaths {}, Kernel::Equal);
    std::vector<double> delivered(demands.size());
    for (std::size_t i = 0; i < demands.size(); ++i) {
      const nc::NodeId s[1] = {demands[i].first};
      const nc::NodeId t[1] = {demands[i].second};
      nc::FlowPolicy policy(fg, s, t, cfg, el.get(), &cache);
      delivered[i] = nc::as_float(policy.place(nc::to_units(offers[i], kResolution)).call_delivered, kResolution);
    }
    return delivered;
  }
};

inline double grid_resolution() { return kResolution; }

// ShortestPathCache series: lookups (roots {r}, Reverse, no eligibility) for every root; cold: a fresh
// cache per series (misses: one support search each); warm: the graph's cache, prepopulated by the
// first (untimed) iteration (hits).
inline std::uint64_t dag_cache_series(const Graph& g, const std::vector<int>& roots, bool warm) {
  std::shared_ptr<nc::ShortestPathCache> cache;
  if (warm) {
    auto& w = const_cast<Graph&>(g).warm_cache;
    if (!w) w = std::make_shared<nc::ShortestPathCache>(g.ref());
    cache = w;
  } else {
    cache = std::make_shared<nc::ShortestPathCache>(g.ref());
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

// One operation per row, 10,000 calls per iteration. Construction and preparation
// are untimed. The immutable config getter is the control for the gated getters.
enum class TinyOperation { PlaceZero, Delivered, Program, Config };
struct TinyCalls {
  std::shared_ptr<const nc::Graph> graph;
  nc::FlowLedger ledger;
  nc::FlowPolicy policy;
  TinyCalls()
      : graph(std::make_shared<const nc::Graph>(2, std::vector<nc::NodeId>{0}, std::vector<nc::NodeId>{1},
            std::vector<nc::Amount>{16}, std::vector<nc::Cost>{3}, 1.0)),
        ledger(graph), policy(ledger, std::vector<nc::NodeId>{0}, std::vector<nc::NodeId>{1}, {}) {
    const auto r = policy.replace(4);
    if (r.delivered != 4 || !r.program || r.by_source.size() != 1 || r.by_source[0].delivered != 4 ||
        ledger.residual() != std::vector<nc::Amount>{12}) throw std::runtime_error("tiny fixture preparation failed");
    const auto zero = policy.place(0);
    if (zero.call_delivered != 0 || zero.delivered != 4 || zero.requested != nc::Counter128{0, 4} ||
        zero.routes_changed || zero.weights_changed || zero.reservations_changed || zero.program != r.program ||
        zero.edge_flows || zero.residual || zero.plan) throw std::runtime_error("tiny zero-call result differs");
  }
  std::uint64_t run(TinyOperation operation) {
    std::uint64_t sum = 0;
    for (int i = 0; i < 10000; ++i) {
      switch (operation) {
        case TinyOperation::PlaceZero: sum += static_cast<std::uint64_t>(policy.place(0).delivered); break;
        case TinyOperation::Delivered: sum += static_cast<std::uint64_t>(policy.delivered()); break;
        case TinyOperation::Program: sum += policy.program() != nullptr ? 1U : 0U; break;
        case TinyOperation::Config: sum += policy.config().admission == nc::Admission::Lossless ? 1U : 0U; break;
      }
    }
    return sum;
  }
};

// Fixed layered graph, 32 heads, 32 destinations, 12 layers, eight outgoing links
// per nonterminal node. Offer 32 sends one unit per head; no link can receive more
// than 32, below every capacity 64, so every call must deliver all 32 units.
struct ParallelLedgers {
  std::shared_ptr<const nc::Graph> graph;
  std::array<std::unique_ptr<nc::FlowLedger>, 2> ledgers;
  std::array<std::unique_ptr<nc::FlowPolicy>, 2> policies;
  ParallelLedgers() {
    std::vector<nc::NodeId> source, destination;
    std::vector<nc::Amount> capacity;
    std::vector<nc::Cost> cost;
    constexpr std::array<int, 8> offsets {0, 1, 3, 7, 11, 17, 23, 29};
    for (int layer = 0; layer < 11; ++layer) for (int head = 0; head < 32; ++head) for (int offset : offsets) {
      source.push_back(layer * 32 + head);
      destination.push_back((layer + 1) * 32 + (head + offset) % 32);
      capacity.push_back(64);
      cost.push_back(1 + (head * 17 + offset * 13 + layer) % 7);
    }
    graph = std::make_shared<const nc::Graph>(384, source, destination, capacity, cost, 1.0);
    std::vector<nc::NodeId> heads, targets;
    for (int i = 0; i < 32; ++i) { heads.push_back(i); targets.push_back(352 + i); }
    for (std::size_t i = 0; i < 2; ++i) {
      ledgers[i] = std::make_unique<nc::FlowLedger>(graph);
      policies[i] = std::make_unique<nc::FlowPolicy>(*ledgers[i], heads, targets, nc::FlowPolicyConfig{});
      const auto r = policies[i]->replace(32);
      if (r.delivered != 32 || r.by_source.size() != 32) throw std::runtime_error("parallel fixture preparation failed");
      for (const auto& row : r.by_source) if (row.delivered != 1) throw std::runtime_error("parallel fixture source delivery differs");
    }
  }
  std::uint64_t run(bool parallel) {
    std::array<std::uint64_t, 2> sums {};
    std::array<std::exception_ptr, 2> errors {};
    auto worker = [&](std::size_t i) noexcept {
      try {
        for (int call = 0; call < 32; ++call) {
          const auto r = policies[i]->evaluate(32);
          if (r.delivered != 32 || r.call_delivered != 32 || r.delivery_delta != 0 || r.by_source.size() != 32) {
            throw std::runtime_error("parallel evaluation differs");
          }
          for (const auto& row : r.by_source) if (row.delivered != 1) throw std::runtime_error("parallel source delivery differs");
          sums[i] += static_cast<std::uint64_t>(r.delivered);
        }
      } catch (...) { errors[i] = std::current_exception(); }
    };
    std::array<std::thread, 2> threads;
    try {
      for (std::size_t i = 0; i < 2; ++i) {
        threads[i] = std::thread(worker, i);
        if (!parallel) threads[i].join();
      }
    } catch (...) {
      for (auto& thread : threads) if (thread.joinable()) thread.join();
      throw;
    }
    for (auto& thread : threads) if (thread.joinable()) thread.join();
    for (const auto& error : errors) if (error) std::rethrow_exception(error);
    return sums[0] + sums[1];
  }
};

// K shortest paths: the FEWER_HOPS tie_break, every parallel edge.
inline std::uint64_t ksp(const Graph& g, int src, int dst, int k) {
  const nc::KspOptions opts(k);
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
// N, E, max_edge_cost, edge_cost_gcd; then for the 16 nodes u_k = (k * 7919) mod N: the out-degree and the
// in-degree; for a non-empty out-row the input positions of its first and last edge in
// traversal order and the cost and physical capacity of the first; for a non-empty in-row the
// input positions of its first and last slot in reverse order.
inline std::uint64_t construct_sample(const nc::Graph& g) {
  Digest d;
  const std::int64_t n = g.num_nodes();
  d.add_i(n);
  d.add_i(g.num_edges());
  d.add_i(g.max_edge_cost());
  d.add_i(g.edge_cost_gcd());
  const auto row = g.outgoing_offsets();
  const auto in_row = g.incoming_offsets();
  const auto in_adj = g.incoming_edges();
  const auto input = g.input_positions();
  for (std::int64_t k = 0; k < 16 && n > 0; ++k) {
    const auto u = static_cast<std::size_t>((k * 7919) % n);
    const std::int32_t b = row[u], e = row[u + 1], ib = in_row[u], ie = in_row[u + 1];
    d.add_i(e - b);
    d.add_i(ie - ib);
    if (e > b) {
      d.add_i(input[static_cast<std::size_t>(b)]);
      d.add_i(input[static_cast<std::size_t>(e - 1)]);
      d.add_i(g.edge_costs()[static_cast<std::size_t>(b)]);
      d.add_d(nc::as_float(g.edge_capacity_units()[static_cast<std::size_t>(b)], kResolution));
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
  const auto dst = g.ref().edge_destinations();
  const auto input = g.ref().input_positions();
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

// Timed: one full search (default route choice: multipath, every parallel member) and its digest.
inline std::uint64_t search_run(const Graph& g, int root, bool reverse) {
  const nc::NodeId roots[1] = {root};
  return search_sample(g, run_spf(g, roots, reverse ? nc::Direction::Reverse : nc::Direction::Forward, std::nullopt,
                                  Sel {}, nullptr, nullptr));
}

}  // namespace abaapi
