// Portable CPU harness: bit-exact parity of queue/workspace variants against the
// production shortest_paths, plus sequential single-thread timing.
//
//   harness correctness                      -> parity matrix, exit 1 on mismatch
//   harness time <variant> [samples] [maxW]   -> CSV timing rows on stdout
//   harness sweep <variant> [samples] [maxW]  -> cost-range sweep rows
//
// variant = ref | heap-fresh | heap-full | heap-sparse | bucket-fresh |
//           bucket-full | bucket-sparse | bsort-sparse | bsort-fresh
#include "spf_variants.hpp"

#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <memory>
#include <random>
#include <string>
#include <sys/resource.h>

using namespace netgraph::core;
using spfx::INF;
using Result = std::pair<std::vector<Cost>, PredDAG>;
using Clock = std::chrono::steady_clock;

static double ms(Clock::time_point start) {
  return std::chrono::duration<double, std::milli>(Clock::now() - start).count();
}

struct Case {
  std::string name;
  StrictMultiDiGraph graph;
  std::unique_ptr<bool[]> node_mask, edge_mask;
  std::vector<Cap> residual;
  NodeId fixed_src = -1, fixed_dst = -1;
  std::span<const bool> nm() const { return {node_mask.get(), size_t(graph.num_nodes())}; }
  std::span<const bool> em() const { return {edge_mask.get(), size_t(graph.num_edges())}; }
};

// Verbatim port of the research generator (same RNG draws), plus new kinds:
//   uniform  : random topology, every cost == 31
//   wideK    : random topology, costs 1..K (K via max_cost argument)
//   pseudo   : 2-tier fabric with a zero-cost pseudo source/sink attachment, as
//              NetGraph's max-flow builds; fixed src/dst = the pseudo nodes
//   zerocyc  : random graph with a share of zero-cost edges (zero-cost cycles)
static Case make_case(std::string kind, int n, int degree = 8, bool masked = false,
                      bool large = false, unsigned seed = 1729, Cost max_cost = 31) {
  std::vector<NodeId> src, dst;
  std::vector<Cost> costs;
  std::vector<Cap> caps;
  std::mt19937 rng(seed);
  auto add = [&](int u, int v, Cost c) {
    src.push_back(u); dst.push_back(v);
    costs.push_back(c + (large ? (Cost(1) << 34) : 0));
    caps.push_back(1.0 + double(rng() % 100) / 7.0);
  };
  NodeId fsrc = -1, fdst = -1;
  if (kind == "chain") {
    for (int i = 0; i < n - 1; ++i) add(i, i + 1, 1 + rng() % 31);
  } else if (kind == "grid") {
    int side = int(std::sqrt(n));
    if (side * side != n) throw std::runtime_error("grid must be square");
    for (int u = 0; u < n; ++u) {
      if (u % side + 1 < side) { add(u, u + 1, 1); add(u + 1, u, 1); }
      if (u + side < n) { add(u, u + side, 1); add(u + side, u, 1); }
    }
  } else if (kind == "fabric" || kind == "clos") {
    int spines = kind == "clos" ? n / 2 : 32;
    for (int u = 0; u < n - spines; ++u) for (int v = n - spines; v < n; ++v) {
      add(u, v, 1); add(v, u, 1);
    }
  } else if (kind == "uniform") {
    for (int u = 0; u < n; ++u) {
      add(u, (u + 1) % n, 31);
      for (int k = 1; k < degree; ++k) add(u, rng() % n, 31);
    }
  } else if (kind == "wide") {
    for (int u = 0; u < n; ++u) {
      add(u, (u + 1) % n, 1 + Cost(rng() % max_cost));
      for (int k = 1; k < degree; ++k) add(u, rng() % n, 1 + Cost(rng() % max_cost));
    }
  } else if (kind == "zerocyc") {
    for (int u = 0; u < n; ++u) {
      add(u, (u + 1) % n, rng() % 4 == 0 ? 0 : 1 + rng() % 31);
      for (int k = 1; k < degree; ++k) add(u, rng() % n, rng() % 4 == 0 ? 0 : 1 + rng() % 31);
    }
  } else if (kind == "pseudo") {
    // n real nodes: leaves [0, n/2), spines [n/2, n); then pseudo src n, pseudo dst n+1.
    int leaves = n / 2;
    for (int u = 0; u < leaves; ++u) for (int v = leaves; v < n; ++v) { add(u, v, 1); add(v, u, 1); }
    int k = std::max(1, leaves / 4);
    for (int u = 0; u < k; ++u) add(n, u, 0);
    for (int u = leaves - k; u < leaves; ++u) add(u, n + 1, 0);
    fsrc = n; fdst = n + 1; n += 2;
  } else {
    for (int u = 0; u < n; ++u) {
      add(u, (u + 1) % n, 1 + rng() % 31);
      for (int k = 1; k < degree; ++k) add(u, rng() % n, 1 + rng() % 31);
    }
  }
  Case c{kind + (masked ? "_masked" : "") + (large ? "_int64" : ""),
         StrictMultiDiGraph::from_arrays(n, src, dst, caps, costs)};
  c.fixed_src = fsrc; c.fixed_dst = fdst;
  int e = c.graph.num_edges();
  c.node_mask = std::make_unique<bool[]>(n);
  c.edge_mask = std::make_unique<bool[]>(e);
  c.residual.assign(c.graph.capacity_view().begin(), c.graph.capacity_view().end());
  for (int i = 0; i < n; ++i) c.node_mask[i] = !masked || rng() % 23 != 0;
  for (int i = 0; i < e; ++i) {
    c.edge_mask[i] = !masked || rng() % 11 != 0;
    if (masked && rng() % 13 == 0) c.residual[i] = 0;
  }
  if (fsrc >= 0) { c.node_mask[fsrc] = true; c.node_mask[fdst] = true; }
  return c;
}

struct Query {
  NodeId src; std::optional<NodeId> dst; bool multipath; EdgeSelection sel;
  bool use_residual; bool use_masks; bool reverse = false;
};

static Result run_ref(const Case& c, const Query& q) {
  auto res = q.use_residual ? std::span<const Cap>(c.residual) : std::span<const Cap>{};
  auto nm = q.use_masks ? c.nm() : std::span<const bool>{};
  auto em = q.use_masks ? c.em() : std::span<const bool>{};
  if (q.reverse) return shortest_paths_to(c.graph, q.src, q.multipath, q.sel, res, nm, em, {});
  return shortest_paths(c.graph, q.src, q.dst, q.multipath, q.sel, res, nm, em);
}

struct Runner {
  std::string variant;
  spfx::HeapQueue heap; spfx::BucketHeapQueue bucket; spfx::BucketPow2Queue bucket2; spfx::BucketSortQueue bsort;
  spfx::Workspace ws;
  spfx::CostProfile prof;
  std::size_t max_width = 1u << 20;
  bool fallback = false;
  void prepare(const Case& c) {
    prof = spfx::profile_costs(c.graph);
    fallback = static_cast<std::size_t>(prof.max_cost / prof.stride) + 1 > max_width;
    ws = spfx::Workspace{};
    if (variant.find("fresh") != std::string::npos) ws.mode = spfx::WsMode::Fresh;
    else if (variant.find("full") != std::string::npos) ws.mode = spfx::WsMode::ReuseFull;
    else if (variant.find("hybrid") != std::string::npos) ws.mode = spfx::WsMode::ReuseHybrid;
    else ws.mode = spfx::WsMode::ReuseSparse;
  }
  Result run(const Case& c, const Query& q) {
    if (variant == "ref" || q.reverse) return run_ref(c, q);   // reverse: production only
    spfx::Workspace fresh; fresh.mode = spfx::WsMode::Fresh;
    spfx::Workspace& w = (ws.mode == spfx::WsMode::Fresh) ? fresh : ws;
    auto res = q.use_residual ? std::span<const Cap>(c.residual) : std::span<const Cap>{};
    auto nm = q.use_masks ? c.nm() : std::span<const bool>{};
    auto em = q.use_masks ? c.em() : std::span<const bool>{};
    const bool use_heap = fallback || variant.rfind("heap", 0) == 0;
    if (use_heap) return spfx::spf_variant(c.graph, q.src, q.dst, q.multipath, q.sel, res, nm, em, heap, w, prof.max_cost, prof.stride);
    if (variant.rfind("bsort", 0) == 0) return spfx::spf_variant(c.graph, q.src, q.dst, q.multipath, q.sel, res, nm, em, bsort, w, prof.max_cost, prof.stride);
    if (variant.rfind("bucket2", 0) == 0) return spfx::spf_variant(c.graph, q.src, q.dst, q.multipath, q.sel, res, nm, em, bucket2, w, prof.max_cost, prof.stride);
    return spfx::spf_variant(c.graph, q.src, q.dst, q.multipath, q.sel, res, nm, em, bucket, w, prof.max_cost, prof.stride);
  }
};

static bool same(const Result& a, const Result& b) {
  return a.first == b.first && a.second.parent_offsets == b.second.parent_offsets &&
         a.second.parents == b.second.parents && a.second.via_edges == b.second.via_edges;
}

static NodeId far_node(const std::vector<Cost>& d) {
  NodeId best = -1; Cost bd = -1;
  for (std::size_t v = 0; v < d.size(); ++v) if (d[v] != INF && d[v] > bd) { bd = d[v]; best = NodeId(v); }
  return best;
}

static int correctness() {
  const char* variants[] = {"heap-fresh", "heap-full", "heap-sparse", "heap-hybrid", "bucket-fresh",
                            "bucket-full", "bucket-sparse", "bucket-hybrid", "bucket2-fresh", "bucket2-sparse",
                            "bucket2-hybrid", "bsort-fresh", "bsort-sparse"};
  std::vector<Case> cases;
  for (unsigned seed = 1; seed <= 12; ++seed) cases.push_back(make_case("random", 37 + seed * 3, 6, true, false, seed));
  for (unsigned seed = 1; seed <= 6; ++seed) cases.push_back(make_case("zerocyc", 40 + seed * 5, 5, seed % 2, false, seed));
  cases.push_back(make_case("random", 1024));
  cases.push_back(make_case("random", 512, 8, true, true, 7));
  cases.push_back(make_case("grid", 1024));
  cases.push_back(make_case("uniform", 500));
  cases.push_back(make_case("wide", 600, 8, true, false, 3, 100000));
  cases.push_back(make_case("pseudo", 64));
  cases.push_back(make_case("pseudo", 200, 8, true, false, 5));
  cases.push_back(make_case("chain", 300));
  cases.push_back(make_case("fabric", 100));
  // Dense parallel edges and self-loops: a tiny multigraph.
  {
    std::vector<NodeId> s, d; std::vector<Cost> c; std::vector<Cap> k; std::mt19937 r(99);
    for (int i = 0; i < 400; ++i) { s.push_back(r() % 12); d.push_back(r() % 12); c.push_back(r() % 3); k.push_back(double(r() % 5)); }
    Case mc{"multi", StrictMultiDiGraph::from_arrays(12, s, d, k, c)};
    mc.node_mask = std::make_unique<bool[]>(12); mc.edge_mask = std::make_unique<bool[]>(400);
    mc.residual.assign(mc.graph.capacity_view().begin(), mc.graph.capacity_view().end());
    for (int i = 0; i < 12; ++i) mc.node_mask[i] = r() % 9 != 0;
    for (int i = 0; i < 400; ++i) { mc.edge_mask[i] = r() % 7 != 0; if (r() % 5 == 0) mc.residual[i] = kMinCap * (r() % 3); }
    cases.push_back(std::move(mc));
  }
  std::size_t checked = 0, failed = 0;
  for (auto& c : cases) {
    const int n = c.graph.num_nodes();
    std::vector<NodeId> sources;
    if (c.fixed_src >= 0) sources = {c.fixed_src, 0, 3};
    else for (int i = 0; i < 6; ++i) sources.push_back((i * 997 + 1) % n);
    for (const char* vname : variants) {
      Runner r; r.variant = vname; r.prepare(c);
      for (NodeId s : sources) {
        for (int mp = 0; mp < 2; ++mp) for (int me = 0; me < 2; ++me) for (int tb = 0; tb < 2; ++tb)
        for (int rc = 0; rc < 2; ++rc) for (int ur = 0; ur < 2; ++ur) for (int um = 0; um < 2; ++um) {
          Query q; q.src = s; q.multipath = mp; q.sel.multi_edge = me;
          q.sel.tie_break = tb ? EdgeTieBreak::PreferHigherResidual : EdgeTieBreak::Deterministic;
          q.sel.require_capacity = rc; q.use_residual = ur; q.use_masks = um;
          std::vector<std::optional<NodeId>> dsts = {std::nullopt};
          auto full = run_ref(c, q);
          if (c.fixed_dst >= 0) dsts.push_back(c.fixed_dst);
          auto row = c.graph.row_offsets_view(); auto col = c.graph.col_indices_view();
          if (row[s + 1] > row[s]) dsts.push_back(col[row[s]]);
          NodeId far = far_node(full.first); if (far >= 0) dsts.push_back(far);
          dsts.push_back(s); dsts.push_back((s + n / 2) % n);
          for (auto& d : dsts) {
            q.dst = d;
            auto expected = run_ref(c, q);
            for (int rep = 0; rep < 2; ++rep) {   // second run exercises reuse/reset paths
              auto actual = r.run(c, q);
              ++checked;
              if (!same(expected, actual)) {
                ++failed;
                if (failed <= 20)
                  std::cerr << "MISMATCH case=" << c.name << " n=" << n << " variant=" << vname
                            << " src=" << s << " dst=" << (d ? std::to_string(*d) : "none")
                            << " mp=" << mp << " me=" << me << " tb=" << tb << " rc=" << rc
                            << " ur=" << ur << " um=" << um << " rep=" << rep << '\n';
              }
            }
          }
        }
      }
    }
  }
  std::cout << "checked=" << checked << " failed=" << failed << '\n';
  return failed ? 1 : 0;
}

struct Workload { std::string label; Case c; int queries; std::string dst_mode; };

static void time_workloads(std::vector<Workload>& wls, const std::string& variant, int samples, std::size_t maxw) {
  std::cout << "case,n,e,queries,dst,variant,fallback,W,median_ms,samples,wall_samples_ms\n";
  for (auto& w : wls) {
    auto& c = w.c; const int n = c.graph.num_nodes();
    std::vector<Query> qs;
    for (int i = 0; i < w.queries; ++i) {
      Query q; q.multipath = true; q.sel = EdgeSelection{}; q.use_residual = true; q.use_masks = true;
      q.src = c.fixed_src >= 0 ? c.fixed_src : (i * 997) % n;
      if (w.dst_mode == "reverse") q.reverse = true;
      else if (w.dst_mode == "fixed") q.dst = c.fixed_dst;
      else if (w.dst_mode == "near") {
        auto row = c.graph.row_offsets_view(); auto col = c.graph.col_indices_view();
        q.dst = col[row[q.src]];
      } else if (w.dst_mode == "far") {
        Query f = q; f.dst = std::nullopt; q.dst = far_node(run_ref(c, f).first);
      }
      qs.push_back(q);
    }
    std::vector<Result> expected; for (auto& q : qs) expected.push_back(run_ref(c, q));
    Runner r; r.variant = variant; r.max_width = maxw; r.prepare(c);
    std::vector<double> times;
    for (int rep = -2; rep < samples; ++rep) {
      std::vector<Result> actual(qs.size());
      auto start = Clock::now();
      for (std::size_t i = 0; i < qs.size(); ++i) actual[i] = r.run(c, qs[i]);
      double el = ms(start);
      for (std::size_t i = 0; i < qs.size(); ++i) if (!same(expected[i], actual[i])) throw std::runtime_error("timed result mismatch in " + w.label);
      if (rep >= 0) times.push_back(el);
    }
    std::vector<double> sorted = times; std::sort(sorted.begin(), sorted.end());
    std::cout << w.label << ',' << n << ',' << c.graph.num_edges() << ',' << w.queries << ',' << w.dst_mode << ','
              << variant << ',' << (r.fallback ? 1 : 0) << ',' << (r.prof.max_cost / r.prof.stride + 1) << ','
              << std::fixed << sorted[sorted.size() / 2] << ',' << samples << ',';
    for (std::size_t j = 0; j < times.size(); ++j) std::cout << (j ? ";" : "") << times[j];
    std::cout << '\n' << std::flush;
  }
}

int main(int argc, char** argv) {
  try {
    std::string mode = argc > 1 ? argv[1] : "correctness";
    if (mode == "correctness") return correctness();
    std::string variant = argc > 2 ? argv[2] : "ref";
    int samples = argc > 3 ? std::stoi(argv[3]) : 7;
    std::size_t maxw = argc > 4 ? std::stoull(argv[4]) : (1u << 20);
    std::vector<Workload> wls;
    if (mode == "time") {
      wls.push_back({"random", make_case("random", 1024), 1, "none"});
      wls.push_back({"random", make_case("random", 4096), 256, "none"});
      wls.push_back({"random", make_case("random", 16384), 64, "none"});
      wls.push_back({"random", make_case("random", 65536), 1, "none"});
      wls.push_back({"random", make_case("random", 4096, 32), 64, "none"});
      wls.push_back({"clos", make_case("clos", 800), 64, "none"});
      wls.push_back({"grid", make_case("grid", 10000), 1, "none"});
      wls.push_back({"chain", make_case("chain", 2048), 1, "none"});
      wls.push_back({"uniform", make_case("uniform", 16384), 64, "none"});
      wls.push_back({"random_masked", make_case("random", 16384, 8, true), 64, "none"});
      wls.push_back({"pseudo", make_case("pseudo", 800), 64, "fixed"});
      wls.push_back({"random", make_case("random", 1024), 64, "far"});
      wls.push_back({"random", make_case("random", 16384), 64, "far"});
      wls.push_back({"random", make_case("random", 65536), 16, "far"});
      wls.push_back({"random", make_case("random", 1024), 64, "near"});
      wls.push_back({"random", make_case("random", 16384), 64, "near"});
      wls.push_back({"random", make_case("random", 65536), 64, "near"});
      wls.push_back({"random", make_case("random", 262144), 16, "near"});
      wls.push_back({"random", make_case("random", 262144), 1, "none"});
      wls.push_back({"random", make_case("random", 4096), 64, "reverse"});
      wls.push_back({"random", make_case("random", 16384), 16, "reverse"});
      wls.push_back({"clos", make_case("clos", 800), 16, "reverse"});
      wls.push_back({"chain", make_case("chain", 2048), 1, "reverse"});
    } else if (mode == "short") {
      wls.push_back({"random", make_case("random", 16384), 64, "none"});
      wls.push_back({"random", make_case("random", 65536), 1, "none"});
      wls.push_back({"clos", make_case("clos", 800), 64, "none"});
      wls.push_back({"chain", make_case("chain", 2048), 1, "none"});
      wls.push_back({"random", make_case("random", 16384), 64, "near"});
      wls.push_back({"random", make_case("random", 65536), 64, "near"});
      wls.push_back({"random", make_case("random", 262144), 16, "near"});
      wls.push_back({"random", make_case("random", 262144), 1, "none"});
      wls.push_back({"random", make_case("random", 4096), 64, "reverse"});
      wls.push_back({"clos", make_case("clos", 800), 16, "reverse"});
    } else if (mode == "faults") {
      auto c = make_case("random", 262144);
      Runner r; r.variant = variant; r.prepare(c);
      std::vector<Query> qs;
      for (int i = 0; i < 16; ++i) {
        Query q; q.multipath = true; q.sel = EdgeSelection{}; q.use_residual = true; q.use_masks = true;
        q.src = (i * 997) % c.graph.num_nodes();
        auto row = c.graph.row_offsets_view(); auto col = c.graph.col_indices_view();
        q.dst = col[row[q.src]]; qs.push_back(q);
      }
      std::cout << "variant,sample,ms,minor_faults\n";
      for (int rep = 0; rep < samples; ++rep) {
        struct rusage a{}, b{}; getrusage(RUSAGE_SELF, &a);
        auto start = Clock::now();
        for (auto& q : qs) { auto res = r.run(c, q); (void)res; }
        double el = ms(start); getrusage(RUSAGE_SELF, &b);
        std::cout << variant << ',' << rep << ',' << el << ',' << (b.ru_minflt - a.ru_minflt) << '\n';
      }
      return 0;
    } else if (mode == "sweep") {
      for (Cost k : {Cost(1), Cost(31), Cost(1023), Cost(32767), Cost(1048575), Cost(33554431)})
        wls.push_back({"wide" + std::to_string(k), make_case("wide", 16384, 8, false, false, 11, k), 16, "none"});
    } else { std::cerr << "unknown mode\n"; return 2; }
    time_workloads(wls, variant, samples, maxw);
    return 0;
  } catch (const std::exception& e) { std::cerr << "FAIL: " << e.what() << '\n'; return 1; }
}
