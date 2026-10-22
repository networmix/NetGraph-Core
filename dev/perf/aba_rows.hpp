// The A/B/A rows of the Core timing harness (README.md, "Protocol"). Each row prepares its
// inputs (untimed), then times a fixed unit of semantic work per iteration through the adapter
// in aba_api.hpp. One process runs one block:
//
//   harness aba-list                         -> one JSON object per row
//   harness aba-run <row> <iterations>       -> one JSON object: per-iteration time stats,
//                                               digest, peak RSS, setup time
//   harness aba-calibrate <row> <seconds>    -> runs the row for about <seconds>, reports the
//                                               median iteration time (the warm-up block)
//
// harness_main, at the end of this file, reads the command line and gives the exit statuses.
// Every iteration returns a digest of its results; all iterations of a block must agree (the
// work is deterministic), and the orchestrator compares digests across sides.
//
// The destination rows read frozen query lists (queries/*.json, compiled into the block below)
// rather than a target picked through the library's edge layout, so every build searches for the
// same targets.
//
// The storage rows (README "Rows"): fourteen graph-construction shapes and six full searches for
// the graph's storage and construction, every input from a generator with fixed seeds, explicit
// integral capacities and costs, frozen with its full-content digest and shape (the second
// frozen block below); their timed digests read only API-independent values.
#pragma once

#include "aba_api.hpp"

#include <algorithm>
#include <array>
#include <cerrno>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <functional>
#include <numeric>
#include <random>
#include <set>
#include <stdexcept>
#include <string>
#include <sys/resource.h>
#include <type_traits>
#include <utility>
#include <vector>

namespace abarows {

using abaapi::Arrays;
using abaapi::Digest;
using abaapi::Elig;
using abaapi::Graph;
using abaapi::Kernel;
using abaapi::Sel;
using Iter = std::function<std::uint64_t()>;

// ---------------------------------------------------------------- frozen query lists
// BEGIN FROZEN QUERY LISTS (the values of queries/*.json; test_tool.py checks that they agree)
inline constexpr std::array<std::pair<int, int>, 64> kNear65k{{
    {1, 2}, {998, 999}, {1995, 1996}, {2992, 2993}, {3989, 3990}, {4986, 116}, {5983, 5984}, {6980, 5695},
    {7977, 4460}, {8974, 3025}, {9971, 1894}, {10968, 10969}, {11965, 4249}, {12962, 9444}, {13959, 2698}, {14956, 3611},
    {15953, 1555}, {16950, 6755}, {17947, 15211}, {18944, 2177}, {19941, 2870}, {20938, 4285}, {21935, 17837}, {22932, 10793},
    {23929, 23930}, {24926, 11028}, {25923, 3885}, {26920, 4095}, {27917, 5336}, {28914, 3084}, {29911, 5525}, {30908, 8194},
    {31905, 13485}, {32902, 2444}, {33899, 3487}, {34896, 3662}, {35893, 9132}, {36890, 7109}, {37887, 10077}, {38884, 7934},
    {39881, 4570}, {40878, 20916}, {41875, 2116}, {42872, 5858}, {43869, 14903}, {44866, 24}, {45863, 16319}, {46860, 5995},
    {47857, 26465}, {48854, 17016}, {49851, 9574}, {50848, 5966}, {51845, 14215}, {52842, 22580}, {53839, 5643}, {54836, 5009},
    {55833, 580}, {56830, 1068}, {57827, 1275}, {58824, 4854}, {59821, 11006}, {60818, 12565}, {61815, 7166}, {62812, 35168}}};
inline constexpr std::uint64_t kNear65kInput = 3859296219677894378ULL;
inline constexpr std::array<std::pair<int, int>, 8> kFar16k{{
    {1, 5378}, {998, 9128}, {1995, 11844}, {2992, 10152}, {3989, 10152}, {4986, 7529}, {5983, 16342}, {6980, 16342}}};
inline constexpr std::uint64_t kFar16kInput = 5806366537376710903ULL;
inline constexpr std::array<std::pair<int, int>, 8> kCanonical16k{{
    {1, 1339}, {998, 1339}, {1995, 2641}, {2992, 1339}, {3989, 10530}, {4986, 1339}, {5983, 1339}, {6980, 11931}}};
inline constexpr std::uint64_t kCanonical16kInput = 11698079179288285235ULL;
inline constexpr std::array<int, 50> kTargetSetDestinations{150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197, 198, 199};
inline constexpr std::uint64_t kTargetSetInput = 4484954733111131726ULL;
inline constexpr std::array<std::pair<const char*, const char*>, 4> kQueryLists{{
    {"spf.dag.near65k", "queries/near65k.json sha256:3f47de99a2b521fafe6d23b29805cdc19aca1306b6af64d13133745dbd8aecd4"},
    {"spf.dag.far16k", "queries/far16k.json sha256:bb84f77497bd2ffc9ba16bf9fdc87d41e848328632c8cb33d7057522062f14ad"},
    {"spf.single.canonical16k", "queries/canonical16k.json sha256:2e842fb1f95908a90603da81904c79d594d306084f125dfcd91942e8051236cd"},
    {"spf.to.targetset_clos", "queries/targetset_clos.json sha256:280f599ac837d892ef6e35a6e2161d657b09e1fea25d89ae38a6668895e76210"}}};
// END FROZEN QUERY LISTS

// The digest of an input (node count, then src, dst and cost of every edge in input order):
// a row refuses a frozen list whose input digest differs from its generator's arrays.
inline std::uint64_t input_digest(const Arrays& a) {
  Digest d;
  d.add_i(a.n);
  for (int e = 0; e < a.num_edges(); ++e) {
    d.add_i(a.src[std::size_t(e)]);
    d.add_i(a.dst[std::size_t(e)]);
    d.add_i(a.cost[std::size_t(e)]);
  }
  return d.h;
}

inline void require_frozen(const Arrays& a, std::uint64_t frozen_input, const char* list) {
  if (input_digest(a) != frozen_input)
    throw std::runtime_error(std::string("the frozen list ") + list + " was not frozen from this row's arrays");
}

// ---------------------------------------------------------------- generators (API-free)

// The harness "random" generator: ring plus random out-edges, costs 1..max_edge_cost.
inline Arrays gen_random(int n, int degree, unsigned seed, std::int64_t max_edge_cost = 31) {
  Arrays a;
  a.n = n;
  std::mt19937 rng(seed);
  auto add = [&](int u, int v, std::int64_t c) {
    a.src.push_back(u);
    a.dst.push_back(v);
    a.cost.push_back(c);
    a.cap.push_back(1.0 + double(rng() % 100) / 7.0);
  };
  for (int u = 0; u < n; ++u) {
    add(u, (u + 1) % n, 1 + std::int64_t(rng() % std::uint64_t(max_edge_cost)));
    for (int k = 1; k < degree; ++k) {
      const int v = int(rng() % unsigned(n));  // one draw per statement: the draw order is fixed
      add(u, v, 1 + std::int64_t(rng() % std::uint64_t(max_edge_cost)));
    }
  }
  return a;
}

// T2's mesh: a bidirectional ring over n nodes plus directed chords up to 4n distinct links in
// ascending (u, v), full capacities of 2^14, 2^15, 2^16 or 2^17 units (physical at kResolution),
// each with 5/8 to 8/8 of it free, costs 1..10; heads 0..h-1 toward node n - 1. `free` is in input order, in units.
struct T2Mesh {
  Arrays a;
  std::vector<std::int64_t> free;
  std::vector<int> heads;
};
inline T2Mesh gen_t2_mesh(int n, int h, unsigned seed) {
  std::mt19937 rng(seed);
  std::set<std::pair<int, int>> links;
  for (int i = 0; i < n; ++i) {
    links.insert({i, (i + 1) % n});
    links.insert({(i + 1) % n, i});
  }
  while (static_cast<int>(links.size()) < 4 * n) {
    const int u = int(rng() % unsigned(n));  // one draw per statement: the draw order is fixed
    const int v = int(rng() % unsigned(n));
    if (u != v) links.insert({u, v});
  }
  T2Mesh m;
  m.a.n = n;
  for (auto [u, v] : links) {
    const std::int64_t full = std::int64_t {1} << (14 + rng() % 4);
    const std::int64_t eighths = 5 + std::int64_t(rng() % 4);
    const std::int64_t cost = 1 + std::int64_t(rng() % 10);
    m.a.src.push_back(u);
    m.a.dst.push_back(v);
    m.a.cost.push_back(cost);
    m.a.cap.push_back(double(full) * abaapi::kResolution);
    m.free.push_back(full * eighths / 8);
  }
  for (int i = 0; i < h; ++i) m.heads.push_back(i);
  return m;
}

// Bidirectional link: input edges 2l (u -> v) and 2l + 1 (v -> u), same cost and capacity.
inline void add_link(Arrays& a, int u, int v, std::int64_t cost, double cap) {
  a.src.push_back(u); a.dst.push_back(v); a.cost.push_back(cost); a.cap.push_back(cap);
  a.src.push_back(v); a.dst.push_back(u); a.cost.push_back(cost); a.cap.push_back(cap);
}

inline double link_cap(std::mt19937& rng) { return 10.0 * double(1 + rng() % 8); }

// Two-tier Clos: leaves [0, L), spines [L, L + S), unit costs, capacities 10..80.
inline Arrays gen_clos(int leaves, int spines, unsigned seed) {
  Arrays a;
  a.n = leaves + spines;
  std::mt19937 rng(seed);
  for (int u = 0; u < leaves; ++u)
    for (int s = 0; s < spines; ++s) add_link(a, u, leaves + s, 1, link_cap(rng));
  return a;
}

// Square grid of bidirectional links; unit costs or costs 1..20.
inline Arrays gen_grid(int side, unsigned seed, bool unit) {
  Arrays a;
  a.n = side * side;
  std::mt19937 rng(seed);
  // One draw per statement (the cost before the capacity): the draw order is fixed.
  auto link = [&](int u, int v) {
    const std::int64_t cost = unit ? 1 : 1 + std::int64_t(rng() % 20);
    add_link(a, u, v, cost, link_cap(rng));
  };
  for (int u = 0; u < a.n; ++u) {
    if (u % side + 1 < side) link(u, u + 1);
    if (u + side < a.n) link(u, u + side);
  }
  return a;
}

// Three-tier Clos: supers [0, X); pod p: leaves then spines; leaf-spine mesh per pod, every
// spine to every super. Leaf id of (pod, i) is returned by fabric3_leaf.
struct Fabric3 { Arrays a; int leaves, spines, supers; };
inline int fabric3_leaf(const Fabric3& f, int pod, int i) { return f.supers + pod * (f.leaves + f.spines) + i; }
inline Fabric3 gen_fabric3(int pods, int leaves, int spines, int supers, unsigned seed) {
  Fabric3 f{{}, leaves, spines, supers};
  f.a.n = supers + pods * (leaves + spines);
  std::mt19937 rng(seed);
  for (int p = 0; p < pods; ++p) {
    const int base = supers + p * (leaves + spines);
    for (int i = 0; i < leaves; ++i)
      for (int j = 0; j < spines; ++j) add_link(f.a, base + i, base + leaves + j, 1, link_cap(rng));
    for (int j = 0; j < spines; ++j)
      for (int k = 0; k < supers; ++k) add_link(f.a, base + leaves + j, k, 1, 4 * link_cap(rng));
  }
  return f;
}

// F distinct failure states of k links each (bidirectional generators: both directions).
inline std::vector<std::vector<int>> failure_states(int num_links, int states, int k, unsigned seed) {
  std::mt19937 rng(seed);
  std::set<std::vector<int>> seen;
  std::vector<std::vector<int>> out;
  while (int(out.size()) < states) {
    std::set<int> links;
    while (int(links.size()) < k) links.insert(int(rng() % unsigned(num_links)));
    std::vector<int> edges;
    for (int l : links) { edges.push_back(2 * l); edges.push_back(2 * l + 1); }
    if (seen.insert(edges).second) out.push_back(edges);
  }
  return out;
}

// ---------------------------------------------------------------- rows

struct Row {
  std::string name, group, desc;
  std::function<Iter()> setup;
};

// The static exclusions and the residual of the masked and residual rows (input order): one
// seeded draw per node, then two per edge (an excluded edge, a zero residual).
struct MasksResidual {
  std::vector<int> xn, xe;
  std::vector<double> res_in;  // input-edge order: the capacity, or 0 for a gated edge
};
inline MasksResidual masks_residual(const Arrays& a) {
  MasksResidual m;
  std::mt19937 rng(17);
  m.res_in = a.cap;
  for (int v = 0; v < a.n; ++v) if (rng() % 23 == 0) m.xn.push_back(v);
  for (int e = 0; e < a.num_edges(); ++e) {
    if (rng() % 11 == 0) m.xe.push_back(e);
    if (rng() % 13 == 0) m.res_in[std::size_t(e)] = 0.0;
  }
  return m;
}

// The full forward searches from `queries` sources, plain or with the static exclusions and the
// residual of masks_residual.
inline Iter spf_full_row(Arrays a, int queries, bool masked_residual) {
  auto g = std::make_shared<Graph>(abaapi::build(a));
  auto el = std::make_shared<Elig>();
  auto res = std::make_shared<std::vector<double>>();
  if (masked_residual) {
    const auto m = masks_residual(a);
    *el = abaapi::eligibility(*g, m.xn, m.xe);
    *res = abaapi::edge_values(*g, m.res_in);
  }
  std::vector<int> srcs;
  for (int i = 0; i < queries; ++i) srcs.push_back(int((std::int64_t(i) * 997 + 1) % a.n));
  return [g, el, res, srcs, masked_residual] {
    Digest d;
    for (int s : srcs)
      d.add(abaapi::spf(*g, s, masked_residual ? el.get() : nullptr, masked_residual ? res.get() : nullptr));
    return d.h;
  };
}

// ---------------------------------------------------------------- storage rows (graph storage and construction)
// README "Rows". Every input comes from a generator below with fixed seeds, one random draw per
// statement (the draw order is part of the input), explicit costs and positive integral
// capacities (exact as doubles and in the adapter's units of 1/4096). Each input is frozen
// in the frozen block by its full-content digest and its shape; a row refuses an input whose
// generator no longer reproduces both.

// The full-content digest of an input: n and E, then src, dst, cost and the capacity's bits of
// every edge in input order (input_digest above omits the capacities).
inline std::uint64_t input_digest_full(const Arrays& a) {
  Digest d;
  d.add_i(a.n);
  d.add_i(a.num_edges());
  for (std::size_t e = 0; e < a.src.size(); ++e) {
    d.add_i(a.src[e]);
    d.add_i(a.dst[e]);
    d.add_i(a.cost[e]);
    d.add_d(a.cap[e]);
  }
  return d.h;
}

// The shape of an input, as frozen: N, E, the number K of distinct costs, the largest out- and
// in-degree, the cost GCD and maximum (0 without edges), the nodes without any edge, the number
// of input positions where src decreases, and whether the input is in ascending (source, cost, destination)
// order (key_sorted; equal keys keep their input order, so this is the B' id order) or in
// descending order (key_reversed). digest: input_digest_full.
struct InputShape {
  std::int64_t n = 0, e = 0, k = 0, max_out = 0, max_in = 0, gcd = 0, max_edge_cost = 0, isolated = 0, src_descents = 0;
  bool key_sorted = false, key_reversed = false;
  std::uint64_t digest = 0;
  bool operator==(const InputShape&) const = default;
};

inline InputShape shape_of(const Arrays& a) {
  InputShape s;
  s.n = a.n;
  s.e = a.num_edges();
  std::vector<std::int64_t> costs(a.cost.begin(), a.cost.end());
  std::sort(costs.begin(), costs.end());
  s.k = std::unique(costs.begin(), costs.end()) - costs.begin();
  std::vector<std::int64_t> out(std::size_t(a.n), 0), in(std::size_t(a.n), 0);
  for (std::size_t e = 0; e < a.src.size(); ++e) {
    ++out[std::size_t(a.src[e])];
    ++in[std::size_t(a.dst[e])];
    s.gcd = std::gcd(s.gcd, a.cost[e]);
    s.max_edge_cost = std::max(s.max_edge_cost, a.cost[e]);
  }
  for (int v = 0; v < a.n; ++v) {
    s.max_out = std::max(s.max_out, out[std::size_t(v)]);
    s.max_in = std::max(s.max_in, in[std::size_t(v)]);
    if (out[std::size_t(v)] == 0 && in[std::size_t(v)] == 0) ++s.isolated;
  }
  s.key_sorted = s.key_reversed = true;
  for (std::size_t e = 0; e + 1 < a.src.size(); ++e) {
    const auto x = std::make_tuple(a.src[e], a.cost[e], a.dst[e]);
    const auto y = std::make_tuple(a.src[e + 1], a.cost[e + 1], a.dst[e + 1]);
    if (y < x) s.key_sorted = false;
    if (x < y) s.key_reversed = false;
    if (a.src[e + 1] < a.src[e]) ++s.src_descents;
  }
  s.digest = input_digest_full(a);
  return s;
}

inline std::string json_shape(const InputShape& s) {
  return "{\"n\":" + std::to_string(s.n) + ",\"e\":" + std::to_string(s.e) + ",\"k\":" + std::to_string(s.k) +
         ",\"max_out\":" + std::to_string(s.max_out) + ",\"max_in\":" + std::to_string(s.max_in) + ",\"gcd\":" +
         std::to_string(s.gcd) + ",\"max_cost\":" + std::to_string(s.max_edge_cost) + ",\"isolated\":" +
         std::to_string(s.isolated) + ",\"src_descents\":" + std::to_string(s.src_descents) + ",\"key_sorted\":" +
         (s.key_sorted ? "true" : "false") + ",\"key_reversed\":" + (s.key_reversed ? "true" : "false") +
         ",\"digest\":\"" + std::to_string(s.digest) + "\"}";
}

inline void push_edge(Arrays& a, int u, int v, std::int64_t cost, double cap) {
  a.src.push_back(u);
  a.dst.push_back(v);
  a.cost.push_back(cost);
  a.cap.push_back(cap);
}

// Fisher-Yates with mt19937_64 (fully specified by the standard, unlike std::shuffle).
inline std::vector<std::size_t> permutation(std::size_t m, std::uint64_t seed) {
  std::vector<std::size_t> p(m);
  std::iota(p.begin(), p.end(), std::size_t {0});
  std::mt19937_64 rng(seed);
  for (std::size_t i = m; i > 1; --i) {
    const std::size_t j = std::size_t(rng() % std::uint64_t(i));
    std::swap(p[i - 1], p[j]);
  }
  return p;
}

// The input permuted: position i takes the edge at position p[i] (every property together).
inline Arrays permuted(const Arrays& a, const std::vector<std::size_t>& p) {
  Arrays s;
  s.n = a.n;
  for (std::size_t i = 0; i < p.size(); ++i) push_edge(s, a.src[p[i]], a.dst[p[i]], a.cost[p[i]], a.cap[p[i]]);
  return s;
}

// A 256 x 256 weighted grid, both directions of every link with independently drawn costs 1..31.
// Each edge draws its cost first, then its integral capacity in 1..100.
inline Arrays gen_wgrid256() {
  Arrays a;
  const int side = 256;
  a.n = side * side;
  std::mt19937 rng(1729);
  auto edge = [&](int u, int v) {
    const std::int64_t c = 1 + std::int64_t(rng() % 31);
    const double cap = double(1 + rng() % 100);
    push_edge(a, u, v, c, cap);
  };
  for (int u = 0; u < a.n; ++u) {
    if (u % side + 1 < side) { edge(u, u + 1); edge(u + 1, u); }
    if (u + side < a.n) { edge(u, u + side); edge(u + side, u); }
  }
  return a;
}

// A Clos with source and sink hubs at positive costs (the graph rejects zero costs): Clos 400 leaves x
// 400 spines in both directions at cost 1, a source hub (node 800) with cost-2 edges to leaves
// 0..99 and a sink hub (node 801) fed by leaves 300..399 at cost 2.
inline Arrays gen_fanout_clos() {
  Arrays a;
  a.n = 802;
  std::mt19937 rng(1729);
  auto edge = [&](int u, int v, std::int64_t c) {
    const double cap = double(1 + rng() % 100);
    push_edge(a, u, v, c, cap);
  };
  for (int u = 0; u < 400; ++u)
    for (int v = 400; v < 800; ++v) { edge(u, v, 1); edge(v, u, 1); }
  for (int u = 0; u < 100; ++u) edge(800, u, 2);
  for (int u = 300; u < 400; ++u) edge(u, 801, 2);
  return a;
}

// A ring plus random out-edges (gen_random's shape) with explicit draws: per node the
// ring edge u -> u + 1 (cost draw, capacity draw), then degree - 1 edges (target draw, cost draw,
// capacity draw). cost_of(rng) draws one cost.
template <class CostOf>
inline Arrays gen_random_drawn(int n, int degree, unsigned seed, CostOf cost_of) {
  Arrays a;
  a.n = n;
  std::mt19937 rng(seed);
  for (int u = 0; u < n; ++u) {
    {
      const std::int64_t c = cost_of(rng);
      const double cap = double(1 + rng() % 100);
      push_edge(a, u, (u + 1) % n, c, cap);
    }
    for (int k = 1; k < degree; ++k) {
      const int v = int(rng() % unsigned(n));
      const std::int64_t c = cost_of(rng);
      const double cap = double(1 + rng() % 100);
      push_edge(a, u, v, c, cap);
    }
  }
  return a;
}

inline std::int64_t cost_1_to_31(std::mt19937& rng) { return 1 + std::int64_t(rng() % 31); }

// Random graph: 131,072 nodes x 8 = 1,048,576 edges, costs 1..31, generator order
// (random1m_nat) or the same edges in a uniform random order (random1m_scat; Fisher-Yates over
// mt19937_64(4242), every property permuted together).
inline Arrays gen_random1m(bool scattered) {
  auto a = gen_random_drawn(131072, 8, 1729, cost_1_to_31);
  return scattered ? permuted(a, permutation(a.src.size(), 4242)) : a;
}

// 65,536 nodes x 8, costs drawn from 1..2^30 (far more than 65,536 distinct values).
inline Arrays gen_wide65k() {
  return gen_random_drawn(65536, 8, 1031, [](std::mt19937& rng) { return 1 + std::int64_t(rng() % (1U << 30)); });
}

// 65,536 nodes x 8 with exactly k distinct costs 1..k: the topology and capacities from
// mt19937(1032) (one topology for both k); the first k edges take a random permutation of 1..k
// (every value present; permutation() with seed 2032), every later edge a uniform draw from 1..k
// (mt19937_64(2033)).
inline Arrays gen_kcosts(std::int64_t k) {
  auto a = gen_random_drawn(65536, 8, 1032, [](std::mt19937&) { return std::int64_t {1}; });
  const auto perm = permutation(std::size_t(k), 2032);
  std::mt19937_64 rng(2033);
  for (std::size_t e = 0; e < a.cost.size(); ++e)
    a.cost[e] = e < std::size_t(k) ? 1 + std::int64_t(perm[e]) : 1 + std::int64_t(rng() % std::uint64_t(k));
  return a;
}

// One node (0) with 1,048,576 out-edges to targets drawn from 131,072 nodes (parallel edges and
// self-loops of node 0 included): per edge the target draw, the cost draw, the capacity draw.
inline Arrays gen_star1m() {
  Arrays a;
  a.n = 131072;
  std::mt19937 rng(1034);
  for (int e = 0; e < (1 << 20); ++e) {
    const int v = int(rng() % unsigned(a.n));
    const std::int64_t c = 1 + std::int64_t(rng() % (1U << 20));
    const double cap = double(1 + rng() % 100);
    push_edge(a, 0, v, c, cap);
  }
  return a;
}

// The same single long out-row with exactly 31 costs (1..31) and distinct capacities 1 + the
// input position (equal (source, cost, destination) keys are told apart by their capacity); reversed: every
// edge turned around, one long in-row of node 0.
inline Arrays gen_star1m_smallk(bool reversed) {
  Arrays a;
  a.n = 131072;
  std::mt19937 rng(1036);
  for (int e = 0; e < (1 << 20); ++e) {
    const int v = int(rng() % unsigned(a.n));
    const std::int64_t c = 1 + std::int64_t(rng() % 31);
    if (reversed) push_edge(a, v, 0, c, double(1 + e));
    else push_edge(a, 0, v, c, double(1 + e));
  }
  return a;
}

// 1,048,576 nodes, 4,096 edges with both endpoints drawn over the whole node range, costs 1..31.
inline Arrays gen_sparse1m() {
  Arrays a;
  a.n = 1 << 20;
  std::mt19937 rng(1035);
  for (int e = 0; e < 4096; ++e) {
    const int u = int(rng() % unsigned(a.n));
    const int v = int(rng() % unsigned(a.n));
    const std::int64_t c = 1 + std::int64_t(rng() % 31);
    const double cap = double(1 + rng() % 100);
    push_edge(a, u, v, c, cap);
  }
  return a;
}

inline Arrays gen_empty1m() {
  Arrays a;
  a.n = 1 << 20;
  return a;
}

// construct.random65k's topology and costs (gen_random(65536, 8, 11)) with integral capacities
// recovered from the same draws (1 + k for the generator's 1 + k / 7), stably sorted by (src, cost,
// dst): the B' id order; reversed: the whole sequence of that input in reverse.
inline Arrays gen_random65k_bprime(bool reversed) {
  auto g = gen_random(65536, 8, 11);
  for (auto& c : g.cap) c = double(1 + std::llround((c - 1.0) * 7.0));
  const auto& [gn, tails, heads, caps, weights] = g;  // the Arrays fields, in declaration order
  std::vector<std::size_t> p(tails.size());
  std::iota(p.begin(), p.end(), std::size_t {0});
  std::stable_sort(p.begin(), p.end(), [&](std::size_t x, std::size_t y) {
    return std::make_tuple(tails[x], weights[x], heads[x]) < std::make_tuple(tails[y], weights[y], heads[y]);
  });
  if (reversed) std::reverse(p.begin(), p.end());
  return permuted(g, p);
}

// The reverse-search heap row: 16,384 nodes x 8, costs 1..2^20 (W far above 65,536: the heap).
inline Arrays gen_wide_fallback16k() {
  return gen_random_drawn(16384, 8, 1037, [](std::mt19937& rng) { return 1 + std::int64_t(rng() % (1U << 20)); });
}

// The bucket-width boundary: one topology (mt19937(1038), 16,384 nodes x 8), each edge's cost
// base + one drawn bit: base 65,534 gives costs {65534, 65535} (gcd 1, W = 65,536: the buckets),
// base 65,535 gives {65535, 65536} (W = 65,537: the heap).
inline Arrays gen_width(std::int64_t base) {
  return gen_random_drawn(16384, 8, 1038, [base](std::mt19937& rng) { return base + std::int64_t(rng() % 2); });
}

struct StorageInput {
  const char* name;
  Arrays (*make)();
};

inline const std::vector<StorageInput>& storage_inputs() {
  static const std::vector<StorageInput> inputs = {
      {"wgrid256", [] { return gen_wgrid256(); }},
      {"fanout_clos", [] { return gen_fanout_clos(); }},
      {"random1m_nat", [] { return gen_random1m(false); }},
      {"random1m_scat", [] { return gen_random1m(true); }},
      {"wide65k", [] { return gen_wide65k(); }},
      {"k65536", [] { return gen_kcosts(65536); }},
      {"k65537", [] { return gen_kcosts(65537); }},
      {"star1m", [] { return gen_star1m(); }},
      {"sparse1m", [] { return gen_sparse1m(); }},
      {"empty1m", [] { return gen_empty1m(); }},
      {"random65k_bprime", [] { return gen_random65k_bprime(false); }},
      {"random65k_reversed", [] { return gen_random65k_bprime(true); }},
      {"star1m_smallk", [] { return gen_star1m_smallk(false); }},
      {"instar1m_smallk", [] { return gen_star1m_smallk(true); }},
      {"wide_fallback16k", [] { return gen_wide_fallback16k(); }},
      {"width65536", [] { return gen_width(65534); }},
      {"width65537", [] { return gen_width(65535); }},
  };
  return inputs;
}

inline Arrays make_storage_input(const std::string& name) {
  for (const auto& i : storage_inputs())
    if (name == i.name) return i.make();
  throw std::runtime_error("no storage input " + name);
}

struct FrozenShape {
  const char* input;
  InputShape shape;
};

// The frozen inputs (shape and digest per input), the search roots and the expected queues.
// BEGIN FROZEN STORAGE INPUTS (queries/storage/inputs.json; test_tool.py checks that they agree)
inline const std::vector<FrozenShape> kStorageInputs {
    {"wgrid256", {65536, 261120, 31, 4, 4, 1, 31, 0, 130305, false, false, 14812220864904460480ULL}},
    {"fanout_clos", {802, 320200, 2, 401, 401, 1, 2, 0, 160000, false, false, 5230821760560145505ULL}},
    {"random1m_nat", {131072, 1048576, 31, 8, 24, 1, 31, 0, 0, false, false, 16852305707948996038ULL}},
    {"random1m_scat", {131072, 1048576, 31, 8, 24, 1, 31, 0, 524810, false, false, 10807403047943366003ULL}},
    {"wide65k", {65536, 524288, 524146, 8, 21, 1, 1073738528, 0, 0, false, false, 13538344444856380895ULL}},
    {"k65536", {65536, 524288, 65536, 8, 21, 1, 65536, 0, 0, false, false, 10929368686095504329ULL}},
    {"k65537", {65536, 524288, 65537, 8, 21, 1, 65537, 0, 0, false, false, 758208004199587604ULL}},
    {"star1m", {131072, 1048576, 662682, 1048576, 27, 1, 1048575, 41, 0, false, false, 15197849622511790248ULL}},
    {"sparse1m", {1048576, 4096, 31, 2, 2, 1, 31, 1040409, 2050, false, false, 12990156200638721186ULL}},
    {"empty1m", {1048576, 0, 0, 0, 0, 0, 0, 1048576, 0, true, true, 18109181668318218688ULL}},
    {"random65k_bprime", {65536, 524288, 31, 8, 21, 1, 31, 0, 0, true, false, 12833760064757678648ULL}},
    {"random65k_reversed", {65536, 524288, 31, 8, 21, 1, 31, 0, 65535, false, true, 11592727304513109687ULL}},
    {"star1m_smallk", {131072, 1048576, 31, 1048576, 23, 1, 31, 46, 0, false, false, 8802474158293546069ULL}},
    {"instar1m_smallk", {131072, 1048576, 31, 23, 1048576, 1, 31, 46, 524217, false, false, 10839275545656040073ULL}},
    {"wide_fallback16k", {16384, 131072, 123230, 8, 23, 1, 1048557, 0, 0, false, false, 1535512841938229597ULL}},
    {"width65536", {16384, 131072, 2, 8, 20, 1, 65535, 0, 0, false, false, 6236947617789042840ULL}},
    {"width65537", {16384, 131072, 2, 8, 20, 1, 65536, 0, 0, false, false, 9851477194174659493ULL}},
};
inline const std::vector<std::pair<const char*, std::vector<int>>> kStorageRoots {
    {"spf.to.wgrid256", {1, 998, 1995, 2992, 3989, 4986, 5983, 6980}},
    {"spf.to.random1m", {1, 998, 1995, 2992, 3989, 4986, 5983, 6980}},
    {"spf.to.wide_fallback16k", {1, 998, 1995, 2992, 3989, 4986, 5983, 6980}},
    {"spf.to.width65536", {1, 998, 1995, 2992, 3989, 4986, 5983, 6980}},
    {"spf.to.width65537", {1, 998, 1995, 2992, 3989, 4986, 5983, 6980}},
    {"spf.full.random1m", {1, 998, 1995, 2992, 3989, 4986, 5983, 6980}},
};
// END FROZEN STORAGE INPUTS

// The roots of a storage search row: 8 nodes (i * 997 + 1) mod N, frozen in the block above.
inline std::vector<int> rule_roots(std::int64_t n) {
  std::vector<int> r;
  for (std::int64_t i = 0; i < 8; ++i) r.push_back(int((i * 997 + 1) % n));
  return r;
}

// The input of a storage row, checked against the frozen shape and digest.
inline Arrays checked_storage_input(const std::string& name, Arrays a) {
  const auto got = shape_of(a);
  for (const auto& f : kStorageInputs) {
    if (name != f.input) continue;
    if (!(got == f.shape))
      throw std::runtime_error("storage input " + name + " differs from its frozen shape: " + json_shape(got) + " frozen " +
                               json_shape(f.shape));
    return a;
  }
  throw std::runtime_error("storage input " + name + " is not frozen");
}

inline Arrays frozen_storage_input(const std::string& name) { return checked_storage_input(name, make_storage_input(name)); }

inline std::vector<int> frozen_roots(const std::string& row, std::int64_t n) {
  for (const auto& [name, roots] : kStorageRoots)
    if (row == name) {
      if (roots != rule_roots(n)) throw std::runtime_error("the frozen roots of " + row + " differ from the rule");
      return roots;
    }
  throw std::runtime_error("no frozen roots for " + row);
}

// Checks the queue choice of a search row: the adapter's (the library's rule: the width W and
// buckets or heap) and the width read from the input arrays (max cost / cost GCD + 1) must both
// equal the frozen expectation.
inline void require_queue(const Graph& g, const Arrays& a, std::int64_t want_width, bool want_buckets) {
  const auto got = abaapi::queue_choice(g);
  std::int64_t gcd = 0, mx = 0;
  for (auto c : a.cost) { gcd = std::gcd(gcd, c); mx = std::max(mx, c); }
  const std::int64_t w_in = mx / std::max<std::int64_t>(1, gcd) + 1;
  if (got.width != want_width || got.buckets != want_buckets || w_in != want_width)
    throw std::runtime_error("queue choice: the adapter reports W = " + std::to_string(got.width) +
                             (got.buckets ? " (buckets)" : " (heap)") + ", the input gives W = " + std::to_string(w_in) +
                             ", the row expects W = " + std::to_string(want_width) + (want_buckets ? " (buckets)" : " (heap)"));
}

struct StorageSearch {
  const char* row;
  const char* input;
  bool reverse;
  std::int64_t width;
  bool buckets;
  const char* desc;
};

inline void add_storage_rows(std::vector<Row>& r) {
  struct C { const char* row; const char* input; const char* desc; };
  const C constructs[] = {
      {"construct.wgrid256", "wgrid256", "construction, 256x256 grid both directions, independent costs 1..31 per direction (many short rows, K = 31)"},
      {"construct.fanout_clos", "fanout_clos", "construction, Clos 400x400 at cost 1 plus a cost-2 source hub to 100 leaves and a sink hub from 100 leaves (K = 2; rows of about 400, two 100-edge hub rows)"},
      {"construct.random1m_nat", "random1m_nat", "construction, random 131,072 x 8 = 1,048,576 edges, costs 1..31, generator (source-major) order"},
      {"construct.random1m_scat", "random1m_scat", "construction, the same 1,048,576 edges in a uniform random input order"},
      {"construct.wide65k", "wide65k", "construction, random 65,536 x 8, costs 1..2^30 (K far above 65,536: the sort fallback)"},
      {"construct.k65536", "k65536", "construction, random 65,536 x 8 with exactly 65,536 distinct costs (the hash limit)"},
      {"construct.k65537", "k65537", "construction, the same topology with exactly 65,537 distinct costs (one above the hash limit)"},
      {"construct.star1m", "star1m", "construction, one node with 1,048,576 out-edges to 131,072 nodes, costs 1..2^20 (one very long row, large K)"},
      {"construct.sparse1m", "sparse1m", "construction, 1,048,576 nodes and 4,096 edges spread over the node range, K = 31 (mostly isolated nodes)"},
      {"construct.empty1m", "empty1m", "construction, 1,048,576 nodes and no edge (empty offsets, maximum and GCD 0)"},
      {"construct.random65k_bprime", "random65k_bprime", "construction, construct.random65k's topology and costs sorted by (src, cost, dst, input): already in B' order"},
      {"construct.random65k_reversed", "random65k_reversed", "construction, the preceding input in reverse order, every property array included"},
      {"construct.star1m_smallk", "star1m_smallk", "construction, one node with 1,048,576 out-edges to 131,072 nodes, exactly 31 costs, distinct capacities"},
      {"construct.instar1m_smallk", "instar1m_smallk", "construction, the preceding edges turned around: one node with 1,048,576 in-edges, exactly 31 costs"},
  };
  for (const auto& c : constructs) {
    const std::string row = c.row, input = c.input;
    Row x {row, "construct", c.desc,
           [input]() -> Iter {
             auto a = std::make_shared<Arrays>(frozen_storage_input(input));
             auto in = std::make_shared<abaapi::ConstructInput>(abaapi::construct_prepare(*a));
             return Iter([a, in] { return abaapi::construct_run(*in); });
           }};
    r.push_back(std::move(x));
  }
  const StorageSearch searches[] = {
      {"spf.to.wgrid256", "wgrid256", true, 32, true, "8 full reverse searches on construct.wgrid256's graph, frozen roots (bucket queue, W = 32)"},
      {"spf.to.random1m", "random1m_scat", true, 32, true, "8 full reverse searches on construct.random1m_scat's graph, frozen roots (large scattered graph, bucket queue)"},
      {"spf.to.wide_fallback16k", "wide_fallback16k", true, 1048558, false, "8 full reverse searches, random 16k x8, costs 1..2^20, gcd 1: the heap (W above 65,536)"},
      {"spf.to.width65536", "width65536", true, 65536, true, "8 full reverse searches, random 16k x8, costs 65,534 and 65,535: the buckets at W = 65,536"},
      {"spf.to.width65537", "width65537", true, 65537, false, "8 full reverse searches, the same topology and roots, costs 65,535 and 65,536: the heap at W = 65,537"},
      {"spf.full.random1m", "random1m_scat", false, 32, true, "8 full forward searches on construct.random1m_scat's graph, frozen roots (forward only: the reverse arrays unread)"},
  };
  for (const auto& sp : searches) {
    const std::string row = sp.row, input = sp.input;
    const bool reverse = sp.reverse, buckets = sp.buckets;
    const std::int64_t width = sp.width;
    Row x {row, "spf", sp.desc,
           [row, input, reverse, width, buckets]() -> Iter {
             const auto a = frozen_storage_input(input);
             const auto roots = frozen_roots(row, a.n);
             auto g = std::make_shared<Graph>(abaapi::build(a));
             require_queue(*g, a, width, buckets);
             return Iter([g, roots, reverse] {
               Digest d;
               for (int v : roots) d.add(abaapi::search_run(*g, v, reverse));
               return d.h;
             });
           }};
    r.push_back(std::move(x));
  }
}

inline std::vector<Row> rows() {
  std::vector<Row> r;
  // ---- SPF classes
  r.push_back({"spf.full.random16k", "spf",
               "full forward SPF, random 16k nodes x8 out-edges, costs 1..31, 8 sources",
               [] { return spf_full_row(gen_random(16384, 8, 1), 8, false); }});
  r.push_back({"spf.full.clos", "spf",
               "full forward SPF, two-tier Clos 400x400 (unit costs, 320k edges), 8 sources",
               [] { return spf_full_row(gen_clos(400, 400, 2), 8, false); }});
  r.push_back({"spf.full.grid10k", "spf", "full forward SPF, 100x100 unit grid, 4 sources",
               [] { return spf_full_row(gen_grid(100, 3, true), 4, false); }});
  r.push_back({"spf.full.elig_residual16k", "spf",
               "full forward SPF with static exclusions and a capacity-gated residual, random 16k, 8 sources",
               [] { return spf_full_row(gen_random(16384, 8, 4), 8, true); }});
  r.push_back({"spf.full.wide_fallback16k", "spf",
               "full forward SPF, costs 1..2^20 (bucket width above the limit: heap route), random 16k, 8 sources",
               [] { return spf_full_row(gen_random(16384, 8, 6, std::int64_t(1) << 20), 8, false); }});
  // The destination rows read frozen (source, target) lists (queries/near65k.json, far16k.json),
  // so every build searches for the same targets. The timed call returns the DAG of the shortest
  // walks that end in the target.
  auto dag_row = [](int n, unsigned seed, const auto& pairs, std::uint64_t frozen, const char* list) {
    auto a = gen_random(n, 8, seed);
    require_frozen(a, frozen, list);
    auto g = std::make_shared<Graph>(abaapi::build(a));
    std::vector<std::pair<int, int>> q(pairs.begin(), pairs.end());
    return Iter([g, q] { Digest d; for (auto [s, t] : q) d.add(abaapi::spf_dag(*g, s, t)); return d.h; });
  };
  r.push_back({"spf.dag.near65k", "spf",
               "destination SPF to a near target (frozen list: the lowest-numbered other out-neighbour), random 65k, "
               "64 queries",
               [=] { return dag_row(65536, 7, kNear65k, kNear65kInput, "near65k"); }});
  r.push_back({"spf.dag.far16k", "spf",
               "destination SPF to the farthest node (frozen list), random 16k, 8 queries",
               [=] { return dag_row(16384, 8, kFar16k, kFar16kInput, "far16k"); }});
  r.push_back({"spf.to.random16k", "spf", "reverse SPF (spf_to), random 16k, 8 destinations",
               [] {
                 auto g = std::make_shared<Graph>(abaapi::build(gen_random(16384, 8, 9)));
                 return Iter([g] { Digest d; for (int i = 0; i < 8; ++i) d.add(abaapi::spf_to(*g, (i * 997 + 1) % 16384)); return d.h; });
               }});
  // The endpoint-set reverse search: the input is the 300-node Clos gen_clos(200, 100, 10) and the
  // destination set, leaves 150..199 (queries/targetset_clos.json); the timed call is a full
  // reverse spf rooted at the destination set with targets absent.
  auto targetset_dests = [] { return std::vector<int>(kTargetSetDestinations.begin(), kTargetSetDestinations.end()); };
  r.push_back({"spf.to.targetset_clos", "spf",
               "endpoint-set reverse SPF rooted at leaves 150..199 of a Clos 200x100, 16 calls (spf(roots, Reverse), "
               "no targets)",
               [=]() -> Iter {
                 auto a = gen_clos(200, 100, 10);
                 require_frozen(a, kTargetSetInput, "targetset_clos");
                 auto ts = std::make_shared<abaapi::TargetSet>(abaapi::targetset_prepare(a, targetset_dests()));
                 return Iter([ts] { Digest d; for (int i = 0; i < 16; ++i) d.add(abaapi::targetset_run(*ts)); return d.h; });
               }});
  // ---- the canonical single path
  // The canonical single path (multipath = false: the complete key under HIGHER_AVAILABLE_CAPACITY
  // on the residual, one member per pair, residual gating) from 8 fixed sources to frozen targets
  // (queries/canonical16k.json: the farthest node over the edges the gating keeps).
  r.push_back({"spf.single.canonical16k", "spf",
               "canonical single path to a fixed target, HIGHER_AVAILABLE_CAPACITY on the residual, one member, residual "
               "gating, random 16k, 8 frozen pairs",
               []() -> Iter {
                 auto a = gen_random(16384, 8, 5);
                 require_frozen(a, kCanonical16kInput, "canonical16k");
                 auto g = std::make_shared<Graph>(abaapi::build(a));
                 auto res = std::make_shared<std::vector<double>>(abaapi::edge_values(*g, masks_residual(a).res_in));
                 Sel s;
                 s.multipath = false;
                 s.parallel_edges = false;
                 s.higher_available_capacity = true;
                 std::vector<std::pair<int, int>> q(kCanonical16k.begin(), kCanonical16k.end());
                 return Iter([g, res, q, s] { Digest d; for (auto [u, v] : q) d.add(abaapi::spf_path(*g, u, v, s, res.get())); return d.h; });
               }});
  // ---- construction
  r.push_back({"construct.random65k", "construct", "graph construction from arrays, random 65k x8 (524k edges)",
               [] { auto a = std::make_shared<Arrays>(gen_random(65536, 8, 11)); return Iter([a] { return abaapi::construct(*a); }); }});
  r.push_back({"construct.clos", "construct", "graph construction from arrays, Clos 400x400 (320k edges)",
               [] { auto a = std::make_shared<Arrays>(gen_clos(400, 400, 12)); return Iter([a] { return abaapi::construct(*a); }); }});
  // ---- placement kernels on fixed DAGs
  auto kernel_row = [](Arrays a, std::vector<std::pair<int, int>> pairs, double amount, Kernel k) {
    auto g = std::make_shared<Graph>(abaapi::build(a));
    auto dags = std::make_shared<std::vector<abaapi::Dag>>();
    for (auto [s, t] : pairs) dags->push_back(abaapi::base_dag(*g, s, t));
    return Iter([g, dags, pairs, amount, k] {
      abaapi::Ledger led(*g);
      Digest d;
      for (std::size_t i = 0; i < pairs.size(); ++i) d.add_d(led.place(pairs[i].first, pairs[i].second, (*dags)[i], amount, k));
      return d.h;
    });
  };
  auto clos_pairs = [] { std::vector<std::pair<int, int>> p; for (int i = 0; i < 64; ++i) p.push_back({i, (i + 32) % 64}); return p; };
  auto grid_pairs = [] {
    std::vector<std::pair<int, int>> p;
    const int side = 24;
    for (int i = 0; i < 32; ++i) { int r0 = i % side, c0 = (i * 7) % side; p.push_back({r0 * side + c0, (side - 1 - r0) * side + (side - 1 - c0)}); }
    return p;
  };
  r.push_back({"kernel.proportional.clos", "kernel",
               "Optimized placement on fixed ShortestPaths DAGs, Clos 64x32, 64 pairs of 50 on one ledger",
               [=] { return kernel_row(gen_clos(64, 32, 13), clos_pairs(), 50.0, Kernel::Optimized); }});
  r.push_back({"kernel.eb.clos", "kernel",
               "Equal placement on fixed ShortestPaths DAGs, Clos 64x32, 64 pairs of 50 on one ledger",
               [=] { return kernel_row(gen_clos(64, 32, 13), clos_pairs(), 50.0, Kernel::Equal); }});
  r.push_back({"kernel.lossy.clos", "kernel",
               "Equal placement with Lossy admission on fixed ShortestPaths DAGs, Clos 64x32, 64 pairs of 50 on one ledger",
               [=] { return kernel_row(gen_clos(64, 32, 13), clos_pairs(), 50.0, Kernel::Lossy); }});
  r.push_back({"kernel.proportional.grid", "kernel",
               "Optimized placement on fixed ShortestPaths DAGs, 24x24 unit grid, 32 corner pairs of 20",
               [=] { return kernel_row(gen_grid(24, 14, true), grid_pairs(), 20.0, Kernel::Optimized); }});
  r.push_back({"kernel.eb.grid", "kernel",
               "Equal placement on fixed ShortestPaths DAGs, 24x24 unit grid, 32 corner pairs of 20",
               [=] { return kernel_row(gen_grid(24, 14, true), grid_pairs(), 20.0, Kernel::Equal); }});
  r.push_back({"kernel.weighted.clos", "kernel",
               "Weighted placement on fixed ShortestPaths DAGs, Clos 64x32, 64 pairs of 50 on one ledger",
               [=] { return kernel_row(gen_clos(64, 32, 13), clos_pairs(), 50.0, Kernel::Weighted); }});
  // ---- scalar ShortestPaths capacity across unique failure states and pairs
  auto capacity_row = [](Kernel k) {
    auto a = gen_clos(32, 16, 15);
    auto g = std::make_shared<Graph>(abaapi::build(a));
    auto states = std::make_shared<std::vector<Elig>>();
    for (auto& fs : failure_states(a.num_edges() / 2, 32, 3, 16)) states->push_back(abaapi::eligibility(*g, {}, fs));
    std::vector<std::pair<int, int>> pairs;
    for (int i = 0; i < 8; ++i) pairs.push_back({i, i + 16});
    return Iter([g, states, pairs, k] {
      Digest d;
      for (auto& el : *states) for (auto [s, t] : pairs) d.add_d(abaapi::shortest_paths_capacity(*g, s, t, k, el));
      return d.h;
    });
  };
  r.push_back({"capacity.costonly_eb.failures", "capacity",
               "scalar ShortestPaths Equal capacity, Clos 32x16, 32 unique failure states (3 links) x 8 pairs",
               [=] { return capacity_row(Kernel::Equal); }});
  r.push_back({"capacity.costonly_prop.failures", "capacity",
               "scalar ShortestPaths Optimized capacity, Clos 32x16, 32 unique failure states (3 links) x 8 pairs",
               [=] { return capacity_row(Kernel::Optimized); }});
  // ---- exact max-flow, totals only versus detailed
  auto exact_row = [](bool detailed) {
    auto a = gen_random(1024, 6, 17);
    auto g = std::make_shared<Graph>(abaapi::build(a));
    auto states = std::make_shared<std::vector<Elig>>();
    std::mt19937 rng(18);
    for (int f = 0; f < 8; ++f) {
      std::vector<int> xe;
      for (int k = 0; k < 16; ++k) xe.push_back(int(rng() % unsigned(a.num_edges())));
      states->push_back(abaapi::eligibility(*g, {}, xe));
    }
    std::vector<std::pair<int, int>> pairs = {{1, 513}, {100, 900}, {300, 7}, {777, 222}};
    return Iter([g, states, pairs, detailed] {
      Digest d;
      for (auto& el : *states) for (auto [s, t] : pairs) d.add(abaapi::maxflow_exact(*g, s, t, el, detailed));
      return d.h;
    });
  };
  r.push_back({"maxflow.exact.totals", "maxflow",
               "exact max-flow totals only, random 1024 x6 weighted, 8 failure states x 4 pairs",
               [=] { return exact_row(false); }});
  r.push_back({"maxflow.exact.detailed", "maxflow",
               "exact max-flow with edge flows, residuals and bottleneck edges, same inputs as totals",
               [=] { return exact_row(true); }});
  // ---- the failure study (evaluate_edge_failures)
  for (int workers : {1, 2}) {
    r.push_back({workers == 1 ? "sens.exact" : "sens.exact.w2", "sensitivity",
                 "exact-model sensitivity (candidates: the positive-flow edges), random 256 x6, 4 pairs",
                 [=] {
                   auto g = std::make_shared<Graph>(abaapi::build(gen_random(256, 6, 19)));
                   auto el = std::make_shared<Elig>(abaapi::eligibility(*g, {}, {}));
                   std::vector<std::pair<int, int>> pairs = {{1, 129}, {30, 200}, {77, 5}, {250, 60}};
                   return Iter([g, el, pairs, workers] { Digest d; for (auto [s, t] : pairs) d.add(abaapi::edge_failures_exact(*g, s, t, *el, workers)); return d.h; });
                 }});
  }
  // Fixed cut/repair populations. Every timed iteration includes the full
  // baseline and complete SensitivityResult, including all classification and copies.
  auto edge_failure_population = [](int kind, int workers = 1) {
    Arrays a;
    const int blocks = 32;
    auto add = [&](int u, int v, double capacity, std::int64_t cost) {
      a.src.push_back(u); a.dst.push_back(v); a.cap.push_back(capacity); a.cost.push_back(cost);
    };
    if (kind == 0 || kind == 3) {
      a.n = 65;
      for (int u = 0; u < 64; ++u) add(u, u + 1, 3, 1 + u % 3);
    } else if (kind == 1 || kind == 4) {
      a.n = blocks + 2;
      for (int i = 0; i < blocks; ++i) {
        add(0, i + 1, 3, 1);
        if (kind == 4) {
          add(0, i + 1, 3, 3);
          add(i + 1, a.n - 1, 3, 2);
        } else {
          add(i + 1, a.n - 1, 3, 2);
          add(i + 1, a.n - 1, 1, 5);
        }
      }
    } else {
      a.n = 3 * blocks + 2;
      for (int i = 0; i < blocks; ++i) {
        const int u = 3 * i + 1, v = u + 1, w = u + 2;
        add(0, u, 3, 1); add(0, v, 3, 3);
        add(u, w, 3, 1); add(v, w, 3, 1); add(w, a.n - 1, 3, 1);
      }
    }
    auto g = std::make_shared<Graph>(abaapi::build(a));
    return Iter([g, destination = a.n - 1, kind, workers] {
      return abaapi::edge_failures_unmasked(*g, destination, kind == 3 ? std::optional<double>(1) : std::nullopt, workers);
    });
  };
  r.push_back({"sens.exact.cut", "sensitivity", "unrestricted optimized edge-failure evaluation, 64-edge nonuniform-cost saturated chain; full baseline/results", [=] { return edge_failure_population(0); }});
  r.push_back({"sens.exact.partial", "sensitivity", "unrestricted optimized edge-failure evaluation, 32 branches with one of three deleted units reroutable; full baseline/results", [=] { return edge_failure_population(1); }});
  r.push_back({"sens.exact.branch", "sensitivity", "32 branches; upstream deletion reroutes all 3 units through a second branch (first repair, zero loss), shared downstream is a cut; full baseline/results", [=] { return edge_failure_population(2); }});
  r.push_back({"sens.exact.reroute", "sensitivity", "32 branches with a parallel upstream member (first repair, zero loss) and a cut downstream link; full baseline/results", [=] { return edge_failure_population(4); }});
  r.push_back({"sens.exact.capped", "sensitivity", "same 64-edge saturated chain at binding offer 1, unchanged fresh-solve fallback; full baseline/results", [=] { return edge_failure_population(3); }});
  r.push_back({"sens.exact.cut.w2", "sensitivity", "unrestricted optimized edge-failure evaluation, 64-edge nonuniform-cost saturated chain; full baseline/results", [=] { return edge_failure_population(0, 2); }});
  r.push_back({"sens.exact.partial.w2", "sensitivity", "unrestricted optimized edge-failure evaluation, 32 branches with one of three deleted units reroutable; full baseline/results", [=] { return edge_failure_population(1, 2); }});
  r.push_back({"sens.exact.branch.w2", "sensitivity", "32 branches; upstream deletion reroutes all 3 units through a second branch (first repair, zero loss), shared downstream is a cut; full baseline/results", [=] { return edge_failure_population(2, 2); }});
  r.push_back({"sens.exact.reroute.w2", "sensitivity", "32 branches with a parallel upstream member (first repair, zero loss) and a cut downstream link; full baseline/results", [=] { return edge_failure_population(4, 2); }});
  r.push_back({"sens.exact.capped.w2", "sensitivity", "same 64-edge saturated chain at binding offer 1, unchanged fresh-solve fallback; full baseline/results", [=] { return edge_failure_population(3, 2); }});
  for (int workers : {1, 2}) {
    r.push_back({workers == 1 ? "sens.policy_all_edges" : "sens.policy_all_edges.w2", "sensitivity",
                 "policy-model sensitivity over every eligible edge, ShortestPaths EB, Clos 16x8 (256 edges), 2 pairs",
                 [=] {
                   auto g = std::make_shared<Graph>(abaapi::build(gen_clos(16, 8, 20)));
                   auto el = std::make_shared<Elig>(abaapi::eligibility(*g, {}, {}));
                   return Iter([g, el, workers] {
                     Digest d;
                     d.add(abaapi::edge_failures_policy_all_edges(*g, 0, 8, *el, workers));
                     d.add(abaapi::edge_failures_policy_all_edges(*g, 3, 12, *el, workers));
                     return d.h;
                   });
                 }});
  }
  // ---- pooled sources and the equal source allocation (several sources of one demand)
  // Eight demands of eight source leaves each to one destination leaf on a Clos 64x32 (unit costs,
  // links of 10..80): demand j has the leaves 8j..8j+7 and the leaf (8j + 36) mod 64. Every source
  // splits over its 32 spines and the eight sources share each spine's one link into the
  // destination: the pooled program, its attribution and (Lossless) its bisection. One ledger.
  auto pooled_row = [](double offer, bool lossy) {
    auto g = std::make_shared<Graph>(abaapi::build(gen_clos(64, 32, 22)));
    std::vector<std::pair<std::vector<int>, int>> demands;
    for (int j = 0; j < 8; ++j) {
      std::vector<int> srcs;
      for (int i = 0; i < 8; ++i) srcs.push_back(8 * j + i);
      demands.push_back({srcs, (8 * j + 36) % 64});
    }
    return Iter([g, demands, offer, lossy] { return abaapi::pooled_place(*g, demands, offer, lossy); });
  };
  r.push_back({"pooled.ecmp.clos", "pooled",
               "pooled ShortestPaths Equal Lossless placement, Clos 64x32, 8 demands of 8 source leaves to one leaf, "
               "400 each on one ledger (the program, attribution and bisection over 8 x 32 + 32 links per demand)",
               [=] { return pooled_row(400.0, false); }});
  r.push_back({"pooled.lossy.clos", "pooled",
               "pooled ShortestPaths Equal Lossy placement, the same Clos and demands, 4,000 each (above capacity: "
               "drops attributed by highest averages on every shared spine link)",
               [=] { return pooled_row(4000.0, true); }});
  // The same demands forwarded by Pinned over from_dag programs of their reverse shortest-path DAGs
  // (built once, before timing: the one pooled builder ShortestPaths Equal uses), and the building
  // of those programs alone.
  auto pooled_dag_row = [](double offer) {
    auto g = std::make_shared<Graph>(abaapi::build(gen_clos(64, 32, 22)));
    std::vector<std::pair<std::vector<int>, int>> demands;
    for (int j = 0; j < 8; ++j) {
      std::vector<int> srcs;
      for (int i = 0; i < 8; ++i) srcs.push_back(8 * j + i);
      demands.push_back({srcs, (8 * j + 36) % 64});
    }
    auto programs = std::make_shared<const abaapi::PooledPrograms>(abaapi::pooled_dag_programs(*g, demands));
    return Iter([g, demands, programs, offer] { return abaapi::pooled_pinned_place(*g, demands, *programs, offer, false); });
  };
  r.push_back({"pooled.from_dag.clos", "pooled",
               "Pinned over from_dag programs (built before timing) of the same Clos 64x32 demands, Equal Lossless, "
               "400 each on one ledger: the forwarding and bisection of pooled.ecmp.clos without the route search",
               [=] { return pooled_dag_row(400.0); }});
  r.push_back({"pooled.from_dag_build.clos", "pooled",
               "building the from_dag programs of the same 8 demands from their reverse shortest-path DAGs (the "
               "search, trimming and the pooled plan; reported, no budget)",
               [] {
                 auto g = std::make_shared<Graph>(abaapi::build(gen_clos(64, 32, 22)));
                 std::vector<std::pair<std::vector<int>, int>> demands;
                 for (int j = 0; j < 8; ++j) {
                   std::vector<int> srcs;
                   for (int i = 0; i < 8; ++i) srcs.push_back(8 * j + i);
                   demands.push_back({srcs, (8 * j + 36) % 64});
                 }
                 return Iter([g, demands] { return abaapi::pooled_dag_build(*g, demands); });
               }});
  // ---- several demands at once (apply_actions): the eight pooled Lossy demands above given 4,000
  // each, as eight scalar replace() calls (the reference), as one ORDERED list (the same work: an
  // equal digest, budget +3 %) and as one SIMULTANEOUS list (every pooled demand forwarded by the
  // joint evaluator at once, budget 2x the ORDERED list; the demands share no link, so its digest
  // equals too).
  auto batch_row = [](abaapi::Batch mode) {
    auto g = std::make_shared<Graph>(abaapi::build(gen_clos(64, 32, 22)));
    std::vector<std::pair<std::vector<int>, int>> demands;
    for (int j = 0; j < 8; ++j) {
      std::vector<int> srcs;
      for (int i = 0; i < 8; ++i) srcs.push_back(8 * j + i);
      demands.push_back({srcs, (8 * j + 36) % 64});
    }
    return Iter([g, demands, mode] { return abaapi::batch_place(*g, demands, 4000.0, mode); });
  };
  r.push_back({"batch.scalar.clos", "batch",
               "the eight pooled Lossy demands of pooled.lossy.clos on a fresh ledger, then eight scalar replace(4000) "
               "calls (the reference of the list rows)",
               [=] { return batch_row(abaapi::Batch::Scalar); }});
  r.push_back({"batch.ordered.clos", "batch",
               "the same construction, then one ORDERED apply_actions of the eight Replace(4000): the scalar calls' work "
               "on a working state, one publication (digest equal to batch.scalar.clos)",
               [=] { return batch_row(abaapi::Batch::Ordered); }});
  r.push_back({"batch.simultaneous.clos", "batch",
               "the same construction, then one SIMULTANEOUS apply_actions of the eight Replace(4000): every program "
               "prepared on one view and forwarded at once by the joint loss evaluator (the demands share no link)",
               [=] { return batch_row(abaapi::Batch::Simultaneous); }});
  // The ring (abaapi::ring_fixture): 64 barriers in one cyclic group settling at round 6, as one
  // program (the reference) and as a SIMULTANEOUS list of its 64 tunnels (equal traffic amounts and
  // digests; the list has 64 publishing policies rather than one, with a 2x operating budget).
  r.push_back({"program.ring.lossy", "program",
               "one Pinned Lossy replace(192) of a 64-head program on a ring of 64 links of 3 units, each tunnel over two "
               "consecutive links: 64 barriers in one cyclic group, 6 rounds",
               [] {
                 auto ring = abaapi::ring_fixture();
                 return Iter([ring] { return abaapi::ring_program(*ring); });
               }});
  r.push_back({"batch.simultaneous.ring", "batch",
               "the same ring as 64 one-tunnel Pinned Lossy policies, one SIMULTANEOUS apply_actions of Replace(3) each: "
               "the same 64 barriers and 6 rounds across 64 demands (digest equal to program.ring.lossy)",
               [] {
                 auto ring = abaapi::ring_fixture();
                 return Iter([ring] { return abaapi::ring_list(*ring); });
               }});
  r.push_back({"global.equal.random", "global",
               "evaluate_demand() of CentralizedTE with the equal source allocation, random 1024 x6, 2 demands of 4 sources (the bisection's "
               "probe solves and the final solve)",
               [] {
                 auto g = std::make_shared<Graph>(abaapi::build(gen_random(1024, 6, 23)));
                 const std::vector<std::pair<std::vector<int>, int>> demands = {{{1, 257, 513, 769}, 600},
                                                                                {{50, 300, 550, 800}, 999}};
                 return Iter([g, demands] {
                   Digest d;
                   for (const auto& [srcs, dst] : demands) d.add(abaapi::global_equal_capacity(*g, srcs, dst));
                   return d.h;
                 });
               }});
  r.push_back({"sens.equal_global", "sensitivity",
               "evaluate_edge_failures() of CentralizedTE with the equal source allocation (every edge a candidate, each solved afresh with its "
               "bisection), random 64 x4 (256 edges), 3 sources to one destination",
               [] {
                 auto g = std::make_shared<Graph>(abaapi::build(gen_random(64, 4, 24)));
                 const std::vector<int> srcs = {1, 17, 33};
                 return Iter([g, srcs] { return abaapi::edge_failures_equal_global(*g, srcs, 50); });
               }});
  // ---- forwarding programs: no older Core has them; they guard later changes.
  // One demand of eight source leaves to one leaf on a Clos 32x16 (links of 10..80): its optimized
  // CentralizedTE program (paths with solved amounts, captured source weights), solved for 2,000.
  auto program_demand = [] {
    std::vector<int> srcs;
    for (int i = 0; i < 8; ++i) srcs.push_back(i);
    return srcs;
  };
  r.push_back({"program.forward.global", "program",
               "forward() of an installed optimized CentralizedTE program, Clos 32x16, 8 sources to one leaf, 32 offers "
               "of 100..3,200 (the injection by captured weights and the per-tunnel split, no solve)",
               [=] {
                 auto g = std::make_shared<Graph>(abaapi::build(gen_clos(32, 16, 25)));
                 auto inst = std::make_shared<abaapi::Installed>(*g, program_demand(), 20, 2000.0);
                 std::vector<double> offers;
                 for (int i = 1; i <= 32; ++i) offers.push_back(100.0 * i);
                 return Iter([g, inst, offers] { return inst->forward_all(offers); });
               }});
  auto failure_row = [=](bool stale) {
    auto a = gen_clos(32, 16, 25);
    auto g = std::make_shared<Graph>(abaapi::build(a));
    auto prog = std::make_shared<abaapi::Program>(abaapi::global_program(*g, program_demand(), 20, 2000.0));
    auto states = std::make_shared<abaapi::Masks>(abaapi::masks(*g, failure_states(a.num_edges() / 2, 32, 3, 26)));
    return Iter([g, prog, states, stale] { return abaapi::pinned_failures(*g, *prog, *states, 2400.0, stale); });
  };
  r.push_back({"program.replay.lossy", "program",
               "Pinned Lossy replace() of that CentralizedTE program under 32 failure states (3 links) as the mask, each on a "
               "fresh ledger: failed tunnels, the heads' shares over the survivors, the loss rule",
               [=] { return failure_row(false); }});
  r.push_back({"program.stale.lossy", "program",
               "the same with the failure states as black holes (stale forwarding: the program unchanged, traffic "
               "dropped at the failures)",
               [=] { return failure_row(true); }});
  r.push_back({"program.replace.costonly", "program",
               "ShortestPaths Equal replace() of 8 source leaves to one leaf, Clos 32x16, under 16 changing loads of "
               "an outside demand (static preparation retained across replacements, the released view and retiring commit)",
               [=] {
                 auto g = std::make_shared<Graph>(abaapi::build(gen_clos(32, 16, 25)));
                 std::vector<double> loads;
                 for (int i = 0; i < 16; ++i) loads.push_back(50.0 * (i % 8));
                 return Iter([g, loads, srcs = program_demand()] {
                   return abaapi::replace_under_load(*g, srcs, 20, {9, 10, 11, 12}, 21, loads, 800.0);
                 });
               }});
  // ---- T1 tunnels (the fabric of 4 pods, 12 demands of 120: one policy per demand, one ledger)
  auto t1_pairs = [](const Fabric3& f) {
    std::vector<std::pair<int, int>> pairs;
    for (int p = 0; p < 4; ++p) for (int i = 0; i < 3; ++i) pairs.push_back({fabric3_leaf(f, p, i), fabric3_leaf(f, (p + 1) % 4, 7 - i)});
    return pairs;
  };
  r.push_back({"t1.setup.fabric", "t1",
               "T1 DistributedTE(TunnelCount(16)) replace(120, design 120) per demand, fabric 4 pods, 12 demands on one "
               "ledger: CSPF setup on the unreserved view, the holds' commit, Lossless forwarding",
               [=] {
                 auto f = gen_fabric3(4, 8, 4, 8, 21);
                 auto g = std::make_shared<Graph>(abaapi::build(f.a));
                 return Iter([g, pairs = t1_pairs(f)] { return abaapi::t1_setup(*g, pairs, 120.0, 16); });
               }});
  r.push_back({"t1.forward.fabric", "t1",
               "the same 12 T1 programs (Lossy) set up once, then forward() of each under 32 offers (no setup, the holds "
               "kept)",
               [=] {
                 auto f = gen_fabric3(4, 8, 4, 8, 21);
                 auto g = std::make_shared<Graph>(abaapi::build(f.a));
                 auto fabric = abaapi::t1_installed(*g, t1_pairs(f), 120.0, 16);
                 std::vector<double> offers;
                 for (int i = 0; i < 32; ++i) offers.push_back(40.0 + 10.0 * i);
                 return Iter([g, fabric, offers] { return abaapi::t1_forward(*fabric, offers); });
               }});
  // ---- T2 tunnels (the same fabric and demands; and a mesh: a bidirectional ring plus directed
  // chords, 4n links, full capacities 2^14..2^17 units, 5/8 to 8/8 of each free, costs 1..10, heads
  // 0..h-1, destination n-1)
  auto t2_row = [=](int count, bool up_to, bool headroom) {
    auto f = gen_fabric3(4, 8, 4, 8, 21);
    auto g = std::make_shared<Graph>(abaapi::build(f.a));
    return Iter([g, pairs = t1_pairs(f), count, up_to, headroom] {
      return abaapi::t2_fabric(*g, pairs, 120.0, count, up_to, headroom);
    });
  };
  r.push_back({"t2.exact16.headroom", "t2",
               "T2 CentralizedTE(tunnel_count=TunnelCount(16), prefer_headroom=True), replace(120, design 120) per demand, fabric 4 pods, 12 demands "
               "on one ledger: sizing, the canonical augmentation with Phi, the lossless bound, Lossless forwarding",
               [=] { return t2_row(16, false, true); }});
  r.push_back({"t2.exact16.cost", "t2",
               "the same without prefer_headroom: sizing and the canonical augmentation without Phi",
               [=] { return t2_row(16, false, false); }});
  r.push_back({"t2.upto256.headroom", "t2",
               "T2 CentralizedTE(tunnel_count=TunnelCount(256, up_to=True), prefer_headroom=True), the same demands: count selection over 1..256 with "
               "the proved bounds, then the chosen count's layout",
               [=] { return t2_row(256, true, true); }});
  r.push_back({"t2.size.mesh", "t2",
               "T2 TunnelCount(16) without prefer_headroom, replace(+inf) of 10 heads on a 500-node mesh (2000 links, 5/8 to 8/8 free): "
               "sizing-dominated",
               [] {
                 auto m = gen_t2_mesh(500, 10, 1744);
                 auto g = std::make_shared<Graph>(abaapi::build(m.a));
                 return Iter([g, m] { return abaapi::t2_mesh(*g, m.free, m.heads, 499, 16, false, false); });
               }});
  // ---- MSD loop
  r.push_back({"msd.loop", "msd",
               "MSD bracket and bisection over ShortestPaths EB placement with a per-series DAG cache, Clos 32x16, 64 demands",
               [] {
                 auto g = std::make_shared<Graph>(abaapi::build(gen_clos(32, 16, 22)));
                 auto el = std::make_shared<Elig>(abaapi::eligibility(*g, {}, {}));
                 std::mt19937 rng(23);
                 std::vector<std::pair<int, int>> dem;
                 std::vector<double> base;
                 while (dem.size() < 64) {
                   int s = int(rng() % 32), t = int(rng() % 32);
                   if (s != t) { dem.push_back({s, t}); base.push_back(double(1 + rng() % 5)); }
                 }
                 return Iter([g, el, dem, base] {
                   abaapi::ShortestPathsSeries series(*g, *el);
                   int probes = 0;
                   auto feasible = [&](double alpha) {
                     ++probes;
                     std::vector<double> v(base.size());
                     for (std::size_t i = 0; i < v.size(); ++i) v[i] = base[i] * alpha;
                     auto delivered = series.place_all(dem, v);
                     for (std::size_t i = 0; i < v.size(); ++i) if (v[i] - delivered[i] > abaapi::grid_resolution()) return false;
                     return true;
                   };
                   // NetGraph MaximumSupportedDemand: start 1, growth 2, resolution 0.01.
                   double lo = 0, hi = 0, alpha = 1;
                   if (feasible(alpha)) { lo = alpha; for (int k = 0; k < 32; ++k) { alpha *= 2; if (!feasible(alpha)) { hi = alpha; break; } lo = alpha; } }
                   else { hi = alpha; for (int k = 0; k < 32; ++k) { alpha /= 2; if (feasible(alpha)) { lo = alpha; break; } hi = alpha; } }
                   for (int k = 0; k < 32 && hi - lo > 0.01; ++k) { double mid = (lo + hi) / 2; (feasible(mid) ? lo : hi) = mid; }
                   Digest d;
                   d.add_d(lo);
                   d.add_i(probes);
                   return d.h;
                 });
               }});
  r.push_back({"msd.loop.epochs", "msd",
               "MSD bracket and bisection over ShortestPaths EB placement, one ledger and 64 policies per series, remove all before each probe, Clos 32x16",
               [] {
                 auto g = std::make_shared<Graph>(abaapi::build(gen_clos(32, 16, 22)));
                 auto el = std::make_shared<Elig>(abaapi::eligibility(*g, {}, {}));
                 std::mt19937 rng(23);
                 std::vector<std::pair<int, int>> dem;
                 std::vector<double> base;
                 while (dem.size() < 64) {
                   int s = int(rng() % 32), t = int(rng() % 32);
                   if (s != t) { dem.push_back({s, t}); base.push_back(double(1 + rng() % 5)); }
                 }
                 return Iter([g, el, dem, base] {
                   abaapi::ShortestPathsEpochSeries series(*g, *el, dem);
                   int probes = 0;
                   auto feasible = [&](double alpha) {
                     ++probes;
                     std::vector<double> v(base.size());
                     for (std::size_t i = 0; i < v.size(); ++i) v[i] = base[i] * alpha;
                     auto delivered = series.place_all(v);
                     for (std::size_t i = 0; i < v.size(); ++i) if (v[i] - delivered[i] > abaapi::grid_resolution()) return false;
                     return true;
                   };
                   // NetGraph MaximumSupportedDemand: start 1, growth 2, resolution 0.01.
                   double lo = 0, hi = 0, alpha = 1;
                   if (feasible(alpha)) { lo = alpha; for (int k = 0; k < 32; ++k) { alpha *= 2; if (!feasible(alpha)) { hi = alpha; break; } lo = alpha; } }
                   else { hi = alpha; for (int k = 0; k < 32; ++k) { alpha /= 2; if (feasible(alpha)) { lo = alpha; break; } hi = alpha; } }
                   for (int k = 0; k < 32 && hi - lo > 0.01; ++k) { double mid = (lo + hi) / 2; (feasible(mid) ? lo : hi) = mid; }
                   Digest d;
                   d.add_d(lo);
                   d.add_i(probes);
                   return d.h;
                 });
               }});
  // ---- k shortest paths, the DAG cache, eligibility
  auto prefix_population = [](int stem, int width, bool ties, bool parallels) {
    Arrays a;
    a.n = stem + width + 2;
    auto add = [&](int u, int v, double cap, std::int64_t cost) {
      a.src.push_back(u); a.dst.push_back(v); a.cap.push_back(cap); a.cost.push_back(cost);
    };
    for (int v = 0; v < stem; ++v) {
      add(v, v + 1, 100, 1);
      if (parallels) add(v, v + 1, 50, 1);
    }
    for (int i = 0; i < width; ++i) {
      add(stem, stem + i + 1, 1 + i % 9, 1);
      add(stem + i + 1, a.n - 1, 1 + i % 9, ties ? 1 : 1 + i % 9);
    }
    return a;
  };
  r.push_back({"ksp.longprefix_k64", "ksp",
               "96-edge common prefix then 64 two-edge branches, costs 2..10, k=64, FEWER_HOPS Expand; complete ordered paths",
               [=] {
                 auto g = std::make_shared<Graph>(abaapi::build(prefix_population(96, 64, false, false)));
                 return Iter([g] { return abaapi::ksp(*g, 0, 161, 64); });
               }});
  r.push_back({"ksp.weighted_k32", "ksp",
               "random 128 x6 seed64, costs1..31, k=32, four endpoint pairs, FEWER_HOPS Expand",
               [] {
                 auto g = std::make_shared<Graph>(abaapi::build(gen_random(128, 6, 64)));
                 return Iter([g] {
                   Digest d;
                   for (int s : {0, 16, 32, 48}) d.add(abaapi::ksp(*g, s, 127 - s, 32));
                   return d.h;
                 });
               }});
  r.push_back({"ksp.collapse", "ksp",
               "12-edge common prefix with two parallel members per hop, 32 branches, k=32, HIGHER_CAPACITY Collapse",
               [=] {
                 auto g = std::make_shared<Graph>(abaapi::build(prefix_population(12, 32, false, true)));
                 return Iter([g] { return abaapi::ksp(*g, 0, 45, 32, false, abaapi::KspTie::Capacity); });
               }});
  r.push_back({"ksp.ties", "ksp",
               "three-edge prefix and 64 equal-cost branches with varied capacities, k=64, each of the four ties, Expand",
               [=] {
                 auto g = std::make_shared<Graph>(abaapi::build(prefix_population(3, 64, true, false)));
                 return Iter([g] {
                   Digest d;
                   for (auto tie : {abaapi::KspTie::Fewer, abaapi::KspTie::Capacity, abaapi::KspTie::Available, abaapi::KspTie::Headroom})
                     d.add(abaapi::ksp(*g, 0, 68, 64, true, tie));
                   return d.h;
                 });
               }});
  r.push_back({"ksp.k1", "ksp",
               "same 96-edge prefix/64 branches as ksp.longprefix_k64, k=1, unchanged first-search-only work",
               [=] {
                 auto g = std::make_shared<Graph>(abaapi::build(prefix_population(96, 64, false, false)));
                 return Iter([g] { return abaapi::ksp(*g, 0, 161, 1); });
               }});
  r.push_back({"ksp.fabric3_k8", "ksp",
               "K shortest paths, k = 8, three-tier fabric 8 pods x 16 leaves, 8 pairs",
               [] {
                 auto f = gen_fabric3(8, 16, 8, 16, 24);
                 auto g = std::make_shared<Graph>(abaapi::build(f.a));
                 std::vector<std::pair<int, int>> pairs;
                 for (int p = 1; p < 8; ++p) pairs.push_back({fabric3_leaf(f, 0, p % 4), fabric3_leaf(f, p, 0)});
                 pairs.push_back({fabric3_leaf(f, 0, 0), fabric3_leaf(f, 0, 9)});
                 return Iter([g, pairs] { Digest d; for (auto [s, t] : pairs) d.add(abaapi::ksp(*g, s, t, 8)); return d.h; });
               }});
  auto cache_row = [](bool warm) {
    auto g = std::make_shared<Graph>(abaapi::build(gen_clos(64, 32, 25)));
    std::vector<int> roots;
    for (int i = 0; i < 64; ++i) roots.push_back(i);
    abaapi::dag_cache_series(*g, {}, warm);  // setup: creates the warm row's cache
    return Iter([g, roots, warm] { return abaapi::dag_cache_series(*g, roots, warm); });
  };
  r.push_back({"cache.cold", "cache", "ShortestPathCache series from a cleared cache: 64 roots, Clos 64x32",
               [=] { return cache_row(false); }});
  r.push_back({"cache.warm", "cache", "ShortestPathCache series on a prepopulated cache: 64 roots, Clos 64x32",
               [=] { return cache_row(true); }});
  r.push_back({"elig.build_intersect", "eligibility",
               "Eligibility build for 64 failure states and intersect with a common state, Clos 64x32",
               [] {
                 auto a = gen_clos(64, 32, 26);
                 auto g = std::make_shared<Graph>(abaapi::build(a));
                 auto states = failure_states(a.num_edges() / 2, 64, 4, 27);
                 std::vector<int> common = {0, 1, 2, 3};
                 abaapi::eligibility_build_intersect(*g, common, {});
                 return Iter([g, states, common] { return abaapi::eligibility_build_intersect(*g, common, states); });
               }});
  for (const auto& [name, operation] : std::array<std::pair<const char*, abaapi::TinyOperation>, 4>{{
      {"concurrency.place_zero", abaapi::TinyOperation::PlaceZero},
      {"concurrency.delivered", abaapi::TinyOperation::Delivered},
      {"concurrency.program", abaapi::TinyOperation::Program},
      {"concurrency.config", abaapi::TinyOperation::Config}}}) {
    r.push_back({name, "concurrency", "10,000 calls of one operation on a prepared one-link LOSSLESS policy (4 of 16 units)",
      [operation] {
        auto fixture = std::make_shared<abaapi::TinyCalls>();
        return Iter([fixture, operation] { return fixture->run(operation); });
      }});
  }
  r.push_back({"cache.frozen.reads", "cache", "64 direct frozen-cache hits, Clos 64x32, same roots as cache.warm",
    [] {
      auto g = std::make_shared<Graph>(abaapi::build(gen_clos(64, 32, 25)));
      std::vector<int> roots;
      for (int i = 0; i < 64; ++i) roots.push_back(i);
      (void)abaapi::dag_cache_series(*g, roots, true);
      g->warm_cache->freeze();
      return Iter([g, roots] { return abaapi::dag_cache_series(*g, roots, true); });
    }});
  for (bool parallel : {false, true}) {
    r.push_back({parallel ? "concurrency.two_ledgers_parallel" : "concurrency.two_ledgers_serial", "concurrency",
      "Two independent ledgers, 32 evaluate(32) calls each using static preparation installed in setup; two thread launches/joins in both modes; layered graph 384 nodes/2816 edges",
      [parallel] {
        auto fixture = std::make_shared<abaapi::ParallelLedgers>();
        return Iter([fixture, parallel] { return fixture->run(parallel); });
      }});
  }
  add_storage_rows(r);  // the storage rows come last in the list
  return r;
}

// ---------------------------------------------------------------- driver

inline double peak_rss_bytes(int who) {
  struct rusage ru{};
  getrusage(who, &ru);
#if defined(__APPLE__)
  return double(ru.ru_maxrss);  // bytes on macOS
#else
  return double(ru.ru_maxrss) * 1024.0;  // KiB on Linux
#endif
}

inline std::string json_escape(const std::string& s) {
  std::string o;
  for (char c : s) { if (c == '"' || c == '\\') o += '\\'; o += c; }
  return o;
}

using Clock = std::chrono::steady_clock;
static_assert(Clock::is_steady && std::is_integral_v<Clock::rep>, "durations are whole clock ticks, never negative");

// The durations of a block's timed iterations, in clock ticks, in memory that does not grow with
// the iteration count, so every block's peak RSS holds the same bookkeeping whatever its count
// (README "Memory"): a count per duration below kShort ticks (65.5 us with the nanosecond clocks
// of macOS and Linux), and the longer durations listed. Each listed duration lasted at least
// kShort ticks, so the list's kLong slots fill only after 4.3 s of such iterations; beyond them it
// grows. Both are allocated and written before the row's setup. The statistics are those of the
// sorted list of every duration, value for value.
class IterationTimes {
 public:
  static constexpr std::size_t kShort = std::size_t{1} << 16;
  static constexpr std::size_t kLong = std::size_t{1} << 16;

  IterationTimes() : short_(kShort, 0), long_(kLong, 0) {
    // Volatile stores, one per page: the pages are resident from here on, not from the first
    // duration that lands on them.
    for (std::vector<std::uint64_t>* v : {&short_, &long_})
      for (std::size_t i = 0; i < v->size(); i += 512) static_cast<volatile std::uint64_t*>(v->data())[i] = 0;
  }

  void add(Clock::duration d) {
    const auto t = static_cast<std::uint64_t>(d.count());
    if (t < kShort) {
      ++short_[t];
    } else {
      if (n_long_ < long_.size()) long_[n_long_] = t;
      else long_.push_back(t);
      ++n_long_;
    }
    ++n_;
  }

  std::size_t size() const { return n_; }

  // Sorts the long durations; once, after the last add.
  void finish() { std::sort(long_.begin(), long_.begin() + std::ptrdiff_t(n_long_)); }

  // The k-th smallest duration (k < size(), from 0) in seconds, converted as one duration is.
  double at(std::size_t k) const {
    std::size_t seen = 0;
    for (std::size_t t = 0; t < kShort; ++t) {
      seen += short_[t];
      if (k < seen) return seconds(t);
    }
    return seconds(long_[k - seen]);
  }

  // The sum of every duration in seconds, the smallest added first.
  double sum() const {
    double s = 0;
    for (std::size_t t = 0; t < kShort; ++t)
      for (std::uint64_t c = short_[t]; c > 0; --c) s += seconds(t);
    for (std::size_t i = 0; i < n_long_; ++i) s += seconds(long_[i]);
    return s;
  }

 private:
  static double seconds(std::uint64_t t) {
    return std::chrono::duration<double>(Clock::duration(Clock::rep(t))).count();
  }

  std::vector<std::uint64_t> short_;  // the count of every duration below kShort ticks
  std::vector<std::uint64_t> long_;   // the longer durations, the first n_long_ entries
  std::size_t n_long_ = 0;
  std::size_t n_ = 0;
};

inline const char* query_list_of(const std::string& row) {
  for (const auto& [name, list] : kQueryLists)
    if (row == name) return list;
  return nullptr;
}

// The harness's command line over a row list (harness.cpp passes rows()):
//
//   aba-list | aba-run <row> <iterations> | aba-calibrate <row> <seconds>
//
// Exit status: 0 success; 1 a failure inside an iteration, the untimed first one included (the
// message on stderr after "FAIL: ");
// 2 a usage error (an unknown mode, a missing argument, an iteration count that is not a positive
// integer, seconds that are not a positive finite number) or an unknown row; 3 a digest that
// changed between iterations; 4 a row whose input or setup check failed.
inline int harness_main(int argc, char** argv, const std::vector<Row>& all) {
  const std::string mode = argc > 1 ? argv[1] : "";
  if (mode == "aba-list" && argc == 2) {
    for (const auto& row : all) {
      const char* list = query_list_of(row.name);
      std::printf("{\"row\":\"%s\",\"group\":\"%s\",\"desc\":\"%s\"%s%s%s}\n", row.name.c_str(),
                  row.group.c_str(), json_escape(row.desc).c_str(), list ? ",\"queries\":\"" : "", list ? list : "",
                  list ? "\"" : "");
    }
    return 0;
  }
  const bool calibrate = mode == "aba-calibrate";
  long long iterations = 0;
  double seconds = 0.0;
  bool usage = (mode != "aba-run" && !calibrate) || argc != 4;
  if (!usage) {
    char* end = nullptr;
    errno = 0;
    if (calibrate) {
      seconds = std::strtod(argv[3], &end);
      usage = *argv[3] == '\0' || *end != '\0' || errno != 0 || !std::isfinite(seconds) || seconds <= 0.0;
    } else {
      iterations = std::strtoll(argv[3], &end, 10);
      usage = *argv[3] == '\0' || *end != '\0' || errno != 0 || iterations < 1;
    }
  }
  if (usage) {
    std::fprintf(stderr, "usage: harness aba-list | aba-run <row> <iterations> | aba-calibrate <row> <seconds>\n");
    return 2;
  }
  const std::string name = argv[2];
  auto it = std::find_if(all.begin(), all.end(), [&](const Row& r) { return r.name == name; });
  if (it == all.end()) { std::fprintf(stderr, "unknown row %s\n", name.c_str()); return 2; }
  IterationTimes times;  // before the setup: the same memory in every block's peak
  auto t0 = Clock::now();
  Iter fn;
  try {
    fn = it->setup();
  } catch (const std::exception& e) {
    // A row whose input or setup check fails (a frozen digest, a queue assertion) exits 4.
    std::fprintf(stderr, "row %s setup failed: %s\n", name.c_str(), e.what());
    return 4;
  }
  std::uint64_t digest = 0;
  double setup_s = 0.0;  // the setup and the untimed first iteration
  auto start = Clock::now();
  try {
    digest = fn();  // untimed first iteration: the reference digest and warm state
    setup_s = std::chrono::duration<double>(Clock::now() - t0).count();
    start = Clock::now();
    for (long long i = 0; calibrate ? std::chrono::duration<double>(Clock::now() - start).count() < seconds : i < iterations; ++i) {
      auto a = Clock::now();
      const std::uint64_t got = fn();
      times.add(Clock::now() - a);
      if (got != digest) { std::fprintf(stderr, "nondeterministic digest in %s at iteration %lld\n", name.c_str(), i); return 3; }
    }
  } catch (const std::exception& e) {
    std::fprintf(stderr, "FAIL: %s\n", e.what());
    return 1;
  }
  const double block_s = std::chrono::duration<double>(Clock::now() - start).count();
  times.finish();
  const std::size_t n = times.size();
  auto q = [&](double p) { return n == 0 ? 0.0 : times.at(std::min(n - 1, std::size_t(p * double(n - 1) + 0.5))); };
  const double mean = n == 0 ? 0.0 : times.sum() / double(n);
  const double median = n == 0 ? 0.0 : (n % 2 ? times.at(n / 2) : 0.5 * (times.at(n / 2 - 1) + times.at(n / 2)));
  std::printf("{\"row\":\"%s\",\"iterations\":%zu,\"median_s\":%.9g,\"mean_s\":%.9g,"
              "\"min_s\":%.9g,\"p10_s\":%.9g,\"p90_s\":%.9g,\"max_s\":%.9g,\"block_s\":%.6f,\"setup_s\":%.6f,"
              "\"digest\":\"%016llx\",\"rss_self_bytes\":%.0f,\"rss_children_bytes\":%.0f}\n",
              name.c_str(), n, median, mean, n == 0 ? 0.0 : times.at(0), q(0.1), q(0.9),
              n == 0 ? 0.0 : times.at(n - 1), block_s, setup_s, static_cast<unsigned long long>(digest),
              peak_rss_bytes(RUSAGE_SELF), peak_rss_bytes(RUSAGE_CHILDREN));
  return 0;
}

}  // namespace abarows
