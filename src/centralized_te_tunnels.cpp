/*
  T2 layout: a controller gives every live head the same number k of
  single-path tunnels with one common whole-unit planned bandwidth.

  The network: the allowed edges (an edge and both endpoints allowed) that do not leave a
  destination (a path stops at the first destination it meets), on V, the free capacity with this
  policy's own traffic released. A head is live when it reaches a destination over edges with
  V >= 1 (one reverse traversal); the other roster heads are dead: they keep their place in the
  roster (their positions j * H + i are holes) and get no tunnel. a live heads, K = a * k tunnels.

  Sizing. Tunnels of a common size s > 0 exist exactly when K path units flow from a virtual root
  (an arc of capacity k to each live head) to the destinations with edge capacities
  min(K, floor(V / s)): the quota flow, a maximum flow of the flow core with costs ignored
  (max_flow_subset: its value capped at K, each live head's virtual arc capped at k; a
  destination's virtual arc never binds, the value being at most K). The quota network keeps
  only the edges with V >= 1 on a walk from a live head to a destination. U, the largest integer size that passes, by integer bisection below the bounds
  every feasible size obeys: max V, floor(sum of a live head's out-edges' V / k), floor(sum of the
  destinations' in-edges' V / K) (and, for "up to", U(k - 1)). The exact fluid optimum s* (a
  fraction V[e] / j with j <= K) by bisection over dyadic sizes N / 2^t, t = 2 * bit_width(K) bits
  (two such fractions differ by at least 1 / K^2 > 2^-t), from [U, U + 1) (or [2^-bit_width(K), 1)
  when U = 0, where every edge with V >= 1 has capacity K), then the one breakpoint in the final
  interval. Capacities are exact: V * 2^t < 2^123, and min(K, .) is decided by a product below 2^92
  before a division whose quotient fits 64 bits.

  Working size (finite design D): u = min(U, max(1, max over live heads of ceil(RR(D, H, i) / k)))
  with i the head's roster position (dead heads' shares stay in the round robin); unlimited: U.

  Routes at u: count vectors x with x[e] <= b[e] = min(K, floor(V[e] / u)); minimum
  J = sum cost * x, then (prefer_headroom) minimum Phi = sum over slots of floor(2^32 (B - V + j u) / B)
  (B: static capacity; j u <= V, so a slot price is an integer in [0, 2^32]), then the canonical
  unit augmentation: K units from zero flow, each on the residual root-to-sink path of minimum
  (J, Phi) price, then fewest arcs, then the smallest sequence of arc keys ((0, EdgeId, direction)
  for links, (1, head, direction) for root arcs, (2, destination, direction) for sink arcs). One
  Dijkstra per unit with pair-valued potentials (signed 128-bit pairs, compared
  lexicographically): it settles every vertex up to the sink's distance; the arcs whose reduced
  price closes a distance exactly are tight; a BFS from the sink over reversed tight arcs counts
  hops; the walk from the root takes, at each vertex, the smallest-key tight arc that lowers the
  hop count. Settled vertices add their distance to their potential, the others the sink's
  (reduced prices stay non-negative: a negative one is an internal error). Per link at most two
  residual arcs: forward at slot x + 1 (x < b), reverse at slot x (x > 0); reverse root and sink
  arcs never lie on a simple root-to-sink path and are not kept. J < K * (cost sum) < 2^92 and
  Phi <= 2^32 * sum x < 2^94 (sum x counts every link of every tunnel, so Phi can pass 2^64): 128-bit,
  each link's slot prices summed in closed form (t2_objectives).

  Program: x cut into unit paths per live head in ascending NodeId, each the walk along the
  smallest positive EdgeId (decompose(), the controllers' cut, its paths expanded into units);
  each head's k paths sorted by EdgeId sequence; tunnels in round order (round j, then roster
  position i), index j * H + i.

  Up to L: k = 1..L upwards; U(k) never increases with k (a feasible (k + 1)-layout minus one
  tunnel per head is feasible), so U(k) = 0 ends the scan (every later score is 0). With M = H k
  roster positions the score G(min(D, F'(k))) (G: the live heads' part of a round-robin total;
  F'(k): the largest total whose live shares fit on the program at U(k), dead shares ignored, in
  [M U, M (U + 1)); unlimited D: G(F'(k))) lies between G(min(D, M U)) and G(min(D, M (U + 1) - 1)):
  a lower bound equal to G(D) selects k without its program, an upper bound at most the
  incumbent's skips k; nothing else is skipped. The smallest k reaching G(D) wins, else the
  largest score, then the smallest k; it is then laid out at its own u. F'(k), M (U + 1) and the
  bounds reach 2^92 (dead heads make M larger than K): 128-bit.
*/
#include <algorithm>
#include <bit>
#include <numeric>
#include <optional>
#include <queue>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>

#include "netgraph/core/detail/min_cost_flow.hpp"
#include "netgraph/core/detail/int128.hpp"
#include "netgraph/core/detail/tunnel_setup.hpp"

namespace netgraph::core::detail {

namespace {

[[noreturn]] void internal(const char* what) { throw std::logic_error(std::string("T2: ") + what); }

U128 u128(std::uint64_t x) noexcept { return U128 {0, x}; }
U128 half(U128 v) noexcept { return U128 {v.hi >> 1, (v.lo >> 1) | (v.hi << 63)}; }
// floor(z / d) and z mod d for d >= 1 and a quotient below 2^64.
std::pair<std::uint64_t, std::uint64_t> divmod(U128 z, std::uint64_t d) noexcept {
  const std::uint64_t q = div_u128(z, d);
  return {q, sub(z, mul_u64(q, d)).lo};
}

// The edges T2 may use and the endpoints, on the free capacity V.
struct Base {
  const Graph& g;
  std::span<const NodeId> roster;
  std::span<const NodeId> dests;
  std::span<const Amount> V;
  std::vector<std::uint8_t> usable {};  // per edge: allowed and not leaving a destination
  std::vector<std::uint8_t> back {};    // per node: reaches a destination over usable edges with V >= 1
};

Base make_base(const Graph& g, const Eligibility* elig, std::span<const NodeId> roster, std::span<const NodeId> dests,
               std::span<const Amount> V) {
  Base b {g, roster, dests, V};
  const auto n = static_cast<std::size_t>(g.num_nodes());
  const auto m = static_cast<std::size_t>(g.num_edges());
  const auto src = g.edge_sources();
  const auto dst = g.edge_destinations();
  std::vector<std::uint8_t> is_dst(n, 0);
  for (NodeId d : dests) is_dst[static_cast<std::size_t>(d)] = 1;
  b.usable.assign(m, 0);
  for (std::size_t e = 0; e < m; ++e) {
    const auto u = src[e];
    const auto v = dst[e];
    if (is_dst[static_cast<std::size_t>(u)] != 0) continue;
    if (elig != nullptr && (!elig->allows_edge(static_cast<EdgeId>(e)) || !elig->allows_node(u) || !elig->allows_node(v))) continue;
    b.usable[e] = 1;
  }
  // The reverse traversal from the destinations over usable edges with V >= 1.
  b.back.assign(n, 0);
  std::vector<NodeId> stack(dests.begin(), dests.end());
  for (NodeId d : dests) b.back[static_cast<std::size_t>(d)] = 1;
  const auto in_row = g.incoming_offsets();
  const auto in_edge = g.incoming_edges();
  while (!stack.empty()) {
    const auto v = static_cast<std::size_t>(stack.back());
    stack.pop_back();
    for (auto s = in_row[v]; s < in_row[v + 1]; ++s) {
      const auto e = static_cast<std::size_t>(in_edge[static_cast<std::size_t>(s)]);
      if (b.usable[e] == 0 || V[e] < 1) continue;
      const auto u = static_cast<std::size_t>(src[e]);
      if (b.back[u] == 0) {
        b.back[u] = 1;
        stack.push_back(static_cast<NodeId>(u));
      }
    }
  }
  return b;
}

// The quota network of the live heads: the usable edges with V >= 1 on a walk from a live head to
// a destination (ascending EdgeId), and the heads' roster positions.
struct Quota {
  std::vector<NodeId> live {};
  std::vector<std::int64_t> live_pos {};  // roster position of each live head
  std::vector<EdgeId> edges {};
  std::vector<Amount> cap_v {};           // V of each edge
  Amount max_v {0};
};

Quota make_quota(const Base& b) {
  Quota q;
  for (std::size_t i = 0; i < b.roster.size(); ++i) {
    if (b.back[static_cast<std::size_t>(b.roster[i])] == 0) continue;
    q.live.push_back(b.roster[i]);
    q.live_pos.push_back(static_cast<std::int64_t>(i));
  }
  const auto& g = b.g;
  const auto n = static_cast<std::size_t>(g.num_nodes());
  const auto row = g.outgoing_offsets();
  const auto dst = g.edge_destinations();
  std::vector<std::uint8_t> fwd(n, 0);
  std::vector<NodeId> stack(q.live.begin(), q.live.end());
  for (NodeId h : q.live) fwd[static_cast<std::size_t>(h)] = 1;
  while (!stack.empty()) {
    const auto u = static_cast<std::size_t>(stack.back());
    stack.pop_back();
    for (auto e = static_cast<std::size_t>(row[u]); e < static_cast<std::size_t>(row[u + 1]); ++e) {
      if (b.usable[e] == 0 || b.V[e] < 1) continue;
      const auto v = static_cast<std::size_t>(dst[e]);
      if (b.back[v] == 0 || fwd[v] != 0) continue;
      fwd[v] = 1;
      stack.push_back(static_cast<NodeId>(v));
    }
  }
  const auto src = g.edge_sources();
  for (std::size_t e = 0; e < static_cast<std::size_t>(g.num_edges()); ++e) {
    if (b.usable[e] == 0 || b.V[e] < 1) continue;
    if (fwd[static_cast<std::size_t>(src[e])] == 0 || b.back[static_cast<std::size_t>(dst[e])] == 0) continue;
    q.edges.push_back(static_cast<EdgeId>(e));
    q.cap_v.push_back(b.V[e]);
    q.max_v = std::max(q.max_v, b.V[e]);
  }
  return q;
}

// min(K, floor(V / s)) for the size s = N / 2^t (N >= 1, t <= 61, V < 2^62): V * 2^t < 2^123; a
// quotient that may exceed K is first compared with K * N (below 2^92 when N < 2^64; otherwise
// the quotient is below 2^59).
Amount slots(Amount v, U128 N, int t, Amount K) noexcept {
  const U128 num = shl(u128(static_cast<std::uint64_t>(v)), t);
  if (N.hi == 0) {
    if (!(num < mul_u64(static_cast<std::uint64_t>(K), N.lo))) return K;
    return static_cast<Amount>(div_u128(num, N.lo));
  }
  return std::min<Amount>(K, static_cast<Amount>(div_u128(num, N)));
}

// One count k: the quota flows of its sizes.
struct Sizing {
  const Base& b;
  const Quota& q;
  std::int64_t k;
  Amount K;
  std::vector<Amount> caps {};
  std::vector<Amount> head_caps {};

  Sizing(const Base& base, const Quota& quota, std::int64_t count)
      : b(base), q(quota), k(count), K(static_cast<Amount>(quota.live.size()) * count), caps(quota.edges.size()),
        head_caps(quota.live.size(), count) {}

  // Do K unit paths fit at the size N / 2^t?
  bool feasible(U128 N, int t) {
    for (std::size_t i = 0; i < q.edges.size(); ++i) caps[i] = slots(q.cap_v[i], N, t, K);
    return max_flow_subset(b.g, q.edges, caps, q.live, b.dests, K, head_caps) == K;
  }

  // U: the largest integer size that fits (0 when none), at most `limit`.
  Amount common(Amount limit) {
    const auto src = b.g.edge_sources();
    const auto dst = b.g.edge_destinations();
    std::vector<Amount> out(q.live.size(), 0);
    Amount in = 0;
    std::vector<std::uint8_t> is_dst(static_cast<std::size_t>(b.g.num_nodes()), 0);
    for (NodeId d : b.dests) is_dst[static_cast<std::size_t>(d)] = 1;
    for (std::size_t i = 0; i < q.edges.size(); ++i) {
      const auto e = static_cast<std::size_t>(q.edges[i]);
      const auto h = std::lower_bound(q.live.begin(), q.live.end(), src[e]);
      if (h != q.live.end() && *h == src[e]) out[static_cast<std::size_t>(h - q.live.begin())] += q.cap_v[i];
      if (is_dst[static_cast<std::size_t>(dst[e])] != 0) in += q.cap_v[i];
    }
    // Every size that fits obeys these: k tunnels of size s leave each live head on edges with
    // floor(V / s) >= 1 summing to at least k, and K enter the destinations.
    Amount hi = std::min(limit, q.max_v);
    for (Amount x : out) hi = std::min(hi, x / k);
    hi = std::min(hi, in / K);
    if (hi < 1 || !feasible(u128(1), 0)) return 0;
    if (feasible(u128(static_cast<std::uint64_t>(hi)), 0)) return hi;
    Amount lo = 1;  // fits; hi does not
    while (hi - lo > 1) {
      const Amount mid = lo + (hi - lo) / 2;
      if (feasible(u128(static_cast<std::uint64_t>(mid)), 0)) {
        lo = mid;
      } else {
        hi = mid;
      }
    }
    return lo;
  }

  // s* as a reduced fraction (numerator, denominator), given U (live heads exist).
  std::pair<Amount, Amount> fluid(Amount U) {
    const int w = std::bit_width(static_cast<std::uint64_t>(K));  // K < 2^w
    const int t = 2 * w;                                          // 2^-t < 1 / K^2
    const U128 base = shl(u128(static_cast<std::uint64_t>(U)), t);
    // Fractional numerators over 2^t: base + lo fits, base + hi does not.
    U128 lo = U >= 1 ? U128 {} : shl(u128(1), t - w);
    U128 hi = shl(u128(1), t);
    while (sub(hi, lo) > u128(1)) {
      const U128 mid = add(lo, half(sub(hi, lo)));
      if (feasible(add(base, mid), t)) {
        lo = mid;
      } else {
        hi = mid;
      }
    }
    // The one breakpoint V[e] / j in [base + lo, base + hi) / 2^t.
    Amount num = 0;
    Amount den = 0;
    for (std::size_t i = 0; i < q.edges.size(); ++i) {
      const Amount v = q.cap_v[i];
      const Amount jmax = slots(v, add(base, lo), t, K);
      const Amount jmin = slots(v, add(base, hi), t, K) + 1;
      if (jmin > jmax) continue;
      if (jmin != jmax) internal("two breakpoints closer than 1 / K^2");
      if (den == 0) {
        num = v;
        den = jmax;
      } else if (!(mul_u64(static_cast<std::uint64_t>(num), static_cast<std::uint64_t>(jmax)) ==
                   mul_u64(static_cast<std::uint64_t>(v), static_cast<std::uint64_t>(den)))) {
        internal("two fluid sizes in one interval");
      }
    }
    if (den == 0) internal("no fluid breakpoint");
    const Amount d = std::gcd(num, den);
    return {num / d, den / d};
  }
};

// ---------------------------------------------------------------- the canonical augmentation

// A (J, Phi) price, compared lexicographically.
struct Pair {
  I128 j {};
  I128 p {};
  friend bool operator==(const Pair&, const Pair&) = default;
  friend bool operator<(const Pair& a, const Pair& b) noexcept { return a.j != b.j ? a.j < b.j : a.p < b.p; }
  friend Pair operator+(const Pair& a, const Pair& b) noexcept { return {a.j + b.j, a.p + b.p}; }
  friend Pair operator-(const Pair& a, const Pair& b) noexcept { return {a.j - b.j, a.p - b.p}; }
};

// The headroom price of slot j (j u <= V <= B): floor(2^32 (B - V + j u) / B), at most 2^32.
Amount slot_price(Amount B, Amount V, Amount j, Amount u) noexcept {
  const auto occupied = static_cast<std::uint64_t>(B - V + j * u);
  return static_cast<Amount>(div_u128(mul_u64(occupied, std::uint64_t {1} << 32), static_cast<std::uint64_t>(B)));
}

// The low 128 bits of v q (exact when the product is below 2^128).
U128 mul_lo(U128 v, std::uint64_t q) noexcept {
  const U128 low = mul_u64(v.lo, q);
  return U128 {low.hi + v.hi * q, low.lo};
}

// n (n - 1) / 2 for n >= 1.
U128 pairs_below(std::uint64_t n) noexcept { return half(mul_u64(n, n - 1)); }

// The sum of slot_price(B, V, j, u) over j = 1..x (x u <= V <= B, u >= 1) in closed form: with
// a = 2^32 u and b = 2^32 (B - V + u), the floor sum of (a i + b) / B over i = 0..x-1, by Euclid's
// reduction (O(log B) steps, so a link's x never needs a loop). a / B and b / B are at most 2^32
// (u <= V <= B and B - V + u <= B); after them every value is below 2^62 except a n + b, below
// 2^124, whose quotient by m is below n + 1. Every term is part of the sum, at most 2^32 x < 2^94.
U128 slot_price_sum(Amount B, Amount V, Amount u, Amount x) noexcept {
  U128 sum {};
  if (x <= 0) return sum;
  auto n = static_cast<std::uint64_t>(x);
  auto m = static_cast<std::uint64_t>(B);
  const auto [qa, a0] = divmod(shl(u128(static_cast<std::uint64_t>(u)), 32), m);
  const auto [qb, b0] = divmod(shl(u128(static_cast<std::uint64_t>(B - V + u)), 32), m);
  sum = add(mul_lo(pairs_below(n), qa), mul_u64(n, qb));
  std::uint64_t a = a0;
  std::uint64_t c = b0;
  while (true) {
    if (a >= m) {
      sum = add(sum, mul_lo(pairs_below(n), a / m));
      a %= m;
    }
    if (c >= m) {
      sum = add(sum, mul_u64(n, c / m));
      c %= m;
    }
    const U128 y = add(mul_u64(a, n), u128(c));
    if (y < u128(m)) return sum;
    std::tie(n, c) = divmod(y, m);
    std::swap(m, a);
  }
}

// An arc's key: (type, id, direction), compared lexicographically.
struct Key {
  int type {0};
  std::int64_t id {0};
  int dir {0};
  friend auto operator<=>(const Key&, const Key&) = default;
};

// The count vector of the canonical augmentation (aligned with q.edges) for k tunnels per live
// head at size u.
std::vector<Amount> canonical(const Base& b, const Quota& q, std::int64_t k, Amount u, bool headroom) {
  const auto& g = b.g;
  const Amount K = static_cast<Amount>(q.live.size()) * k;
  const auto src = g.edge_sources();
  const auto dst = g.edge_destinations();
  const auto cost = g.edge_costs();
  const auto cap = g.edge_capacity_units();
  // The links with a slot at u, their local endpoints; root = L, sink = L + 1.
  std::vector<std::int32_t> local(static_cast<std::size_t>(g.num_nodes()), -1);
  std::vector<NodeId> node_of;
  auto add_node = [&](NodeId v) {
    auto& s = local[static_cast<std::size_t>(v)];
    if (s < 0) {
      s = static_cast<std::int32_t>(node_of.size());
      node_of.push_back(v);
    }
    return s;
  };
  struct Link {
    std::size_t q {0};  // position in q.edges
    EdgeId e {0};
    std::int32_t tail {0};
    std::int32_t head {0};
    Amount slots {0};
    Amount x {0};
    Pair fwd {};  // the price of slot x + 1
    Pair rev {};  // minus the price of slot x
  };
  std::vector<Link> links;
  for (NodeId h : q.live) add_node(h);
  for (std::size_t i = 0; i < q.edges.size(); ++i) {
    const Amount s = std::min(K, q.cap_v[i] / u);
    if (s == 0) continue;
    const auto e = static_cast<std::size_t>(q.edges[i]);
    links.push_back(Link {i, q.edges[i], add_node(src[e]), add_node(dst[e]), s});
  }
  const auto L = static_cast<std::int32_t>(node_of.size());
  const std::int32_t root = L;
  const std::int32_t sink = L + 1;
  const auto N = static_cast<std::size_t>(L + 2);
  std::vector<std::uint8_t> is_dst(N, 0);
  for (NodeId d : b.dests) {
    const auto s = local[static_cast<std::size_t>(d)];
    if (s >= 0) is_dst[static_cast<std::size_t>(s)] = 1;
  }
  auto price = [&](const Link& l, Amount j) -> Pair {
    const auto ei = static_cast<std::size_t>(l.e);
    const Amount B = cap[ei];
    return {I128(cost[ei]), I128(headroom && B > 0 ? slot_price(B, b.V[ei], j, u) : 0)};
  };
  auto reprice = [&](Link& l) {
    if (l.x < l.slots) l.fwd = price(l, l.x + 1);
    if (l.x > 0) l.rev = Pair {} - price(l, l.x);
  };
  for (auto& l : links) reprice(l);
  // Per local node its forward links (by tail) and reverse links (by head), CSR.
  std::vector<std::int32_t> out_first(N + 1, 0);
  std::vector<std::int32_t> in_first(N + 1, 0);
  for (const auto& l : links) {
    ++out_first[static_cast<std::size_t>(l.tail) + 1];
    ++in_first[static_cast<std::size_t>(l.head) + 1];
  }
  for (std::size_t v = 0; v < N; ++v) {
    out_first[v + 1] += out_first[v];
    in_first[v + 1] += in_first[v];
  }
  std::vector<std::int32_t> out_links(links.size());
  std::vector<std::int32_t> in_links(links.size());
  {
    std::vector<std::int32_t> o(out_first.begin(), out_first.end() - 1);
    std::vector<std::int32_t> n(in_first.begin(), in_first.end() - 1);
    for (std::size_t i = 0; i < links.size(); ++i) {
      out_links[static_cast<std::size_t>(o[static_cast<std::size_t>(links[i].tail)]++)] = static_cast<std::int32_t>(i);
      in_links[static_cast<std::size_t>(n[static_cast<std::size_t>(links[i].head)]++)] = static_cast<std::int32_t>(i);
    }
  }
  std::vector<std::int64_t> root_flow(q.live.size(), 0);
  std::vector<Pair> pi(N);
  std::vector<Pair> dist(N);
  std::vector<std::uint8_t> reached(N, 0);
  std::vector<std::uint8_t> settled(N, 0);
  std::vector<std::int32_t> order;  // the settled vertices
  // A residual arc out of a vertex: its head, reduced price, key and what an augmentation does
  // (kind 0 a link forward, 1 a link reverse, 2 a root arc, 3 a sink arc; index).
  struct Arc {
    std::int32_t to {0};
    Pair rp {};
    Key key {};
    int kind {0};
    std::int32_t index {0};
  };
  auto arcs_of = [&](std::int32_t v, auto&& visit) {
    if (v == root) {
      for (std::size_t h = 0; h < q.live.size(); ++h) {
        if (root_flow[h] >= k) continue;
        const auto to = static_cast<std::int32_t>(h);  // live heads are the first local nodes
        visit(Arc {to, pi[static_cast<std::size_t>(root)] - pi[static_cast<std::size_t>(to)], Key {1, q.live[h], 0}, 2,
                   static_cast<std::int32_t>(h)});
      }
      return;
    }
    const auto vi = static_cast<std::size_t>(v);
    for (auto s = out_first[vi]; s < out_first[vi + 1]; ++s) {
      const auto i = out_links[static_cast<std::size_t>(s)];
      const auto& l = links[static_cast<std::size_t>(i)];
      if (l.x >= l.slots) continue;
      visit(Arc {l.head, l.fwd + pi[vi] - pi[static_cast<std::size_t>(l.head)], Key {0, l.e, 0}, 0, i});
    }
    for (auto s = in_first[vi]; s < in_first[vi + 1]; ++s) {
      const auto i = in_links[static_cast<std::size_t>(s)];
      const auto& l = links[static_cast<std::size_t>(i)];
      if (l.x <= 0) continue;
      visit(Arc {l.tail, l.rev + pi[vi] - pi[static_cast<std::size_t>(l.tail)], Key {0, l.e, 1}, 1, i});
    }
    if (is_dst[vi] != 0) {
      visit(Arc {sink, pi[vi] - pi[static_cast<std::size_t>(sink)], Key {2, node_of[vi], 0}, 3, v});
    }
  };
  // The residual arcs into a vertex: visit(tail, reduced price) (none into the root).
  std::vector<std::int32_t> dest_local;
  for (std::size_t v = 0; v < N; ++v) {
    if (is_dst[v] != 0) dest_local.push_back(static_cast<std::int32_t>(v));
  }
  auto in_arcs_of = [&](std::int32_t w, auto&& visit) {
    const auto wi = static_cast<std::size_t>(w);
    if (w == sink) {
      for (std::int32_t d : dest_local) visit(d, pi[static_cast<std::size_t>(d)] - pi[wi]);
      return;
    }
    if (w == root) return;
    for (auto s = in_first[wi]; s < in_first[wi + 1]; ++s) {
      const auto& l = links[static_cast<std::size_t>(in_links[static_cast<std::size_t>(s)])];
      if (l.x < l.slots) visit(l.tail, l.fwd + pi[static_cast<std::size_t>(l.tail)] - pi[wi]);
    }
    for (auto s = out_first[wi]; s < out_first[wi + 1]; ++s) {
      const auto& l = links[static_cast<std::size_t>(out_links[static_cast<std::size_t>(s)])];
      if (l.x > 0) visit(l.head, l.rev + pi[static_cast<std::size_t>(l.head)] - pi[wi]);
    }
    if (wi < q.live.size() && root_flow[wi] < k) visit(root, pi[static_cast<std::size_t>(root)] - pi[wi]);
  };
  using Entry = std::pair<Pair, std::int32_t>;
  auto later = [](const Entry& a, const Entry& b) { return b.first < a.first; };
  std::priority_queue<Entry, std::vector<Entry>, decltype(later)> heap(later);
  std::vector<std::int32_t> hop(N, -1);
  std::vector<std::int32_t> frontier;
  for (Amount unit = 0; unit < K; ++unit) {
    std::fill(reached.begin(), reached.end(), 0);
    std::fill(settled.begin(), settled.end(), 0);
    order.clear();
    dist[static_cast<std::size_t>(root)] = Pair {};
    reached[static_cast<std::size_t>(root)] = 1;
    heap.push({Pair {}, root});
    while (!heap.empty()) {
      const auto [d, v] = heap.top();
      const auto vi = static_cast<std::size_t>(v);
      if (settled[vi] != 0 || dist[vi] < d) {
        heap.pop();
        continue;
      }
      // Every vertex at the sink's distance or nearer is settled; the rest is not on a shortest path.
      if (settled[static_cast<std::size_t>(sink)] != 0 && dist[static_cast<std::size_t>(sink)] < d) break;
      heap.pop();
      settled[vi] = 1;
      order.push_back(v);
      if (v == sink) continue;
      arcs_of(v, [&](const Arc& a) {
        if (a.rp < Pair {}) internal("a negative reduced price");
        const auto w = static_cast<std::size_t>(a.to);
        const Pair nd = d + a.rp;
        if (settled[w] != 0) return;
        if (reached[w] == 0 || nd < dist[w]) {
          reached[w] = 1;
          dist[w] = nd;
          heap.push({nd, a.to});
        }
      });
    }
    while (!heap.empty()) heap.pop();
    if (settled[static_cast<std::size_t>(sink)] == 0) internal("no augmenting path below K units");
    const Pair dt = dist[static_cast<std::size_t>(sink)];
    // The hop counts back from the sink over the tight arcs (settled at both ends, the reduced
    // price closing the distances exactly): a breadth-first search over the arcs into each vertex,
    // so only vertices on shortest paths are reached.
    for (std::int32_t v : order) hop[static_cast<std::size_t>(v)] = -1;
    hop[static_cast<std::size_t>(sink)] = 0;
    frontier.assign(1, sink);
    for (std::size_t f = 0; f < frontier.size(); ++f) {
      const std::int32_t w = frontier[f];
      const auto wi = static_cast<std::size_t>(w);
      in_arcs_of(w, [&](std::int32_t a, const Pair& rp) {
        const auto ai = static_cast<std::size_t>(a);
        if (settled[ai] == 0 || hop[ai] >= 0 || !(dist[ai] + rp == dist[wi])) return;
        hop[ai] = hop[wi] + 1;
        frontier.push_back(a);
      });
    }
    if (hop[static_cast<std::size_t>(root)] < 0) internal("the sink is not reached over tight arcs");
    // The walk: at each vertex the smallest-key tight arc one hop nearer the sink; one unit along it.
    std::int32_t v = root;
    while (v != sink) {
      const auto vi = static_cast<std::size_t>(v);
      std::optional<Arc> best;
      arcs_of(v, [&](const Arc& arc) {
        const auto w = static_cast<std::size_t>(arc.to);
        if (settled[w] == 0 || hop[w] != hop[vi] - 1 || !(dist[vi] + arc.rp == dist[w])) return;
        if (!best || arc.key < best->key) best = arc;
      });
      if (!best) internal("the walk found no tight arc");
      if (best->kind == 0) {
        ++links[static_cast<std::size_t>(best->index)].x;
        reprice(links[static_cast<std::size_t>(best->index)]);
      } else if (best->kind == 1) {
        --links[static_cast<std::size_t>(best->index)].x;
        reprice(links[static_cast<std::size_t>(best->index)]);
      } else if (best->kind == 2) {
        ++root_flow[static_cast<std::size_t>(best->index)];
      }
      v = best->to;
    }
    for (std::size_t w = 0; w < N; ++w) pi[w] = pi[w] + (settled[w] != 0 ? dist[w] : dt);
  }
  for (std::int64_t f : root_flow) {
    if (f != k) internal("a live head without its k tunnels");
  }
  std::vector<Amount> x(q.edges.size(), 0);
  for (const auto& l : links) x[l.q] = l.x;
  return x;
}

// ---------------------------------------------------------------- programs and budgets

// A layout at one count and size: the tunnels in round order and the objectives.
struct Laid {
  std::vector<TunnelSetupResult> tunnels {};
  Counter128 route_cost {};
  Counter128 headroom_penalty {};
};

Laid lay(const Base& b, const Quota& q, std::int64_t k, Amount u, bool headroom) {
  const auto& g = b.g;
  const auto x = canonical(b, q, k, u, headroom);
  std::vector<Amount> dense(static_cast<std::size_t>(g.num_edges()), 0);
  for (std::size_t i = 0; i < q.edges.size(); ++i) dense[static_cast<std::size_t>(q.edges[i])] = x[i];
  const auto objectives = t2_objectives(g, q.edges, x, b.V, u, headroom);
  const auto cost = g.edge_costs();
  // Unit paths per live head (the controllers' cut, each path expanded), sorted by EdgeId sequence.
  auto dec = decompose(g, {}, std::move(dense), q.live, b.dests);
  std::vector<std::vector<std::vector<EdgeId>>> per_head(q.live.size());
  for (auto& p : dec.paths) {
    auto& rows = per_head[static_cast<std::size_t>(p.roster_index)];
    for (Amount n = 0; n < p.weight; ++n) rows.push_back(p.edges);
  }
  for (auto& rows : per_head) {
    if (static_cast<std::int64_t>(rows.size()) != k) internal("a live head's paths are not k");
    std::sort(rows.begin(), rows.end());
  }
  Laid out;
  out.route_cost = objectives.route_cost;
  out.headroom_penalty = objectives.headroom_penalty;
  out.tunnels.reserve(q.live.size() * static_cast<std::size_t>(k));
  for (std::int64_t j = 0; j < k; ++j) {
    for (std::size_t h = 0; h < q.live.size(); ++h) {
      TunnelSetupResult t;
      t.roster_index = static_cast<std::int32_t>(q.live_pos[h]);
      t.local_index = static_cast<std::int32_t>(j);
      t.planned = u;
      t.edges = per_head[h][static_cast<std::size_t>(j)];
      for (EdgeId e : t.edges) t.path_cost += cost[static_cast<std::size_t>(e)];
      out.tunnels.push_back(std::move(t));
    }
  }
  return out;
}

// G(z): the live heads' part of a round-robin total z over the H roster positions.
U128 live_part(const Quota& q, std::int64_t H, U128 z) noexcept {
  const auto [per, rest] = divmod(z, static_cast<std::uint64_t>(H));
  U128 g = mul_u64(per, q.live.size());
  for (std::int64_t i : q.live_pos) {
    if (static_cast<std::uint64_t>(i) < rest) g = add(g, u128(1));
  }
  return g;
}

// F'(k): the largest total z whose live shares fit on the tunnels of the program at U (each tunnel
// at position j H + i receives RR(z, M, j H + i); dead positions ignored), in [M U, M (U + 1)).
U128 live_budget(const Base& b, const std::vector<TunnelSetupResult>& tunnels, std::int64_t H, std::int64_t k, Amount U) {
  const auto M = static_cast<std::uint64_t>(H * k);
  std::vector<Amount> load(static_cast<std::size_t>(b.g.num_edges()), 0);
  std::vector<EdgeId> touched;
  auto fits = [&](U128 z) {
    const auto [per, rest] = divmod(z, M);
    for (const auto& t : tunnels) {
      const auto pos = static_cast<std::uint64_t>(t.local_index * H + t.roster_index);  // below H k <= 2^30
      const auto amount = static_cast<Amount>(per) + (pos < rest ? 1 : 0);
      for (EdgeId e : t.edges) {
        auto& x = load[static_cast<std::size_t>(e)];
        if (x == 0) touched.push_back(e);
        x += amount;
      }
    }
    bool ok = true;
    for (EdgeId e : touched) {
      ok = ok && load[static_cast<std::size_t>(e)] <= b.V[static_cast<std::size_t>(e)];
      load[static_cast<std::size_t>(e)] = 0;
    }
    touched.clear();
    return ok;
  };
  U128 lo = mul_u64(M, static_cast<std::uint64_t>(U));  // fits: every tunnel gets U
  U128 hi = add(lo, u128(M));                           // does not: every tunnel would get U + 1
  while (sub(hi, lo) > u128(1)) {
    const U128 mid = add(lo, half(sub(hi, lo)));
    if (fits(mid)) {
      lo = mid;
    } else {
      hi = mid;
    }
  }
  return lo;
}

// The working size of a finite design D at count k: the largest live design share per tunnel,
// at least 1, at most U.
Amount working_size(const Quota& q, std::int64_t H, std::int64_t k, Amount U, std::optional<Amount> design) {
  if (!design.has_value()) return U;
  Amount need = 1;
  for (std::int64_t i : q.live_pos) need = std::max(need, (rr(*design, H, i) + k - 1) / k);
  return std::min(U, need);
}

} // namespace

T2Layout t2_layout(const Graph& g, const Eligibility* elig, std::span<const NodeId> roster,
                   std::span<const NodeId> dests, const TunnelCount& tunnels, bool headroom, std::optional<Amount> design,
                   std::span<const Amount> V) {
  const Base b = make_base(g, elig, roster, dests, V);
  const Quota q = make_quota(b);
  const auto H = static_cast<std::int64_t>(roster.size());
  T2Layout out;
  out.live = q.live;
  out.count = tunnels.up_to ? 1 : tunnels.per_head;
  if (q.live.empty()) return out;
  std::int64_t k = tunnels.per_head;
  Amount U = 0;
  std::optional<Laid> cached;  // the chosen count's program at U(k), when "up to" built it
  if (!tunnels.up_to) {
    U = Sizing(b, q, k).common(q.max_v);
  } else {
    // Scores are live parts of totals (below 2^63); `target` is G(D) for a finite design.
    const std::optional<U128> D = design.has_value() ? std::optional<U128>(u128(static_cast<std::uint64_t>(*design)))
                                                     : std::nullopt;
    const std::optional<U128> target = D.has_value() ? std::optional<U128>(live_part(q, H, *D)) : std::nullopt;
    auto capped = [&](U128 z) { return D.has_value() && *D < z ? *D : z; };
    k = 1;
    U128 best {};
    Amount limit = q.max_v;
    for (std::int64_t c = 1; c <= tunnels.per_head; ++c) {
      Sizing sz(b, q, c);
      const Amount Uc = sz.common(limit);
      if (Uc == 0) break;  // every later count fits no unit either: score 0
      limit = Uc;
      if (c == 1) U = Uc;  // the incumbent's size until a better count
      const auto M = static_cast<std::uint64_t>(H * c);
      const U128 lower = live_part(q, H, capped(mul_u64(M, static_cast<std::uint64_t>(Uc))));
      const U128 upper = live_part(q, H, capped(sub(mul_u64(M, static_cast<std::uint64_t>(Uc) + 1), u128(1))));
      if (target.has_value() && lower == *target) {
        k = c;
        U = Uc;
        cached.reset();
        break;
      }
      if (!(best < upper)) continue;
      Laid at_u = lay(b, q, c, Uc, headroom);
      const U128 score = live_part(q, H, capped(live_budget(b, at_u.tunnels, H, c, Uc)));
      if (best < score) {
        best = score;
        k = c;
        U = Uc;
        cached = std::move(at_u);
      }
      if (target.has_value() && score == *target) break;
    }
  }
  out.count = k;
  out.common = U;
  if (U == 0) {
    // No whole unit per tunnel fits (live heads exist): no tunnels; the fluid size names the cause.
    std::tie(out.fluid_num, out.fluid_den) = Sizing(b, q, k).fluid(0);
    return out;
  }
  std::tie(out.fluid_num, out.fluid_den) = Sizing(b, q, k).fluid(U);
  out.planned = working_size(q, H, k, U, design);
  Laid laid = cached.has_value() && out.planned == U ? std::move(*cached) : lay(b, q, k, out.planned, headroom);
  out.tunnels = std::move(laid.tunnels);
  out.route_cost = laid.route_cost;
  out.headroom_penalty = laid.headroom_penalty;
  return out;
}

T2Objectives t2_objectives(const Graph& g, std::span<const EdgeId> edges, std::span<const Amount> x,
                           std::span<const Amount> V, Amount u, bool headroom) {
  const auto cost = g.edge_costs();
  const auto cap = g.edge_capacity_units();
  U128 J {};
  U128 Phi {};
  for (std::size_t i = 0; i < edges.size(); ++i) {
    if (x[i] == 0) continue;
    const auto e = static_cast<std::size_t>(edges[i]);
    J = add(J, mul_u64(static_cast<std::uint64_t>(cost[e]), static_cast<std::uint64_t>(x[i])));
    if (!headroom || cap[e] == 0) continue;
    if (V[e] > cap[e] || x[i] > V[e] / u) internal("a link's slots exceed its free capacity");
    Phi = add(Phi, slot_price_sum(cap[e], V[e], u, x[i]));
  }
  return {Counter128 {J.hi, J.lo}, Counter128 {Phi.hi, Phi.lo}};
}

} // namespace netgraph::core::detail
