// PredDAG and Path: from_csr's complete validation in the documented order, the
// canonical predecessor order, structural equality, Path validation in its order, cost and
// to_dag().
#include <gtest/gtest.h>

#include <algorithm>
#include <random>
#include <vector>

#include "netgraph/core/detail/pred_dag_builder.hpp"
#include "netgraph/core/pred_dag.hpp"
#include "test_utils.hpp"

using namespace netgraph::core;
using namespace netgraph::core::test;

namespace {

// 0 -> 1 (two parallel members), 1 -> 2, 0 -> 2 and 2 -> 0; ids follow (src, cost, dst, input).
Graph diamondish() {
  return make_graph(3, {{0, 1, 5, 1}, {0, 1, 5, 1}, {1, 2, 5, 1}, {0, 2, 5, 3}, {2, 0, 5, 1}});
}

using V = std::vector<std::int32_t>;

} // namespace

TEST(PredDAG, FromCsrRoundTripAndAccessors) {
  auto g = diamondish();
  // ids: 0: 0->1, 1: 0->1, 2: 0->2 (cost 3), 3: 1->2, 4: 2->0
  const V off = {0, 0, 2, 4};
  const V par = {0, 0, 0, 1};
  const V via = {0, 1, 2, 3};
  auto dag = PredDAG::from_csr(g, off, par, via);
  EXPECT_EQ(dag.graph_uid(), g.uid());
  EXPECT_EQ(dag.parent_offsets(), off);
  EXPECT_EQ(V(dag.parents().begin(), dag.parents().end()), par);
  EXPECT_EQ(V(dag.via_edges().begin(), dag.via_edges().end()), via);
  EXPECT_EQ(ascending_edges(dag), (std::vector<EdgeId>{0, 1, 2, 3}));
  // Structural equality: the builder's value from the same edge set is equal.
  const std::vector<EdgeId> edges = {3, 1, 0, 2};
  EXPECT_TRUE(detail::PredDAGBuilder::of_edges(g, edges) == dag);
  auto empty = PredDAG::from_csr(g, V {0, 0, 0, 0}, V {}, V {});
  EXPECT_FALSE(empty == dag);
  EXPECT_TRUE(ascending_edges(empty).empty());
  // Equality includes the graph: the same CSR on a twin graph (the same arrays, another uid) is
  // another value, from from_csr and from the builder alike.
  const auto twin = diamondish();
  ASSERT_NE(twin.uid(), g.uid());
  const auto there = PredDAG::from_csr(twin, off, par, via);
  EXPECT_EQ(there.parent_offsets(), off);
  EXPECT_EQ(V(there.via_edges().begin(), there.via_edges().end()), via);
  EXPECT_FALSE(there == dag);
  EXPECT_FALSE(detail::PredDAGBuilder::of_edges(twin, edges) == dag);
  EXPECT_TRUE(detail::PredDAGBuilder::of_edges(twin, edges) == there);
}

// The builder's entries for any edge subset in any input order: children ascending, via
// ascending within a child, the parent each via's tail, and the same value as from_csr. The
// builder groups a subset of fewer than N / 8 edges by sorting (child, via) pairs and a larger one
// by a counting sort over the nodes: small graphs take mostly the counting sort, and graphs of 512
// to 767 nodes with four to seven children of five to nine in-edges each (at most 63 edges, often
// more than a sort keeps in insertion order) take only the pair sort. Each population is guarded by
// its DAGs with a child of several entries, whose order is what differs between the two.
TEST(PredDAG, BuilderGroupsAnyEdgeSubsetCanonically) {
  std::mt19937 rng(5);
  // Checks the builder's value for one subset; true when a child has several entries.
  const auto check = [&](const Graph& g, std::vector<EdgeId> subset, int trial) {
    const int n = g.num_nodes();
    std::shuffle(subset.begin(), subset.end(), rng);
    const auto dag = detail::PredDAGBuilder::of_edges(g, subset);
    V off(static_cast<std::size_t>(n) + 1, 0);
    std::vector<std::pair<NodeId, EdgeId>> entries;
    for (EdgeId e : subset) entries.emplace_back(g.edge_dst()[static_cast<std::size_t>(e)], e);
    std::sort(entries.begin(), entries.end());
    V par;
    V via;
    for (auto [child, e] : entries) {
      ++off[static_cast<std::size_t>(child) + 1];
      par.push_back(g.edge_src()[static_cast<std::size_t>(e)]);
      via.push_back(e);
    }
    for (std::size_t i = 1; i < off.size(); ++i) off[i] += off[i - 1];
    EXPECT_EQ(dag.parent_offsets(), off) << trial;
    EXPECT_EQ(V(dag.parents().begin(), dag.parents().end()), par) << trial;
    EXPECT_EQ(V(dag.via_edges().begin(), dag.via_edges().end()), via) << trial;
    EXPECT_TRUE(PredDAG::from_csr(g, off, par, via) == dag) << trial;
    std::sort(subset.begin(), subset.end());
    EXPECT_EQ(ascending_edges(dag), subset) << trial;
    const auto o = dag.parent_offsets();
    bool several = false;
    for (std::size_t i = 1; i < o.size(); ++i) several = several || o[i] - o[i - 1] >= 2;
    return several;
  };
  int grouped = 0;  // DAGs with a child of several entries (their order is checked)
  for (int trial = 0; trial < 200; ++trial) {
    const int n = 2 + static_cast<int>(rng() % 9);
    std::vector<InEdge> es;
    for (int k = 0, m = 1 + static_cast<int>(rng() % 30); k < m; ++k) {  // u < v: every subset is a DAG
      const auto u = static_cast<NodeId>(rng() % static_cast<unsigned>(n - 1));
      const auto v = static_cast<NodeId>(u + 1 + static_cast<NodeId>(rng() % static_cast<unsigned>(n - 1 - u)));
      es.push_back({u, v, 1, static_cast<Cost>(1 + rng() % 3)});
    }
    const auto g = make_graph(n, es);
    std::vector<EdgeId> subset;
    for (EdgeId e = 0; e < g.num_edges(); ++e) {
      if (rng() % 2 == 0) subset.push_back(e);
    }
    grouped += check(g, subset, trial) ? 1 : 0;
  }
  EXPECT_GT(grouped, 50);
  int sparse_grouped = 0;  // pair-sorted DAGs with a child of several entries
  for (int trial = 0; trial < 200; ++trial) {
    const int n = 512 + static_cast<int>(rng() % 256);
    std::vector<InEdge> es;
    for (int c = 0, children = 4 + static_cast<int>(rng() % 4); c < children; ++c) {
      const auto v = static_cast<NodeId>(1 + rng() % static_cast<unsigned>(n - 1));
      for (int k = 0, d = 5 + static_cast<int>(rng() % 5); k < d; ++k) {  // u < v; parallel members possible
        es.push_back({static_cast<NodeId>(rng() % static_cast<unsigned>(v)), v, 1, static_cast<Cost>(1 + rng() % 3)});
      }
    }
    const auto g = make_graph(n, es);
    std::vector<EdgeId> subset;
    for (EdgeId e = 0; e < g.num_edges(); ++e) {
      if (rng() % 8 != 0) subset.push_back(e);
    }
    ASSERT_LT(subset.size() * 8, static_cast<std::size_t>(n));  // the pair sort
    sparse_grouped += check(g, subset, 1000 + trial) ? 1 : 0;
  }
  EXPECT_GT(sparse_grouped, 100);
}

TEST(PredDAG, FromCsrShapeChecksInOrder) {
  auto g = diamondish();
  // offsets of length N + 1
  EXPECT_TRUE(contains(error_of([&] { (void)PredDAG::from_csr(g, V {0, 0, 0}, V {}, V {}); }), "offsets has length 3"));
  // offsets[0] == 0 (also before monotone)
  EXPECT_TRUE(contains(error_of([&] { (void)PredDAG::from_csr(g, V {1, 0, 2, 1}, V {}, V {}); }), "offsets[0] = 1"));
  // monotone
  EXPECT_TRUE(contains(error_of([&] { (void)PredDAG::from_csr(g, V {0, 2, 1, 2}, V {0, 0}, V {0, 1}); }),
                       "not monotone at 2"));
  // offsets[N] == len(parents) == len(via)
  EXPECT_TRUE(contains(error_of([&] { (void)PredDAG::from_csr(g, V {0, 0, 2, 2}, V {0}, V {0, 1}); }),
                       "offsets[N] = 2 but parents has 1"));
  EXPECT_TRUE(contains(error_of([&] { (void)PredDAG::from_csr(g, V {0, 0, 2, 2}, V {0, 0}, V {0}); }),
                       "and via has 1"));
  // Shape before entries: a bad parent and a shape defect report the shape.
  EXPECT_TRUE(contains(error_of([&] { (void)PredDAG::from_csr(g, V {0, 0, 2, 3}, V {9, 0}, V {0, 1}); }),
                       "offsets[N] = 3"));
  // Every shape defect at once, then each reported one removed: each step is reported before
  // all the later ones (and before the bad entry 9, 9).
  EXPECT_EQ(error_of([&] { (void)PredDAG::from_csr(g, V {1, 0, 0}, V {9}, V {9}); }),
            "PredDAG.from_csr: offsets has length 3, expected N + 1 = 4");
  EXPECT_EQ(error_of([&] { (void)PredDAG::from_csr(g, V {1, 0, 2, 5}, V {9}, V {9}); }),
            "PredDAG.from_csr: offsets[0] = 1, expected 0");
  EXPECT_EQ(error_of([&] { (void)PredDAG::from_csr(g, V {0, 2, 1, 5}, V {9}, V {9}); }),
            "PredDAG.from_csr: offsets are not monotone at 2 (2 > 1)");
  EXPECT_EQ(error_of([&] { (void)PredDAG::from_csr(g, V {0, 0, 2, 5}, V {9}, V {9}); }),
            "PredDAG.from_csr: offsets[N] = 5 but parents has 1 and via has 1 entries");
}

TEST(PredDAG, FromCsrEntryChecksInOrder) {
  auto g = diamondish();
  const V off = {0, 0, 2, 2};
  // parents in range, then via in range (entry by entry)
  EXPECT_EQ(error_of([&] { (void)PredDAG::from_csr(g, off, V {3, 0}, V {0, 1}); }),
            "PredDAG.from_csr: parents[0] = 3 out of node range [0, 3)");
  EXPECT_TRUE(contains(error_of([&] { (void)PredDAG::from_csr(g, off, V {0, 0}, V {0, 5}); }), "via[1] = 5"));
  EXPECT_EQ(error_of([&] { (void)PredDAG::from_csr(g, off, V {-1, 0}, V {-1, 1}); }),
            "PredDAG.from_csr: parents[0] = -1 out of node range [0, 3)");
  // via endpoints: wrong parent, wrong child
  EXPECT_TRUE(contains(error_of([&] { (void)PredDAG::from_csr(g, off, V {1, 0}, V {0, 1}); }), "joins (0, 1), not (1, 1)"));
  EXPECT_TRUE(contains(error_of([&] { (void)PredDAG::from_csr(g, V {0, 0, 0, 1}, V {0}, V {0}); }),
                       "joins (0, 1), not (0, 2)"));
  // duplicate (child, via) before the order check
  EXPECT_TRUE(contains(error_of([&] { (void)PredDAG::from_csr(g, off, V {0, 0}, V {1, 1}); }), "duplicate entry"));
  // a cycle (0 -> 1 -> 2 -> 0) before the order check: node 1's two entries descend (via 1, then 0)
  EXPECT_TRUE(contains(error_of([&] { (void)PredDAG::from_csr(g, V {0, 1, 3, 4}, V {2, 0, 0, 1}, V {4, 1, 0, 3}); }),
                       "cycle"));
  // canonical order: ascending via within a node
  EXPECT_TRUE(contains(error_of([&] { (void)PredDAG::from_csr(g, off, V {0, 0}, V {1, 0}); }),
                       "not in ascending via EdgeId"));
  // Every entry defect at once, then each reported one removed. Node 0 has (2, via 4), node 1
  // (0, via 1), (0, via 0), (0, via 0): a cycle 0 -> 1 -> 2 -> 0, a duplicate and a descending
  // pair; node 2 has a wrong endpoint (entry 4) before an out-of-range via (entry 5), so the
  // ranges of every entry come before any endpoint.
  const auto first = [&](V o, V p, V w) { return error_of([&] { (void)PredDAG::from_csr(g, o, p, w); }); };
  EXPECT_EQ(first(V {0, 1, 4, 6}, V {2, 0, 0, 0, 0, 1}, V {4, 1, 0, 0, 0, 9}),
            "PredDAG.from_csr: via[5] = 9 out of edge range [0, 5)");
  EXPECT_EQ(first(V {0, 1, 4, 6}, V {2, 0, 0, 0, 0, 1}, V {4, 1, 0, 0, 0, 3}),
            "PredDAG.from_csr: entry 4 of node 2: via edge 0 joins (0, 1), not (0, 2)");
  EXPECT_EQ(first(V {0, 1, 4, 6}, V {2, 0, 0, 0, 0, 1}, V {4, 1, 0, 0, 2, 3}),
            "PredDAG.from_csr: duplicate entry (child 1, via 0)");
  EXPECT_EQ(first(V {0, 1, 3, 5}, V {2, 0, 0, 0, 1}, V {4, 1, 0, 2, 3}), "PredDAG.from_csr: the entries contain a cycle");
  EXPECT_EQ(first(V {0, 0, 2, 4}, V {0, 0, 0, 1}, V {1, 0, 2, 3}),
            "PredDAG.from_csr: the entries of node 1 are not in ascending via EdgeId (1 then 0)");
  EXPECT_EQ(first(V {0, 0, 2, 4}, V {0, 0, 0, 1}, V {0, 1, 2, 3}), "");
}

TEST(PredDAG, SelfLoopIsACycle) {
  auto g = make_graph(2, {{0, 0, 1, 1}, {0, 1, 1, 1}});
  // id 0: 0->0, id 1: 0->1
  EXPECT_TRUE(contains(error_of([&] { (void)PredDAG::from_csr(g, V {0, 1, 1}, V {0}, V {0}); }), "cycle"));
}

TEST(PredDAG, EmptyGraph) {
  auto g = make_graph(0, {});
  auto dag = PredDAG::from_csr(g, V {0}, V {}, V {});
  EXPECT_EQ(dag.parent_offsets(), V {0});
  EXPECT_EQ(error_of([&] { (void)PredDAG::from_csr(g, V {}, V {}, V {}); }),
            "PredDAG.from_csr: offsets has length 0, expected N + 1 = 1");
}

TEST(Path, CostSrcDstAndToDag) {
  auto g = diamondish();
  // Grouped first hop (members 1 and 0, caller order kept), then 1 -> 2.
  Path p(g, {{1, 0}, {3}});
  EXPECT_EQ(p.src(), 0);
  EXPECT_EQ(p.dst(), 2);
  EXPECT_EQ(p.cost(), 2);
  EXPECT_EQ(p.hops(), (std::vector<std::vector<EdgeId>>{{1, 0}, {3}}));
  EXPECT_EQ(p.graph_uid(), g.uid());
  auto dag = p.to_dag();
  EXPECT_EQ(dag.parent_offsets(), (V {0, 0, 2, 3}));
  EXPECT_EQ(V(dag.via_edges().begin(), dag.via_edges().end()), (V {0, 1, 3}));
  EXPECT_EQ(V(dag.parents().begin(), dag.parents().end()), (V {0, 0, 1}));
  // to_dag gives a canonical value: from_csr accepts it unchanged.
  EXPECT_TRUE(PredDAG::from_csr(g, dag.parent_offsets(), dag.parents(), dag.via_edges()) == dag);
  // Cost is the minimum member cost per hop.
  auto g2 = make_graph(2, {{0, 1, 1, 4}, {0, 1, 1, 2}});
  EXPECT_EQ(Path(g2, {{0, 1}}).cost(), 2);
  EXPECT_TRUE(Path(g2, {{0, 1}}) == Path(g2, {{0, 1}}));
  EXPECT_FALSE(Path(g2, {{0, 1}}) == Path(g2, {{1, 0}}));
  // Equality includes the graph: the same hops on a twin graph (the same arrays, another uid).
  auto twin = make_graph(2, {{0, 1, 1, 4}, {0, 1, 1, 2}});
  ASSERT_NE(twin.uid(), g2.uid());
  EXPECT_EQ(Path(twin, {{0, 1}}).hops(), Path(g2, {{0, 1}}).hops());
  EXPECT_FALSE(Path(twin, {{0, 1}}) == Path(g2, {{0, 1}}));
}

TEST(Path, ValidationInOrder) {
  auto g = diamondish();
  // range before everything (an empty hop elsewhere is reported later)
  EXPECT_TRUE(contains(error_of([&] { Path(g, {{}, {7}}); }), "member 7 out of edge range"));
  EXPECT_TRUE(contains(error_of([&] { Path(g, {{-1}}); }), "member -1 out of edge range"));
  // duplicate member before an empty hop
  EXPECT_TRUE(contains(error_of([&] { Path(g, {{0, 0}, {}}); }), "duplicate member 0"));
  EXPECT_EQ(error_of([&] { Path(g, {}); }), "Path: the hop list is empty");
  EXPECT_TRUE(contains(error_of([&] { Path(g, {{0}, {}}); }), "hop 1 is empty"));
  // members of a hop share (u, v)
  EXPECT_TRUE(contains(error_of([&] { Path(g, {{0, 2}}); }), "do not share (u, v)"));
  // contiguous hops
  EXPECT_TRUE(contains(error_of([&] { Path(g, {{0}, {2}}); }), "not contiguous"));
  // no repeated node: 0 -> 2 -> 0
  EXPECT_TRUE(contains(error_of([&] { Path(g, {{2}, {4}}); }), "node 0 is repeated"));
  // Several defects at once, then each reported one removed (an empty hop list has no member, so
  // it is never reported with another defect). The later defect sits in an earlier hop where it
  // can: the order is by step, not by position.
  EXPECT_EQ(error_of([&] { Path(g, {{0, 0, 2}, {}, {7}}); }), "Path: hop 2 member 7 out of edge range [0, 5)");
  EXPECT_EQ(error_of([&] { Path(g, {{0, 0, 2}, {}, {4}}); }), "Path: duplicate member 0");
  EXPECT_EQ(error_of([&] { Path(g, {{0, 2}, {}, {4}}); }), "Path: hop 1 is empty");
  // hop 2 mixes 0 -> 2 and 0 -> 1; hops 0 and 1 do not meet (1 then 2); node 0 repeats
  EXPECT_EQ(error_of([&] { Path(g, {{1}, {4}, {2, 0}}); }),
            "Path: the members of hop 2 do not share (u, v): edge 0 joins (0, 1), edge 2 joins (0, 2)");
  EXPECT_EQ(error_of([&] { Path(g, {{1}, {4}, {2}}); }), "Path: hops 0 and 1 are not contiguous (1 then 2)");
  EXPECT_EQ(error_of([&] { Path(g, {{2}, {4}, {0}}); }), "Path: node 0 is repeated");
}
