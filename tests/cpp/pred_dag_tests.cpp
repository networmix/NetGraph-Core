// PredDAG and Path (APPROACH 5): from_csr's complete validation in the documented order, the
// canonical predecessor order, structural equality, Path validation in its order, cost and
// to_dag().
#include <gtest/gtest.h>

#include <algorithm>
#include <functional>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

#include "netgraph/core/detail/pred_dag_builder.hpp"
#include "netgraph/core/pred_dag.hpp"
#include "test_utils.hpp"

using namespace netgraph::core;
using namespace netgraph::core::test;

namespace {

std::string error_of(const std::function<void()>& fn) {
  try {
    fn();
  } catch (const std::invalid_argument& e) {
    return e.what();
  }
  return "";
}

bool contains(const std::string& s, const std::string& part) { return s.find(part) != std::string::npos; }

// 0 -> 1 (two parallel members), 1 -> 2, 0 -> 2; ids follow (src, cost, dst, input).
Graph diamondish() {
  return make_graph(3, {{0, 1, 5, 1}, {0, 1, 5, 1}, {1, 2, 5, 1}, {0, 2, 5, 3}, {2, 0, 5, 1}});
}

using V = std::vector<std::int32_t>;

// The DAG's edges in ascending EdgeId (the builder's view of the storage).
std::vector<EdgeId> ascending(const PredDAG& dag) {
  const auto e = detail::PredDAGBuilder::ascending(dag);
  return {e.begin(), e.end()};
}

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
  EXPECT_EQ(ascending(dag), (std::vector<EdgeId>{0, 1, 2, 3}));
  // Structural equality: the builder's value from the same edge set is equal.
  const std::vector<EdgeId> edges = {3, 1, 0, 2};
  EXPECT_TRUE(detail::PredDAGBuilder::of_edges(g, edges) == dag);
  auto empty = PredDAG::from_csr(g, V {0, 0, 0, 0}, V {}, V {});
  EXPECT_FALSE(empty == dag);
  EXPECT_TRUE(ascending(empty).empty());
}

// The builder's entries for any edge subset in any input order: children ascending, via
// ascending within a child, the parent each via's tail, and the same value as from_csr.
TEST(PredDAG, BuilderGroupsAnyEdgeSubsetCanonically) {
  std::mt19937 rng(5);
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
    EXPECT_EQ(ascending(dag), subset) << trial;
  }
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
}

TEST(PredDAG, FromCsrEntryChecksInOrder) {
  auto g = diamondish();
  const V off = {0, 0, 2, 2};
  // parents in range, then via in range (entry by entry)
  EXPECT_TRUE(contains(error_of([&] { (void)PredDAG::from_csr(g, off, V {3, 0}, V {0, 1}); }), "parents[0] = 3"));
  EXPECT_TRUE(contains(error_of([&] { (void)PredDAG::from_csr(g, off, V {0, 0}, V {0, 5}); }), "via[1] = 5"));
  EXPECT_TRUE(contains(error_of([&] { (void)PredDAG::from_csr(g, off, V {-1, 0}, V {-1, 1}); }), "parents[0] = -1"));
  // via endpoints: wrong parent, wrong child
  EXPECT_TRUE(contains(error_of([&] { (void)PredDAG::from_csr(g, off, V {1, 0}, V {0, 1}); }), "joins (0, 1), not (1, 1)"));
  EXPECT_TRUE(contains(error_of([&] { (void)PredDAG::from_csr(g, V {0, 0, 0, 1}, V {0}, V {0}); }),
                       "joins (0, 1), not (0, 2)"));
  // duplicate (child, via) before the order check
  EXPECT_TRUE(contains(error_of([&] { (void)PredDAG::from_csr(g, off, V {0, 0}, V {1, 1}); }), "duplicate entry"));
  // a cycle (0 -> 1 -> 2 -> 0) before the order check
  EXPECT_TRUE(contains(error_of([&] { (void)PredDAG::from_csr(g, V {0, 1, 2, 3}, V {2, 0, 1}, V {4, 0, 3}); }),
                       "cycle"));
  // canonical order: ascending via within a node
  EXPECT_TRUE(contains(error_of([&] { (void)PredDAG::from_csr(g, off, V {0, 0}, V {1, 0}); }),
                       "not in ascending via EdgeId"));
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
  EXPECT_FALSE(error_of([&] { (void)PredDAG::from_csr(g, V {}, V {}, V {}); }).empty());
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
}

TEST(Path, ValidationInOrder) {
  auto g = diamondish();
  // range before everything (an empty hop elsewhere is reported later)
  EXPECT_TRUE(contains(error_of([&] { Path(g, {{}, {7}}); }), "member 7 out of edge range"));
  EXPECT_TRUE(contains(error_of([&] { Path(g, {{-1}}); }), "member -1 out of edge range"));
  // duplicate member before an empty hop
  EXPECT_TRUE(contains(error_of([&] { Path(g, {{0, 0}, {}}); }), "duplicate member 0"));
  EXPECT_TRUE(contains(error_of([&] { Path(g, {}); }), "hop list is empty"));
  EXPECT_TRUE(contains(error_of([&] { Path(g, {{0}, {}}); }), "hop 1 is empty"));
  // members of a hop share (u, v)
  EXPECT_TRUE(contains(error_of([&] { Path(g, {{0, 2}}); }), "do not share (u, v)"));
  // contiguous hops
  EXPECT_TRUE(contains(error_of([&] { Path(g, {{0}, {2}}); }), "not contiguous"));
  // no repeated node: 0 -> 2 -> 0
  EXPECT_TRUE(contains(error_of([&] { Path(g, {{2}, {4}}); }), "node 0 is repeated"));
}
