// This translation unit deliberately sees only the public ledger definition.
#include <gtest/gtest.h>

#include "netgraph/core/flow_ledger.hpp"

namespace netgraph::core::test {
namespace {
template<class T> concept Complete = requires { sizeof(T); };
static_assert(!Complete<detail::LedgerEntry>);
}

TEST(FlowLedgerHeader, IncompleteEntrySupportsPublicLifecycle) {
  const auto graph = std::make_shared<const Graph>(2, std::vector<NodeId> {0}, std::vector<NodeId> {1},
      std::vector<Amount> {10}, std::vector<Cost> {3}, 1.0);
  const std::vector<Amount> service {7}, reservation {5};
  FlowLedger first(graph), second(graph, service), third(graph, service, reservation);
  FlowLedger copied(third), moved(std::move(copied));
  EXPECT_EQ(first.residual(), (std::vector<Amount> {10}));
  EXPECT_EQ(second.residual(), service);
  EXPECT_EQ(moved.residual(), service);
  EXPECT_EQ(moved.unreserved(), reservation);
  EXPECT_TRUE(moved.records().empty());
  EXPECT_THROW((void)copied.uid(), std::logic_error);
  moved.reset();
  EXPECT_EQ(moved.residual(), (std::vector<Amount> {10}));
  std::optional<FlowLedger> local;
  local.emplace(graph); local.reset();
}

} // namespace netgraph::core::test
