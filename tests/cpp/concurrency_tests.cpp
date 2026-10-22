// Native synchronization, with explicit rendezvous rather than scheduling delays.
#include <gtest/gtest.h>

#include <algorithm>
#include <atomic>
#include <chrono>
#include <condition_variable>
#include <exception>
#include <future>
#include <functional>
#include <memory>
#include <mutex>
#include <optional>
#include <stdexcept>
#include <thread>
#include <type_traits>
#include <vector>

#include "netgraph/core/detail/cache_access.hpp"
#include "netgraph/core/detail/concurrency.hpp"
#include "netgraph/core/detail/ledger_access.hpp"
#include "netgraph/core/detail/policy_engine.hpp"
#include "netgraph/core/evaluators.hpp"
#include "netgraph/core/flow_policy.hpp"

namespace netgraph::core::test {
namespace {

using detail::CacheGate;
using detail::ConcurrentStage;
using detail::LedgerAccess;
using namespace std::chrono_literals;

std::shared_ptr<const Graph> one_link() {
  return std::make_shared<const Graph>(2, std::vector<NodeId> {0}, std::vector<NodeId> {1},
                                       std::vector<Amount> {10}, std::vector<Cost> {3}, 1.0);
}

#if NETGRAPH_CORE_PRIVATE_HOOKS
struct Rendezvous {
  std::mutex mutex;
  std::condition_variable changed;
  std::vector<ConcurrentStage> stages;
  std::optional<ConcurrentStage> fail;
  std::optional<ConcurrentStage> pause;
  bool released {false};
  // Reserve before the hook is armed: probe bookkeeping must not allocate while locked.
  Rendezvous() { stages.reserve(32); }
  static void notify(void* context, ConcurrentStage stage, const void*) {
    auto& self = *static_cast<Rendezvous*>(context);
    {
      std::unique_lock lock(self.mutex);
      self.stages.push_back(stage);
      self.changed.notify_all();
      if (self.pause == stage && !self.changed.wait_for(lock, 5s, [&] { return self.released; })) {
        throw std::runtime_error("native rendezvous timed out");
      }
    }
    if (self.fail == stage) throw std::runtime_error("injected admission failure");
  }
  bool wait(ConcurrentStage stage, std::chrono::milliseconds timeout = 5000ms) {
    std::unique_lock lock(mutex);
    return changed.wait_for(lock, timeout, [&] {
      return std::find(stages.begin(), stages.end(), stage) != stages.end();
    });
  }
  void release() {
    const std::lock_guard lock(mutex);
    released = true;
    changed.notify_all();
  }
};

struct ProbeScope {
  detail::ConcurrentProbe probe;
  const detail::ConcurrentProbe* previous;
  explicit ProbeScope(Rendezvous& context)
      : probe {Rendezvous::notify, &context}, previous(detail::current_concurrent_probe()) {
    detail::set_concurrent_probe(&probe);
  }
  ~ProbeScope() { detail::set_concurrent_probe(previous); }
};
#endif

TEST(Concurrency, LedgerReadsWaitAtTheSameGate) {
#if NETGRAPH_CORE_PRIVATE_HOOKS
  FlowLedger ledger(one_link());
  FlowPolicy policy(ledger, std::vector<NodeId>{0}, std::vector<NodeId>{1}, {});
  (void)policy.replace(3);
  const auto uid = ledger.uid();
  const auto graph = ledger.graph();
  const auto record = policy.records().front();
  const auto program = policy.program();
  const std::vector<std::pair<const char*, std::function<void()>>> readers {
    {"ledger.uid", [&] { EXPECT_EQ(ledger.uid(), uid); }},
    {"ledger.graph", [&] { EXPECT_EQ(ledger.graph(), graph); }},
    {"ledger.edge_capacity_units", [&] { EXPECT_EQ(ledger.edge_capacity_units().front(), 10); }},
    {"ledger.residual", [&] { EXPECT_EQ(ledger.residual(), (std::vector<Amount>{7})); }},
    {"ledger.carried", [&] { EXPECT_EQ(ledger.carried(), (std::vector<Amount>{3})); }},
    {"ledger.reserved", [&] { EXPECT_EQ(ledger.reserved(), (std::vector<Amount>{0})); }},
    {"ledger.unreserved", [&] { EXPECT_EQ(ledger.unreserved(), (std::vector<Amount>{10})); }},
    {"ledger.records", [&] { const auto rows = ledger.records(); ASSERT_EQ(rows.size(), 1); EXPECT_EQ(rows.front().flow_id, record.flow_id); EXPECT_EQ(rows.front().delivered, 3); }},
    {"ledger.record", [&] { EXPECT_EQ(ledger.record(record.flow_id).delivered, 3); }},
    {"ledger.saturated_boundary_edges", [&] { EXPECT_TRUE(ledger.saturated_boundary_edges(std::vector<NodeId>{0}).empty()); }},
    {"policy.delivered", [&] { EXPECT_EQ(policy.delivered(), 3); }},
    {"policy.records", [&] { const auto rows = policy.records(); ASSERT_EQ(rows.size(), 1); EXPECT_EQ(rows.front().flow_id, record.flow_id); EXPECT_EQ(rows.front().delivered, 3); }},
    {"policy.record", [&] { EXPECT_EQ(policy.record(record.flow_id).delivered, 3); }},
    {"policy.program", [&] { EXPECT_EQ(policy.program(), program); }},
  };
  for (const auto& [name, read] : readers) {
    SCOPED_TRACE(name);
    auto gate = LedgerAccess::lock(ledger);
    Rendezvous rendezvous;
    auto result = std::async(std::launch::async, [&] {
      const ProbeScope probe(rendezvous);
      read();
    });
    EXPECT_TRUE(rendezvous.wait(ConcurrentStage::LedgerWait));
    EXPECT_FALSE(rendezvous.wait(ConcurrentStage::LedgerAcquired, 100ms));
    EXPECT_EQ(result.wait_for(0s), std::future_status::timeout);
    gate.unlock();
    result.get();
  }
  // Immutable policy configuration is explicitly outside the gate's reader set.
  auto gate = LedgerAccess::lock(ledger);
  auto immutable = std::async(std::launch::async, [&] { return policy.config().admission; });
  const auto status = immutable.wait_for(5s);
  EXPECT_EQ(status, std::future_status::ready);
  gate.unlock();
  EXPECT_EQ(immutable.get(), Admission::Lossless);
#else
  GTEST_SKIP() << "requires native rendezvous hooks";
#endif
}

TEST(Concurrency, PolicyReadinessFollowsEngineInstallation) {
#if NETGRAPH_CORE_PRIVATE_HOOKS
  FlowLedger ledger(one_link());
  auto policy = detail::PolicyAccess::shell();
  auto prepared = detail::PolicyAccess::prepare(ledger, std::vector<NodeId>{0}, std::vector<NodeId>{1}, {}, nullptr, nullptr);
  Rendezvous rendezvous;
  rendezvous.pause = ConcurrentStage::PolicyEngineInstalled;
  auto construction = std::async(std::launch::async, [&] {
    const ProbeScope probe(rendezvous);
    detail::PolicyAccess::construct(*policy, std::move(prepared), nullptr);
  });
  const bool paused = rendezvous.wait(ConcurrentStage::PolicyEngineInstalled);
  EXPECT_TRUE(paused);
  // Read only the atomic readiness bit: a wrong order fails safely even when the
  // engine has not yet been installed. Do not dereference an initializing shell.
  EXPECT_FALSE(detail::PolicyAccess::constructed(*policy));
  EXPECT_EQ(construction.wait_for(0s), std::future_status::timeout);
  rendezvous.release();
  construction.get();
  ASSERT_TRUE(detail::PolicyAccess::constructed(*policy));
  EXPECT_EQ(policy->delivered(), 0);
  EXPECT_EQ(policy->place(3).delivered, 3);
#else
  GTEST_SKIP() << "requires native rendezvous hooks";
#endif
}

// Both frozen readers and building callers enroll behind the lifecycle request. The
// queue marker proves order, independent of when each OS thread happens to run next.
TEST(Concurrency, CacheLifecycleCannotBeOvertaken) {
#if NETGRAPH_CORE_PRIVATE_HOOKS
  for (bool frozen : {false, true}) {
    CacheGate gate;
    {
      const CacheGate::Lease lease(gate, true);
      gate.freeze(frozen);
    }
    std::optional<CacheGate::Lease> active;
    active.emplace(gate);
    Rendezvous writer_ready;
    Rendezvous later_ready;
    std::vector<int> order;
    std::mutex order_mutex;
    auto writer = std::async(std::launch::async, [&] {
      const ProbeScope probe(writer_ready);
      const CacheGate::Lease lease(gate, true);
      const std::lock_guard lock(order_mutex);
      order.push_back(1);
    });
    EXPECT_TRUE(writer_ready.wait(ConcurrentStage::CacheQueued));
    auto later = std::async(std::launch::async, [&] {
      const ProbeScope probe(later_ready);
      const CacheGate::Lease lease(gate);
      const std::lock_guard lock(order_mutex);
      order.push_back(2);
    });
    EXPECT_TRUE(later_ready.wait(ConcurrentStage::CacheQueued));
    EXPECT_EQ(writer.wait_for(0s), std::future_status::timeout);
    EXPECT_EQ(later.wait_for(0s), std::future_status::timeout);
    active.reset();
    writer.get();
    later.get();
    EXPECT_EQ(order, (std::vector<int> {1, 2}));
  }
#else
  GTEST_SKIP() << "requires native rendezvous hooks";
#endif
}

TEST(Concurrency, FailedAdmissionLeavesNoQueuedNodeOrLease) {
#if NETGRAPH_CORE_PRIVATE_HOOKS
  for (ConcurrentStage point : {ConcurrentStage::CacheQueued, ConcurrentStage::CacheAcquired}) {
    CacheGate gate;
    Rendezvous failure;
    failure.fail = point;
    {
      const ProbeScope probe(failure);
      EXPECT_THROW((CacheGate::Lease(gate)), std::runtime_error);
    }
    ASSERT_FALSE(gate.has_waiters());
    // Both a following lifecycle operation and lookup must terminate. CacheQueued
    // fails while the node is enrolled; CacheAcquired fails after ownership changes.
    {
      const CacheGate::Lease lifecycle(gate, true);
      gate.freeze(true);
    }
    const CacheGate::Lease reader(gate);
    EXPECT_TRUE(gate.frozen());
  }
#else
  GTEST_SKIP() << "requires native rendezvous hooks";
#endif
}

TEST(Concurrency, ImmutableGraphAnchorSurvivesMove) {
  auto graph = one_link();
  std::weak_ptr<const Graph> weak = graph;
  std::optional<FlowLedger> source;
  source.emplace(graph);
  const auto& anchor = source->graph();
  graph.reset();
  {
    FlowLedger target(std::move(*source));
    EXPECT_EQ(anchor.get(), target.graph().get());
    EXPECT_THROW((void)source->graph(), std::logic_error);
  }
  ASSERT_FALSE(weak.expired());
  EXPECT_EQ(anchor->num_edges(), 1);
  source.reset();
  EXPECT_TRUE(weak.expired());
  // A copied ledger owns its immutable anchor after every original owner dies.
  // Check liveness before dereferencing so a borrowing-copy plant fails safely.
  std::optional<FlowLedger> copy;
  {
    auto original_graph = one_link();
    weak = original_graph;
    FlowLedger original(original_graph);
    copy.emplace(original);
  }
  ASSERT_FALSE(weak.expired());
  EXPECT_EQ(copy->graph()->num_edges(), 1);
  copy.reset();
  EXPECT_TRUE(weak.expired());
}

TEST(Concurrency, SameLedgerCallsHaveASerialHistory) {
  FlowLedger ledger(one_link());
  FlowPolicy a(ledger, std::vector<NodeId>{0}, std::vector<NodeId>{1}, {});
  FlowPolicy b(ledger, std::vector<NodeId>{0}, std::vector<NodeId>{1}, {});
  std::vector<std::future<Amount>> work;
  for (int i = 0; i < 4; ++i) {
    work.push_back(std::async(std::launch::async, [&, i] {
      Amount delivered = 0;
      for (int n = 0; n < 10; ++n) delivered += (i % 2 == 0 ? a : b).place(1).call_delivered;
      return delivered;
    }));
  }
  Amount delivered = 0;
  for (auto& result : work) delivered += result.get();
  EXPECT_EQ(delivered, 10);
  EXPECT_EQ(a.delivered() + b.delivered(), 10);
  EXPECT_EQ(ledger.residual(), std::vector<Amount>{0});
}

TEST(Concurrency, IndependentLedgersDoNotShareAGate) {
  const auto graph = one_link();
  FlowLedger a(graph), b(graph);
  FlowPolicy policy(b, std::vector<NodeId>{0}, std::vector<NodeId>{1}, {});
  const auto gate = LedgerAccess::lock(a);
  auto result = std::async(std::launch::async, [&] { return policy.place(3).call_delivered; });
  EXPECT_EQ(result.wait_for(5s), std::future_status::ready);
  EXPECT_EQ(result.get(), 3);
}

TEST(Concurrency, WeightedSourceCallsSerializeAndDoNotShareAGateAcrossLedgers) {
  const auto graph=std::make_shared<const Graph>(3,std::vector<NodeId>{0,1},std::vector<NodeId>{2,2},
      std::vector<Amount>{10,10},std::vector<Cost>{1,3},1.0);
  CentralizedTE model(SourceAllocation::Weighted);
  model.source_weights=std::make_shared<const std::vector<Amount>>(std::vector<Amount>{2,1,0});
  FlowLedger a(graph),b(graph);
  FlowPolicy p(a,std::vector<NodeId>{0,1},std::vector<NodeId>{2},FlowPolicyConfig{model});
  FlowPolicy q(b,std::vector<NodeId>{0,1},std::vector<NodeId>{2},FlowPolicyConfig{model});
  {
    const auto gate=LedgerAccess::lock(a);
    auto result=std::async(std::launch::async,[&]{return q.place(3).call_delivered;});
    EXPECT_EQ(result.wait_for(5s),std::future_status::ready);EXPECT_EQ(result.get(),3);
  }
  std::vector<std::future<Amount>> tasks;
  for(int i=0;i<4;++i)tasks.push_back(std::async(std::launch::async,[&]{Amount sum=0;for(int j=0;j<10;++j)sum+=p.place(1).call_delivered;return sum;}));
  Amount sum=0;for(auto& t:tasks)sum+=t.get();
  EXPECT_EQ(sum,15);EXPECT_EQ(p.delivered(),15);
  const auto r=p.replace(15);
  EXPECT_EQ(r.by_source[0].delivered,10);EXPECT_EQ(r.by_source[1].delivered,5);
}

TEST(Concurrency, ScalarAndListsDoNotLeaseTheConstructionCache) {
  const auto graph = one_link();
  ShortestPathCache cache(*graph);
  FlowLedger ledger(graph);
  FlowPolicy a(ledger, std::vector<NodeId>{0}, std::vector<NodeId>{1}, {}, nullptr, &cache);
  FlowPolicy b(ledger, std::vector<NodeId>{0}, std::vector<NodeId>{1}, {}, nullptr, &cache);
  const auto before = cache.stats();
  const CacheGate::Lease lease(detail::DagCacheAccess::gate(cache), true);
  auto result = std::async(std::launch::async, [&] {
    const auto scalar = a.place(2);
    const std::vector<Action> actions {Place{b, 3}};
    return scalar.call_delivered + apply_actions(ledger, actions)[0].call_delivered;
  });
  EXPECT_EQ(result.wait_for(5s), std::future_status::ready);
  EXPECT_EQ(result.get(), 5);
  EXPECT_EQ(cache.stats().hits, before.hits);
  EXPECT_EQ(cache.stats().misses, before.misses);
}

TEST(Concurrency, ArgumentChecksPrecedeEveryStateCheck) {
  FlowLedger ledger(one_link());
  FlowPolicy a(ledger, std::vector<NodeId>{0}, std::vector<NodeId>{1}, {});
  FlowPolicy b(ledger, std::vector<NodeId>{0}, std::vector<NodeId>{1}, {});
  auto text = [](auto&& body) {
    try { body(); } catch (const std::invalid_argument& error) { return std::string(error.what()); }
    return std::string("no exception");
  };
  EXPECT_EQ(text([&] { (void)a.forward(-1); }), "offer must be an integer in [0, 2^62) or permitted +inf");
  EXPECT_EQ(text([&] { (void)a.forward(1); }), "forward requires an installed program");
  const std::vector<Action> actions {Forward{a, 1}, Replace{b, -1}};
  EXPECT_EQ(text([&] { (void)apply_actions(ledger, actions); }), "offer must be an integer in [0, 2^62) or permitted +inf");
  EXPECT_TRUE(ledger.records().empty());
  EXPECT_FALSE(a.program());
  EXPECT_FALSE(b.program());
  FlowLedger source(one_link());
  FlowLedger moved(std::move(source));
  EXPECT_EQ(text([&] { (void)apply_actions(source, {}, static_cast<BatchMode>(-1)); }),
            "mode: -1 is not a value of BatchMode");
}

TEST(Concurrency, GateCoversExchangeAndEveryEpoch) {
#if NETGRAPH_CORE_PRIVATE_HOOKS
  FlowLedger ledger(one_link());
  FlowPolicy a(ledger, std::vector<NodeId>{0}, std::vector<NodeId>{1}, {});
  FlowPolicy b(ledger, std::vector<NodeId>{0}, std::vector<NodeId>{1}, {});
  Rendezvous writer_stage, reader_stage;
  writer_stage.pause = ConcurrentStage::AfterLedgerExchange;
  auto writer = std::async(std::launch::async, [&] {
    const ProbeScope probe(writer_stage);
    const std::vector<Action> actions {Replace{a, 2}, Replace{b, 3}};
    return apply_actions(ledger, actions);
  });
  EXPECT_TRUE(writer_stage.wait(ConcurrentStage::AfterLedgerExchange));
  auto reader = std::async(std::launch::async, [&] {
    const ProbeScope probe(reader_stage);
    const auto gate = LedgerAccess::lock(ledger);
    const auto& ae = detail::PolicyAccess::engine(a);
    const auto& be = detail::PolicyAccess::engine(b);
    return ae.epoch.program != nullptr && be.epoch.program != nullptr && LedgerAccess::carried(ledger)[0] == 5;
  });
  EXPECT_TRUE(reader_stage.wait(ConcurrentStage::LedgerWait));
  EXPECT_EQ(reader.wait_for(0s), std::future_status::timeout);
  writer_stage.release();
  EXPECT_EQ(writer.get().size(), 2U);
  EXPECT_TRUE(reader.get());
#else
  GTEST_SKIP() << "requires native rendezvous hooks";
#endif
}

TEST(Concurrency, FailureBeforePublicationReleasesTheGate) {
#if NETGRAPH_CORE_PRIVATE_HOOKS
  FlowLedger ledger(one_link());
  FlowPolicy policy(ledger, std::vector<NodeId>{0}, std::vector<NodeId>{1}, {});
  Rendezvous failure;
  failure.fail = ConcurrentStage::BeforePublication;
  {
    const ProbeScope probe(failure);
    EXPECT_THROW((void)policy.place(3), std::runtime_error);
  }
  EXPECT_TRUE(ledger.records().empty());
  EXPECT_FALSE(policy.program());
  EXPECT_EQ(policy.place(3).call_delivered, 3);
#else
  GTEST_SKIP() << "requires native rendezvous hooks";
#endif
}

TEST(Concurrency, LowestWorkerFailureWins) {
#if NETGRAPH_CORE_PRIVATE_HOOKS
  struct Failures {
    std::mutex mutex;
    std::condition_variable changed;
    bool later_failed {false};
    static void notify(void* context, ConcurrentStage stage, const void* object) {
      if (stage != ConcurrentStage::WorkerBegin) return;
      auto& self = *static_cast<Failures*>(context);
      const auto index = *static_cast<const std::size_t*>(object);
      std::unique_lock lock(self.mutex);
      if (index == 2) {
        self.later_failed = true;
        self.changed.notify_all();
        throw std::runtime_error("task two failed first");
      }
      if (index == 0) {
        if (!self.changed.wait_for(lock, 5s, [&] { return self.later_failed; })) throw std::runtime_error("worker rendezvous timed out");
        throw std::runtime_error("task zero failed last");
      }
    }
  } failures;
  const auto graph = one_link();
  const detail::ConcurrentProbe probe {Failures::notify, &failures};
  const detail::ScopedConcurrentProbe installed(&probe);
  const std::vector<std::pair<std::vector<NodeId>, std::vector<NodeId>>> pairs(3, {{0}, {1}});
  try {
    (void)evaluate_pairs(*graph, pairs, {}, {}, 3);
    FAIL() << "both failing workers were ignored";
  } catch (const std::runtime_error& error) {
    EXPECT_STREQ(error.what(), "task zero failed last");
  }
#else
  GTEST_SKIP() << "requires native rendezvous hooks";
#endif
}

TEST(Concurrency, PartialThreadLaunchJoinsAndReleasesLease) {
#if NETGRAPH_CORE_PRIVATE_HOOKS
  struct Failure {
    std::size_t after;
    static void notify(void* context, ConcurrentStage stage, const void* object) {
      if (stage == ConcurrentStage::ThreadLaunch && *static_cast<const std::size_t*>(object) == static_cast<Failure*>(context)->after) {
        throw std::runtime_error("injected thread launch failure");
      }
    }
  } failure;
  const auto graph = one_link();
  ShortestPathCache cache(*graph);
  (void)cache.lookup(*graph, std::vector<NodeId>{1}, Direction::Reverse, nullptr);
  cache.freeze();
  const std::vector<std::pair<std::vector<NodeId>, std::vector<NodeId>>> pairs(3, {{0}, {1}});
  EvaluationOptions options;
  options.cache = &cache;
  for (std::size_t after : {0U, 1U, 2U}) {
    failure.after = after;
    const auto before = cache.stats();
    const detail::ConcurrentProbe probe {Failure::notify, &failure};
    {
      const detail::ScopedConcurrentProbe installed(&probe);
      EXPECT_THROW((void)evaluate_pairs(*graph, pairs, options, {}, 3), std::runtime_error);
    }
    EXPECT_EQ(cache.stats().hits, before.hits);
    cache.thaw();
    cache.freeze();
    EXPECT_EQ(evaluate_pairs(*graph, pairs, options, {}, 3).size(), 3U);
  }
#else
  GTEST_SKIP() << "requires native rendezvous hooks";
#endif
}

TEST(Concurrency, NativeArgumentsAreCopiedBeforeWaiting) {
#if NETGRAPH_CORE_PRIVATE_HOOKS
  // All scalar operations and both list modes have their own public copy point.
  for (int entry = 0; entry < 6; ++entry) {
  const auto graph = one_link();
  FlowLedger ledger(graph);
  FlowPolicy policy(ledger, std::vector<NodeId>{0}, std::vector<NodeId>{1}, {.admission = Admission::Lossy});
  (void)policy.replace(2);
  OutputOptions outputs;
  outputs.residual = true;
  auto gate = LedgerAccess::lock(ledger);
  Rendezvous waiting;
  auto result = std::async(std::launch::async, [&] {
    const ProbeScope probe(waiting);
    switch (entry) {
      case 0: return policy.place(3, outputs);
      case 1: return policy.evaluate(3, {}, outputs);
      case 2: return policy.replace(3, {}, outputs);
      case 3: return policy.forward(3, nullptr, outputs);
      default: return apply_actions(ledger, std::vector<Action>{Replace{policy, 3}},
          entry == 4 ? BatchMode::Ordered : BatchMode::Simultaneous, outputs)[0];
    }
  });
  EXPECT_TRUE(waiting.wait(ConcurrentStage::LedgerWait));
  outputs.residual = false;  // after the call's copy, before it is admitted
  gate.unlock();
  const auto value = result.get();
  EXPECT_TRUE(value.residual.has_value()) << entry;
  if (value.residual) EXPECT_EQ(*value.residual, std::vector<Amount>{entry == 0 ? 5 : 7}) << entry;
  EXPECT_EQ(ledger.residual(), std::vector<Amount>{entry == 1 ? 8 : (entry == 0 ? 5 : 7)}) << entry;
  }
#else
  GTEST_SKIP() << "requires native rendezvous hooks";
#endif
}

TEST(Concurrency, NativeEvaluatorInputsAreOwnedWhileQueued) {
#if NETGRAPH_CORE_PRIVATE_HOOKS
  for (int entry = 0; entry < 3; ++entry) {
  const auto graph = one_link();
  ShortestPathCache cache(*graph);
  (void)cache.lookup(*graph, std::vector<NodeId>{1}, Direction::Reverse, nullptr);
  std::optional<CacheGate::Lease> active;
  active.emplace(detail::DagCacheAccess::gate(cache));
  EvaluationOptions options;
  options.cache = &cache;
  options.offer = std::optional<Amount>{3};
  std::vector<std::pair<std::vector<NodeId>, std::vector<NodeId>>> pairs {{{0}, {1}}};
  std::vector<NodeId> sources {0}, destinations {1};
  std::vector<const Eligibility*> masks {nullptr};
  const auto foreign = one_link();
  const Eligibility wrong_mask = foreign->eligibility({}, {});
  Rendezvous waiting;
  auto result = std::async(std::launch::async, [&] {
    const ProbeScope probe(waiting);
    if (entry == 0) return evaluate_pairs(*graph, pairs, options, masks, 1)[0].delivered;
    if (entry == 1) return evaluate_demand(*graph, sources, destinations, options).delivered;
    return evaluate_edge_failures(*graph, sources, destinations, options, 1).baseline.delivered;
  });
  EXPECT_TRUE(waiting.wait(ConcurrentStage::CacheQueued));
  pairs[0].second[0] = 99;
  options.offer = std::optional<Amount>{7};
  sources[0] = 99;
  destinations[0] = 99;
  masks[0] = &wrong_mask;
  active.reset();
  EXPECT_EQ(result.get(), 3) << entry;
  EXPECT_GE(cache.stats().hits, 1) << entry;
  }
#else
  GTEST_SKIP() << "requires native rendezvous hooks";
#endif
}

TEST(Concurrency, FrozenTallyPublishesInsideTheLease) {
#if NETGRAPH_CORE_PRIVATE_HOOKS
  const auto graph = one_link();
  ShortestPathCache cache(*graph);
  (void)cache.lookup(*graph, std::vector<NodeId>{1}, Direction::Reverse, nullptr);
  cache.freeze();
  const auto before = cache.stats().hits;
  EvaluationOptions options;
  options.cache = &cache;
  options.offer = std::optional<Amount>{3};
  const std::vector<std::pair<std::vector<NodeId>, std::vector<NodeId>>> pairs(3, {{0}, {1}});
  Rendezvous finishing, waiting;
  finishing.pause = ConcurrentStage::CacheTallyComplete;
  auto result = std::async(std::launch::async, [&] {
    const ProbeScope probe(finishing);
    return evaluate_pairs(*graph, pairs, options, {}, 3);
  });
  EXPECT_TRUE(finishing.wait(ConcurrentStage::CacheTallyComplete));
  EXPECT_EQ(cache.stats().hits, before + 3);
  auto lifecycle = std::async(std::launch::async, [&] {
    const ProbeScope probe(waiting);
    cache.thaw();
  });
  EXPECT_TRUE(waiting.wait(ConcurrentStage::CacheQueued));
  EXPECT_EQ(lifecycle.wait_for(0s), std::future_status::timeout);
  finishing.release();
  EXPECT_EQ(result.get().size(), 3U);
  lifecycle.get();
  EXPECT_EQ(cache.stats().hits, before + 3);
#else
  GTEST_SKIP() << "requires native rendezvous hooks";
#endif
}

TEST(Concurrency, NativeCacheRootsAreOwnedWhileQueued) {
#if NETGRAPH_CORE_PRIVATE_HOOKS
  const auto graph = one_link();
  ShortestPathCache cache(*graph);
  const auto retained = cache.lookup(*graph, std::vector<NodeId>{1}, Direction::Reverse, nullptr);
  std::optional<CacheGate::Lease> active;
  active.emplace(detail::DagCacheAccess::gate(cache));
  std::vector<NodeId> roots {1};
  Rendezvous waiting;
  auto result = std::async(std::launch::async, [&] {
    const ProbeScope probe(waiting);
    return cache.lookup(*graph, roots, Direction::Reverse, nullptr);
  });
  EXPECT_TRUE(waiting.wait(ConcurrentStage::CacheQueued));
  roots[0] = 0;
  active.reset();
  EXPECT_EQ(result.get(), retained);
  EXPECT_EQ(cache.stats().entries, 1);
  EXPECT_EQ(cache.stats().hits, 1);
#else
  GTEST_SKIP() << "requires native rendezvous hooks";
#endif
}

static_assert(std::is_same_v<decltype(std::declval<const FlowLedger&>().graph()),
                             const std::shared_ptr<const Graph>&>);
static_assert(noexcept(std::declval<const ShortestPathCache&>().stats()));
static_assert(!noexcept(std::declval<ShortestPathCache&>().freeze()));
static_assert(!noexcept(std::declval<ShortestPathCache&>().thaw()));

TEST(Concurrency, TunnelRosterChecksRunInParallelBeforeAdmission) {
#if NETGRAPH_CORE_PRIVATE_HOOKS
  struct Overlap {
    std::mutex mutex;
    std::condition_variable changed;
    std::size_t searches {0};
    bool admitted {false};
    static void notify(void* context, ConcurrentStage stage, const void*) {
      auto& self = *static_cast<Overlap*>(context);
      std::unique_lock lock(self.mutex);
      if (stage == ConcurrentStage::LedgerWait || stage == ConcurrentStage::CacheQueued) self.admitted = true;
      if (stage != ConcurrentStage::TunnelRosterSearch) return;
      EXPECT_FALSE(self.admitted);
      ++self.searches;
      self.changed.notify_all();
      if (self.searches == 1 && !self.changed.wait_for(lock, 2s, [&] { return self.searches >= 2; })) {
        throw std::runtime_error("roster checks did not overlap");
      }
    }
  };
  const auto graph = one_link();
  const std::vector<std::pair<std::vector<NodeId>, std::vector<NodeId>>> pairs(8, {{0}, {1}});
  CentralizedTE central(SourceAllocation::Equal);
  central.tunnel_count = TunnelCount{1};
  for (const Routing& routing : {Routing{DistributedTE{TunnelCount{1}}}, Routing{central}}) {
    Overlap overlap;
    EvaluationOptions options;
    options.config.routing = routing;
    options.offer = Amount{1};
    options.design_demand = Amount{1};
    const detail::ConcurrentProbe probe {Overlap::notify, &overlap};
    const detail::ScopedConcurrentProbe installed(&probe);
    const auto result = evaluate_pairs(*graph, pairs, options, {}, 2);
    ASSERT_EQ(result.size(), 8U);
    for (const auto& r : result) EXPECT_EQ(r.delivered, 1);
    EXPECT_EQ(overlap.searches, pairs.size());
    EXPECT_TRUE(overlap.admitted);
  }
#else
  GTEST_SKIP() << "requires native rendezvous hooks";
#endif
}

TEST(Concurrency, TunnelPreflightReportsTheLowestErrorBeforeLaterChecks) {
#if NETGRAPH_CORE_PRIVATE_HOOKS
  struct Failures {
    std::mutex mutex;
    std::condition_variable changed;
    bool later {false};
    static void notify(void* context, ConcurrentStage stage, const void* object) {
      if (stage != ConcurrentStage::WorkerBegin) return;
      auto& self = *static_cast<Failures*>(context);
      const auto index = *static_cast<const std::size_t*>(object);
      std::unique_lock lock(self.mutex);
      if (index == 0) {
        if (!self.changed.wait_for(lock, 2s, [&] { return self.later; })) throw std::runtime_error("no later worker");
        throw std::runtime_error("first pair failed last");
      }
      self.later = true;
      self.changed.notify_all();
      throw std::runtime_error("second pair failed first");
    }
  } failures;
  const auto graph = one_link();
  EvaluationOptions options;
  options.config.routing = DistributedTE{TunnelCount{1}};
  options.offer = Amount{1};
  options.design_demand = Amount{1};
  const std::vector<std::pair<std::vector<NodeId>, std::vector<NodeId>>> pairs {{{0}, {1}}, {{0}, {1}}, {{}, {1}}};
  const detail::ConcurrentProbe probe {Failures::notify, &failures};
  const detail::ScopedConcurrentProbe installed(&probe);
  try {
    (void)evaluate_pairs(*graph, pairs, options, {}, 2);
    FAIL() << "the first-pass errors were ignored";
  } catch (const std::runtime_error& error) {
    EXPECT_STREQ(error.what(), "first pair failed last");
  }
#else
  GTEST_SKIP() << "requires native rendezvous hooks";
#endif
}

TEST(Concurrency, TunnelPreflightJoinsAfterPartialLaunchFailure) {
#if NETGRAPH_CORE_PRIVATE_HOOKS
  struct Failure {
    std::size_t after {0};
    static void notify(void* context, ConcurrentStage stage, const void* object) {
      if (stage == ConcurrentStage::ThreadLaunch &&
          *static_cast<const std::size_t*>(object) == static_cast<Failure*>(context)->after) {
        throw std::runtime_error("injected preflight launch failure");
      }
    }
  } failure;
  const auto graph = one_link();
  const std::vector<std::pair<std::vector<NodeId>, std::vector<NodeId>>> pairs(8, {{0}, {1}});
  EvaluationOptions options;
  options.config.routing = DistributedTE{TunnelCount{1}};
  options.offer = Amount{1};
  options.design_demand = Amount{1};
  for (std::size_t after : {0U, 1U, 2U}) {
    failure.after = after;
    const detail::ConcurrentProbe probe {Failure::notify, &failure};
    {
      const detail::ScopedConcurrentProbe installed(&probe);
      EXPECT_THROW((void)evaluate_pairs(*graph, pairs, options, {}, 3), std::runtime_error);
    }
    const auto retry = evaluate_pairs(*graph, pairs, options, {}, 3);
    ASSERT_EQ(retry.size(), pairs.size());
    for (const auto& result : retry) EXPECT_EQ(result.delivered, 1);
  }
#else
  GTEST_SKIP() << "requires native rendezvous hooks";
#endif
}

static_assert(std::is_nothrow_destructible_v<FlowPolicy>);
static_assert(!noexcept(std::declval<const FlowPolicy&>().program()));
static_assert(std::is_nothrow_move_assignable_v<detail::PolicyEpoch>);
static_assert(std::is_nothrow_move_assignable_v<PlacementResult>);

} // namespace
} // namespace netgraph::core::test
