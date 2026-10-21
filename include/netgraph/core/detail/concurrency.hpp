/* Native-only rendezvous and failure probes of a hooks build. The caller owns the
   probe and its context; callbacks must never use Python or allocate at publication. */
#pragma once

#include "netgraph/core/detail/hooks_flag.hpp"

namespace netgraph::core::detail {

enum class ConcurrentStage {
  LedgerWait,
  LedgerAcquired,
  CacheQueued,
  CacheAcquired,
  BeforePublication,
  AfterLedgerExchange,
  NativeComplete,
  WorkerBegin,
  ThreadLaunch,
  NativeBegin,
  CacheTallyComplete,
  PolicyEngineInstalled,
  TunnelRosterSearch,
  RouteDefinitionCheck,
  RouteWalk,
  ArgumentBegin,
  Count,
};

#if NETGRAPH_CORE_PRIVATE_HOOKS
struct ConcurrentProbe {
  void (*notify)(void*, ConcurrentStage, const void*);
  void* context;
};
[[nodiscard]] const ConcurrentProbe* current_concurrent_probe() noexcept;
void set_concurrent_probe(const ConcurrentProbe* probe) noexcept;
struct ScopedConcurrentProbe {
  const ConcurrentProbe* previous;
  explicit ScopedConcurrentProbe(const ConcurrentProbe* probe) noexcept : previous(current_concurrent_probe()) { set_concurrent_probe(probe); }
  ~ScopedConcurrentProbe() { set_concurrent_probe(previous); }
};
void concurrent_event(ConcurrentStage stage, const void* object);
#else
inline void concurrent_event(ConcurrentStage, const void*) noexcept {}
#endif

// Publication probes observe a no-fail tail. A throwing test observer is a broken
// observer, never an exception the public operation can translate after committing.
inline void publication_event(ConcurrentStage stage, const void* object) noexcept {
  concurrent_event(stage, object);
}

} // namespace netgraph::core::detail
