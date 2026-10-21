// One probe slot shared by the native library and its Python bindings. Keeping the
// TLS behind functions also makes its identity independent of TU visibility/LTO.
#include "netgraph/core/detail/concurrency.hpp"
#if NETGRAPH_CORE_PRIVATE_HOOKS
namespace netgraph::core::detail {
namespace {
thread_local const ConcurrentProbe* probe_slot = nullptr;
}
const ConcurrentProbe* current_concurrent_probe() noexcept { return probe_slot; }
void set_concurrent_probe(const ConcurrentProbe* probe) noexcept { probe_slot = probe; }
void concurrent_event(ConcurrentStage stage, const void* object) {
  if (probe_slot != nullptr) probe_slot->notify(probe_slot->context, stage, object);
}
} // namespace netgraph::core::detail
#endif
