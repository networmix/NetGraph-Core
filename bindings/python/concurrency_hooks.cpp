// Native rendezvous controls for progress/atomicity tests. Python is called only
// by run(), outside every native lock; the probes themselves touch native state.
#include "registration.hpp"
#if NETGRAPH_CORE_PRIVATE_HOOKS
#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <condition_variable>
#include <cstdio>
#include <cstdlib>
#include <memory>
#include <mutex>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

#include "netgraph/core/detail/concurrency.hpp"
#include "netgraph/core/detail/concurrency_models.hpp"
#include "netgraph/core/flow_policy.hpp"
#include "netgraph/core/flow_ledger.hpp"
#include "netgraph/core/shortest_path_cache.hpp"
#include "native_call.hpp"
#include "values.hpp"

namespace py = pybind11;
namespace netgraph::core::bindings {
namespace {
using Stage = detail::ConcurrentStage;
constexpr std::array<const char*, static_cast<std::size_t>(Stage::Count)> names {
  "ledger_wait", "ledger_acquired", "cache_queued", "cache_acquired",
  "before_publication", "after_exchange", "native_complete", "worker_begin", "thread_launch", "native_begin",
  "cache_tally_complete", "policy_engine_installed", "tunnel_roster_search", "route_definition_check", "route_walk", "argument_begin"
};
Stage stage_of(const std::string& name) {
  for (std::size_t i = 0; i < names.size(); ++i) if (name == names[i]) return static_cast<Stage>(i);
  throw py::value_error("unknown native probe stage: " + name);
}
struct NativeProbe {
  std::mutex mutex;
  std::condition_variable changed;
  std::array<std::size_t, names.size() + 2> counts {};
  std::optional<Stage> pause, fail;
  std::vector<std::size_t> failing_workers;
  bool reverse_failures {false};
  bool later_failed {false};
  bool released {false};
  int auto_release_ms {0};
  detail::ConcurrentProbe probe {notify, this};

  NativeProbe(const std::string& paused, const std::string& failed, std::vector<std::size_t> workers, bool reverse,
              int automatic)
      : failing_workers(std::move(workers)), reverse_failures(reverse), auto_release_ms(automatic) {
    if (automatic < 0 || automatic > 10000) throw py::value_error("auto_release_ms must be in [0, 10000]");
    if (!paused.empty()) pause = stage_of(paused);
    if (!failed.empty()) fail = stage_of(failed);
  }
  static void notify(void* context, Stage stage, const void* object) {
    auto& self = *static_cast<NativeProbe*>(context);
    std::unique_lock lock(self.mutex);
    ++self.counts[static_cast<std::size_t>(stage)];
    self.changed.notify_all();
    if (stage == Stage::WorkerBegin && !self.failing_workers.empty()) {
      const auto index = *static_cast<const std::size_t*>(object);
      if (std::find(self.failing_workers.begin(), self.failing_workers.end(), index) != self.failing_workers.end()) {
        if (self.reverse_failures && index == self.failing_workers.front()) {
          if (!self.changed.wait_for(lock, std::chrono::seconds(10), [&] { return self.later_failed; })) {
            throw std::runtime_error("worker rendezvous timed out");
          }
        } else {
          self.later_failed = true;
          self.changed.notify_all();
        }
        throw std::runtime_error("injected worker failure " + std::to_string(index));
      }
    }
    if (self.pause == stage && !self.released) {
      ++self.counts[names.size()];
      const auto duration = std::chrono::milliseconds(self.auto_release_ms == 0 ? 10000 : self.auto_release_ms);
      const bool signalled = self.changed.wait_for(lock, duration, [&] { return self.released; });
      ++self.counts[names.size() + 1];
      if (!signalled && self.auto_release_ms == 0) throw std::runtime_error("native rendezvous timed out");
      self.released = true;
    }
    if (self.fail == stage) throw std::runtime_error("injected native failure: " + std::string(names[static_cast<std::size_t>(stage)]));
  }
  bool wait(Stage stage, double seconds, std::size_t count) {
    std::unique_lock lock(mutex);
    return changed.wait_for(lock, std::chrono::duration<double>(seconds), [&] { return counts[static_cast<std::size_t>(stage)] >= count; });
  }
  void release() {
    const std::lock_guard lock(mutex);
    released = true;
    changed.notify_all();
  }
};

// Py_AtExit runs after Python teardown. Keep only native owners here, so this
// fixture cannot accidentally retain a Python global through the daemon frame.
struct ShutdownWaiter {
  std::shared_ptr<FlowLedger> ledger;
  std::shared_ptr<ShortestPathCache> cache;
  std::shared_ptr<NativeProbe> probe;
  std::unique_ptr<FlowPolicy> policy;
};
std::atomic<ShutdownWaiter*> shutdown_waiter {nullptr};
void native_shutdown() noexcept {
  std::unique_ptr<ShutdownWaiter> state(shutdown_waiter.exchange(nullptr));
  if (!state) return;
  if (Py_IsInitialized() != 0) std::_Exit(91);
  {
    const std::lock_guard lock(state->probe->mutex);
    if (state->probe->counts[names.size()] <= state->probe->counts[names.size() + 1]) std::_Exit(92);
  }
  std::fputs("CHILD: native finalizer waiting after Python teardown\n", stdout);
  std::fflush(stdout);
  try {
    state->policy.reset();
    if (state->cache) state->cache->thaw();
    state.reset();
  } catch (...) {
    std::_Exit(93);
  }
  std::fputs("CHILD: native finalizer completed\n", stdout);
  std::fflush(stdout);
}
} // namespace
void register_concurrency_hooks(py::module_& m) {
  py::class_<NativeProbe, std::shared_ptr<NativeProbe>>(m, "_ConcurrencyProbe")
    .def(py::init<const std::string&, const std::string&, std::vector<std::size_t>, bool, int>(),
         py::arg("pause") = "", py::arg("fail") = "", py::arg("workers") = std::vector<std::size_t>{}, py::arg("reverse") = false,
         py::arg("auto_release_ms") = 0)
    .def("run", [](NativeProbe& self, const py::function& callback, const py::args& args) {
      const detail::ScopedConcurrentProbe scope(&self.probe);
      return callback(*args);
    }, py::arg("callable"))
    .def("wait", [](NativeProbe& self, const std::string& stage, double seconds, std::size_t count) {
      const auto point = stage_of(stage);
      // Do not run native_call's entry probe inside the probe's own wait.
      py::gil_scoped_release detached;
      return self.wait(point, seconds, count);
    }, py::arg("stage"), py::arg("seconds") = 5.0, py::arg("count") = 1)
    .def("release", &NativeProbe::release, py::call_guard<py::gil_scoped_release>())
    .def("counts", [](NativeProbe& self) {
      std::array<std::size_t, names.size() + 2> snapshot;
      {
        py::gil_scoped_release detached;
        const std::lock_guard lock(self.mutex);
        snapshot = self.counts;
      }
      py::dict out;
      for (std::size_t i = 0; i < names.size(); ++i) out[names[i]] = snapshot[i];
      out["pause_begin"] = snapshot[names.size()];
      out["pause_end"] = snapshot[names.size() + 1];
      return out;
    });
  m.def("_at_native_exit", [](const std::shared_ptr<FlowLedger>& ledger,
      const std::shared_ptr<ShortestPathCache>& cache, const std::shared_ptr<NativeProbe>& probe) {
    if ((!ledger && !cache) || (ledger && cache) || !probe) throw py::value_error("one native exit resource and a probe are required");
    auto state = native_call([&] {
      auto s = std::make_unique<ShutdownWaiter>();
      s->ledger = ledger;
      s->cache = cache;
      s->probe = probe;
      if (ledger) s->policy = std::make_unique<FlowPolicy>(*ledger, std::vector<NodeId>{0}, std::vector<NodeId>{1}, FlowPolicyConfig{});
      return s;
    });
    ShutdownWaiter* empty = nullptr;
    if (!shutdown_waiter.compare_exchange_strong(empty, state.get())) {
      native_call([&] { state.reset(); });
      throw py::value_error("a native exit waiter is already registered");
    }
    if (Py_AtExit(native_shutdown) != 0) {
      shutdown_waiter.store(nullptr);
      native_call([&] { state.reset(); });
      throw std::runtime_error("native exit registration failed");
    }
    (void)state.release();
  });
  m.def("_concurrency_model", [](const std::string& name) {
    py::dict out;
    if (name == "publication-locked" || name == "publication-unlocked") {
      const auto result = detail::publication_model(name == "publication-locked");
      out["states"] = result.states;
      out["mixed"] = result.mixed;
    } else {
      const auto result = detail::phase_model(name);
      out["states"] = result.states;
      out["partial"] = result.partial;
      out["python_under_lock"] = result.python_under_lock;
      out["late_hits"] = result.late_hits;
    }
    return out;
  }, py::arg("protocol"));
}
} // namespace netgraph::core::bindings
#endif
