// The Python/native boundary. Only native values may cross this scope: all Python
// conversion and ownership work is done by the caller before or after it.
#pragma once

#include <pybind11/pybind11.h>
#include <chrono>
#include <thread>
#include <utility>
#include <type_traits>
#include "netgraph/core/detail/concurrency.hpp"

namespace netgraph::core::bindings {
// The GIL released for native work: every native_call, and a FlowPolicy's destructor (its holder's
// deleter, which waits for the ledger gate); made with the GIL held. A thread whose call began
// while Python was initialized and that ends it after finalization began (a daemon thread) never
// retakes the GIL: CPython before 3.13.8 would end the thread with pthread_exit, which on glibc
// unwinds into this noexcept destructor and terminates the process. It parks forever instead, as
// later CPython does. A call that began during finalization (the finalizing thread destroying
// module globals) retakes the GIL as usual. Finalization starting between the check and the
// retake surfaces as a forced unwind from PyEval_RestoreThread, which parks the thread the same way.
class Detached {
 public:
  Detached() : initialized_(Py_IsInitialized() != 0), state_(PyEval_SaveThread()) {}
  ~Detached() {
    if (initialized_ && Py_IsInitialized() == 0) park();
    attach(state_);
  }
  Detached(const Detached&) = delete;
  Detached& operator=(const Detached&) = delete;

 private:
  [[noreturn]] static void park() noexcept {
    for (;;) std::this_thread::sleep_for(std::chrono::hours(1));
  }
  static void attach(PyThreadState* state) {
    try {
      PyEval_RestoreThread(state);
    } catch (...) {
      park();
    }
  }
  bool initialized_;
  PyThreadState* state_;
};

// Native argument work can precede computation. It uses the same detach scope,
// without announcing that the computation or its publication has begun.
template<bool Computation = true, class Body>
decltype(auto) native_call(Body&& body) {
  const Detached detached;
  if constexpr (Computation) detail::concurrent_event(detail::ConcurrentStage::NativeBegin, nullptr);
  else detail::concurrent_event(detail::ConcurrentStage::ArgumentBegin, nullptr);
  if constexpr (std::is_void_v<std::invoke_result_t<Body>>) {
    std::forward<Body>(body)();
    if constexpr (Computation) detail::publication_event(detail::ConcurrentStage::NativeComplete, nullptr);
  } else {
    auto result = std::forward<Body>(body)();
    if constexpr (Computation) detail::publication_event(detail::ConcurrentStage::NativeComplete, nullptr);
    return result;
  }
}
// Run only the graph-size part of a first pass detached, preserving the
// attached argument conversions and the surrounding check order.
inline void detached_argument_work(void (*body)(void*), void* context) {
  native_call<false>([&] { body(context); });
}
} // namespace netgraph::core::bindings
