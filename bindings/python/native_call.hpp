// The Python/native boundary. Only native values may cross this scope: all Python
// conversion and ownership work is done by the caller before or after it.
#pragma once

#include <pybind11/pybind11.h>
#include <utility>
#include <type_traits>
#include "netgraph/core/detail/concurrency.hpp"

namespace netgraph::core::bindings {
// Native argument work can precede computation. It uses the same detach scope,
// without announcing that the computation or its publication has begun.
template<bool Computation = true, class Body>
decltype(auto) native_call(Body&& body) {
  pybind11::gil_scoped_release detached;
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
