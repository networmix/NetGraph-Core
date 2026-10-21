// Throwing NumPy conversions, with a private allocation-test window around the cast itself.
#pragma once

#include <pybind11/numpy.h>

#include "netgraph/core/detail/hooks_flag.hpp"

namespace netgraph::core::bindings {

#if NETGRAPH_CORE_PRIVATE_HOOKS
inline thread_local bool in_numpy_conversion = false;

class NumpyConversionScope {
 public:
  NumpyConversionScope() noexcept : previous_(in_numpy_conversion) { in_numpy_conversion = true; }
  ~NumpyConversionScope() { in_numpy_conversion = previous_; }
  NumpyConversionScope(const NumpyConversionScope&) = delete;
  NumpyConversionScope& operator=(const NumpyConversionScope&) = delete;

 private:
  bool previous_;
};
#endif

template <class T, int Flags = pybind11::array::forcecast>
pybind11::array_t<T, Flags> numpy_conversion(const pybind11::object& input) {
#if NETGRAPH_CORE_PRIVATE_HOOKS
  const NumpyConversionScope scope;
#endif
  return pybind11::array_t<T, Flags>(input);
}

} // namespace netgraph::core::bindings
