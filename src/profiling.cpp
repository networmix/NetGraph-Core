/*
  Profiling implementation - singleton definitions.

  These are defined in a .cpp file (not inline in the header) to avoid ODR
  violations when linking a static library into a shared library (Python module).
  With inline definitions, separate static variable instances can exist in each
  translation unit, causing profiling data to be recorded in one instance while
  being read from another.
*/
#include "netgraph/core/profiling.hpp"

#include <cstdlib>

namespace netgraph::core {

bool profiling_enabled() noexcept {
    static const bool enabled = [] {
        const char* env = std::getenv("NGRAPH_CORE_PROFILE");
        return env && env[0] == '1';
    }();
    return enabled;
}

ProfilingStats& ProfilingStats::instance() {
    static ProfilingStats inst;
    return inst;
}

} // namespace netgraph::core
