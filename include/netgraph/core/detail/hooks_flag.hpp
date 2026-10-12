/* The build flag of the private hooks (internal): NETGRAPH_CORE_PRIVATE_HOOKS, 1 when the test
   controls of private_hooks.hpp are compiled in. CMake defines it for the library and, through
   the library target, for the extension and the tests (ON for the editable development build
   and the test suites, OFF for an installed wheel); a native caller compiling without that
   definition gets OFF. Code that only asks whether the hooks exist includes this header, not the
   hook declarations, which reach the policy and the trace. */
#pragma once

#ifndef NETGRAPH_CORE_PRIVATE_HOOKS
#define NETGRAPH_CORE_PRIVATE_HOOKS 0
#endif
