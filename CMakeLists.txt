cmake_minimum_required(VERSION 3.23)
project(netgraph_core LANGUAGES CXX)

# =============================================================================
# Options. Every build states them on its command line (pyproject.toml for pip builds,
# dev/cpp-tests.sh for the C++ test builds), so a value cached by an earlier configure of the same
# build directory never decides a later build.
# =============================================================================

option(NETGRAPH_CORE_BUILD_TESTS "Build C++ unit tests" OFF)
# The private hooks (test controls, and diagnostics such as the placement trace that the tests
# compare with the reference under tests/reference/policy/; neither is public API):
# ON for the editable development build (pyproject.toml) and the hooks-on test suites
# (dev/cpp-tests.sh), OFF for a wheel and for the hooks-off test suite.
option(NETGRAPH_CORE_PRIVATE_HOOKS "Build the private test hooks into the library and the extension" OFF)
option(NETGRAPH_CORE_COVERAGE "Enable C++ coverage instrumentation" OFF)
option(NETGRAPH_CORE_SANITIZE "Enable Address/Undefined sanitizers" OFF)
option(NETGRAPH_CORE_TSAN "Enable ThreadSanitizer" OFF)
option(NETGRAPH_CORE_NATIVE "Optimize for current CPU (not portable)" OFF)
if(NETGRAPH_CORE_SANITIZE AND NETGRAPH_CORE_TSAN)
  # ThreadSanitizer cannot be combined with AddressSanitizer in one binary.
  message(FATAL_ERROR "NETGRAPH_CORE_SANITIZE and NETGRAPH_CORE_TSAN are mutually exclusive")
endif()

# Settings for every target: CMake applies them to targets created later, so they come first.
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
# Silence pybind11 FindPython compatibility warning by selecting the modern mode
set(PYBIND11_FINDPYTHON NEW CACHE STRING "Use modern FindPython in pybind11")

# macOS: enforce consistent deployment target and architecture for wheels/builds
if(APPLE)
  # Allow environment to override, otherwise choose a sensible default
  if(DEFINED ENV{MACOSX_DEPLOYMENT_TARGET})
    set(CMAKE_OSX_DEPLOYMENT_TARGET "$ENV{MACOSX_DEPLOYMENT_TARGET}" CACHE STRING "" FORCE)
  elseif(NOT DEFINED CMAKE_OSX_DEPLOYMENT_TARGET)
    set(CMAKE_OSX_DEPLOYMENT_TARGET "15.0" CACHE STRING "" FORCE)
  endif()
  if(NOT DEFINED CMAKE_OSX_ARCHITECTURES)
    # Default to arm64 on Apple Silicon hosts
    if(CMAKE_SYSTEM_PROCESSOR MATCHES "arm64|aarch64")
      set(CMAKE_OSX_ARCHITECTURES "arm64" CACHE STRING "" FORCE)
    endif()
  endif()
endif()

# Dependencies. pybind11 is exactly the version pinned in pyproject.toml's build requirements,
# whose licence notice the wheel ships (third_party/pybind11/LICENSE).
find_package(pybind11 3.1.0 EXACT CONFIG QUIET)
if(NOT pybind11_FOUND)
  message(STATUS "pybind11 3.1.0 not found via config; using FetchContent")
  include(FetchContent)
  FetchContent_Declare(
    pybind11
    URL https://github.com/pybind/pybind11/archive/refs/tags/v3.1.0.tar.gz
    URL_HASH SHA256=ef712655692a2e9bf7bb7874c022564a45f91d847ddee987e720cd9e28849665
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE
  )
  FetchContent_MakeAvailable(pybind11)
endif()

# =============================================================================
# Targets
# =============================================================================

add_library(netgraph_core STATIC
  src/concurrency.cpp
  src/graph.cpp
  src/eligibility.cpp
  src/min_cost_flow.cpp
  src/route_structures.cpp
  src/search.cpp
  src/route_kernels.cpp
  src/flow_ledger.cpp
  src/shortest_path_cache.cpp
  src/flow_policy.cpp
  src/evaluators.cpp
  src/forwarding_program.cpp
  src/forwarding.cpp
  src/distributed_te.cpp
  src/centralized_te_tunnels.cpp
)
target_include_directories(netgraph_core PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
)
# PUBLIC: the library, the extension and the tests see one value (0 or 1).
target_compile_definitions(netgraph_core PUBLIC NETGRAPH_CORE_PRIVATE_HOOKS=$<BOOL:${NETGRAPH_CORE_PRIVATE_HOOKS}>)
target_compile_features(netgraph_core PUBLIC cxx_std_20)

# Python extension: the module entry and its registration files (bindings/python/), in the
# order the entry calls them; the shared helpers are headers there.
pybind11_add_module(_netgraph_core
  bindings/python/module.cpp
  bindings/python/options.cpp
  bindings/python/units.cpp
  bindings/python/eligibility.cpp
  bindings/python/graph.cpp
  bindings/python/flow_solver_hook.cpp
  bindings/python/predecessor_dag.cpp
  bindings/python/path.cpp
  bindings/python/flow_record.cpp
  bindings/python/flow_ledger.cpp
  bindings/python/forwarding_program.cpp
  bindings/python/routing_models.cpp
  bindings/python/flow_policy_config.cpp
  bindings/python/placement_result.cpp
  bindings/python/shortest_path_cache.cpp
  bindings/python/flow_policy.cpp
  bindings/python/apply_actions.cpp
  bindings/python/evaluators.cpp
  bindings/python/path_enumeration.cpp
  bindings/python/search.cpp
  bindings/python/hooks.cpp
  bindings/python/concurrency_hooks.cpp
)
# The library's PUBLIC usage requirements (C++20, the hooks flag, the libc++ settings below) reach
# the extension through this link.
target_link_libraries(_netgraph_core PRIVATE netgraph_core)
if(APPLE)
  # Prefer maximum runtime compatibility for libc++ on macOS
  target_compile_definitions(netgraph_core PUBLIC _LIBCPP_DISABLE_AVAILABILITY=1 _LIBCPP_ABI_VERSION=1)
endif()

# The targets that get the project's warnings and either the link-time optimisation (where the
# compiler supports it) or, in a coverage or sanitizer build, the instrumentation: the library,
# the extension and, in a test build, the two test executables (added below). The header
# self-containment check below is not one of them. The extra code-generation flags go to the
# library and the extension only.
set(NETGRAPH_CORE_PROJECT_TARGETS netgraph_core _netgraph_core)

if(NETGRAPH_CORE_BUILD_TESTS)
  include(FetchContent)
  FetchContent_Declare(
    googletest
    URL https://github.com/google/googletest/archive/refs/tags/v1.14.0.zip
    URL_HASH SHA256=1f357c27ca988c3f7c6b4bf68a9395005ac6761f034046e9dde0896e3aba00e4
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE
  )
  # MSVC: link googletest against the same (shared) C runtime as the project.
  set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
  FetchContent_MakeAvailable(googletest)
  # googletest itself stays at C++17 (under C++20 its char8_t printer warns with Clang 21); the
  # project's tests below are C++20 like the library.
  set_target_properties(gtest gtest_main gmock gmock_main PROPERTIES CXX_STANDARD 17)
  enable_testing()
  add_executable(netgraph_core_tests
    tests/cpp/graph_tests.cpp
    tests/cpp/graph_layout_tests.cpp
    tests/cpp/graph_concurrency_tests.cpp
    tests/cpp/concurrency_tests.cpp
    tests/cpp/min_cost_flow_tests.cpp
    tests/cpp/min_cost_flow_deletion_tests.cpp
    tests/cpp/eligibility_tests.cpp
    tests/cpp/route_structures_tests.cpp
    tests/cpp/search_tests.cpp
    tests/cpp/ksp_tests.cpp
    tests/cpp/shortest_path_cache_tests.cpp
    tests/cpp/route_kernels_tests.cpp
    tests/cpp/shortest_paths_creation_tests.cpp
    tests/cpp/flow_ledger_tests.cpp
    tests/cpp/flow_policy_fixture_tests.cpp
    tests/cpp/flow_policy_tests.cpp
    tests/cpp/epoch_static_tests.cpp
    tests/cpp/weighted_source_tests.cpp
    tests/cpp/zero_capacity_tests.cpp
    tests/cpp/pooled_source_allocation_tests.cpp
    tests/cpp/program_tests.cpp
    tests/cpp/joint_forwarding_tests.cpp
    tests/cpp/apply_actions_tests.cpp
    tests/cpp/distributed_te_tests.cpp
    tests/cpp/centralized_te_tunnels_tests.cpp
    tests/cpp/contract_examples_tests.cpp
    tests/cpp/public_names_tests.cpp
  )
  target_compile_definitions(netgraph_core_tests PRIVATE NETGRAPH_CORE_TEST_DATA="${CMAKE_CURRENT_SOURCE_DIR}/tests/data")
  # Allocation-failure tests replace the global operator new: a separate executable.
  add_executable(netgraph_core_fault_tests
    tests/cpp/fault_injection.cpp
    tests/cpp/flow_ledger_fault_tests.cpp
    tests/cpp/flow_policy_fault_tests.cpp
    tests/cpp/apply_actions_fault_tests.cpp
    tests/cpp/min_cost_flow_deletion_fault_tests.cpp
    tests/cpp/search_fault_tests.cpp
    tests/cpp/graph_fault_tests.cpp
  )
  foreach(tgt netgraph_core_tests netgraph_core_fault_tests)
    target_link_libraries(${tgt} PRIVATE netgraph_core GTest::gtest_main)
    target_include_directories(${tgt} PRIVATE tests/cpp)
  endforeach()
  # Every public header compiles on its own: one generated source file per header, holding only
  # its include, built with the test executables. No public header names a private detail/
  # header in an #include (the configure stops on one).
  file(GLOB NETGRAPH_CORE_PUBLIC_HEADERS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}/include
       CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/include/netgraph/core/*.hpp)
  set(NETGRAPH_CORE_HEADER_CHECKS)
  foreach(header ${NETGRAPH_CORE_PUBLIC_HEADERS})
    file(STRINGS ${CMAKE_CURRENT_SOURCE_DIR}/include/${header} private_includes
         REGEX "^[ \t]*#[ \t]*include[ \t]*[\"<][^\">]*detail/")
    if(private_includes)
      message(FATAL_ERROR "public header ${header} includes a private detail/ header: ${private_includes}")
    endif()
    string(MAKE_C_IDENTIFIER ${header} id)
    file(GENERATE OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/header_checks/${id}.cpp CONTENT "#include \"${header}\"\n")
    list(APPEND NETGRAPH_CORE_HEADER_CHECKS ${CMAKE_CURRENT_BINARY_DIR}/header_checks/${id}.cpp)
  endforeach()
  add_library(netgraph_core_header_checks OBJECT ${NETGRAPH_CORE_HEADER_CHECKS})
  target_link_libraries(netgraph_core_header_checks PRIVATE netgraph_core)
  add_dependencies(netgraph_core_tests netgraph_core_header_checks)
  list(APPEND NETGRAPH_CORE_PROJECT_TARGETS netgraph_core_tests netgraph_core_fault_tests)
  include(GoogleTest)
  gtest_discover_tests(netgraph_core_tests)
  # These sweeps make about 19,000 and 13,000 injection attempts, each with a clean retry;
  # evaluate_edge_failures also creates and joins worker threads. Debug sanitizer runs on small runners
  # can exceed the runner's 120-second default. Keep every attempt and allow ten minutes only
  # for these two tests; all other tests retain the runner's timeout.
  set(NETGRAPH_CORE_LONG_FAULT_TESTS
      "PolicyFaults.EvaluatorsLeaveNoPoisonedScratch:PolicyFaults.SensitivityWorkersFailWithoutTerminating")
  gtest_discover_tests(netgraph_core_fault_tests TEST_FILTER "-${NETGRAPH_CORE_LONG_FAULT_TESTS}")
  gtest_discover_tests(netgraph_core_fault_tests TEST_FILTER "${NETGRAPH_CORE_LONG_FAULT_TESTS}"
                      PROPERTIES TIMEOUT 600)
endif()

# Warnings: the project's set on the library, the extension and the test executables (the list
# above); third-party code keeps its own.
if(MSVC)
  set(NETGRAPH_CORE_WARNINGS /W4 /permissive-)
else()
  set(NETGRAPH_CORE_WARNINGS -Wall -Wextra -Wpedantic)
endif()
foreach(tgt ${NETGRAPH_CORE_PROJECT_TARGETS})
  target_compile_options(${tgt} PRIVATE ${NETGRAPH_CORE_WARNINGS})
endforeach()

# =============================================================================
# Performance optimizations (disabled when coverage/sanitizers are enabled)
# =============================================================================

# Apply the extra optimization settings when coverage and sanitizers are off (whatever the build
# type: a Debug build gets them too).
if(NOT NETGRAPH_CORE_COVERAGE AND NOT NETGRAPH_CORE_SANITIZE AND NOT NETGRAPH_CORE_TSAN)
  # LTO (Link Time Optimization) - enables cross-module inlining. LANGUAGES CXX: googletest
  # enables C as well, and the check would otherwise require LTO support for C too and turn LTO
  # off in every test build. The test executables are LTO-linked like the extension, so the C++
  # tests exercise the code generation that ships.
  include(CheckIPOSupported)
  check_ipo_supported(RESULT ipo_supported OUTPUT ipo_error LANGUAGES CXX)
  if(ipo_supported)
    message(STATUS "Enabling LTO (Link Time Optimization)")
    set_target_properties(${NETGRAPH_CORE_PROJECT_TARGETS} PROPERTIES INTERPROCEDURAL_OPTIMIZATION TRUE)
  else()
    message(STATUS "LTO not supported: ${ipo_error}")
  endif()

  # MSVC gets no extra flags: its build type's optimisation (/O2 in Release) and CMake's LTO (/GL
  # with the matching /LTCG link option) suffice.
  if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    # GCC-only: allow aggressive inlining by promising symbols won't be interposed
    if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
      target_compile_options(netgraph_core PRIVATE -fno-semantic-interposition)
      target_compile_options(_netgraph_core PRIVATE -fno-semantic-interposition)
    endif()

    # Safe floating-point optimizations (we don't use errno, FP exceptions, or change rounding mode)
    target_compile_options(netgraph_core PRIVATE -fno-math-errno -fno-trapping-math)
    target_compile_options(_netgraph_core PRIVATE -fno-math-errno -fno-trapping-math)

    # Aggressive loop optimizations
    target_compile_options(netgraph_core PRIVATE -funroll-loops)
    target_compile_options(_netgraph_core PRIVATE -funroll-loops)

    # Native CPU instructions (AVX2, etc.) - only for local builds
    if(NETGRAPH_CORE_NATIVE)
      message(STATUS "Enabling native CPU optimizations (-march=native)")
      target_compile_options(netgraph_core PRIVATE -march=native)
      target_compile_options(_netgraph_core PRIVATE -march=native)
    endif()
  endif()
endif()

# =============================================================================
# Coverage and sanitizer options (mutually exclusive with the extra optimizations above). A request
# the compiler cannot honour stops the configure: a build that ran its tests without the requested
# checker would otherwise report success.
# =============================================================================

if(NETGRAPH_CORE_COVERAGE)
  message(STATUS "Enabling coverage instrumentation (disables optimizations)")
  if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    # Atomic counters (GCC and Clang): the evaluators' worker threads corrupt the default
    # non-atomic profile counters under concurrent updates (gcov then emits impossible values
    # such as "branch 4 taken -1").
    set(COV_FLAGS -O0 -g --coverage -fprofile-update=atomic)
    foreach(tgt ${NETGRAPH_CORE_PROJECT_TARGETS})
      target_compile_options(${tgt} PRIVATE ${COV_FLAGS})
      target_link_options(${tgt} PRIVATE --coverage)
    endforeach()
  else()
    message(FATAL_ERROR "NETGRAPH_CORE_COVERAGE requires GCC or Clang")
  endif()
endif()

if(NETGRAPH_CORE_SANITIZE)
  message(STATUS "Enabling sanitizers (without the extra optimizations)")
  if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    # Must be a CMake list: a single quoted string reaches the compiler as one
    # argument and clang rejects "-fsanitize=address,undefined -fno-omit-frame-pointer".
    #
    # -fno-sanitize-recover=undefined makes UBSan findings FATAL. Without it UBSan
    # prints a "runtime error:" diagnostic and lets the program continue, so the test
    # still exits 0 and ctest reports success -- undefined behaviour would be detected
    # and then silently ignored, which is indistinguishable from not checking at all.
    set(SAN_FLAGS -fsanitize=address,undefined -fno-omit-frame-pointer
                  -fno-sanitize-recover=undefined)
    # The test executables link the sanitized static library, so they are built and linked with
    # the same flags to pull in the sanitizer runtime.
    foreach(tgt ${NETGRAPH_CORE_PROJECT_TARGETS})
      target_compile_options(${tgt} PRIVATE ${SAN_FLAGS})
      target_link_options(${tgt} PRIVATE ${SAN_FLAGS})
    endforeach()
  else()
    message(FATAL_ERROR "NETGRAPH_CORE_SANITIZE requires GCC or Clang")
  endif()
endif()

if(NETGRAPH_CORE_TSAN)
  message(STATUS "Enabling ThreadSanitizer (without the extra optimizations)")
  if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    # The data-race check of the evaluators' worker threads and the shared ShortestPathCache. A CMake
    # list, as SAN_FLAGS above. A report makes the process exit non-zero (TSan's exit code 66,
    # or an abort under the halt_on_error=1 of dev/cpp-tests.sh), so ctest fails on any race. -O1
    # as the TSan documentation advises (the Debug build type adds -g); the instrumented
    # suites at -O0 would be several times slower for no extra coverage.
    set(TSAN_FLAGS -fsanitize=thread -fno-omit-frame-pointer -O1)
    # The test executables link the instrumented static library: same compile and link flags, so
    # they carry the TSan runtime and their own accesses are instrumented too.
    foreach(tgt ${NETGRAPH_CORE_PROJECT_TARGETS})
      target_compile_options(${tgt} PRIVATE ${TSAN_FLAGS})
      target_link_options(${tgt} PRIVATE ${TSAN_FLAGS})
    endforeach()
  else()
    message(FATAL_ERROR "NETGRAPH_CORE_TSAN requires GCC or Clang")
  endif()
endif()

# Install
install(TARGETS _netgraph_core
  LIBRARY DESTINATION .
  RUNTIME DESTINATION .
)
