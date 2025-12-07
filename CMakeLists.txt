cmake_minimum_required(VERSION 3.23)
project(netgraph_core LANGUAGES CXX)
# Silence pybind11 FindPython compatibility warning by selecting the modern mode
set(PYBIND11_FINDPYTHON NEW CACHE STRING "Use modern FindPython in pybind11")
# Optional C++ tests
option(NETGRAPH_CORE_BUILD_TESTS "Build C++ unit tests" OFF)
if(NETGRAPH_CORE_BUILD_TESTS)
  include(FetchContent)
  FetchContent_Declare(
    googletest
    URL https://github.com/google/googletest/archive/refs/tags/v1.14.0.zip
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE
  )
  FetchContent_MakeAvailable(googletest)
  enable_testing()
add_executable(netgraph_core_tests
  tests/cpp/smoke_test.cpp
  tests/cpp/flow_policy_tests.cpp
  tests/cpp/strict_multidigraph_tests.cpp
  tests/cpp/shortest_paths_tests.cpp
  tests/cpp/flow_state_tests.cpp
  tests/cpp/flow_graph_tests.cpp
  tests/cpp/max_flow_tests.cpp
  tests/cpp/k_shortest_paths_tests.cpp
  tests/cpp/masking_tests.cpp
)
  target_link_libraries(netgraph_core_tests PRIVATE netgraph_core GTest::gtest_main)
  target_include_directories(netgraph_core_tests PRIVATE tests/cpp)
  if(NETGRAPH_CORE_COVERAGE AND (CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang"))
    target_compile_options(netgraph_core_tests PRIVATE -O0 -g --coverage)
    target_link_options(netgraph_core_tests PRIVATE --coverage)
  endif()
  include(GoogleTest)
  gtest_discover_tests(netgraph_core_tests)
endif()

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

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

# Dependencies
find_package(pybind11 3 CONFIG QUIET)
if(NOT pybind11_FOUND)
  message(STATUS "pybind11 not found via config; using FetchContent")
  include(FetchContent)
  FetchContent_Declare(
    pybind11
    GIT_REPOSITORY https://github.com/pybind/pybind11.git
    GIT_TAG v3.0.0
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE
  )
  FetchContent_MakeAvailable(pybind11)
endif()


add_library(netgraph_core STATIC
  src/strict_multidigraph.cpp
  src/shortest_paths.cpp
  src/k_shortest_paths.cpp
  src/max_flow.cpp
  src/flow_state.cpp
  src/flow_graph.cpp
  src/flow_policy.cpp
  src/cpu_backend.cpp
)
target_include_directories(netgraph_core PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
)
target_compile_definitions(netgraph_core PRIVATE _USE_MATH_DEFINES)
## C++ standard is set globally above.
target_compile_features(netgraph_core PUBLIC cxx_std_20)

# Warnings
if(MSVC)
  target_compile_options(netgraph_core PRIVATE /W4 /permissive-)
else()
  target_compile_options(netgraph_core PRIVATE -Wall -Wextra -Wpedantic)
endif()

# Python extension
pybind11_add_module(_netgraph_core bindings/python/module.cpp)
target_link_libraries(_netgraph_core PRIVATE netgraph_core)
target_compile_features(_netgraph_core PRIVATE cxx_std_20)
if(APPLE)
  # Prefer maximum runtime compatibility for libc++ on macOS
  target_compile_definitions(netgraph_core PUBLIC _LIBCPP_DISABLE_AVAILABILITY=1 _LIBCPP_ABI_VERSION=1)
  target_compile_definitions(_netgraph_core PRIVATE _LIBCPP_DISABLE_AVAILABILITY=1 _LIBCPP_ABI_VERSION=1)
endif()

# =============================================================================
# Performance optimizations (disabled when coverage/sanitizers are enabled)
# =============================================================================

option(NETGRAPH_CORE_COVERAGE "Enable C++ coverage instrumentation" OFF)
option(NETGRAPH_CORE_SANITIZE "Enable Address/Undefined sanitizers" OFF)
option(NETGRAPH_CORE_NATIVE "Optimize for current CPU (not portable)" OFF)
option(NETGRAPH_CORE_FAST_MATH "Enable fast math (may affect precision)" OFF)
option(NETGRAPH_CORE_PGO_GENERATE "Instrument for PGO profile collection" OFF)
option(NETGRAPH_CORE_PGO_USE "Use PGO profile data (set NETGRAPH_CORE_PGO_DIR)" OFF)
set(NETGRAPH_CORE_PGO_DIR "${CMAKE_BINARY_DIR}/pgo" CACHE PATH "Directory for PGO profile data")

# Only apply optimizations when not in debug/instrumentation mode
if(NOT NETGRAPH_CORE_COVERAGE AND NOT NETGRAPH_CORE_SANITIZE)
  # LTO (Link Time Optimization) - enables cross-module inlining
  include(CheckIPOSupported)
  check_ipo_supported(RESULT ipo_supported OUTPUT ipo_error)
  if(ipo_supported)
    message(STATUS "Enabling LTO (Link Time Optimization)")
    set_target_properties(netgraph_core PROPERTIES INTERPROCEDURAL_OPTIMIZATION TRUE)
    set_target_properties(_netgraph_core PROPERTIES INTERPROCEDURAL_OPTIMIZATION TRUE)
  else()
    message(STATUS "LTO not supported: ${ipo_error}")
  endif()

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

    # Fast math - trades IEEE compliance for speed
    if(NETGRAPH_CORE_FAST_MATH)
      message(STATUS "Enabling fast math optimizations")
      target_compile_options(netgraph_core PRIVATE -ffast-math)
      target_compile_options(_netgraph_core PRIVATE -ffast-math)
    endif()

    # PGO (Profile-Guided Optimization) - uses runtime profile to guide inlining decisions
    # GCC: -fprofile-generate/-fprofile-use, Clang: -fprofile-instr-*
    if(NETGRAPH_CORE_PGO_GENERATE)
      message(STATUS "PGO: Instrumenting for profile generation -> ${NETGRAPH_CORE_PGO_DIR}")
      file(MAKE_DIRECTORY ${NETGRAPH_CORE_PGO_DIR})
      if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
        target_compile_options(netgraph_core PRIVATE -fprofile-generate=${NETGRAPH_CORE_PGO_DIR})
        target_compile_options(_netgraph_core PRIVATE -fprofile-generate=${NETGRAPH_CORE_PGO_DIR})
        target_link_options(netgraph_core PRIVATE -fprofile-generate=${NETGRAPH_CORE_PGO_DIR})
        target_link_options(_netgraph_core PRIVATE -fprofile-generate=${NETGRAPH_CORE_PGO_DIR})
      else()
        # Clang/AppleClang: use LLVM instrumentation PGO
        target_compile_options(netgraph_core PRIVATE -fprofile-instr-generate=${NETGRAPH_CORE_PGO_DIR}/default.profraw)
        target_compile_options(_netgraph_core PRIVATE -fprofile-instr-generate=${NETGRAPH_CORE_PGO_DIR}/default.profraw)
        target_link_options(netgraph_core PRIVATE -fprofile-instr-generate=${NETGRAPH_CORE_PGO_DIR}/default.profraw)
        target_link_options(_netgraph_core PRIVATE -fprofile-instr-generate=${NETGRAPH_CORE_PGO_DIR}/default.profraw)
      endif()
    elseif(NETGRAPH_CORE_PGO_USE)
      message(STATUS "PGO: Using profile data from ${NETGRAPH_CORE_PGO_DIR}")
      if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
        target_compile_options(netgraph_core PRIVATE -fprofile-use=${NETGRAPH_CORE_PGO_DIR} -fprofile-correction)
        target_compile_options(_netgraph_core PRIVATE -fprofile-use=${NETGRAPH_CORE_PGO_DIR} -fprofile-correction)
      else()
        # Clang/AppleClang: use merged .profdata file
        target_compile_options(netgraph_core PRIVATE -fprofile-instr-use=${NETGRAPH_CORE_PGO_DIR}/default.profdata)
        target_compile_options(_netgraph_core PRIVATE -fprofile-instr-use=${NETGRAPH_CORE_PGO_DIR}/default.profdata)
      endif()
    endif()
  elseif(MSVC)
    # /O2: full optimization, /GL: whole program optimization (MSVC's LTO)
    target_compile_options(netgraph_core PRIVATE /O2 /GL)
    target_compile_options(_netgraph_core PRIVATE /O2 /GL)
    if(NETGRAPH_CORE_FAST_MATH)
      target_compile_options(netgraph_core PRIVATE /fp:fast)
      target_compile_options(_netgraph_core PRIVATE /fp:fast)
    endif()
  endif()
endif()

# =============================================================================
# Debug/instrumentation options (mutually exclusive with optimizations above)
# =============================================================================

if(NETGRAPH_CORE_COVERAGE)
  message(STATUS "Enabling coverage instrumentation (disables optimizations)")
  if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    foreach(tgt netgraph_core _netgraph_core)
      target_compile_options(${tgt} PRIVATE -O0 -g --coverage)
      target_link_options(${tgt} PRIVATE --coverage)
    endforeach()
  else()
    message(WARNING "NETGRAPH_CORE_COVERAGE requires GCC or Clang")
  endif()
endif()

if(NETGRAPH_CORE_SANITIZE)
  message(STATUS "Enabling sanitizers (disables optimizations)")
  if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    set(SAN_FLAGS "-fsanitize=address,undefined -fno-omit-frame-pointer")
    foreach(tgt netgraph_core _netgraph_core)
      target_compile_options(${tgt} PRIVATE ${SAN_FLAGS})
      target_link_options(${tgt} PRIVATE ${SAN_FLAGS})
    endforeach()
  else()
    message(WARNING "NETGRAPH_CORE_SANITIZE requires GCC or Clang")
  endif()
endif()

# Install
install(TARGETS _netgraph_core
  LIBRARY DESTINATION .
  RUNTIME DESTINATION .
)
