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

# Optional coverage instrumentation for GCC/Clang
option(NETGRAPH_CORE_COVERAGE "Enable C++ coverage instrumentation" OFF)
if(NETGRAPH_CORE_COVERAGE)
  message(STATUS "Enabling coverage instrumentation for C++ targets")
  if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    foreach(tgt netgraph_core _netgraph_core)
      # Use no optimizations and include debug info for accurate line mapping
      target_compile_options(${tgt} PRIVATE -O0 -g --coverage)
      # Link coverage runtime
      target_link_options(${tgt} PRIVATE --coverage)
    endforeach()
  else()
    message(WARNING "NETGRAPH_CORE_COVERAGE is set but compiler '${CMAKE_CXX_COMPILER_ID}' is not supported")
  endif()
endif()

# Optional sanitizers for debug testing
option(NETGRAPH_CORE_SANITIZE "Enable Address/Undefined sanitizers" OFF)
if(NETGRAPH_CORE_SANITIZE)
  if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    set(SAN_FLAGS "-fsanitize=address,undefined -fno-omit-frame-pointer")
    foreach(tgt netgraph_core _netgraph_core)
      target_compile_options(${tgt} PRIVATE ${SAN_FLAGS})
      target_link_options(${tgt} PRIVATE ${SAN_FLAGS})
    endforeach()
  else()
    message(WARNING "NETGRAPH_CORE_SANITIZE is set but compiler '${CMAKE_CXX_COMPILER_ID}' is not supported")
  endif()
endif()

# Install
install(TARGETS _netgraph_core
  LIBRARY DESTINATION .
  RUNTIME DESTINATION .
)
