#!/usr/bin/env bash
# Propagate configuration, build, and test failures; bound local CPU use.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."
jobs=${CMAKE_BUILD_PARALLEL_LEVEL:-4}
case "$jobs" in ''|*[!0-9]*|0) echo 'CMAKE_BUILD_PARALLEL_LEVEL must be positive.' >&2; exit 2 ;; esac
tsan=OFF
case "${1:-release}" in
    release) build_dir=build/cpp-tests; build_type=Release; sanitize=OFF ;;
    sanitize) build_dir=build/cpp-sanitize; build_type=Debug; sanitize=ON ;;
    # ThreadSanitizer (data races of the worker threads and the shared DagCache), its own build
    # directory; the longest instrumented test takes about 30 s, inside the common timeout.
    tsan) build_dir=build/cpp-tsan; build_type=Debug; sanitize=OFF; tsan=ON ;;
    *) echo 'Expected release, sanitize or tsan.' >&2; exit 2 ;;
esac
if [[ "$(uname -s)" == Darwin ]]; then
    export CC="$(xcrun --find clang)"
    export CXX="$(xcrun --find clang++)"
    export MACOSX_DEPLOYMENT_TARGET="${MACOSX_DEPLOYMENT_TARGET:-15.0}"
    leak_detection=0
else
    leak_detection=1
fi
cmake -S . -B "$build_dir" -DNETGRAPH_CORE_BUILD_TESTS=ON -DNETGRAPH_CORE_PRIVATE_HOOKS=ON \
    -DNETGRAPH_CORE_SANITIZE="$sanitize" -DNETGRAPH_CORE_TSAN="$tsan" -DCMAKE_BUILD_TYPE="$build_type"
cmake --build "$build_dir" --config "$build_type" --parallel "$jobs"
if [[ "$sanitize" == ON ]]; then
    export ASAN_OPTIONS="${ASAN_OPTIONS:-detect_leaks=$leak_detection}"
fi
if [[ "$tsan" == ON ]]; then
    # A report aborts the test (halt_on_error=1; TSan's own default is exit code 66 at exit):
    # ctest fails on any race.
    export TSAN_OPTIONS="${TSAN_OPTIONS:-halt_on_error=1 second_deadlock_stack=1}"
fi
ctest --test-dir "$build_dir" --build-config "$build_type" \
    --output-on-failure --parallel "$jobs" --timeout 120 --no-tests=error
