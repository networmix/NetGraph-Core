#!/usr/bin/env bash
# Propagate configuration, build, and test failures; bound local CPU use.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."
jobs=${CMAKE_BUILD_PARALLEL_LEVEL:-4}
case "$jobs" in ''|*[!0-9]*|0) echo 'CMAKE_BUILD_PARALLEL_LEVEL must be positive.' >&2; exit 2 ;; esac
case "${1:-release}" in
    release) build_dir=build/cpp-tests; build_type=Release; sanitize=OFF ;;
    sanitize) build_dir=build/cpp-sanitize; build_type=Debug; sanitize=ON ;;
    *) echo 'Expected release or sanitize.' >&2; exit 2 ;;
esac
if [[ "$(uname -s)" == Darwin ]]; then
    export CC="$(xcrun --find clang)"
    export CXX="$(xcrun --find clang++)"
    export MACOSX_DEPLOYMENT_TARGET="${MACOSX_DEPLOYMENT_TARGET:-15.0}"
    leak_detection=0
else
    leak_detection=1
fi
cmake -S . -B "$build_dir" -DNETGRAPH_CORE_BUILD_TESTS=ON \
    -DNETGRAPH_CORE_SANITIZE="$sanitize" -DCMAKE_BUILD_TYPE="$build_type"
cmake --build "$build_dir" --config "$build_type" --parallel "$jobs"
if [[ "$sanitize" == ON ]]; then
    export ASAN_OPTIONS="${ASAN_OPTIONS:-detect_leaks=$leak_detection}"
fi
ctest --test-dir "$build_dir" --build-config "$build_type" \
    --output-on-failure --parallel "$jobs" --timeout 120 --no-tests=error
