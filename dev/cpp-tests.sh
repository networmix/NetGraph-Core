#!/usr/bin/env bash
# The C++ test runner (make cpp-test, cpp-test-off, sanitize-test, tsan-test, check-ci, cov).
#
#   bash dev/cpp-tests.sh MODE
#
#   MODE      private hooks  build type             build directory
#   release   ON             Release (LTO)          build/cpp-tests
#   off       OFF            Release (LTO)          build/cpp-tests-off   the shipped configuration
#   sanitize  ON             Debug, ASan + UBSan    build/cpp-sanitize
#   tsan      ON             Debug, TSan            build/cpp-tsan
#   coverage  ON             Debug, --coverage      build/cpp-tests-cov   used by make cov
#
# Each mode has its own build directory and states every project option on the command line, so a
# value cached by an earlier configure never decides this one. One generator (Ninja), parallelism
# bounded by CMAKE_BUILD_PARALLEL_LEVEL (default 4), ctest fails when it finds no test, and every
# configuration, build or test failure fails the script. The skipped tests are checked: none in a
# hooks-on mode, exactly EXPECTED_OFF_SKIPS (the tests that need a hook) in mode off, so a test that
# silently became hooks-only is caught.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."

# The C++ tests that need the private hooks and skip without them (GTEST_SKIP). Inside tests that
# still run, three checks are weaker without them: the helper check_deletion
# (min_cost_flow_deletion_tests.cpp) checks each re-solve's final result, but observes each repair step
# only through the probe hook (see that file's header); PolicyFixtures.EveryFixtureReplaysIdentically
# compares no trace, and PolicyFixtures.ACorruptedExpectationFails corrupts only what that replay
# compares (see flow_policy_fixture_tests.cpp's header).
EXPECTED_OFF_SKIPS="
FlowCoreDeletion.CheckpointsEqualTheIterativeUpdate
FlowCoreDeletion.RerunningTheFailedSearchAddsExactlyOneSearch
ForwardingSize.TheCallChecksEachCountInOrder
GraphLayout.CombinedKeyPassesOnlyBelowNTimesKEqualsE
GraphLayout.HashAndFallbackBoundaryAt65536DistinctCosts
GraphLayout.ProbeBudgetFallbackInCollectionAndGrowth
GraphLayout.ProbeBudgetIsCumulativeAcrossTheLargeTableGrowths
GraphLayout.RadixFallbackDigitsOnBothSidesOfTheCombinedKeyBoundary
GraphLayout.RunsOfEqualCostsSkipTheLookup
Search.QueueDispatchAtTheWidthLimit
"

die() { echo "dev/cpp-tests.sh: $*" >&2; exit 1; }

jobs=${CMAKE_BUILD_PARALLEL_LEVEL:-4}
case "$jobs" in ''|*[!0-9]*|0) echo 'CMAKE_BUILD_PARALLEL_LEVEL must be positive.' >&2; exit 2 ;; esac
mode=${1:-}
hooks=ON; sanitize=OFF; tsan=OFF; coverage=OFF
case "$mode" in
    release) build_dir=build/cpp-tests; build_type=Release ;;
    off) build_dir=build/cpp-tests-off; build_type=Release; hooks=OFF ;;
    sanitize) build_dir=build/cpp-sanitize; build_type=Debug; sanitize=ON ;;
    # ThreadSanitizer (data races of the worker threads and the shared ShortestPathCache); the longest
    # instrumented test takes about 30 s, inside the common timeout.
    tsan) build_dir=build/cpp-tsan; build_type=Debug; tsan=ON ;;
    coverage) build_dir=build/cpp-tests-cov; build_type=Debug; coverage=ON ;;
    *) echo 'Usage: bash dev/cpp-tests.sh release|off|sanitize|tsan|coverage' >&2; exit 2 ;;
esac

command -v cmake >/dev/null 2>&1 || die "cmake not found on PATH (make dev installs it into venv/bin)"
command -v ctest >/dev/null 2>&1 || die "ctest not found on PATH (make dev installs it into venv/bin)"
command -v ninja >/dev/null 2>&1 || die "ninja not found on PATH (make dev installs it into venv/bin)"
cache="$build_dir/CMakeCache.txt"
if [[ -f "$cache" ]] && ! grep -qx 'CMAKE_GENERATOR:INTERNAL=Ninja' "$cache"; then
    die "$build_dir was configured by another generator ($(sed -n 's/^CMAKE_GENERATOR:INTERNAL=//p' "$cache")); remove that directory and rerun"
fi

if [[ "$(uname -s)" == Darwin ]]; then
    CC=$(xcrun --find clang)
    CXX=$(xcrun --find clang++)
    export CC CXX
    export MACOSX_DEPLOYMENT_TARGET="${MACOSX_DEPLOYMENT_TARGET:-15.0}"
    leak_detection=0
else
    leak_detection=1
fi
if [[ "$sanitize" == ON ]]; then
    # A sanitizer report must fail, including during test discovery at build time. Append
    # the failure settings after the caller's options; keep the platform's leak policy.
    export ASAN_OPTIONS="${ASAN_OPTIONS:-detect_leaks=$leak_detection}:halt_on_error=1:exitcode=1"
    export UBSAN_OPTIONS="${UBSAN_OPTIONS:+$UBSAN_OPTIONS:}halt_on_error=1:exitcode=1"
fi
if [[ "$tsan" == ON ]]; then
    export TSAN_OPTIONS="${TSAN_OPTIONS:-second_deadlock_stack=1}:halt_on_error=1:exitcode=1"
fi
cmake -S . -B "$build_dir" -G Ninja -DCMAKE_BUILD_TYPE="$build_type" \
    -DNETGRAPH_CORE_BUILD_TESTS=ON -DNETGRAPH_CORE_PRIVATE_HOOKS="$hooks" \
    -DNETGRAPH_CORE_SANITIZE="$sanitize" -DNETGRAPH_CORE_TSAN="$tsan" \
    -DNETGRAPH_CORE_COVERAGE="$coverage" -DNETGRAPH_CORE_NATIVE=OFF
cmake --build "$build_dir" --parallel "$jobs"
if [[ "$coverage" == ON ]]; then
    # Test discovery runs each binary at build time and writes profiles that are no test's coverage.
    find "$build_dir" -name '*.gcda' -delete
fi
log="$build_dir/ctest.log"
ctest --test-dir "$build_dir" --output-on-failure --parallel "$jobs" --timeout 120 --no-tests=error \
    2>&1 | tee "$log"

# ctest reports a GTEST_SKIP as "***Skipped" (gtest_discover_tests' skip pattern).
skipped=$(sed -n 's/^ *[0-9]*\/[0-9]* Test *#[0-9]*: \([^ ]*\) .*\*\*\*Skipped.*/\1/p' "$log" | LC_ALL=C sort)
if [[ "$hooks" == OFF ]]; then
    expected=$(echo "$EXPECTED_OFF_SKIPS" | sed '/^$/d' | LC_ALL=C sort)
    if [[ "$skipped" != "$expected" ]]; then
        diff <(echo "$expected") <(echo "$skipped") >&2 || true
        die "the skipped tests differ from the expected hook-only tests (< expected, > skipped)"
    fi
    echo "hooks off: exactly the $(echo "$expected" | wc -l | tr -d ' ') expected hook-only tests skipped"
elif [[ -n "$skipped" ]]; then
    die "tests skipped in a hooks-on build: $(echo "$skipped" | tr '\n' ' ')"
fi
