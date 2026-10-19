#!/usr/bin/env bash
# The Python-side check entrypoints. Each one only reads the tree: formatting and fixes are
# `make format` and `make hooks`.
#
#   bash dev/run-checks.sh pytest on [PYTEST_ARGS...]
#       The Python suite against the netgraph_core that $PYTHON imports, which must carry the
#       private hooks (the editable development build): no test may skip. (make check-ci)
#   bash dev/run-checks.sh pytest off [PYTEST_ARGS...]
#       The Python suite against an installed, non-editable netgraph_core without the hooks (a
#       wheel, or a build from the sdist): exactly the tests that need a hook skip, listed in
#       dev/pytest_guard.py, which runs and checks both modes. (release CI)
#   bash dev/run-checks.sh cov
#       C++ and Python coverage (make cov): a coverage venv of its own (build/coverage-venv, pinned
#       by .github/ci-constraints.txt; the development venv is not touched), the extension and the
#       C++ tests built with coverage in their own build directories, every build, test and
#       report failure fails the run, and so does coverage that was not measured: no profile
#       data written in either build directory, no file of the C++ report under one of its
#       filters, or a report that counts no line. Reports in build/coverage/.
#
# PYTHON names the interpreter (default python3); CMAKE_BUILD_PARALLEL_LEVEL bounds the builds,
# ctest and gcovr (default 4; exported, so it also reaches the extension's build under pip); GCOV
# names the gcov command gcovr runs (default: xcrun llvm-cov gcov on macOS, llvm-cov gcov when CXX
# is clang, gcov otherwise).
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."
ROOT=$PWD
PYTHON=${PYTHON:-python3}
JOBS=${CMAKE_BUILD_PARALLEL_LEVEL:-4}
case "$JOBS" in ''|*[!0-9]*|0) echo 'CMAKE_BUILD_PARALLEL_LEVEL must be positive.' >&2; exit 2 ;; esac
export CMAKE_BUILD_PARALLEL_LEVEL=$JOBS

die() { echo "dev/run-checks.sh: $*" >&2; exit 1; }

if [[ "$(uname -s)" == Darwin ]]; then
    CC=$(xcrun --find clang)
    CXX=$(xcrun --find clang++)
    export CC CXX
    export MACOSX_DEPLOYMENT_TARGET="${MACOSX_DEPLOYMENT_TARGET:-15.0}"
fi

# run_pytest on|off [PYTEST_ARGS...]: the suite, then the check of what was imported and skipped
# (dev/pytest_guard.py, which the Windows CI jobs run directly; PYTEST_ADDOPTS and its reports stay
# the caller's).
run_pytest() {
    "$PYTHON" "$ROOT/dev/pytest_guard.py" "$@"
}

coverage() {
    local cvenv=build/coverage-venv ext=build/coverage-ext cpp=build/cpp-tests-cov out=build/coverage
    local gcov dir root filters=(include/ src/ bindings/)
    rm -rf "$cvenv" "$out"
    mkdir -p "$out"
    echo "== coverage venv $cvenv (pins: .github/ci-constraints.txt)"
    "$PYTHON" -m venv "$cvenv"
    local cpy=$ROOT/$cvenv/bin/python
    "$cpy" -m pip install --disable-pip-version-check -q -U -c .github/ci-constraints.txt pip
    "$cpy" -m pip install --disable-pip-version-check -q -c .github/ci-constraints.txt \
        scikit-build-core pybind11 cmake ninja numpy pytest pytest-cov gcovr
    export PATH="$ROOT/$cvenv/bin:$PATH"
    echo "== the extension with coverage (editable, hooks on, Debug) in $ext"
    "$cpy" -m pip install --disable-pip-version-check --no-build-isolation --no-deps -e . \
        -C build-dir="$ext" -C cmake.build-type=Debug \
        -C cmake.define.NETGRAPH_CORE_COVERAGE=ON -C cmake.define.NETGRAPH_CORE_PRIVATE_HOOKS=ON
    # gcovr reads every profile (.gcda file) it finds and reports no data, or old data, without
    # failing: the profiles of earlier runs go first, and each build directory must then hold
    # profiles that this run's tests wrote.
    find "$ext" -name '*.gcda' -delete
    [[ ! -d "$cpp" ]] || find "$cpp" -name '*.gcda' -delete
    echo "== Python tests with coverage"
    PYTHON=$cpy run_pytest on --cov=netgraph_core --cov-report=term-missing \
        --cov-report=xml:"$out/coverage-python.xml"
    echo "== C++ tests with coverage"
    bash dev/cpp-tests.sh coverage
    for dir in "$ext" "$cpp"; do
        [[ -n "$(find "$dir" -name '*.gcda' -print -quit)" ]] ||
            die "no C++ profile data in $dir: the tests that ran wrote no coverage there"
    done
    echo "== C++ coverage report (gcovr: the library, the bindings and the headers)"
    if [[ "$(uname -s)" == Darwin ]]; then
        gcov="xcrun llvm-cov gcov"
    elif [[ "${CXX:-}" == *clang* ]]; then
        gcov="llvm-cov gcov"
    else
        gcov=gcov
    fi
    gcov=${GCOV:-$gcov}
    "$cpy" -m gcovr --root . --gcov-executable "$gcov" -j "$JOBS" \
        "${filters[@]/#/--filter=}" --xml-pretty -o "$out/coverage-cpp.xml" "$ext" "$cpp"
    "$cpy" dev/coverage_summary.py "$out/coverage-python.xml" "$out/coverage-cpp.xml" \
        --html="$out/coverage-combined.html"
    # A filter that matches no file (a moved directory) would leave that part out silently.
    for root in "${filters[@]}"; do
        grep -q "filename=\"$root" "$out/coverage-cpp.xml" ||
            die "the C++ coverage report has no file under $root"
    done
    echo "Coverage ready in $out/: coverage-python.xml, coverage-cpp.xml, coverage-combined.html"
}

case "${1:-}" in
    pytest)
        case "${2:-}" in on|off) ;; *) die "usage: dev/run-checks.sh pytest on|off [PYTEST_ARGS...]" ;; esac
        hooks=$2; shift 2
        run_pytest "$hooks" "$@"
        ;;
    cov) coverage ;;
    *) echo 'Usage: bash dev/run-checks.sh pytest on|off [PYTEST_ARGS...] | cov' >&2; exit 2 ;;
esac
