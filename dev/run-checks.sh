#!/usr/bin/env bash
# The Python-side check entrypoints. Each one only reads the tree: formatting and fixes are
# `make format` and `make hooks`.
#
#   bash dev/run-checks.sh pytest on [PYTEST_ARGS...]
#       The Python suite against the netgraph_core that $PYTHON imports, which must carry the
#       private hooks (the editable development build): no test may skip. (make check-ci)
#   bash dev/run-checks.sh pytest off [PYTEST_ARGS...]
#       The Python suite against an installed, non-editable netgraph_core without the hooks (a
#       wheel, or a build from the sdist): exactly the tests that need a hook skip. (release CI)
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

# The Python tests that need the private hooks and skip without them, by pytest node id, one per
# line (an id can contain spaces); every replay of the fixtures produced by the reference under
# tests/reference/policy/ (the parametrizations of test_fixture) skips as well.
EXPECTED_OFF_SKIPS="
tests/py/test_apportionment.py::test_apportionment_entries
tests/py/test_dag_cache.py::test_frozen_statistics_survive_python_allocation_failures
tests/py/test_dag_cache.py::test_the_frozen_sweep_rejects_a_wrong_call
tests/py/test_flow_ledger.py::test_a_flow_id_read_passes_an_allocation_failure_through
tests/py/test_flow_policy.py::test_evaluator_results_survive_python_allocation_failures
tests/py/test_flow_policy.py::test_owner_hook_is_the_owner_of_the_records
tests/py/test_flow_policy.py::test_python_publication_is_all_or_nothing[False-costonly]
tests/py/test_flow_policy.py::test_python_publication_is_all_or_nothing[False-global]
tests/py/test_flow_policy.py::test_python_publication_is_all_or_nothing[True-costonly]
tests/py/test_flow_policy.py::test_python_publication_is_all_or_nothing[True-global]
tests/py/test_flow_policy_fixtures.py::test_a_corrupted_expectation_fails[ledger]
tests/py/test_flow_policy_fixtures.py::test_a_corrupted_expectation_fails[policy]
tests/py/test_flow_policy_fixtures.py::test_a_corrupted_expectation_fails[result]
tests/py/test_flow_policy_fixtures.py::test_a_corrupted_expectation_fails[trace]
tests/py/test_flow_policy_fixtures.py::test_every_other_replay_clause_reports_its_operation
tests/py/test_flow_policy_fixtures.py::test_an_expected_value_of_another_type_fails[bool as integer]
tests/py/test_flow_policy_fixtures.py::test_an_expected_value_of_another_type_fails[integer as float]
tests/py/test_ksp.py::test_python_result_lists_are_all_or_nothing
tests/py/test_native_allocation.py::test_native_allocation_failure_raises_and_changes_nothing[FlowLedger.reset]
tests/py/test_native_allocation.py::test_native_allocation_failure_raises_and_changes_nothing[batch_capacity.global]
tests/py/test_native_allocation.py::test_native_allocation_failure_raises_and_changes_nothing[capacity.frozen_cache]
tests/py/test_native_allocation.py::test_native_allocation_failure_raises_and_changes_nothing[capacity.global]
tests/py/test_native_allocation.py::test_native_allocation_failure_raises_and_changes_nothing[ksp]
tests/py/test_native_allocation.py::test_native_allocation_failure_raises_and_changes_nothing[place.costonly]
tests/py/test_native_allocation.py::test_native_allocation_failure_raises_and_changes_nothing[place.global]
tests/py/test_native_allocation.py::test_native_allocation_failure_raises_and_changes_nothing[place_traced.global]
tests/py/test_native_allocation.py::test_native_allocation_failure_raises_and_changes_nothing[resolve_to_paths]
tests/py/test_native_allocation.py::test_native_allocation_failure_raises_and_changes_nothing[sensitivity.global]
tests/py/test_native_allocation.py::test_native_allocation_failure_raises_and_changes_nothing[spf.full]
tests/py/test_native_allocation.py::test_native_allocation_failure_raises_and_changes_nothing[spf.single_path]
tests/py/test_native_allocation.py::test_native_allocation_failure_raises_and_changes_nothing[spf.targets_residual]
tests/py/test_native_allocation.py::test_the_row_checks_reject_a_wrong_call
tests/py/test_search.py::test_python_result_construction_is_all_or_nothing
tests/py/test_search.py::test_the_python_allocation_sweep_rejects_a_wrong_call
"
EXPECTED_OFF_SKIPPED_PREFIX='tests/py/test_flow_policy_fixtures.py::test_fixture['


die() { echo "dev/run-checks.sh: $*" >&2; exit 1; }
TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT

if [[ "$(uname -s)" == Darwin ]]; then
    CC=$(xcrun --find clang)
    CXX=$(xcrun --find clang++)
    export CC CXX
    export MACOSX_DEPLOYMENT_TARGET="${MACOSX_DEPLOYMENT_TARGET:-15.0}"
fi

# run_pytest on|off [PYTEST_ARGS...]: the suite, then the check of what was imported and skipped.
# A small plugin records the skipped node ids and the imported modules (PYTEST_ADDOPTS and its
# reports stay the caller's).
run_pytest() {
    local hooks=$1 tmp=$TMP/pytest; shift
    rm -rf "$tmp"; mkdir -p "$tmp"
    cat > "$tmp/ngc_skip_record.py" <<'PY'
import json
import os
import sys

_seen, _skipped = set(), {}


def pytest_runtest_logreport(report):
    _seen.add(report.nodeid)
    if report.skipped:
        _skipped[report.nodeid] = str(report.longrepr[-1] if isinstance(report.longrepr, tuple) else report.longrepr)


def pytest_sessionfinish(session, exitstatus):
    ext = sys.modules.get("_netgraph_core")
    pkg = sys.modules.get("netgraph_core")
    record = {
        "seen": sorted(_seen),
        "skipped": _skipped,
        "hooks": hasattr(ext, "_max_flows_per_call"),
        "modules": [getattr(ext, "__file__", None), getattr(pkg, "__file__", None)],
    }
    with open(os.environ["NGC_SKIP_RECORD"], "w") as f:
        json.dump(record, f)
PY
    NGC_SKIP_RECORD="$tmp/record.json" PYTHONPATH="$tmp${PYTHONPATH:+:$PYTHONPATH}" \
        "$PYTHON" -m pytest -p ngc_skip_record "$@"
    [[ -s "$tmp/record.json" ]] || die "pytest wrote no record of its run"
    echo "$EXPECTED_OFF_SKIPS" | sed '/^$/d' > "$tmp/expected.txt"
    "$PYTHON" - "$hooks" "$tmp/record.json" "$tmp/expected.txt" "$EXPECTED_OFF_SKIPPED_PREFIX" "$ROOT" <<'PY'
import json
import os
import sys

hooks, record_path, expected_path, prefix, root = sys.argv[1:]
record = json.load(open(record_path))
# One node id per line: an id can hold spaces (a parametrization id such as "[integer as
# float]"), so the list is never split on other whitespace.
with open(expected_path) as f:
    expected = {line.strip() for line in f if line.strip()}
skipped = record["skipped"]
problems = []
if None in record["modules"]:
    problems.append(f"netgraph_core or _netgraph_core was not imported: {record['modules']}")
if hooks == "on":
    if not record["hooks"]:
        problems.append("the extension lacks the private hooks (expected the editable build)")
    if skipped:
        problems.append(f"{len(skipped)} tests skipped in a hooks-on run: {sorted(skipped)}")
else:
    if record["hooks"]:
        problems.append("the extension carries the private hooks (expected a wheel or sdist build)")
    inside = [m for m in record["modules"] if m and os.path.realpath(m).startswith(os.path.realpath(root) + os.sep)]
    if inside:
        problems.append(f"imported from the source tree, not an installation: {inside}")
    replays = {n for n in record["seen"] if n.startswith(prefix)}
    if not replays:
        problems.append(f"no test {prefix}...] ran")
    want = expected | replays
    got = set(skipped)
    for name in sorted(want - got):
        problems.append(f"expected to skip, did not: {name}")
    for name in sorted(got - want):
        problems.append(f"skipped, not a known hook-only test: {name}")
    for name, reason in sorted(skipped.items()):
        if "needs the private hooks" not in reason:
            problems.append(f"skipped for another reason: {name}: {reason}")
for p in problems:
    print(f"dev/run-checks.sh: {p}", file=sys.stderr)
if problems:
    sys.exit(1)
if hooks == "off":
    print(f"hooks off: exactly the {len(skipped)} expected hook-only tests skipped; imported {record['modules'][0]}")
else:
    print(f"hooks on: no test skipped; imported {record['modules'][0]}")
PY
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
