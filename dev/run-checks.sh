#!/bin/bash
# NetGraph-Core: Comprehensive checks and tests
# - Pre-commit (auto-fix, then verify)
# - C++ tests via ctest (optional, if cmake is available)
# - Python tests (with coverage if pytest-cov is available)

set -e
set -u
set -o pipefail

# Prefer provided PYTHON env (passed from Makefile), fallback to python3
PYTHON=${PYTHON:-python3}
# Prepend the directory of the configured Python (usually venv/bin) to PATH
PY_BIN_DIR=$(dirname "$PYTHON")
export PATH="$PY_BIN_DIR:$PATH"
# If running within project venv, also prepend its bin to PATH so cmake/ninja resolve from venv
if [ -n "${VIRTUAL_ENV:-}" ] && [ -x "$VIRTUAL_ENV/bin/python" ]; then
    export PATH="$VIRTUAL_ENV/bin:$PATH"
fi

# Prefer Apple Command Line Tools compilers to avoid Homebrew libc++ ABI mismatches
APPLE_CLANG=$(xcrun --find clang 2>/dev/null || true)
APPLE_CLANGXX=$(xcrun --find clang++ 2>/dev/null || true)
DEFAULT_MACOSX=15.0
if [ -n "$APPLE_CLANG" ] && [ -n "$APPLE_CLANGXX" ]; then
    export CC="$APPLE_CLANG"
    export CXX="$APPLE_CLANGXX"
    export CMAKE_ARGS="${CMAKE_ARGS:-} -DCMAKE_C_COMPILER=$CC -DCMAKE_CXX_COMPILER=$CXX"
fi
export MACOSX_DEPLOYMENT_TARGET="${MACOSX_DEPLOYMENT_TARGET:-$DEFAULT_MACOSX}"

# Ensure required tools are present
if ! "$PYTHON" -m pre_commit --version &> /dev/null; then
    echo "❌ pre-commit is not installed. Please run 'make dev' first."
    exit 1
fi

if ! "$PYTHON" -m pytest --version &> /dev/null; then
    echo "❌ pytest is not installed. Please run 'make dev' first."
    exit 1
fi

# Ensure pre-commit hooks are installed (if in a git repo)
if [ -d .git ] && [ ! -f .git/hooks/pre-commit ]; then
    echo "⚠️  Pre-commit hooks not installed. Installing now..."
    "$PYTHON" -m pre_commit install
    echo ""
fi

# First pass: allow auto-fixes
echo "🏃 Running pre-commit (first pass: apply auto-fixes if needed)..."
set +e
"$PYTHON" -m pre_commit run --all-files
first_pass_status=$?
set -e

if [ $first_pass_status -ne 0 ]; then
    echo "ℹ️  Some hooks modified files or reported issues. Re-running checks..."
fi

# Second pass: verify all checks
echo "🏃 Running pre-commit (second pass: verify all checks)..."
if ! "$PYTHON" -m pre_commit run --all-files; then
    echo ""
    echo "❌ Pre-commit checks failed after applying fixes. Please address the issues above."
    exit 1
fi

autofixed=0
if [ $first_pass_status -ne 0 ]; then
    autofixed=1
fi

echo ""
echo "✅ Pre-commit checks passed!"
echo ""

# Optional: Run C++ tests via ctest if cmake is available (and not skipped)
if [ "${SKIP_CPP_TESTS:-0}" = "1" ]; then
    echo "⏭️  SKIP_CPP_TESTS=1 set. Skipping C++ tests."
else
    echo "🔧 Checking for CMake to run C++ tests..."
    if command -v cmake >/dev/null 2>&1; then
        if [ "${SKIP_CPP_TESTS:-0}" != "1" ]; then
            echo "🧪 Running C++ tests (ctest)..."
            BUILD_DIR="build/cpp-tests"
            mkdir -p "$BUILD_DIR"
            GEN_ARGS=""
            if command -v ninja >/dev/null 2>&1; then
                GEN_ARGS="-G Ninja"
            fi
            cmake -S . -B "$BUILD_DIR" -DNETGRAPH_CORE_BUILD_TESTS=ON -DCMAKE_BUILD_TYPE=Release $GEN_ARGS
            cmake --build "$BUILD_DIR" --config Release -j
            # Determine parallelism
            if command -v sysctl >/dev/null 2>&1; then
                NPROC=$(sysctl -n hw.ncpu 2>/dev/null || echo 2)
            elif command -v nproc >/dev/null 2>&1; then
                NPROC=$(nproc)
            else
                NPROC=2
            fi
            # Add per-test timeout to avoid hangs
            ctest --test-dir "$BUILD_DIR" --output-on-failure -j "$NPROC" --timeout 120
        fi
    else
        echo "⚠️  CMake not found. Skipping C++ tests."
    fi
fi

echo ""

# Run Python tests (with coverage if available)
echo "🧪 Running Python tests..."
# Ensure the extension is built for the active Python and importable
echo "🔧 Installing project in editable mode for current Python..."
"$PYTHON" -m pip install -e . >/dev/null
if "$PYTHON" -c "import pytest_cov" >/dev/null 2>&1; then
    "$PYTHON" -m pytest --cov=netgraph_core --cov-report=term-missing
else
    "$PYTHON" -m pytest
fi

if [ $? -eq 0 ]; then
    echo ""
    if [ $autofixed -eq 1 ]; then
        echo "🎉 All checks and tests passed. Auto-fixes were applied by pre-commit."
    else
        echo "🎉 All checks and tests passed."
    fi
else
    echo ""
    if [ $autofixed -eq 1 ]; then
        echo "❌ Some tests failed. Note: auto-fixes were applied earlier by pre-commit."
    else
        echo "❌ Some tests failed."
    fi
    exit 1
fi
