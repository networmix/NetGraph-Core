# NetGraph-Core Development Makefile

.PHONY: help venv clean-venv dev install install-native install-pgo check check-ci lint format test qt build clean check-dist publish-test publish info hooks check-python cpp-test cov sanitize-test rebuild

.DEFAULT_GOAL := help

# --------------------------------------------------------------------------
# Python interpreter detection
# --------------------------------------------------------------------------
# VENV_BIN: path to local virtualenv bin directory
VENV_BIN := $(PWD)/venv/bin

# PY_BEST: scan for newest supported Python (used when creating new venvs)
# Supports 3.11-3.13 to match requires-python >=3.11
PY_BEST := $(shell for v in 3.13 3.12 3.11; do command -v python$$v >/dev/null 2>&1 && { echo python$$v; exit 0; }; done; command -v python3 2>/dev/null || command -v python 2>/dev/null)

# PY_PATH: active python3/python on PATH (respects CI setup-python and activated venvs)
PY_PATH := $(shell command -v python3 2>/dev/null || command -v python 2>/dev/null)

# PYTHON: interpreter used for all commands
#   1. Use local venv if present
#   2. Otherwise use active python on PATH (important for CI)
#   3. Fall back to best available version
#   4. Final fallback to 'python3' literal for clear error messages
PYTHON ?= $(if $(wildcard $(VENV_BIN)/python),$(VENV_BIN)/python,$(if $(PY_PATH),$(PY_PATH),$(if $(PY_BEST),$(PY_BEST),python3)))

# Derived tool commands (always use -m to ensure correct environment)
PIP := $(PYTHON) -m pip
PYTEST := $(PYTHON) -m pytest
RUFF := $(PYTHON) -m ruff
PRECOMMIT := $(PYTHON) -m pre_commit

# Detect Apple Command Line Tools compilers (prefer system toolchain on macOS)
APPLE_CLANG := $(shell xcrun --find clang 2>/dev/null)
APPLE_CLANGXX := $(shell xcrun --find clang++ 2>/dev/null)
DEFAULT_MACOSX := 15.0

help:
	@echo "🔧 NetGraph-Core Development Commands"
	@echo ""
	@echo "Setup & Installation:"
	@echo "  make venv          - Create a local virtualenv (./venv)"
	@echo "  make dev           - Full development environment (package + dev deps + hooks)"
	@echo "  make install       - Install package (default optimizations: LTO, loop unrolling)"
	@echo "  make install-native - Install with CPU-specific optimizations (faster, not portable)"
	@echo "  make install-pgo   - Profile-guided optimization (two-phase build, experimental)"
	@echo "  make clean-venv    - Remove virtual environment"
	@echo "  make rebuild       - Clean and rebuild (respects CMAKE_ARGS)"
	@echo ""
	@echo "Code Quality & Testing:"
	@echo "  make check         - Run pre-commit (auto-fix) + C++/Python tests, then lint"
	@echo "  make check-ci      - Run non-mutating lint + tests (CI entrypoint)"
	@echo "  make lint          - Run only linting (non-mutating: ruff + pyright)"
	@echo "  make format        - Auto-format code with ruff"
	@echo "  make test          - Run tests with coverage"
	@echo "  make qt            - Run quick tests only (exclude slow/benchmark)"
	@echo "  make cpp-test      - Build and run C++ tests"
	@echo "  make cov           - Coverage summary + XML + single-page combined HTML"
	@echo "  make sanitize-test - Build and run C++ tests with sanitizers"
	@echo ""
	@echo "Build & Package:"
	@echo "  make build         - Build distribution packages"
	@echo "  make clean         - Clean build artifacts and cache files"
	@echo "  make check-dist    - Check distribution packages with twine"
	@echo "  make publish-test  - Publish to Test PyPI"
	@echo "  make publish       - Publish to PyPI"
	@echo ""
	@echo "Utilities:"
	@echo "  make info          - Show project information"
	@echo "  make hooks         - Run pre-commit on all files"
	@echo "  make check-python  - Check if venv Python matches best available"

# Allow callers to pass CMAKE_ARGS and MACOSX_DEPLOYMENT_TARGET consistently
ENV_MACOS := $(if $(MACOSX_DEPLOYMENT_TARGET),MACOSX_DEPLOYMENT_TARGET=$(MACOSX_DEPLOYMENT_TARGET),MACOSX_DEPLOYMENT_TARGET=$(DEFAULT_MACOSX))
ENV_CC := $(if $(APPLE_CLANG),CC=$(APPLE_CLANG),)
ENV_CXX := $(if $(APPLE_CLANGXX),CXX=$(APPLE_CLANGXX),)
ENV_CMAKE := $(if $(APPLE_CLANGXX),CMAKE_ARGS="$(strip $(CMAKE_ARGS) -DCMAKE_C_COMPILER=$(APPLE_CLANG) -DCMAKE_CXX_COMPILER=$(APPLE_CLANGXX))",$(if $(CMAKE_ARGS),CMAKE_ARGS="$(CMAKE_ARGS)",))
DEV_ENV := $(ENV_MACOS) $(ENV_CC) $(ENV_CXX) $(ENV_CMAKE)

dev:
	@echo "🚀 Setting up development environment..."
	@if [ ! -x "$(VENV_BIN)/python" ]; then \
		if [ -z "$(PY_BEST)" ]; then \
			echo "❌ Error: No Python interpreter found (python3 or python)"; \
			exit 1; \
		fi; \
		echo "🐍 Creating virtual environment with $(PY_BEST) ..."; \
		$(PY_BEST) -m venv venv || { echo "❌ Failed to create venv"; exit 1; }; \
		if [ ! -x "$(VENV_BIN)/python" ]; then \
			echo "❌ Error: venv creation failed - $(VENV_BIN)/python not found"; \
			exit 1; \
		fi; \
		$(VENV_BIN)/python -m pip install -U pip setuptools wheel; \
	fi
	@echo "📦 Installing dev dependencies..."
	@$(DEV_ENV) $(VENV_BIN)/python -m pip install -e .'[dev]'
	@echo "🔗 Installing pre-commit hooks..."
	@$(VENV_BIN)/python -m pre_commit install --install-hooks
	@echo "✅ Dev environment ready. Activate with: source venv/bin/activate"
	@$(MAKE) check-python

venv:
	@echo "🐍 Creating virtual environment in ./venv ..."
	@if [ -z "$(PY_BEST)" ]; then \
		echo "❌ Error: No Python interpreter found (python3 or python)"; \
		exit 1; \
	fi
	@$(PY_BEST) -m venv venv || { echo "❌ Failed to create venv"; exit 1; }
	@if [ ! -x "$(VENV_BIN)/python" ]; then \
		echo "❌ Error: venv creation failed - $(VENV_BIN)/python not found"; \
		exit 1; \
	fi
	@$(VENV_BIN)/python -m pip install -U pip setuptools wheel
	@echo "✅ venv ready. Activate with: source venv/bin/activate"

clean-venv:
	@rm -rf venv/

install:
	@echo "📦 Installing package (editable, default optimizations)"
	@$(DEV_ENV) $(PIP) install -e .

install-native:
	@echo "📦 Installing with native CPU optimizations (-march=native)"
	@echo "   Note: Binary will only work on this CPU architecture"
	@$(ENV_MACOS) $(ENV_CC) $(ENV_CXX) CMAKE_ARGS="$(strip $(CMAKE_ARGS) -DNETGRAPH_CORE_NATIVE=ON)" $(PIP) install -e .

# PGO profile stored outside build/ so it survives rebuild
PGO_DIR := $(PWD)/.pgo-profile

install-pgo:
	@echo "📦 PGO Build Phase 1/3: Instrumenting..."
	@rm -rf build/ $(PGO_DIR)
	@mkdir -p $(PGO_DIR)
	@$(ENV_MACOS) $(ENV_CC) $(ENV_CXX) CMAKE_ARGS="$(strip $(CMAKE_ARGS) -DNETGRAPH_CORE_PGO_GENERATE=ON -DNETGRAPH_CORE_PGO_DIR=$(PGO_DIR) -DNETGRAPH_CORE_NATIVE=ON)" $(PIP) install -e .
	@echo "📦 PGO Build Phase 2/3: Collecting profile (running benchmark)..."
	@$(PYTHON) dev/benchmark_profiling_overhead.py --mesh-size 25 --spf-iters 1000 --flow-iters 200 >/dev/null 2>&1
	@# Clang: merge raw profiles into .profdata (profile may be in cwd or PGO_DIR)
	@if ls $(PGO_DIR)/*.profraw >/dev/null 2>&1; then \
		echo "   Merging Clang profile data..."; \
		xcrun llvm-profdata merge -output=$(PGO_DIR)/default.profdata $(PGO_DIR)/*.profraw 2>/dev/null || \
		llvm-profdata merge -output=$(PGO_DIR)/default.profdata $(PGO_DIR)/*.profraw; \
	elif ls *.profraw >/dev/null 2>&1; then \
		mv *.profraw $(PGO_DIR)/; \
		xcrun llvm-profdata merge -output=$(PGO_DIR)/default.profdata $(PGO_DIR)/*.profraw 2>/dev/null || \
		llvm-profdata merge -output=$(PGO_DIR)/default.profdata $(PGO_DIR)/*.profraw; \
	fi
	@echo "📦 PGO Build Phase 3/3: Rebuilding with profile data..."
	@rm -rf build/
	@$(ENV_MACOS) $(ENV_CC) $(ENV_CXX) CMAKE_ARGS="$(strip $(CMAKE_ARGS) -DNETGRAPH_CORE_PGO_USE=ON -DNETGRAPH_CORE_PGO_DIR=$(PGO_DIR) -DNETGRAPH_CORE_NATIVE=ON)" $(PIP) install -e .
	@echo "✅ PGO build complete"

check:
	@PYTHON=$(PYTHON) bash dev/run-checks.sh

check-ci:
	@$(MAKE) lint
	@$(MAKE) cpp-test
	@$(MAKE) test

lint:
	@$(RUFF) format --check .
	@$(RUFF) check .
	@$(PYTHON) -m pyright

format:
	@$(RUFF) format .

test:
	@$(PYTEST)

qt:
	@$(PYTEST) --no-cov -m "not slow and not benchmark"

build:
	@echo "🏗️  Building distribution packages..."
	@if $(PYTHON) -c "import build" >/dev/null 2>&1; then \
		$(PYTHON) -m build; \
	else \
		echo "❌ build module not installed. Install dev dependencies with: make dev"; \
		exit 1; \
	fi

check-dist:
	@echo "🔍 Checking distribution packages..."
	@if $(PYTHON) -c "import twine" >/dev/null 2>&1; then \
		$(PYTHON) -m twine check dist/*; \
	else \
		echo "❌ twine not installed. Install dev dependencies with: make dev"; \
		exit 1; \
	fi

publish-test:
	@echo "📦 Publishing to Test PyPI..."
	@if $(PYTHON) -c "import twine" >/dev/null 2>&1; then \
		$(PYTHON) -m twine upload --repository testpypi dist/*; \
	else \
		echo "❌ twine not installed. Install dev dependencies with: make dev"; \
		exit 1; \
	fi

publish:
	@echo "🚀 Publishing to PyPI..."
	@if $(PYTHON) -c "import twine" >/dev/null 2>&1; then \
		$(PYTHON) -m twine upload dist/*; \
	else \
		echo "❌ twine not installed. Install dev dependencies with: make dev"; \
		exit 1; \
	fi

clean:
	@echo "🧹 Cleaning build artifacts and cache files..."
	@rm -rf build/ dist/ *.egg-info/ || true
	@rm -rf .pytest_cache .ruff_cache .mypy_cache htmlcov htmlcov-python .coverage coverage.xml coverage-*.xml .benchmarks .pytest-benchmark || true
	@rm -rf Testing CTestTestfile.cmake || true
	@find . -path "./venv" -prune -o -type f -name "*.pyc" -delete 2>/dev/null || true
	@find . -path "./venv" -prune -o -type d -name "__pycache__" -exec rm -rf {} + 2>/dev/null || true
	@find . -path "./venv" -prune -o -type f -name "*.pyo" -delete 2>/dev/null || true
	@find . -path "./venv" -prune -o -type f -name "*~" -delete 2>/dev/null || true
	@find . -path "./venv" -prune -o -type f -name "*.orig" -delete 2>/dev/null || true
	@echo "✅ Cleanup complete!"

info:
	@echo "Python (active): $$($(PYTHON) --version)"
	@echo "Python (best):   $$($(PY_BEST) --version 2>/dev/null || echo 'missing')"
	@$(MAKE) check-python
	@echo "Ruff: $$($(RUFF) --version 2>/dev/null || echo 'missing')"
	@echo "Pyright: $$($(PYTHON) -m pyright --version 2>/dev/null | head -1 || echo 'missing')"
	@echo "Pytest: $$($(PYTEST) --version 2>/dev/null || echo 'missing')"
	@echo "CMake: $$(cmake --version 2>/dev/null | head -1 || echo 'missing')"
	@echo "Ninja: $$(ninja --version 2>/dev/null || echo 'missing')"

check-python:
	@if [ -x "$(VENV_BIN)/python" ]; then \
		VENV_VER=$$($(VENV_BIN)/python -c "import sys; print(f'{sys.version_info.major}.{sys.version_info.minor}')" 2>/dev/null || echo "unknown"); \
		BEST_VER=$$($(PY_BEST) -c "import sys; print(f'{sys.version_info.major}.{sys.version_info.minor}')" 2>/dev/null || echo "unknown"); \
		if [ -n "$$VENV_VER" ] && [ -n "$$BEST_VER" ] && [ "$$VENV_VER" != "$$BEST_VER" ]; then \
			echo "⚠️  WARNING: venv Python ($$VENV_VER) != best available Python ($$BEST_VER)"; \
			echo "   Run 'make clean-venv && make dev' to recreate venv if desired"; \
		fi; \
	fi

hooks:
	@$(PRECOMMIT) run --all-files || (echo "Some pre-commit hooks failed. Fix and re-run." && exit 1)

cpp-test:
	@echo "🧪 Building and running C++ tests..."
	@BUILD_DIR="build/cpp-tests"; \
		mkdir -p "$$BUILD_DIR"; \
		GEN_ARGS=""; \
		if command -v ninja >/dev/null 2>&1; then GEN_ARGS="-G Ninja"; fi; \
		cmake -S . -B "$$BUILD_DIR" -DNETGRAPH_CORE_BUILD_TESTS=ON -DCMAKE_BUILD_TYPE=Release $$GEN_ARGS; \
		cmake --build "$$BUILD_DIR" --config Release -j; \
		if command -v sysctl >/dev/null 2>&1; then \
			NPROC=$$(sysctl -n hw.ncpu 2>/dev/null || echo 2); \
		elif command -v nproc >/dev/null 2>&1; then \
			NPROC=$$(nproc); \
		else \
			NPROC=2; \
		fi; \
		ctest --test-dir "$$BUILD_DIR" --output-on-failure -j "$$NPROC" --timeout 120

cov:
	@echo "📦 Reinstalling with C++ coverage instrumentation..."
	@$(PIP) install -U scikit-build-core "pybind11>=2.11"
	@PIP_NO_BUILD_ISOLATION=1 CMAKE_ARGS="-DNETGRAPH_CORE_COVERAGE=ON" $(PIP) install -e .'[dev]'
	@echo "🧪 Running Python tests with coverage..."
	@mkdir -p build/coverage
	@$(PYTEST) --cov=netgraph_core --cov-report=term-missing --cov-report=xml:build/coverage/coverage-python.xml
	@echo "🛠️  Building and running C++ tests with coverage..."
	@BUILD_DIR="build/cpp-tests-cov"; \
		mkdir -p "$$BUILD_DIR"; \
		GEN_ARGS=""; \
		if command -v ninja >/dev/null 2>&1; then GEN_ARGS="-G Ninja"; fi; \
		cmake -S . -B "$$BUILD_DIR" -DNETGRAPH_CORE_BUILD_TESTS=ON -DNETGRAPH_CORE_COVERAGE=ON -DCMAKE_BUILD_TYPE=Debug $$GEN_ARGS; \
		cmake --build "$$BUILD_DIR" --config Debug -j; \
		ctest --test-dir "$$BUILD_DIR" --output-on-failure || echo "⚠️  Some C++ tests failed (continuing for coverage)"
	@echo "📈 Generating C++ coverage XML (gcovr)..."
	@$(PYTHON) -m gcovr --root . \
		--object-directory build \
		--object-directory build/cpp-tests-cov \
		--filter 'include/netgraph' --filter 'src' --exclude 'tests' --exclude 'bindings/.*' --exclude '.*pybind11.*' --exclude '_deps/pybind11-src/.*' \
		--gcov-ignore-errors=all \
		--xml-pretty -o build/coverage/coverage-cpp.xml
	@echo ""
	@echo "================ Python + C++ coverage (summary) ================"
	@$(PYTHON) dev/coverage_summary.py build/coverage/coverage-python.xml build/coverage/coverage-cpp.xml --html=build/coverage/coverage-combined.html
	@echo ""
	@echo "✅ Coverage ready in build/coverage/: coverage-python.xml, coverage-cpp.xml, coverage-combined.html"

sanitize-test:
	@echo "🧪 Building and running C++ tests with sanitizers..."
	@BUILD_DIR="build/cpp-sanitize"; \
		mkdir -p "$$BUILD_DIR"; \
		GEN_ARGS=""; \
		if command -v ninja >/dev/null 2>&1; then GEN_ARGS="-G Ninja"; fi; \
		cmake -S . -B "$$BUILD_DIR" -DNETGRAPH_CORE_BUILD_TESTS=ON -DNETGRAPH_CORE_SANITIZE=ON -DCMAKE_BUILD_TYPE=Debug $$GEN_ARGS; \
		cmake --build "$$BUILD_DIR" --config Debug -j; \
		ASAN_OPTIONS=detect_leaks=1 ctest --test-dir "$$BUILD_DIR" --output-on-failure || echo "⚠️  Some sanitizer tests failed"

# Clean + reinstall in dev mode (respects CMAKE_ARGS and MACOSX_DEPLOYMENT_TARGET)
# Uses active PYTHON (venv or PATH) to avoid environment mismatches
rebuild: clean
	@echo "🔨 Rebuilding with: $(PYTHON)"
	@$(DEV_ENV) $(PIP) install -e .'[dev]'
