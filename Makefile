# NetGraph-Core Development Makefile

.PHONY: help venv clean-venv dev install check check-ci lint format test qt clean build info hooks cov cpp-test rebuild check-python sanitize-test

.DEFAULT_GOAL := help

VENV_BIN := $(PWD)/venv/bin
# Use dynamic (recursive) assignment so a newly created venv is picked up
# Prefer explicit versions first (newest to oldest), then python3 on PATH
PY_FIND := $(shell for v in 3.13 3.12 3.11 3.10; do command -v python$$v >/dev/null 2>&1 && { echo python$$v; exit 0; }; done; command -v python3 2>/dev/null || command -v python 2>/dev/null)
PYTHON ?= $(if $(wildcard $(VENV_BIN)/python),$(VENV_BIN)/python,$(PY_FIND))
PIP = $(PYTHON) -m pip
PYTEST = $(PYTHON) -m pytest
RUFF = $(PYTHON) -m ruff
PRECOMMIT = $(PYTHON) -m pre_commit

# Detect Apple Command Line Tools compilers (prefer system toolchain on macOS)
APPLE_CLANG := $(shell xcrun --find clang 2>/dev/null)
APPLE_CLANGXX := $(shell xcrun --find clang++ 2>/dev/null)
DEFAULT_MACOSX := 15.0

help:
	@echo "🔧 NetGraph-Core Development Commands"
	@echo "  make venv         - Create a local virtualenv (./venv)"
	@echo "  make dev          - Install dev deps and pre-commit"
	@echo "  make install      - Editable install (no dev deps)"
	@echo "  make check        - Pre-commit (auto-fix) + C++/Python tests, then lint"
	@echo "  make check-ci     - Non-mutating lint + tests (CI entrypoint)"
	@echo "  make lint         - Ruff + Pyright"
	@echo "  make format       - Ruff format"
	@echo "  make test         - Run pytest"
	@echo "  make qt           - Quick tests (exclude slow/benchmark if marked)"
	@echo "  make cov          - Coverage summary + XML + single-page combined HTML"
	@echo "  make build        - Build wheel"
	@echo "  make clean        - Clean build artifacts"
	@echo "  make hooks        - Run pre-commit on all files"
	@echo "  make info         - Show tool versions"
	@echo "  make check-python - Check if venv Python matches system Python"
	@echo "  make rebuild      - Clean and rebuild with system Python (respects CMAKE_ARGS)"

# Allow callers to pass CMAKE_ARGS and MACOSX_DEPLOYMENT_TARGET consistently
ENV_MACOS := $(if $(MACOSX_DEPLOYMENT_TARGET),MACOSX_DEPLOYMENT_TARGET=$(MACOSX_DEPLOYMENT_TARGET),MACOSX_DEPLOYMENT_TARGET=$(DEFAULT_MACOSX))
ENV_CC := $(if $(APPLE_CLANG),CC=$(APPLE_CLANG),)
ENV_CXX := $(if $(APPLE_CLANGXX),CXX=$(APPLE_CLANGXX),)
ENV_CMAKE := $(if $(APPLE_CLANGXX),CMAKE_ARGS="$(strip $(CMAKE_ARGS) -DCMAKE_C_COMPILER=$(APPLE_CLANG) -DCMAKE_CXX_COMPILER=$(APPLE_CLANGXX))",$(if $(CMAKE_ARGS),CMAKE_ARGS="$(CMAKE_ARGS)",))
DEV_ENV := $(ENV_MACOS) $(ENV_CC) $(ENV_CXX) $(ENV_CMAKE)

dev:
	@echo "🚀 Setting up development environment..."
	@if [ ! -x "$(VENV_BIN)/python" ]; then \
		if [ -z "$(PY_FIND)" ]; then \
			echo "❌ Error: No Python interpreter found (python3 or python)"; \
			exit 1; \
		fi; \
		echo "🐍 Creating virtual environment with $(PY_FIND) ..."; \
		$(PY_FIND) -m venv venv || { echo "❌ Failed to create venv"; exit 1; }; \
		if [ ! -x "$(VENV_BIN)/python" ]; then \
			echo "❌ Error: venv creation failed - $(VENV_BIN)/python not found"; \
			exit 1; \
		fi; \
		$(VENV_BIN)/python -m pip install -U pip wheel; \
	fi
	@echo "📦 Installing dev dependencies..."
	@$(DEV_ENV) $(VENV_BIN)/python -m pip install -e .'[dev]'
	@echo "🔗 Installing pre-commit hooks..."
	@$(VENV_BIN)/python -m pre_commit install --install-hooks
	@echo "✅ Dev environment ready. Activate with: source venv/bin/activate"
	@$(MAKE) check-python

venv:
	@echo "🐍 Creating virtual environment in ./venv ..."
	@if [ -z "$(PY_FIND)" ]; then \
		echo "❌ Error: No Python interpreter found (python3 or python)"; \
		exit 1; \
	fi
	@$(PY_FIND) -m venv venv || { echo "❌ Failed to create venv"; exit 1; }
	@if [ ! -x "$(VENV_BIN)/python" ]; then \
		echo "❌ Error: venv creation failed - $(VENV_BIN)/python not found"; \
		exit 1; \
	fi
	@$(VENV_BIN)/python -m pip install -U pip wheel
	@echo "✅ venv ready. Activate with: source venv/bin/activate"

clean-venv:
	@rm -rf venv/

install:
	@echo "📦 Installing package (editable)"
	@$(DEV_ENV) $(PIP) install -e .


check:
	@PYTHON=$(PYTHON) bash dev/run-checks.sh
	@$(MAKE) lint

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

clean:
	@echo "🧹 Cleaning build artifacts and cache files..."
	@rm -rf build/ dist/ *.egg-info/ || true
	@find . -type d -name "__pycache__" -exec rm -rf {} + 2>/dev/null || true
	@find . -type f -name "*.pyc" -delete 2>/dev/null || true
	@find . -type f -name "*.pyo" -delete 2>/dev/null || true
	@rm -f .coverage coverage-*.xml coverage-*.html || true
	@rm -rf htmlcov-python || true
	@rm -rf Testing CTestTestfile.cmake || true
	@echo "✅ Cleanup complete!"

info:
	@echo "Python (active): $$($(PYTHON) --version)"
	@echo "Python (system): $$($(PY_FIND) --version 2>/dev/null || echo 'missing')"
	@$(MAKE) check-python
	@echo "Ruff: $$($(RUFF) --version 2>/dev/null || echo 'missing')"
	@echo "Pyright: $$($(PYTHON) -m pyright --version 2>/dev/null | head -1 || echo 'missing')"
	@echo "Pytest: $$($(PYTEST) --version 2>/dev/null || echo 'missing')"
	@echo "CMake: $$(cmake --version 2>/dev/null | head -1 || echo 'missing')"
	@echo "Ninja: $$(ninja --version 2>/dev/null || echo 'missing')"

check-python:
	@if [ -x "$(VENV_BIN)/python" ]; then \
		VENV_VER=$$($(VENV_BIN)/python -c "import sys; print(f'{sys.version_info.major}.{sys.version_info.minor}')" 2>/dev/null || echo "unknown"); \
		SYS_VER=$$($(PY_FIND) -c "import sys; print(f'{sys.version_info.major}.{sys.version_info.minor}')" 2>/dev/null || echo "unknown"); \
		if [ -n "$$VENV_VER" ] && [ -n "$$SYS_VER" ] && [ "$$VENV_VER" != "$$SYS_VER" ]; then \
			echo "⚠️  WARNING: venv Python ($$VENV_VER) != system Python ($$SYS_VER)"; \
			echo "   Run 'make clean-venv && make dev' to recreate venv with system Python"; \
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
		# Determine parallelism
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
	@$(PIP) install -U scikit-build-core "pybind11>=3"
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
		ctest --test-dir "$$BUILD_DIR" --output-on-failure || true
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

.PHONY: sanitize-test
sanitize-test:
	@echo "🧪 Building and running C++ tests with sanitizers..."
	@BUILD_DIR="build/cpp-sanitize"; \
		mkdir -p "$$BUILD_DIR"; \
		GEN_ARGS=""; \
		if command -v ninja >/dev/null 2>&1; then GEN_ARGS="-G Ninja"; fi; \
		cmake -S . -B "$$BUILD_DIR" -DNETGRAPH_CORE_BUILD_TESTS=ON -DNETGRAPH_CORE_SANITIZE=ON -DCMAKE_BUILD_TYPE=Debug $$GEN_ARGS; \
		cmake --build "$$BUILD_DIR" --config Debug -j; \
		ASAN_OPTIONS=detect_leaks=1 ctest --test-dir "$$BUILD_DIR" --output-on-failure || true

# Clean + reinstall in dev mode (respects CMAKE_ARGS and MACOSX_DEPLOYMENT_TARGET)
# Always uses system Python to avoid venv version mismatches
rebuild: clean
	@echo "🔨 Rebuilding with system Python: $(PY_FIND)"
	@$(DEV_ENV) $(PY_FIND) -m pip install -e .'[dev]'
