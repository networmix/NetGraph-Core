# NetGraph-Core Development Makefile

.PHONY: help venv clean-venv dev install install-native check check-ci lint format test build clean check-dist publish-test publish info hooks check-python cpp-test cpp-test-off cov sanitize-test tsan-test stubs-check perf-tool-test rebuild

.DEFAULT_GOAL := help

# --------------------------------------------------------------------------
# Python interpreter detection
# --------------------------------------------------------------------------
# VENV_BIN: the checkout's virtualenv bin directory (CURDIR is the directory make runs in, -C included)
VENV_BIN := $(CURDIR)/venv/bin

# PY_BEST: scan for newest supported Python (used when creating new venvs)
# Supports 3.11-3.14 to match requires-python >=3.11
PY_BEST := $(shell for v in 3.14 3.13 3.12 3.11; do command -v python$$v >/dev/null 2>&1 && { echo python$$v; exit 0; }; done; command -v python3 2>/dev/null || command -v python 2>/dev/null)

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
	@echo "  make install       - Install package (editable, private hooks; LTO, loop unrolling)"
	@echo "  make install-native - Install with CPU-specific optimizations (faster, not portable)"
	@echo "  make clean-venv    - Remove virtual environment"
	@echo "  make rebuild       - Clean and rebuild (respects CMAKE_ARGS)"
	@echo ""
	@echo "Code Quality & Testing:"
	@echo "  make check         - Rebuild the editable extension, then everything check-ci runs"
	@echo "  make check-ci      - Lint, fixture data check, the pytest guard's rules, C++ tests (hooks on and off), Python tests"
	@echo "                       (the CI entrypoint)"
	@echo "                       (check and check-ci only verify: they never edit a file)"
	@echo "  make lint          - Run only linting (non-mutating: ruff + pyright)"
	@echo "  make format        - Auto-format code with ruff (edits files)"
	@echo "  make test          - Run Python tests (pytest)"
	@echo "  make cpp-test      - Build and run C++ tests (private hooks on)"
	@echo "  make cpp-test-off  - Build and run C++ tests without the hooks (the shipped configuration)"
	@echo "  make cov           - Coverage (own venv and builds) + XML + single-page combined HTML"
	@echo "  make sanitize-test - Build and run C++ tests with sanitizers"
	@echo "  make tsan-test     - Build and run C++ tests with ThreadSanitizer"
	@echo "  make stubs-check   - Type stub: drift, wheel/sdist contents, typed clients (own venvs)"
	@echo "                       (regenerate the stub: bash dev/stubs/run.sh generate)"
	@echo "  make perf-tool-test - Tests of the A/B/A timing tool in dev/perf (after changing it; CI: macOS)"
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
	@echo "  make hooks         - Run pre-commit on all files (applies its fixes)"
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

# The project options are stated in pyproject.toml for every build; a variant overrides one of them
# with -C (a config setting), so it cannot stay in the shared build directory's cache.
install-native:
	@echo "📦 Installing with native CPU optimizations (-march=native)"
	@echo "   Note: Binary will only work on this CPU architecture"
	@$(DEV_ENV) $(PIP) install -e . -C cmake.define.NETGRAPH_CORE_NATIVE=ON

# Verification only: neither target edits a file (make format and make hooks apply fixes).
check:
	@$(MAKE) install
	@$(MAKE) check-ci

# The fixture data check: tests/data/policy_fixtures.json still equals what its frozen reference
# computes (tests/reference/policy/README.md; standard library only; -B writes no __pycache__).
check-ci:
	@$(MAKE) lint
	@$(PYTHON) -I -B tests/reference/policy/export_policy_fixtures.py --check
	@$(PYTHON) -I -B dev/test_pytest_guard.py
	@$(MAKE) cpp-test
	@$(MAKE) cpp-test-off
	@PYTHON=$(PYTHON) bash dev/run-checks.sh pytest on

# Pyright resolves imports (NumPy in the type stub, the extension) through the interpreter it
# is given; without --pythonpath it would take whichever python3 comes first on PATH.
lint:
	@$(RUFF) format --check .
	@$(RUFF) check .
	@$(PYTHON) -m pyright --pythonpath "$$($(PYTHON) -c 'import sys; print(sys.executable)')"

format:
	@$(RUFF) format .

test:
	@$(PYTEST)

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
	@rm -rf .pytest_cache .ruff_cache || true
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
	@PATH="$(VENV_BIN):$$PATH" bash dev/cpp-tests.sh release

cpp-test-off:
	@PATH="$(VENV_BIN):$$PATH" bash dev/cpp-tests.sh off

# The summary script's own rules first (its refusals make a failed measurement fail; no build).
cov:
	@$(PYTHON) -I -B dev/test_coverage_summary.py
	@PYTHON=$(PYTHON) bash dev/run-checks.sh cov

sanitize-test:
	@PATH="$(VENV_BIN):$$PATH" bash dev/cpp-tests.sh sanitize

tsan-test:
	@PATH="$(VENV_BIN):$$PATH" bash dev/cpp-tests.sh tsan

stubs-check:
	@PYTHON=$(PYTHON) bash dev/stubs/run.sh check

perf-tool-test:
	@$(PYTHON) -B dev/perf/test_tool.py

# Clean + reinstall in dev mode (respects CMAKE_ARGS and MACOSX_DEPLOYMENT_TARGET)
# Uses active PYTHON (venv or PATH) to avoid environment mismatches
rebuild: clean
	@echo "🔨 Rebuilding with: $(PYTHON)"
	@$(DEV_ENV) $(PIP) install -e .'[dev]'
