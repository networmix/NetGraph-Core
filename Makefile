# NetGraph-Core Development Makefile

.PHONY: help venv clean-venv dev install check check-ci lint format test qt clean build info hooks cov cpp-test

.DEFAULT_GOAL := help

VENV_BIN := $(PWD)/venv/bin
# Use dynamic (recursive) assignment so a newly created venv is picked up
PYTHON = $(if $(wildcard $(VENV_BIN)/python),$(VENV_BIN)/python,python3)
PIP = $(PYTHON) -m pip
PYTEST = $(PYTHON) -m pytest
RUFF = $(PYTHON) -m ruff
PRECOMMIT = $(PYTHON) -m pre_commit

help:
	@echo "🔧 NetGraph-Core Development Commands"
	@echo "  make venv       - Create a local virtualenv (./venv)"
	@echo "  make dev        - Install dev deps and pre-commit"
	@echo "  make install    - Editable install (no dev deps)"
	@echo "  make check      - Pre-commit (auto-fix) + C++/Python tests, then lint"
	@echo "  make check-ci   - Non-mutating lint + tests (CI entrypoint)"
	@echo "  make lint       - Ruff + Pyright"
	@echo "  make format     - Ruff format"
	@echo "  make test       - Run pytest"
	@echo "  make qt         - Quick tests (exclude slow/benchmark if marked)"
	@echo "  make cov        - Coverage summary + XML + single-page combined HTML"
	@echo "  make build      - Build wheel"
	@echo "  make clean      - Clean build artifacts"
	@echo "  make hooks      - Run pre-commit on all files"
	@echo "  make info       - Show tool versions"

dev:
	@echo "🚀 Setting up development environment..."
	@if [ ! -x "$(VENV_BIN)/python" ]; then \
		echo "🐍 Creating virtual environment in ./venv ..."; \
		python3 -m venv venv; \
		$(VENV_BIN)/python -m pip install -U pip wheel; \
	fi
	@echo "📦 Installing dev dependencies..."
	@$(VENV_BIN)/python -m pip install -e .[dev]
	@echo "🔗 Installing pre-commit hooks..."
	@$(VENV_BIN)/python -m pre_commit install --install-hooks
	@echo "✅ Dev environment ready. Activate with: source venv/bin/activate"

venv:
	@echo "🐍 Creating virtual environment in ./venv ..."
	@python3 -m venv venv
	@$(VENV_BIN)/python -m pip install -U pip wheel
	@echo "✅ venv ready. Activate with: source venv/bin/activate"

clean-venv:
	@rm -rf venv/

install:
	@echo "📦 Installing package (editable)"
	@$(PIP) install -e .


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
	@$(PYTEST) -m "not slow and not benchmark"

build:
	@$(PYTHON) -m build

clean:
	@rm -rf build/ dist/ *.egg-info/ **/__pycache__ */__pycache__
	@rm -f .coverage coverage-*.xml coverage-*.html
	@rm -rf htmlcov-python
	@rm -rf Testing CTestTestfile.cmake

info:
	@echo "Python: $$($(PYTHON) --version)"
	@echo "Ruff: $$($(RUFF) --version 2>/dev/null || echo 'missing')"
	@echo "Pyright: $$($(PYTHON) -m pyright --version 2>/dev/null | head -1 || echo 'missing')"
	@echo "Pytest: $$($(PYTEST) --version 2>/dev/null || echo 'missing')"
	@echo "CMake: $$(cmake --version 2>/dev/null | head -1 || echo 'missing')"
	@echo "Ninja: $$(ninja --version 2>/dev/null || echo 'missing')"

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
	@PIP_NO_BUILD_ISOLATION=1 CMAKE_ARGS="-DNETGRAPH_CORE_COVERAGE=ON" $(PIP) install -e .[dev]
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
