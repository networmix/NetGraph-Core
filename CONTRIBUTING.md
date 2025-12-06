# Contributing to NetGraph-Core

This document provides guidelines for setting up your development environment and submitting contributions.

## Development Setup

NetGraph-Core is a hybrid C++/Python project. You will need:

- Python 3.11+
- C++20 compatible compiler (GCC 10+, Clang 12+, MSVC 2019+)
- CMake 3.23+
- Ninja (recommended)

### 1. Clone the repository

```bash
git clone https://github.com/networmix/NetGraph-Core.git
cd NetGraph-Core
```

### 2. Create a virtual environment and install dependencies

We use a `Makefile` to simplify development tasks.

```bash
make dev
source venv/bin/activate
```

This will:

- Create a virtual environment in `./venv`
- Install all development dependencies
- Install pre-commit hooks

## Workflow

### Running Tests

Run all checks (linting + C++ tests + Python tests):

```bash
make check
```

Run specific test suites:

```bash
make cpp-test   # C++ tests (GoogleTest)
make py-test    # Python tests (pytest)
```

### Code Style

- **Python**: We use `ruff` for linting and formatting, and `pyright` for static type checking.
- **C++**: We follow standard C++20 practices.

Auto-fix Python formatting:

```bash
make format
```

## Release Process

Releases are automated via GitHub Actions when a new tag is pushed.

1. Bump version in `pyproject.toml` and `python/netgraph_core/_version.py`.
2. Commit and push.
3. Create and push a tag:

   ```bash
   git tag v0.1.0
   git push origin v0.1.0
   ```

4. The CI pipeline will build wheels, sdist, and publish to PyPI.
