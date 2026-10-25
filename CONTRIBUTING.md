# Contributing to NetGraph-Core

This document provides guidelines for setting up your development environment and submitting contributions.

## Development Setup

NetGraph-Core is a hybrid C++/Python project. You will need:

- Python 3.11+
- C++20 compatible compiler (GCC 10+, Clang 12+, MSVC 2019+)
- CMake 3.23+
- Ninja (the C++ test runner uses it; `make dev` installs CMake and Ninja into the venv)

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

Run all checks: rebuild the editable extension, then lint (ruff, pyright), the check that the
policy fixture data still equals its frozen reference (`tests/reference/policy/README.md`), the
C++ tests with and without the private test hooks, and the Python tests. It only verifies and
never edits a file; `make check-ci` is the same without the rebuild (CI installs the package
itself):

```bash
make check
```

Run specific test suites:

```bash
make cpp-test       # C++ tests (GoogleTest), private hooks on
make cpp-test-off   # C++ tests without the hooks (the shipped configuration)
make test           # Python tests (pytest)
make sanitize-test  # C++ tests under ASan/UBSan (required for C++ changes)
make tsan-test      # C++ tests under ThreadSanitizer (required for C++ changes)
make cov            # Coverage of the C++ and Python suites (its own venv and builds)
```

### Code Style

- **Python**: We use `ruff` for linting and formatting, and `pyright` for static type checking.
- **C++**: We follow standard C++20 practices.

Apply fixes (these edit files):

```bash
make format   # ruff formatting
make hooks    # every pre-commit hook, including ruff --fix and whitespace fixers
```

The commit hook and `make hooks` run Pyright as `make lint` does, from the `make dev` venv: it
needs the built extension to resolve the package's imports.

### Type Stub

`python/netgraph_core/__init__.pyi` is generated from the compiled module. After changing a
binding's signature or docstring, run `bash dev/stubs/run.sh generate` and commit the stub;
`make stubs-check` (also in CI) fails when the stub drifts or the typed client checks fail.

### Behaviour Documentation

[docs/semantics.md](docs/semantics.md) is the maintained statement of what the library computes.
A change to observable behaviour (a result, an order, an outcome, an error, a default) updates it
in the same change, and its examples must still print what the document shows:
`tests/py/test_docs_examples.py` runs them, and README's quick start, in the Python tests.

### Mathematical Proofs

`proofs/` is a Lean project that proves properties of the round-robin rule;
[proofs/README.md](proofs/README.md) gives the commands and
[dev/proofs/README.md](dev/proofs/README.md) the rules for adding or changing a theorem. Ordinary
native and Python development does not need Lean. A theorem is accepted only with a plain
explanation of its statement and its proof, supporting lemmas included, that an engineer without
mathematical training can follow: a worked example and a boundary example, and a failing case
for each assumption that can fail, with every input, step and result explained.

## Release Process

Releases are automated via GitHub Actions when a new tag is pushed.

1. Move the pending notes into a new `## [X.Y.Z] - YYYY-MM-DD` section in `CHANGELOG.md`,
   dated the day you tag. Mark anything that breaks consumers as **BREAKING**.
2. Bump `version` in `pyproject.toml` to the same `X.Y.Z`. Release CI compares the tag
   with this value and refuses to publish on a mismatch.
3. Commit and push, then wait for the `Lean verification` workflow to pass on the pushed
   branch: publication waits for the same workflow on the tag, so a failure there blocks it.
4. Create and push a matching tag:

   ```bash
   git tag vX.Y.Z
   git push origin vX.Y.Z
   ```

5. The CI pipeline builds wheels and the sdist, tests them, and publishes to PyPI via
   Trusted Publishing. `make publish` / `make publish-test` upload whatever is in your
   local `dist/` and are for Test PyPI or emergencies only -- real releases go through
   the tag.
