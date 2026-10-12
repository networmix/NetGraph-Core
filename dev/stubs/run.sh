#!/usr/bin/env bash
# The generated type stub python/netgraph_core/__init__.pyi: regenerate it, or check it.
#
#   bash dev/stubs/run.sh generate
#       Build a wheel of this checkout (the shipped configuration, without the private hooks),
#       install it into a fresh venv and rewrite the stub from it (dev/stubs/generate.py:
#       pybind11-stubgen, the reviewed overrides of dev/stubs/overrides.pyi, ruff format).
#       Run after any change to a binding's signature or docstring, then commit the stub.
#   bash dev/stubs/run.sh check                                         (make stubs-check, CI)
#       0. the tools' rules (no build, a few seconds): dev/stubs/test_generate.py provokes each
#          refusal of generate.py once on a small synthetic stub and expects its message, and
#          runs its --write and --check on it; dev/stubs/test_checks.py provokes each refusal of
#          the contents check (step 1) on synthetic file lists, and each abnormal end of Pyright
#          and each comparison rule of the client check (step 3) on a fake Pyright;
#       1. contents: the sdist holds exactly the intended source files of this checkout (the
#          rule of pyproject.toml's sdist.include, applied to the files git tracks or would
#          track) and the wheel exactly the package, its extension and its metadata; both carry
#          the stub, py.typed and both licence notices, and no _docs.py
#          (dev/stubs/check_contents.py);
#       2. drift: the committed stub equals a fresh generation from the wheel;
#       3. clients, against the installed wheel and against an installation built from the
#          sdist: tests/typing/valid_client.py type-checks without a diagnostic and runs, and
#          every marked line of tests/typing/invalid_client.py gets exactly its expected Pyright
#          diagnostic (dev/stubs/check_clients.py).
#
# Everything happens in STUBS_DIR (default build/stubs, emptied first): its own venvs with the
# pinned tools of pyproject.toml's dev extra (pybind11-stubgen, pyright, ruff, build) and the
# pins of .github/ci-constraints.txt, its own build directory. The development venv and the
# shared build directory are not touched. PYTHON names the interpreter (default python3, 3.11
# or newer); CMAKE_BUILD_PARALLEL_LEVEL bounds the builds.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/../.."
ROOT=$PWD
PYTHON=${PYTHON:-python3}
mode=${1:-}
case "$mode" in
    generate|check) ;;
    *) echo 'Usage: bash dev/stubs/run.sh generate|check' >&2; exit 2 ;;
esac
die() { echo "dev/stubs/run.sh: $*" >&2; exit 1; }

"$PYTHON" -c 'import sys; sys.exit(sys.version_info < (3, 11))' \
    || die "$PYTHON is older than Python 3.11"
if [[ "$mode" == check ]]; then
    echo "== the tools' rules"
    "$PYTHON" -I -B dev/stubs/test_generate.py
    "$PYTHON" -I -B dev/stubs/test_checks.py
fi
if [[ "$(uname -s)" == Darwin ]]; then
    CC=$(xcrun --find clang)
    CXX=$(xcrun --find clang++)
    export CC CXX
    export MACOSX_DEPLOYMENT_TARGET="${MACOSX_DEPLOYMENT_TARGET:-15.0}"
fi
# The pins reach every pip call and every build environment: PIP_CONSTRAINT the direct installs
# and the environments of `python -m build` (which fills them with a direct install),
# PIP_BUILD_CONSTRAINT the environment pip makes itself to build the sdist (pip 26.2 and newer
# apply only build constraints there; a pip older than 25.3 ignores the variable and passes
# PIP_CONSTRAINT on).
export PIP_CONSTRAINT="$ROOT/.github/ci-constraints.txt"
export PIP_BUILD_CONSTRAINT="$ROOT/.github/ci-constraints.txt"
export PIP_DISABLE_PIP_VERSION_CHECK=1

WORK=${STUBS_DIR:-build/stubs}
rm -rf "$WORK"
mkdir -p "$WORK"
WORK=$(cd "$WORK" && pwd)

# The tools as pinned in pyproject.toml's dev extra.
TOOLS=$("$PYTHON" - <<'PY'
import re
import tomllib

with open("pyproject.toml", "rb") as f:
    dev = tomllib.load(f)["project"]["optional-dependencies"]["dev"]
want = {"pybind11-stubgen", "pyright", "ruff", "build"}
found = {re.split(r"[\s=<>!~;\[]", r, maxsplit=1)[0].lower(): r for r in dev}
missing = want - set(found)
if missing:
    raise SystemExit(f"pyproject.toml's dev extra lacks {sorted(missing)}")
# The shell splits the printed list on spaces.
spaced = [found[name] for name in sorted(want) if any(c.isspace() for c in found[name])]
if spaced:
    raise SystemExit(f"pyproject.toml's dev extra writes {spaced} with a space; remove it")
print(" ".join(found[name] for name in sorted(want)))
PY
)

# new_env NAME: a fresh venv with NumPy and the tools; prints its interpreter. It runs inside
# $(...), where bash turns set -e off, so every step returns its own failure: the caller's
# assignment then fails and stops the run.
new_env() {
    "$PYTHON" -m venv "$WORK/$1" || return
    "$WORK/$1/bin/python" -m pip install -q -U pip >&2 || return  # the pinned pip (PIP_CONSTRAINT)
    # shellcheck disable=SC2086  # TOOLS: requirement strings without spaces (checked above)
    "$WORK/$1/bin/python" -m pip install -q numpy $TOOLS >&2 || return
    echo "$WORK/$1/bin/python"
}

echo "== tools: $TOOLS"
WPY=$(new_env wheel-env)
echo "== the wheel of this checkout (hooks off; build directory $WORK/wheel-build)"
"$WPY" -m build --wheel --outdir "$WORK/dist" --config-setting=build-dir="$WORK/wheel-build" . \
    > "$WORK/build-wheel.log" 2>&1 || { tail -40 "$WORK/build-wheel.log" >&2; die "the wheel build failed"; }
wheels=("$WORK"/dist/*.whl)
[[ ${#wheels[@]} -eq 1 && -f ${wheels[0]} ]] || die "expected one wheel in $WORK/dist"
"$WPY" -m pip install -q --no-index --no-deps "${wheels[0]}"

if [[ "$mode" == generate ]]; then
    "$WPY" -I dev/stubs/generate.py --write --raw "$WORK/raw"
    echo "Regenerated; check it with: bash dev/stubs/run.sh check"
    exit 0
fi

echo "== the sdist"
"$WPY" -m build --sdist --outdir "$WORK/dist" . > "$WORK/build-sdist.log" 2>&1 \
    || { tail -40 "$WORK/build-sdist.log" >&2; die "the sdist build failed"; }
sdists=("$WORK"/dist/*.tar.gz)
[[ ${#sdists[@]} -eq 1 && -f ${sdists[0]} ]] || die "expected one sdist in $WORK/dist"

echo "== contents"
"$WPY" -I dev/stubs/check_contents.py "${wheels[0]}" "${sdists[0]}" "$WORK"

echo "== drift"
"$WPY" -I dev/stubs/generate.py --check --raw "$WORK/raw"

echo "== clients against the installed wheel"
"$WPY" -I dev/stubs/check_clients.py wheel "$WORK/clients-wheel"

echo "== clients against an installation built from the sdist"
SPY=$(new_env sdist-env)
"$SPY" -m pip install -q --no-deps "${sdists[0]}" > "$WORK/install-sdist.log" 2>&1 \
    || { tail -40 "$WORK/install-sdist.log" >&2; die "installing the sdist failed"; }
"$SPY" -I dev/stubs/check_clients.py sdist "$WORK/clients-sdist"
echo "Stub, package contents and typed clients: OK"
