#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."
mode=${1:-help}
case "$mode" in
    setup|check) ;;
    *) echo 'Usage: bash .superset/workspace.sh {setup|check}'; exit 2 ;;
esac
if [[ "$mode" == setup && ! -x venv/bin/python ]]; then
    make venv
fi
[[ -x venv/bin/python ]] || { echo 'Run bash .superset/workspace.sh setup first.' >&2; exit 1; }
export VIRTUAL_ENV="$PWD/venv"
export PATH="$VIRTUAL_ENV/bin:$PATH"
unset PYTHONPATH PYTHONHOME
export CMAKE_BUILD_PARALLEL_LEVEL="${CMAKE_BUILD_PARALLEL_LEVEL:-4}"
if [[ "$(uname -s)" == Darwin ]]; then
    export CC="$(xcrun --find clang)"
    export CXX="$(xcrun --find clang++)"
    export MACOSX_DEPLOYMENT_TARGET="${MACOSX_DEPLOYMENT_TARGET:-15.0}"
fi
if [[ "$mode" == setup ]]; then
    python -m pip install -e '.[dev]'
    python -m pip check
    python -c 'import netgraph_core, _netgraph_core; print("Workspace environment ready")'
else
    make install
    make check-ci
fi
