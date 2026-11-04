#!/usr/bin/env bash
# Explicit Lean checks. Normal Python/CMake builds never call this script.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/../.."
proof_root=$PWD
proof_python=${PROOFS_PYTHON:-python3}
proof_python=$("$proof_python" -I -c 'import sys; print(sys.executable)')
proof_mode=${1:-}
proof_output=${PROOFS_OUTPUT:-"$proof_root/build/proofs"}
mkdir -p "$proof_output"
proof_output=$(cd "$proof_output" && pwd)
export PROOFS_PYTHON="$proof_python"

check() {
    "$proof_python" -I -B proofs/Checks/audit.py --preflight
    command -v lake >/dev/null || { echo 'Install the toolchain in proofs/lean-toolchain and put lake on PATH.' >&2; return 1; }
    rm -f "$proof_output/check.json"
    (cd proofs && lake --wfail build && lake exe proof-audit --output "$proof_output/check.json") \
        > "$proof_output/check.log" 2>&1 \
        || { tail -60 "$proof_output/check.log" >&2; return 1; }
    "$proof_python" -I -B proofs/Checks/audit.py --verify-report "$proof_output/check.json"
    echo "Lean model check: $proof_output/check.json"
}

case "$proof_mode" in
    check) check ;;
    replay)
        check
        (cd proofs && lake env leanchecker --fresh Checks) > "$proof_output/replay.log" 2>&1 \
            || { tail -60 "$proof_output/replay.log" >&2; exit 1; }
        "$proof_python" -I -B proofs/Checks/audit.py --verify-report "$proof_output/check.json"
        "$proof_python" -I -B dev/proofs/replay_report.py "$proof_output"
        ;;
    test)
        "$proof_python" -I -B dev/proofs/test_tool.py "${@:2}"
        ;;
    bridge)
        : "${PROOFS_NATIVE_PYTHON:?Set PROOFS_NATIVE_PYTHON to the interpreter with the selected wheel installed}"
        : "${PROOFS_WHEEL:?Set PROOFS_WHEEL to that exact wheel file}"
        "$proof_python" -I -B proofs/Checks/audit.py --verify-report "$proof_output/check.json"
        "$PROOFS_NATIVE_PYTHON" -I -B dev/proofs/bridge.py \
            --proof-report "$proof_output/check.json" --wheel "$PROOFS_WHEEL" \
            --output "$proof_output/bridge.json"
        ;;
    *) echo 'Usage: bash dev/proofs/run.sh check|replay|test [--preflight-only]|bridge' >&2; exit 2 ;;
esac
