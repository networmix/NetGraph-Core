#!/bin/bash
set -euo pipefail
TASK=$(cd "$(dirname "$0")" && pwd)
EFFORT=$(cd "$TASK/../.." && pwd)
SCRATCH="$EFFORT/scratch/split"
export CC CXX MACOSX_DEPLOYMENT_TARGET=15.0 CMAKE_BUILD_PARALLEL_LEVEL=4
CC=$(xcrun --find clang)
CXX=$(xcrun --find clang++)
export PIP_NO_INDEX=1 PYTHONDONTWRITEBYTECODE=1
"$EFFORT/scratch/consumer/ref/venv/bin/python" -B "$TASK/prepare_build.py"
"$EFFORT/scratch/wheels/buildenv/bin/python" -B -m pip wheel \
  --no-deps --no-index --no-build-isolation --wheel-dir "$SCRATCH/wheels" "$SCRATCH/core"
"$SCRATCH/venv/bin/python" -B -m pip install --force-reinstall \
  --no-deps --no-index "$SCRATCH"/wheels/netgraph_core-*.whl
"$SCRATCH/venv/bin/python" -B -c \
  'import _netgraph_core as c,sys,pathlib,hashlib,json; p=pathlib.Path(c.__file__).resolve(); assert p.is_relative_to(pathlib.Path(sys.prefix).resolve()); print(json.dumps({"native":str(p),"sha256":hashlib.sha256(p.read_bytes()).hexdigest(),"snapshot":c._split_snapshot(),"summary":"PASS: scratch wheel import and counters"}))'
