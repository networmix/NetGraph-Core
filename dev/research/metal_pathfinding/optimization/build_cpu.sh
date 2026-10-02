#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/../../.."
flags=(-std=c++20 -O3 -flto -funroll-loops -fno-math-errno -fno-trapping-math)
if [[ ${1:-} == sanitize ]]; then
  flags=(-std=c++20 -O1 -g -fsanitize=address,undefined -fno-sanitize-recover=undefined -fno-omit-frame-pointer)
fi
xcrun clang++ "${flags[@]}" -fobjc-arc -Iinclude \
  dev/research/metal_pathfinding/optimization/cpu_algorithms.mm src/strict_multidigraph.cpp \
  src/shortest_paths.cpp src/profiling.cpp -framework Foundation -framework Metal \
  -o "build/metal_pathfinding/cpu-algorithms${1:+-$1}"
