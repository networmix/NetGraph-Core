#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/../.."
mkdir -p build/metal_pathfinding
flags=(-std=c++20 -O3 -flto -funroll-loops -fno-math-errno -fno-trapping-math)
if [[ ${1:-} == sanitize ]]; then
  flags=(-std=c++20 -O1 -g -fsanitize=address,undefined -fno-sanitize-recover=undefined -fno-omit-frame-pointer)
fi
xcrun clang++ "${flags[@]}" -fobjc-arc -Iinclude \
  dev/research/metal_pathfinding/experiment.mm src/strict_multidigraph.cpp \
  src/shortest_paths.cpp src/profiling.cpp -framework Foundation -framework Metal \
  -o "build/metal_pathfinding/experiment${1:+-$1}"
xcrun clang++ -std=c++20 -O3 -fobjc-arc dev/research/metal_pathfinding/capabilities.mm \
  -framework Foundation -framework Metal -o build/metal_pathfinding/capabilities
