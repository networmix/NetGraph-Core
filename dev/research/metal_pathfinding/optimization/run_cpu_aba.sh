#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/../../.."
out=dev/research/metal_pathfinding/optimization/results/cpu-algorithms-aba
mkdir -p "$out"
unset NGRAPH_CORE_PROFILE
{
  date -u '+%Y-%m-%dT%H:%M:%SZ'
  git rev-parse HEAD
  shasum -a 256 dev/research/metal_pathfinding/optimization/cpu_algorithms.mm \
    build/metal_pathfinding/cpu-algorithms build/metal_pathfinding/optimized
} > "$out/environment.txt"
for phase in a1 b a2; do
  top -l 2 -n 8 -o cpu -stats pid,command,cpu,threads > "$out/load-$phase.txt"
  ioreg -r -c AGXAccelerator -l | sed -n '/"PerformanceStatistics"/p' > "$out/gpu-load-$phase.txt"
  if [[ $phase == b ]]; then
    build/metal_pathfinding/optimized optimized full 7 > "$out/$phase.csv" 2> "$out/$phase-stderr.txt"
  else
    build/metal_pathfinding/cpu-algorithms 7 > "$out/$phase.csv" 2> "$out/$phase-stderr.txt"
  fi
done
