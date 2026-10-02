#!/usr/bin/env bash
# Serial fresh-process CPU / Metal / CPU runs. Do not overlap with builds/tests.
set -euo pipefail
cd "$(dirname "$0")/../.."
out=${1:-dev/research/metal_pathfinding/results}
mkdir -p "$out"
unset NGRAPH_CORE_PROFILE
{
  date -u '+%Y-%m-%dT%H:%M:%SZ'
  git rev-parse HEAD
  sw_vers
  sysctl -n machdep.cpu.brand_string hw.memsize hw.ncpu
  xcrun clang++ --version
  shasum -a 256 dev/research/metal_pathfinding/experiment.mm dev/research/metal_pathfinding/sssp.metal
  pmset -g therm
} > "$out/environment.txt"
for phase in a1 b a2; do
  top -l 2 -n 8 -o cpu -stats pid,command,cpu,threads > "$out/load-$phase.txt"
  ioreg -r -c AGXAccelerator -l | sed -n '/"PerformanceStatistics"/p' > "$out/gpu-load-$phase.txt"
  mode=cpu
  [[ $phase != b ]] || mode=gpu
  build/metal_pathfinding/experiment "$mode" dev/research/metal_pathfinding/sssp.metal full 9 \
    > "$out/$phase.csv" 2> "$out/$phase-stderr.txt"
done
