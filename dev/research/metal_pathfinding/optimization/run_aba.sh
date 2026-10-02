#!/usr/bin/env bash
# Run sequentially on the same machine. Do not overlap with builds or tests.
set -euo pipefail
cd "$(dirname "$0")/../../.."
out=${1:-dev/research/metal_pathfinding/optimization/results/aba}
mkdir -p "$out"
unset NGRAPH_CORE_PROFILE
{
  date -u '+%Y-%m-%dT%H:%M:%SZ'
  git rev-parse HEAD
  sw_vers
  sysctl -n machdep.cpu.brand_string hw.memsize hw.ncpu
  xcrun clang++ --version
  shasum -a 256 dev/research/metal_pathfinding/experiment.mm dev/research/metal_pathfinding/sssp.metal \
    dev/research/metal_pathfinding/optimization/experiment.mm dev/research/metal_pathfinding/optimization/kernels.metal \
    build/metal_pathfinding/optimized
  pmset -g therm
} > "$out/environment.txt"
for phase in cpu-a1 base-a1 optimized-b base-a2 cpu-a2; do
  top -l 2 -n 8 -o cpu -stats pid,command,cpu,threads > "$out/load-$phase.txt"
  ioreg -r -c AGXAccelerator -l | sed -n '/"PerformanceStatistics"/p' > "$out/gpu-load-$phase.txt"
  case "$phase" in
    cpu-*) mode=cpu ;;
    base-*) mode=base ;;
    optimized-b) mode=optimized ;;
  esac
  date -u '+%Y-%m-%dT%H:%M:%SZ' > "$out/time-$phase.txt"
  build/metal_pathfinding/optimized "$mode" full 7 > "$out/$phase.csv" 2> "$out/$phase-stderr.txt"
  date -u '+%Y-%m-%dT%H:%M:%SZ' >> "$out/time-$phase.txt"
done
