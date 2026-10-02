#!/usr/bin/env bash
# Fresh process per variant; production reference first and last (A1 ... A2).
set -euo pipefail
cd "$(dirname "$0")/../../.."
out=${1:-dev/perf/spf_queue/results/first}
samples=${2:-9}
mkdir -p "$out"
bin=build/perf/cpu-queue
{
  echo "date: $(date -u +%Y-%m-%dT%H:%M:%SZ)"; echo "commit: $(git rev-parse HEAD)"; sysctl -n machdep.cpu.brand_string
  shasum -a 256 dev/perf/spf_queue/*.hpp dev/perf/spf_queue/*.cpp src/shortest_paths.cpp src/strict_multidigraph.cpp
  "${CXX:-clang++}" --version | head -1
} > "$out/environment.txt"
for phase in ref:a1 heap-fresh:b heap-full:b heap-sparse:b bucket-fresh:b bucket-full:b bucket-sparse:b bsort-fresh:b bsort-sparse:b ref:a2; do
  v=${phase%%:*}; tag=${phase##*:}
  uptime >> "$out/load.txt"
  "$bin" time "$v" "$samples" > "$out/time-$v-$tag.csv"
done
for v in ref heap-sparse bucket-sparse bsort-sparse; do
  uptime >> "$out/load.txt"
  "$bin" sweep "$v" "$samples" 2097152 > "$out/sweep-$v.csv"
done
uptime >> "$out/load.txt"
echo done
