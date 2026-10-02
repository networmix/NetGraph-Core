#!/usr/bin/env bash
# Core-level old/new/old A/B/A on the harness workloads: the `ref` variant is the
# production shortest_paths of whichever sources the binary was built from.
#   cpu-queue-old: previous sources (git archive of the base commit)
#   cpu-queue:     this worktree; B-heap forces the heap route in the new binary
set -euo pipefail
cd "$(dirname "$0")/../../.."
out=${1:-dev/perf/spf_queue/results/core}
samples=${2:-9}
mkdir -p "$out"
{
  echo "date: $(date -u +%Y-%m-%dT%H:%M:%SZ)"; echo "commit: $(git rev-parse HEAD)"
  shasum -a 256 build/perf/cpu-queue-old build/perf/cpu-queue src/shortest_paths.cpp
} > "$out/environment.txt"
for phase in a1-old b-new b-new-heap a2-old; do
  uptime >> "$out/load.txt"
  case $phase in
    *old) bin=build/perf/cpu-queue-old; envs=() ;;
    b-new) bin=build/perf/cpu-queue; envs=() ;;
    b-new-heap) bin=build/perf/cpu-queue; envs=(NGRAPH_CORE_SPF_QUEUE=heap) ;;
  esac
  env ${envs[@]+"${envs[@]}"} "$bin" time ref "$samples" > "$out/time-$phase.csv"
  env ${envs[@]+"${envs[@]}"} "$bin" sweep ref "$samples" 2097152 > "$out/sweep-$phase.csv"
done
uptime >> "$out/load.txt"
echo done
