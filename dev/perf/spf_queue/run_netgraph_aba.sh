#!/usr/bin/env bash
# Old/new/old A/B/A on realistic NetGraph workloads, one fresh process per phase.
#   A1: main worktree build (fbab119)   B: this worktree   B-heap: this worktree with
#   NGRAPH_CORE_SPF_QUEUE=heap (in-binary control)   A2: main worktree build again.
set -euo pipefail
cd "$(dirname "$0")/../../.."
out=${1:-dev/perf/spf_queue/results/netgraph}
samples=${2:-5}
old=${OLD_CORE:-/Users/networmix/ws/NetGraph-Core-main}
new=$PWD
py=/Users/networmix/ws/NetGraph/venv/bin/python
bench=dev/perf/spf_queue/bench_netgraph.py
mkdir -p "$out"
{
  echo "date: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
  echo "new: $(git rev-parse HEAD) (dirty: $(git status --short | wc -l | tr -d ' ') paths)"
  echo "old: $(git -C "$old" rev-parse HEAD)"
  echo "netgraph: $(git -C /Users/networmix/ws/NetGraph rev-parse HEAD)"
  shasum -a 256 "$bench" "$new"/venv/lib/python3.14/site-packages/_netgraph_core*.so "$old"/venv/lib/python3.14/site-packages/_netgraph_core*.so
  sysctl -n machdep.cpu.brand_string
} > "$out/environment.txt"
run() {  # label core-root extra-env
  local label=$1 root=$2; shift 2
  uptime >> "$out/load.txt"
  # PYTHONHASHSEED pins NetGraph's set-iteration order so Monte Carlo result
  # digests are comparable across processes (verified: with a random seed the
  # digests differ between two runs of the same build).
  env "$@" PYTHONHASHSEED=0 PYTHONPATH="$root/venv/lib/python3.14/site-packages:$root/python" "$py" "$bench" "$label" "$samples" > "$out/$label.csv"
}
run a1-old "$old"
run b-new "$new"
run b-new-heap "$new" NGRAPH_CORE_SPF_QUEUE=heap
run a2-old "$old"
uptime >> "$out/load.txt"
echo done
