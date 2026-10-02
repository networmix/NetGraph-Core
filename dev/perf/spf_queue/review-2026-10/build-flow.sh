#!/bin/bash
set -eu
r=dev/perf/spf_queue/review-2026-10
flags=(-std=c++20 -O3 -funroll-loops -fno-math-errno -fno-trapping-math)
for v in old new fixed; do
 case $v in
  old) inc="$r/old-include"; spf="$r/old-shortest_paths.cpp"; graph="$r/old-strict_multidigraph.cpp"; flow=src/max_flow.cpp; state=src/flow_state.cpp;;
  new) inc=include; spf=src/shortest_paths.cpp; graph=src/strict_multidigraph.cpp; flow=src/max_flow.cpp; state=src/flow_state.cpp;;
  fixed) inc=include; spf=src/shortest_paths.cpp; graph=src/strict_multidigraph.cpp; flow="$r/max_flow-fixed.cpp"; state="$r/flow_state-fixed.cpp";;
 esac
 clang++ "${flags[@]}" -I"$inc" -Iinclude "$r/flow-driver.cpp" "$spf" "$graph" "$flow" "$state" src/profiling.cpp -o "build/perf/review-flow-$v"
done
for phase in old1 new1 fixed new2 old2; do
 case $phase in old*) v=old;;new*)v=new;;*)v=fixed;;esac
 uptime >> "$r/results/load-flow.txt"
 build/perf/review-flow-$v > "$r/results/flow-$phase.csv"
 echo "$phase done"
done
