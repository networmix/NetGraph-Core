#!/bin/bash
set -eu
cd "$(dirname "$0")/../../../../.."
r=dev/perf/spf_queue/review-2026-10
flags=(-std=c++20 -O3 -funroll-loops -fno-math-errno -fno-trapping-math)
mkdir -p "$r/old-include/netgraph/core"
for f in shortest_paths strict_multidigraph; do
  git show "fbab119:include/netgraph/core/$f.hpp" > "$r/old-include/netgraph/core/$f.hpp"
  git show "fbab119:src/$f.cpp" > "$r/old-$f.cpp"
done
clang++ "${flags[@]}" -I"$r/old-include" -Iinclude "$r/driver.cpp" "$r/old-shortest_paths.cpp" "$r/old-strict_multidigraph.cpp" src/profiling.cpp -o build/perf/review-old
for v in baseline settled singleton heapitem notouch inline append frontpop combined; do
  clang++ "${flags[@]}" -Iinclude "$r/driver.cpp" "$r/$v.cpp" src/strict_multidigraph.cpp src/profiling.cpp -o "build/perf/review-$v"
  echo "built $v"
done
