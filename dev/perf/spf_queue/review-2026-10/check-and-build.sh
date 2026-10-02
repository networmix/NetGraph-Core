#!/bin/bash
set -eu
r=dev/perf/spf_queue/review-2026-10
for v in baseline settled inline; do
  build/perf/review-$v correctness > "$r/results/correctness-$v.txt"
  echo "correctness $v passed"
done
clang++ -std=c++20 -O3 -funroll-loops -fno-math-errno -fno-trapping-math -Iinclude "$r/driver.cpp" "$r/append.cpp" src/strict_multidigraph.cpp src/profiling.cpp -o build/perf/review-append
for v in old baseline gcd; do
  case $v in
    old) graph="$r/old-strict_multidigraph.cpp"; inc="$r/old-include";;
    baseline) graph=src/strict_multidigraph.cpp; inc=include;;
    gcd) graph="$r/gcd.cpp"; inc=include;;
  esac
  clang++ -std=c++20 -O3 -funroll-loops -fno-math-errno -fno-trapping-math -I"$inc" -Iinclude "$r/construction.cpp" "$graph" -o build/perf/review-build-$v
done
