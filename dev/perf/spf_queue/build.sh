#!/usr/bin/env bash
# Plain C++ build of the parity/timing harness (no Objective-C++/Metal; any
# clang/gcc with C++20).
#   bash build.sh            -> build/perf/cpu-queue           (working tree sources)
#   bash build.sh sanitize   -> build/perf/cpu-queue-sanitize  (ASan/UBSan)
#   bash build.sh old [REV]  -> build/perf/cpu-queue-old       (sources of REV, default fbab119,
#                               extracted with git archive into build/perf/old-src)
set -euo pipefail
cd "$(dirname "$0")/../../.."
mode=${1:-release}
flags=(-std=c++20 -O3 -funroll-loops -fno-math-errno -fno-trapping-math -Wall -Wextra)
src=.
out=build/perf/cpu-queue
case "$mode" in
  release) ;;
  sanitize)
    flags=(-std=c++20 -O1 -g -fsanitize=address,undefined -fno-sanitize-recover=undefined -fno-omit-frame-pointer)
    out=build/perf/cpu-queue-sanitize ;;
  old)
    rev=${2:-fbab119}
    src=build/perf/old-src
    rm -rf "$src" && mkdir -p "$src"
    git archive "$rev" src include | tar -x -C "$src"
    out=build/perf/cpu-queue-old ;;
  *) echo "usage: build.sh [release|sanitize|old [REV]]" >&2; exit 2 ;;
esac
mkdir -p build/perf
${CXX:-clang++} "${flags[@]}" -I"$src/include" \
  dev/perf/spf_queue/harness.cpp \
  "$src/src/strict_multidigraph.cpp" "$src/src/shortest_paths.cpp" "$src/src/profiling.cpp" -o "$out"
echo "built $out"
