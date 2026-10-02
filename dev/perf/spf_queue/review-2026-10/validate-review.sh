#!/bin/bash
set -eu
r=dev/perf/spf_queue/review-2026-10
out=build/perf/review-validation
mkdir -p "$out"
flags=(-std=c++20 -O1 -g -fsanitize=address,undefined -fno-sanitize-recover=undefined -fno-omit-frame-pointer -Iinclude)
for file in src/strict_multidigraph.cpp "$r/combined.cpp" "$r/max_flow-fixed.cpp" "$r/flow_state-fixed.cpp" src/flow_graph.cpp src/flow_policy.cpp src/k_shortest_paths.cpp src/cpu_backend.cpp src/profiling.cpp; do
 clang++ "${flags[@]}" -c "$file" -o "$out/$(basename "$file" .cpp).o"
done
for spec in validate-flow bitgate; do
 if [[ $spec == bitgate ]]; then file=dev/perf/spf_queue/bitgate/dump.cpp; else file="$r/validate-flow.cpp"; fi
 clang++ "${flags[@]}" "$file" "$out"/*.o -o "$out/$spec"
 ASAN_OPTIONS=detect_leaks=0 NGRAPH_CORE_BATCH_THREADS=1 NGRAPH_CORE_SENSITIVITY_THREADS=1 "$out/$spec" > "$r/results/$spec-variant.txt"
done
NGRAPH_CORE_BATCH_THREADS=1 NGRAPH_CORE_SENSITIVITY_THREADS=1 build/perf/bitgate-old > "$r/results/bitgate-old.txt"
cmp "$r/results/bitgate-old.txt" "$r/results/bitgate-variant.txt"
shasum -a 256 "$r/results/bitgate-old.txt" "$r/results/bitgate-variant.txt" > "$r/results/bitgate-sha256.txt"
echo 'targeted flow certificates, cost factors, and bitgate passed under ASan/UBSan'
