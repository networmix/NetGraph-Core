#include <metal_stdlib>
using namespace metal;

// Research kernels: exact int64, nonnegative costs, one writer per distance.
// Incoming CSR includes a host-computed edge eligibility byte, so double
// capacity comparisons keep their CPU semantics.
constant long INF = 0x7fffffffffffffffL;
struct Params { uint n; uint queries; uint round; uint unused; };

kernel void relax_pull(device const int* row [[buffer(0)]],
                       device const int* col [[buffer(1)]],
                       device const long* cost [[buffer(2)]],
                       device const uchar* allowed [[buffer(3)]],
                       device const long* before [[buffer(4)]],
                       device long* after [[buffer(5)]],
                       device atomic_uint* changed [[buffer(6)]],
                       constant Params& p [[buffer(7)]],
                       uint tid [[thread_position_in_grid]]) {
  if (tid >= p.n * p.queries) return;
  uint v = tid % p.n;
  uint base = tid - v;
  long best = before[tid];
  for (int j = row[v]; j < row[v + 1]; ++j) {
    long d = before[base + col[j]];
    if (allowed[j] && d != INF) best = min(best, d + cost[j]);
  }
  after[tid] = best;
  if (best != before[tid])
    atomic_store_explicit(changed + p.round, 1u, memory_order_relaxed);
}

// A whole SSSP in one threadgroup, with independent queries across groups.
// Each group owns both of its N-element distance arrays. Device-memory barriers
// are local to that group; no unsafe grid-wide barrier or in-place data race.
kernel void solve_group(device const int* row [[buffer(0)]],
                        device const int* col [[buffer(1)]],
                        device const long* cost [[buffer(2)]],
                        device const uchar* allowed [[buffer(3)]],
                        device long* a [[buffer(4)]],
                        device long* b [[buffer(5)]],
                        device uint* rounds [[buffer(6)]],
                        constant Params& p [[buffer(7)]],
                        device const int* sources [[buffer(8)]],
                        uint lane [[thread_index_in_threadgroup]],
                        uint width [[threads_per_threadgroup]],
                        uint query [[threadgroup_position_in_grid]]) {
  uint base = query * p.n;
  for (uint v = lane; v < p.n; v += width)
    a[base + v] = (int(v) == sources[query]) ? 0L : INF;
  threadgroup atomic_uint changed;
  threadgroup_barrier(mem_flags::mem_device);
  uint step = 0;
  for (; step < p.n; ++step) {
    if (lane == 0) atomic_store_explicit(&changed, 0u, memory_order_relaxed);
    threadgroup_barrier(mem_flags::mem_threadgroup);
    for (uint v = lane; v < p.n; v += width) {
      long best = a[base + v];
      for (int j = row[v]; j < row[v + 1]; ++j) {
        long d = a[base + col[j]];
        if (allowed[j] && d != INF) best = min(best, d + cost[j]);
      }
      b[base + v] = best;
      if (best != a[base + v])
        atomic_store_explicit(&changed, 1u, memory_order_relaxed);
    }
    threadgroup_barrier(mem_flags::mem_device | mem_flags::mem_threadgroup);
    uint any_changed = atomic_load_explicit(&changed, memory_order_relaxed);
    // Every lane must observe this round before lane 0 resets the flag for the
    // next one. A barrier *before* the loads alone does not prevent that race.
    threadgroup_barrier(mem_flags::mem_threadgroup);
    if (any_changed == 0) break;
    device long* tmp = a; a = b; b = tmp;
  }
  // Host selects the output buffer by round parity; convergence leaves a == b.
  if (lane == 0) rounds[query] = min(step + 1, p.n);
}

// Small-graph variant: 2 * N * 8 bytes of threadgroup memory, N <= 2048 on
// this device. Removes repeated device-memory accesses and barriers per round.
kernel void solve_shared(device const int* row [[buffer(0)]],
                         device const int* col [[buffer(1)]],
                         device const long* cost [[buffer(2)]],
                         device const uchar* allowed [[buffer(3)]],
                         device long* output [[buffer(4)]],
                         device uint* rounds [[buffer(6)]],
                         constant Params& p [[buffer(7)]],
                         device const int* sources [[buffer(8)]],
                         threadgroup long* scratch [[threadgroup(0)]],
                         uint lane [[thread_index_in_threadgroup]],
                         uint width [[threads_per_threadgroup]],
                         uint query [[threadgroup_position_in_grid]]) {
  threadgroup long* a = scratch;
  threadgroup long* b = scratch + p.n;
  for (uint v = lane; v < p.n; v += width)
    a[v] = (int(v) == sources[query]) ? 0L : INF;
  threadgroup atomic_uint changed;
  threadgroup_barrier(mem_flags::mem_threadgroup);
  uint step = 0;
  for (; step < p.n; ++step) {
    if (lane == 0) atomic_store_explicit(&changed, 0u, memory_order_relaxed);
    threadgroup_barrier(mem_flags::mem_threadgroup);
    bool lane_changed = false;
    for (uint v = lane; v < p.n; v += width) {
      long best = a[v];
      for (int j = row[v]; j < row[v + 1]; ++j) {
        long d = a[col[j]];
        if (allowed[j] && d != INF) best = min(best, d + cost[j]);
      }
      b[v] = best;
      lane_changed |= best != a[v];
    }
    if (simd_any(lane_changed) && (lane % 32 == 0))
      atomic_store_explicit(&changed, 1u, memory_order_relaxed);
    threadgroup_barrier(mem_flags::mem_threadgroup);
    threadgroup long* tmp = a; a = b; b = tmp;
    uint any_changed = atomic_load_explicit(&changed, memory_order_relaxed);
    threadgroup_barrier(mem_flags::mem_threadgroup);
    if (any_changed == 0) break;
  }
  for (uint v = lane; v < p.n; v += width) output[query * p.n + v] = a[v];
  if (lane == 0) rounds[query] = min(step + 1, p.n);
}
