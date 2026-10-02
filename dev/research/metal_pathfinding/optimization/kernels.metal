#include <metal_stdlib>
using namespace metal;

struct P {
  uint n, q, round, layout;
  uint current, epoch, lanes, delta;
};
constant long INF64 = 0x7fffffffffffffffL;
constant uint INF32 = 0xffffffffu;

// Same relaxation as phase one, with optional query-fast layout and 32-bit
// distances only after the host proves the entire finite-distance domain fits.
template<typename T>
void pull(device const int* row, device const int* col,
          device const T* weight, device const uchar* enabled,
          device const T* a, device T* b, device atomic_uint* changed,
          constant P& p, uint tid, uint lane, T inf) {
  uint total = p.n * p.q;
  bool dirty = false;
  if (tid < total) {
    uint v = p.layout ? tid / p.q : tid % p.n;
    uint query = p.layout ? tid % p.q : tid / p.n;
    T best = a[tid];
    for (int j = row[v]; j < row[v + 1]; ++j) {
      uint index = p.layout ? uint(col[j]) * p.q + query : query * p.n + uint(col[j]);
      T d = a[index], w = weight[j];
      if (enabled[j] && d != inf && d <= inf - w) best = min(best, T(d + w));
    }
    b[tid] = best;
    dirty = best != a[tid];
  }
  // One flag operation per SIMD group, instead of per changed vertex.
  if (simd_any(dirty) && lane == 0)
    atomic_store_explicit(changed + p.round, 1u, memory_order_relaxed);
}

#define PULL_ARGS(T) device const int* row [[buffer(0)]], device const int* col [[buffer(1)]], \
 device const T* weight [[buffer(2)]], device const uchar* enabled [[buffer(3)]], \
 device const T* a [[buffer(4)]], device T* b [[buffer(5)]], \
 device atomic_uint* changed [[buffer(6)]], constant P& p [[buffer(7)]], \
 uint tid [[thread_position_in_grid]], uint lane [[thread_index_in_simdgroup]]
kernel void pull64(PULL_ARGS(long)) { pull(row,col,weight,enabled,a,b,changed,p,tid,lane,INF64); }
kernel void pull32(PULL_ARGS(uint)) { pull(row,col,weight,enabled,a,b,changed,p,tid,lane,INF32); }

// Frontier queues contain flattened (query * N + node) IDs. Only 32-bit integer
// atomics are required. All queue consumers run after producer dispatch barriers.
kernel void prepare(device atomic_uint* counts [[buffer(0)]],
                    device uint* indirect [[buffer(1)]],
                    constant P& p [[buffer(7)]]) {
  uint count = atomic_load_explicit(counts + p.current, memory_order_relaxed);
  atomic_store_explicit(counts + (1 - p.current), 0u, memory_order_relaxed);
  indirect[0] = max(1u, (count * p.lanes + 255) / 256);
  indirect[1] = 1; indirect[2] = 1;
}

kernel void expand(device const int* row [[buffer(0)]],
                   device const int* col [[buffer(1)]],
                   device const uint* weight [[buffer(2)]],
                   device const uchar* enabled [[buffer(3)]],
                   device atomic_uint* dist [[buffer(4)]],
                   device const uint* in_queue [[buffer(5)]],
                   device uint* out_queue [[buffer(6)]],
                   constant P& p [[buffer(7)]],
                   device atomic_uint* counts [[buffer(8)]],
                   device atomic_uint* marks [[buffer(9)]],
                   device atomic_uint* work [[buffer(10)]],
                   uint tid [[thread_position_in_grid]]) {
  uint slot = tid / p.lanes, lane = tid % p.lanes;
  uint count = atomic_load_explicit(counts + p.current, memory_order_relaxed);
  if (slot >= count) return;
  uint index = in_queue[slot], u = index % p.n, base = index - u;
  uint d = atomic_load_explicit(dist + index, memory_order_relaxed);
  if (lane == 0) atomic_fetch_add_explicit(work + p.round, uint(row[u+1]-row[u]), memory_order_relaxed);
  for (uint j = uint(row[u]) + lane; j < uint(row[u+1]); j += p.lanes) {
    uint w = weight[j];
    if (!enabled[j] || d == INF32 || d > INF32 - w) continue;
    uint target = base + uint(col[j]);
    uint candidate = d + w;
    uint old = atomic_fetch_min_explicit(dist + target, candidate, memory_order_relaxed);
    if (candidate < old) {
      if (p.delta) {
        // Bucket mode: mark pending. Selection and clearing happen in a later
        // dispatch, so a late improvement cannot be lost behind an early reset.
        atomic_store_explicit(marks + target, 1u, memory_order_relaxed);
      } else if (atomic_exchange_explicit(marks + target, p.epoch, memory_order_relaxed) != p.epoch) {
        uint pos = atomic_fetch_add_explicit(counts + (1-p.current), 1u, memory_order_relaxed);
        out_queue[pos] = target; // Epoch deduplication bounds this queue by N*Q.
      }
    }
  }
}

// A simple bucketed label-correcting scheduler. This deliberately does not claim
// the optimized light/heavy split or asynchronous scheduler of published ADDS.
kernel void bucket_reset(device atomic_uint* counts [[buffer(0)]],
                         device atomic_uint* min_distance [[buffer(1)]]) {
  atomic_store_explicit(counts, 0u, memory_order_relaxed);
  atomic_store_explicit(min_distance, INF32, memory_order_relaxed);
}
kernel void bucket_min(device const atomic_uint* dist [[buffer(0)]],
                       device const atomic_uint* pending [[buffer(1)]],
                       device atomic_uint* minimum [[buffer(2)]],
                       constant P& p [[buffer(7)]], uint tid [[thread_position_in_grid]]) {
  if (tid >= p.n*p.q) return;
  if (atomic_load_explicit(pending+tid, memory_order_relaxed))
    atomic_fetch_min_explicit(minimum, atomic_load_explicit(dist+tid, memory_order_relaxed), memory_order_relaxed);
}
kernel void bucket_select(device const atomic_uint* dist [[buffer(0)]],
                          device atomic_uint* pending [[buffer(1)]],
                          device const atomic_uint* minimum [[buffer(2)]],
                          device uint* queue [[buffer(3)]],
                          device atomic_uint* counts [[buffer(4)]],
                          constant P& p [[buffer(7)]], uint tid [[thread_position_in_grid]]) {
  if (tid >= p.n*p.q) return;
  uint lowest = atomic_load_explicit(minimum, memory_order_relaxed);
  if (lowest == INF32) return;
  ulong upper = (ulong(lowest) / p.delta + 1) * p.delta;
  if (atomic_load_explicit(pending+tid, memory_order_relaxed) &&
      ulong(atomic_load_explicit(dist+tid, memory_order_relaxed)) < upper) {
    atomic_store_explicit(pending+tid, 0u, memory_order_relaxed);
    uint pos = atomic_fetch_add_explicit(counts, 1u, memory_order_relaxed);
    queue[pos] = tid;
  }
}

// Counts and fills are deterministic within each node's incoming CSR row.
// The host prefix-sums N+1 counts per query; it never scans edges for GPU DAGs.
struct DP { uint n, q, edges, layout, reverse, unused0, unused1, unused2; };
template<typename T, bool Fill>
void dag(device const int* row, device const int* col, device const long* costs,
         device const int* ids, device const uchar* enabled, device const T* distances,
         device int* offsets, device int* parents, device int* via,
         constant DP& p, uint tid, T inf) {
  if (tid >= p.n*p.q) return;
  // Match distance layout. Query-major frontiers must not be read with a
  // query-fast thread layout (which also scatters the output writes).
  uint v = p.layout ? tid / p.q : tid % p.n;
  uint query = p.layout ? tid % p.q : tid / p.n;
  uint index = p.layout ? v*p.q+query : query*p.n+v;
  T dv = distances[index];
  uint count = 0;
  for (int j=row[v]; j<row[v+1]; ++j) {
    uint u=uint(col[j]);
    T du=distances[p.layout ? u*p.q+query : query*p.n+u];
    long w=costs[j];
    bool tight=enabled[j] && (p.reverse ? (dv!=inf && long(du)==long(dv)+w) :
                                                          (du!=inf && long(dv)==long(du)+w));
    if (tight) {
      if (Fill) {
        uint pos=query*p.edges + uint(offsets[query*(p.n+1)+v]) + count;
        parents[pos]=int(u); via[pos]=ids[j];
      }
      ++count;
    }
  }
  if (!Fill) offsets[query*(p.n+1)+v+1]=int(count);
}
#define DAG_ARGS(T) device const int* row [[buffer(0)]], device const int* col [[buffer(1)]], \
 device const long* costs [[buffer(2)]], device const int* ids [[buffer(3)]], \
 device const uchar* enabled [[buffer(4)]], device const T* distances [[buffer(5)]], \
 device int* offsets [[buffer(6)]], constant DP& p [[buffer(7)]], \
 device int* parents [[buffer(8)]], device int* via [[buffer(9)]], uint tid [[thread_position_in_grid]]
kernel void dag_count64(DAG_ARGS(long)) {dag<long,false>(row,col,costs,ids,enabled,distances,offsets,parents,via,p,tid,INF64);}
kernel void dag_fill64(DAG_ARGS(long)) {dag<long,true>(row,col,costs,ids,enabled,distances,offsets,parents,via,p,tid,INF64);}
kernel void dag_count32(DAG_ARGS(uint)) {dag<uint,false>(row,col,costs,ids,enabled,distances,offsets,parents,via,p,tid,INF32);}
kernel void dag_fill32(DAG_ARGS(uint)) {dag<uint,true>(row,col,costs,ids,enabled,distances,offsets,parents,via,p,tid,INF32);}
