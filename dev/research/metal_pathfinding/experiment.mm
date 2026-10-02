// Standalone research harness, not a registered NetGraph backend.
#import <Foundation/Foundation.h>
#import <Metal/Metal.h>

#include "netgraph/core/shortest_paths.hpp"
#include "netgraph/core/constants.hpp"
#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <condition_variable>
#include <cstring>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <mutex>
#include <numeric>
#include <queue>
#include <random>
#include <stdexcept>
#include <string>
#include <thread>

using namespace netgraph::core;
using Result = std::pair<std::vector<Cost>, PredDAG>;
using Clock = std::chrono::steady_clock;
constexpr Cost INF = std::numeric_limits<Cost>::max();
double ms(Clock::time_point start) {
  return std::chrono::duration<double, std::milli>(Clock::now() - start).count();
}

struct Case {
  std::string name;
  StrictMultiDiGraph graph;
  std::unique_ptr<bool[]> node_mask, edge_mask;
  std::vector<Cap> residual;
  bool reverse = false;
  std::span<const bool> nm() const { return {node_mask.get(), size_t(graph.num_nodes())}; }
  std::span<const bool> em() const { return {edge_mask.get(), size_t(graph.num_edges())}; }
};

Case make_case(std::string kind, int n, int degree = 8, bool masked = false,
               bool large = false, bool reverse = false, unsigned seed = 1729) {
  std::vector<NodeId> src, dst;
  std::vector<Cost> costs;
  std::vector<Cap> caps;
  std::mt19937 rng(seed);
  auto add = [&](int u, int v, Cost c) {
    src.push_back(u); dst.push_back(v);
    costs.push_back(c + (large ? (Cost(1) << 34) : 0));
    caps.push_back(1.0 + double(rng() % 100) / 7.0);
  };
  if (kind == "chain") {
    for (int i = 0; i < n - 1; ++i) add(i, i + 1, 1 + rng() % 31);
  } else if (kind == "grid") {
    int side = int(std::sqrt(n));
    if (side * side != n) throw std::runtime_error("grid must be square");
    for (int u = 0; u < n; ++u) {
      if (u % side + 1 < side) { add(u, u + 1, 1); add(u + 1, u, 1); }
      if (u + side < n) { add(u, u + side, 1); add(u + side, u, 1); }
    }
  } else if (kind == "fabric" || kind == "clos") {
    int spines = kind == "clos" ? n / 2 : 32;
    for (int u = 0; u < n - spines; ++u) for (int v = n - spines; v < n; ++v) {
      add(u, v, 1); add(v, u, 1);
    }
  } else {
    for (int u = 0; u < n; ++u) {
      add(u, (u + 1) % n, 1 + rng() % 31);
      for (int k = 1; k < degree; ++k) add(u, rng() % n, 1 + rng() % 31);
    }
  }
  Case c{kind + (masked ? "_masked" : "") + (large ? "_int64" : "") +
             (reverse ? "_reverse" : ""),
         StrictMultiDiGraph::from_arrays(n, src, dst, caps, costs)};
  c.reverse = reverse;
  int e = c.graph.num_edges();
  c.node_mask = std::make_unique<bool[]>(n);
  c.edge_mask = std::make_unique<bool[]>(e);
  c.residual.assign(c.graph.capacity_view().begin(), c.graph.capacity_view().end());
  for (int i = 0; i < n; ++i) c.node_mask[i] = !masked || rng() % 23 != 0;
  for (int i = 0; i < e; ++i) {
    c.edge_mask[i] = !masked || rng() % 11 != 0;
    if (masked && rng() % 13 == 0) c.residual[i] = 0;
  }
  return c;
}

Result cpu(const Case& c, NodeId source) {
  if (c.reverse)
    return shortest_paths_to(c.graph, source, true, EdgeSelection{}, c.residual, c.nm(), c.em());
  return shortest_paths(c.graph, source, {}, true, EdgeSelection{}, c.residual, c.nm(), c.em());
}

// Stronger CPU comparator: no DAG allocation and no capacity tie-breaking.
std::vector<Cost> cpu_distance(const Case& c, NodeId source) {
  auto& g = c.graph;
  auto row = c.reverse ? g.in_row_offsets_view() : g.row_offsets_view();
  auto col = c.reverse ? g.in_col_indices_view() : g.col_indices_view();
  auto ids = c.reverse ? g.in_adj_edge_index_view() : g.adj_edge_index_view();
  auto costs = g.cost_view();
  std::vector<Cost> dist(g.num_nodes(), INF);
  if (!c.node_mask[source]) return dist;
  using Item = std::pair<Cost, NodeId>;
  std::priority_queue<Item, std::vector<Item>, std::greater<Item>> queue;
  dist[source] = 0; queue.emplace(0, source);
  while (!queue.empty()) {
    auto [d, u] = queue.top(); queue.pop();
    if (d != dist[u]) continue;
    for (int j = row[u]; j < row[u + 1]; ++j) {
      int v = col[j], e = ids[j];
      if (!c.node_mask[v] || !c.edge_mask[e] || c.residual[e] < kMinCap) continue;
      Cost next = d + costs[e];
      if (next < dist[v]) { dist[v] = next; queue.emplace(next, v); }
    }
  }
  return dist;
}

// Reusable workers: CPU batch baseline does not pay thread creation per sample.
class Pool {
  std::mutex mutex;
  std::condition_variable start, done;
  std::vector<std::thread> workers;
  std::function<void(int)> task;
  std::atomic<int> next{0};
  int count = 0, pending = 0, generation = 0;
  bool stop = false;
public:
  explicit Pool(int n) {
    for (int i = 0; i < n; ++i) workers.emplace_back([this] {
      int seen = 0;
      for (;;) {
        std::unique_lock lock(mutex);
        start.wait(lock, [&] { return stop || generation != seen; });
        if (stop) return;
        seen = generation;
        lock.unlock();
        for (int i; (i = next.fetch_add(1)) < count;) task(i);
        lock.lock();
        if (--pending == 0) done.notify_one();
      }
    });
  }
  void run(int n, std::function<void(int)> f) {
    std::unique_lock lock(mutex);
    task = std::move(f); count = n; next = 0; pending = int(workers.size()); ++generation;
    start.notify_all();
    done.wait(lock, [&] { return pending == 0; });
  }
  ~Pool() {
    { std::lock_guard lock(mutex); stop = true; }
    start.notify_all();
    for (auto& t : workers) t.join();
  }
};

struct Metal {
  id<MTLDevice> device;
  id<MTLCommandQueue> queue;
  id<MTLComputePipelineState> pull, group, shared;
  double setup_ms;
  Metal(const char* path) {
    auto t = Clock::now();
    device = MTLCreateSystemDefaultDevice();
    if (!device) throw std::runtime_error("No Metal device");
    queue = [device newCommandQueue];
    NSError* error = nil;
    NSString* code = [NSString stringWithContentsOfFile:[NSString stringWithUTF8String:path]
                                              encoding:NSUTF8StringEncoding error:&error];
    if (!code) throw std::runtime_error(error.localizedDescription.UTF8String);
    id<MTLLibrary> library = [device newLibraryWithSource:code options:nil error:&error];
    if (!library) throw std::runtime_error(error.localizedDescription.UTF8String);
    pull = [device newComputePipelineStateWithFunction:[library newFunctionWithName:@"relax_pull"] error:&error];
    group = [device newComputePipelineStateWithFunction:[library newFunctionWithName:@"solve_group"] error:&error];
    shared = [device newComputePipelineStateWithFunction:[library newFunctionWithName:@"solve_shared"] error:&error];
    if (!pull || !group || !shared) throw std::runtime_error(error.localizedDescription.UTF8String);
    setup_ms = ms(t);
    std::cerr << "device=" << device.name.UTF8String << " setup_ms=" << setup_ms
              << " execution_width=" << pull.threadExecutionWidth << '\n';
  }
  id<MTLBuffer> buffer(const void* ptr, size_t size) {
    auto b = [device newBufferWithLength:std::max(size_t(4), size) options:MTLResourceStorageModeShared];
    if (!b) throw std::runtime_error("Metal allocation failed");
    if (ptr && size) std::memcpy(b.contents, ptr, size);
    return b;
  }
};

struct GpuGraph {
  const Case& c;
  Metal& m;
  id<MTLBuffer> row, col, cost, allowed, a, b, changed, sources;
  std::vector<uint8_t> enabled;
  std::vector<int> valid_sources;
  size_t n, q;
  double upload_ms = 0, gpu_ms = 0, distance_ms = 0;
  int rounds = 0;
  GpuGraph(const Case& testcase, Metal& metal, const std::vector<NodeId>& query_sources)
      : c(testcase), m(metal), n(c.graph.num_nodes()), q(query_sources.size()) {
    auto t = Clock::now();
    auto& g = c.graph;
    // Reverse SSSP pulls from outgoing CSR. Eligibility uses original edge IDs.
    auto rows = c.reverse ? g.row_offsets_view() : g.in_row_offsets_view();
    auto cols = c.reverse ? g.col_indices_view() : g.in_col_indices_view();
    auto ids = c.reverse ? g.adj_edge_index_view() : g.in_adj_edge_index_view();
    std::vector<Cost> costs(ids.size());
    enabled.resize(ids.size());
    for (size_t j = 0; j < ids.size(); ++j) {
      int e = ids[j];
      enabled[j] = c.node_mask[g.edge_src_view()[e]] && c.node_mask[g.edge_dst_view()[e]] &&
                   c.edge_mask[e] && c.residual[e] >= kMinCap;
      costs[j] = g.cost_view()[e];
    }
    valid_sources = query_sources;
    for (auto& source : valid_sources) if (!c.node_mask[source]) source = -1;
    row = m.buffer(rows.data(), rows.size_bytes());
    col = m.buffer(cols.data(), cols.size_bytes());
    cost = m.buffer(costs.data(), costs.size() * sizeof(Cost));
    allowed = m.buffer(enabled.data(), enabled.size());
    a = m.buffer(nullptr, n * q * sizeof(Cost));
    b = m.buffer(nullptr, n * q * sizeof(Cost));
    changed = m.buffer(nullptr, std::max(size_t(8), q) * sizeof(uint32_t));
    sources = m.buffer(valid_sources.data(), q * sizeof(int));
    upload_ms = ms(t);
  }

  void bind(id<MTLComputeCommandEncoder> encoder, id<MTLBuffer> before, id<MTLBuffer> after) {
    [encoder setBuffer:row offset:0 atIndex:0];
    [encoder setBuffer:col offset:0 atIndex:1];
    [encoder setBuffer:cost offset:0 atIndex:2];
    [encoder setBuffer:allowed offset:0 atIndex:3];
    [encoder setBuffer:before offset:0 atIndex:4];
    [encoder setBuffer:after offset:0 atIndex:5];
    [encoder setBuffer:changed offset:0 atIndex:6];
  }
  void wait(id<MTLCommandBuffer> command) {
    [command commit]; [command waitUntilCompleted];
    if (command.status == MTLCommandBufferStatusError)
      throw std::runtime_error(command.error.localizedDescription.UTF8String);
    gpu_ms += (command.GPUEndTime - command.GPUStartTime) * 1000;
  }
  std::vector<Result> run(bool persistent, bool materialize = true, bool shared = false) {
    auto begin = Clock::now();
    gpu_ms = 0; rounds = 0;
    struct Params { uint32_t n, queries, round, unused; } p{uint32_t(n), uint32_t(q), 0, 0};
    std::vector<const Cost*> outputs(q);
    if (persistent) {
      id<MTLCommandBuffer> command = [m.queue commandBuffer];
      auto encoder = [command computeCommandEncoder];
      [encoder setComputePipelineState:shared ? m.shared : m.group];
      bind(encoder, a, b);
      [encoder setBytes:&p length:sizeof(p) atIndex:7];
      [encoder setBuffer:sources offset:0 atIndex:8];
      if (shared) {
        if (n * 16 + m.shared.staticThreadgroupMemoryLength > m.device.maxThreadgroupMemoryLength)
          throw std::runtime_error("threadgroup scratch too large");
        [encoder setThreadgroupMemoryLength:n * 16 atIndex:0];
      }
      [encoder dispatchThreadgroups:MTLSizeMake(q, 1, 1) threadsPerThreadgroup:MTLSizeMake(256, 1, 1)];
      [encoder endEncoding]; wait(command);
      auto counts = static_cast<const uint32_t*>(changed.contents);
      for (size_t i = 0; i < q; ++i) {
        rounds = std::max(rounds, int(counts[i]));
        outputs[i] = static_cast<const Cost*>((!shared && counts[i] % 2 ? b : a).contents) + i * n;
      }
    } else {
      auto initial = static_cast<Cost*>(a.contents);
      std::fill(initial, initial + n * q, INF);
      for (size_t i = 0; i < q; ++i) if (valid_sources[i] >= 0) initial[i * n + valid_sources[i]] = 0;
      id<MTLBuffer> before = a, after = b;
      while (rounds < int(n)) {
        int steps = std::min(8, int(n) - rounds);
        std::memset(changed.contents, 0, 8 * sizeof(uint32_t));
        id<MTLCommandBuffer> command = [m.queue commandBuffer];
        for (int j = 0; j < steps; ++j) {
          // Separate encoders with tracked resources establish the round dependency.
          auto encoder = [command computeCommandEncoder];
          [encoder setComputePipelineState:m.pull]; bind(encoder, before, after);
          p.round = j;
          [encoder setBytes:&p length:sizeof(p) atIndex:7];
          [encoder dispatchThreads:MTLSizeMake(n * q, 1, 1) threadsPerThreadgroup:MTLSizeMake(256, 1, 1)];
          [encoder endEncoding]; std::swap(before, after);
        }
        wait(command); rounds += steps;
        if (static_cast<const uint32_t*>(changed.contents)[steps - 1] == 0) break;
      }
      for (size_t i = 0; i < q; ++i) outputs[i] = static_cast<const Cost*>(before.contents) + i * n;
    }
    std::vector<Result> results(q);
    for (size_t i = 0; i < q; ++i) results[i].first.assign(outputs[i], outputs[i] + n);
    distance_ms = ms(begin);
    if (!materialize) return results;
    for (size_t query = 0; query < q; ++query) {
      auto& d = results[query].first;
      auto& dag = results[query].second;
      auto& g = c.graph;
      // Positive-cost ECMP reconstruction. Zero costs intentionally require fallback.
      auto esrc = g.edge_src_view(); auto edst = g.edge_dst_view(); auto costs = g.cost_view();
      auto tight = [&](size_t e) {
        int u = esrc[e], v = edst[e];
        return c.node_mask[u] && c.node_mask[v] && c.edge_mask[e] && c.residual[e] >= kMinCap &&
               (c.reverse ? (d[v] != INF && d[u] == d[v] + costs[e]) :
                            (d[u] != INF && d[v] == d[u] + costs[e]));
      };
      dag.parent_offsets.assign(n + 1, 0);
      for (size_t e = 0; e < costs.size(); ++e) if (tight(e)) ++dag.parent_offsets[edst[e] + 1];
      std::partial_sum(dag.parent_offsets.begin(), dag.parent_offsets.end(), dag.parent_offsets.begin());
      dag.parents.resize(dag.parent_offsets.back()); dag.via_edges.resize(dag.parents.size());
      auto cursor = dag.parent_offsets;
      for (size_t e = 0; e < costs.size(); ++e) if (tight(e)) {
        int pos = cursor[edst[e]]++; dag.parents[pos] = esrc[e]; dag.via_edges[pos] = int(e);
      }
    }
    return results;
  }
};

std::vector<std::tuple<int, int, int>> canonical(const PredDAG& dag) {
  std::vector<std::tuple<int, int, int>> entries;
  for (size_t v = 0; v + 1 < dag.parent_offsets.size(); ++v)
    for (int j = dag.parent_offsets[v]; j < dag.parent_offsets[v + 1]; ++j)
      entries.emplace_back(v, dag.parents[j], dag.via_edges[j]);
  std::sort(entries.begin(), entries.end()); return entries;
}

void verify(const std::vector<Result>& expected, const std::vector<Result>& actual) {
  if (expected.size() != actual.size()) throw std::runtime_error("query size mismatch");
  for (size_t q = 0; q < expected.size(); ++q) {
    if (expected[q].first != actual[q].first) {
      auto mismatch = std::mismatch(expected[q].first.begin(), expected[q].first.end(), actual[q].first.begin());
      throw std::runtime_error("distance mismatch at query " + std::to_string(q) + " node " +
          std::to_string(mismatch.first - expected[q].first.begin()) + " expected " +
          std::to_string(*mismatch.first) + " actual " + std::to_string(*mismatch.second));
    }
    if (canonical(expected[q].second) != canonical(actual[q].second))
      throw std::runtime_error("DAG edge set mismatch at query " + std::to_string(q));
  }
}

double median(std::vector<double> x) { std::sort(x.begin(), x.end()); return x[x.size() / 2]; }

void witnesses(Metal& metal) {
  auto fixture = [](int n, std::vector<NodeId> src, std::vector<NodeId> dst,
                    std::vector<Cost> costs, std::vector<Cap> caps) {
    Case c{"witness", StrictMultiDiGraph::from_arrays(n, src, dst, caps, costs)};
    c.node_mask = std::make_unique<bool[]>(n);
    c.edge_mask = std::make_unique<bool[]>(costs.size());
    std::fill_n(c.node_mask.get(), n, true);
    std::fill_n(c.edge_mask.get(), costs.size(), true);
    c.residual.assign(c.graph.capacity_view().begin(), c.graph.capacity_view().end());
    return c;
  };
  auto check = [&](Case& c, NodeId source, const std::vector<Cost>& distances, size_t edges) {
    auto expected = cpu(c, source);
    if (expected.first != distances || expected.second.via_edges.size() != edges)
      throw std::runtime_error("hand-derived witness mismatch");
    GpuGraph gpu(c, metal, {source});
    for (int method = 0; method < 3; ++method)
      verify({expected}, gpu.run(method != 0, true, method == 2));
  };
  auto c = fixture(5, {0, 0, 0, 1, 2, 1}, {1, 1, 2, 3, 3, 1},
                   {1, 1, 1, 1, 1, 1}, {1, 1, 1, 1, 1, 1});
  check(c, 0, {0, 1, 1, 2, INF}, 5);
  std::cout << "PASS diamond: all ECMP alternatives, parallel edges, positive self-loop excluded, isolated node\n";
  c.reverse = true;
  check(c, 3, {2, 1, 1, 0, INF}, 5);
  c.reverse = false;
  std::cout << "PASS reverse diamond: hand-derived distances and forward DAG\n";
  c.node_mask[0] = false;
  check(c, 0, {INF, INF, INF, INF, INF}, 0);
  c.node_mask[0] = true;
  c.node_mask[1] = false;
  check(c, 0, {0, INF, 1, 2, INF}, 2);
  c.node_mask[1] = true;
  std::cout << "PASS masked source and transit node\n";
  // All costs equal, graph order is (src,dst), making edge 0 the first parallel edge.
  c.residual[0] = std::nextafter(kMinCap, 0.0);
  c.residual[1] = kMinCap;
  check(c, 0, {0, 1, 1, 2, INF}, 4);
  std::cout << "PASS double capacity threshold: below rejected, exact threshold admitted\n";
  c.edge_mask[1] = false;
  check(c, 0, {0, INF, 1, 2, INF}, 2);
  std::cout << "PASS edge mask and residual gate combined\n";

  Cost big = Cost(1) << 54;
  auto wide = fixture(5, {0, 0, 1, 2}, {1, 2, 3, 3},
                      {big, big, big + 1, big + 2}, {1, 1, 1, 1});
  check(wide, 0, {0, big, big, 2 * big + 1, INF}, 3);
  std::cout << "PASS exact int64 costs above 2^53 and one-unit path difference\n";

  auto zero = fixture(3, {0, 1, 2}, {1, 2, 1}, {0, 0, 0}, {1, 1, 1});
  auto expected = cpu(zero, 0);
  GpuGraph zg(zero, metal, {0});
  auto actual = zg.run(true);
  if (actual[0].first != expected.first || expected.second.via_edges.size() != 2 ||
      actual[0].second.via_edges.size() != 3)
    throw std::runtime_error("zero-cost gap witness changed");
  std::cout << "PASS expected gap: zero-cost distances agree, tight-edge GPU DAG has 1->2->1 cycle; CPU DAG has 2 edges\n";
  if (float(16777216) != float(16777217) ||
      float(std::nextafter(kMinCap, 0.0)) != float(kMinCap))
    throw std::runtime_error("precision gap witness changed");
  std::cout << "PASS precision gaps: float32 collapses 16777216/16777217 and admits nextafter(kMinCap,0)\n";

  auto tie = fixture(2, {0, 0}, {1, 1}, {1, 1}, {1, 1 + 1e-8});
  EdgeSelection sel;
  sel.multi_edge = false; sel.tie_break = EdgeTieBreak::PreferHigherResidual;
  auto t = shortest_paths(tie.graph, 0, {}, false, sel);
  if (t.second.via_edges != std::vector<EdgeId>{1} || float(1 + 1e-8) != float(1))
    throw std::runtime_error("capacity tie witness changed");
  std::cout << "PASS capacity tie gap: CPU chooses edge 1 at 1+1e-8; float32 collapses both capacities\n";
}

void benchmark(Case& c, int queries, const std::string& mode, Metal* metal, Pool& pool, int samples) {
  std::vector<NodeId> sources(queries);
  for (int i = 0; i < queries; ++i) sources[i] = (i * 997) % c.graph.num_nodes();
  std::vector<Result> expected(queries);
  for (int i = 0; i < queries; ++i) expected[i] = cpu(c, sources[i]);
  std::unique_ptr<GpuGraph> gpu;
  if (metal) gpu = std::make_unique<GpuGraph>(c, *metal, sources);
  std::vector<std::string> methods = mode == "cpu" ?
    std::vector<std::string>{"cpu_full", "cpu_distance", "cpu_pool10", "cpu_pool10_distance"} :
    std::vector<std::string>{"metal_pull", "metal_group"};
  if (metal && size_t(c.graph.num_nodes()) * 16 + metal->shared.staticThreadgroupMemoryLength <= metal->device.maxThreadgroupMemoryLength)
    methods.push_back("metal_shared");
  for (const auto& method : methods) {
    std::vector<double> times, device_times, distance_times;
    int exact_order = 0;
    for (int rep = -2; rep < samples; ++rep) {
      @autoreleasepool {
        auto start = Clock::now();
        std::vector<Result> actual;
        if (method == "cpu_full") {
          actual.resize(queries);
          for (int i = 0; i < queries; ++i) actual[i] = cpu(c, sources[i]);
        } else if (method == "cpu_pool10") {
          actual.resize(queries);
          pool.run(queries, [&](int i) { actual[i] = cpu(c, sources[i]); });
        } else if (method == "cpu_distance" || method == "cpu_pool10_distance") {
          actual.resize(queries);
          auto task = [&](int i) { actual[i].first = cpu_distance(c, sources[i]); };
          if (method == "cpu_distance") {
            for (int i = 0; i < queries; ++i) task(i);
          } else pool.run(queries, task);
        } else actual = gpu->run(method != "metal_pull", true, method == "metal_shared");
        double elapsed = ms(start);
        // Full validation is outside timing, on every sample, not just warmup.
        if (method == "cpu_distance" || method == "cpu_pool10_distance") {
          for (int i = 0; i < queries; ++i) if (actual[i].first != expected[i].first)
            throw std::runtime_error("CPU distance baseline mismatch");
        } else {
          verify(expected, actual);
          exact_order = 0;
          for (int i = 0; i < queries; ++i)
            if (expected[i].second.parents == actual[i].second.parents &&
                expected[i].second.via_edges == actual[i].second.via_edges) ++exact_order;
        }
        if (rep >= 0) {
          times.push_back(elapsed);
          device_times.push_back(gpu ? gpu->gpu_ms : 0);
          distance_times.push_back(gpu ? gpu->distance_ms : elapsed);
        }
      }
    }
    std::cout << c.name << ',' << c.graph.num_nodes() << ',' << c.graph.num_edges() << ',' << queries
              << ',' << method << ',' << std::fixed << std::setprecision(6) << median(times)
              << ',' << median(device_times) << ',' << median(distance_times) << ','
              << (gpu ? gpu->upload_ms : 0) << ',' << (gpu ? gpu->rounds : 0) << ',' << samples
              << ",pass," << exact_order << ',';
    for (size_t i = 0; i < times.size(); ++i) std::cout << (i ? ";" : "") << times[i];
    std::cout << '\n' << std::flush;
  }
}

int main(int argc, char** argv) {
  @autoreleasepool {
    try {
      if (argc < 4) throw std::runtime_error("usage: experiment cpu|gpu kernel-path smoke|full|correctness [samples]");
      std::string mode = argv[1], suite = argv[3];
      if (mode != "cpu" && mode != "gpu") throw std::runtime_error("invalid mode");
      int samples = argc > 4 ? std::stoi(argv[4]) : 7;
      if (samples < 1) throw std::runtime_error("samples must be positive");
      Pool pool(10);
      std::unique_ptr<Metal> metal;
      if (mode == "gpu") metal = std::make_unique<Metal>(argv[2]);
      if (suite == "witnesses") {
        if (!metal) throw std::runtime_error("witnesses require gpu mode");
        witnesses(*metal); return 0;
      }
      std::cout << "case,n,e,queries,method,wall_ms,gpu_ms,distance_wall_ms,upload_ms,rounds,samples,correct,exact_order_queries,wall_samples_ms\n";
      auto run = [&](Case c, std::initializer_list<int> qs) {
        for (int q : qs) benchmark(c, q, mode, metal.get(), pool, samples);
      };
      if (suite == "stress") {
        run(make_case("grid", 1024), {64});
        run(make_case("random", 1024, 8, true, true, true), {64});
      } else if (suite == "correctness") {
        for (unsigned seed = 1; seed <= 12; ++seed) {
          run(make_case("random", 37 + seed * 3, 6, true, seed % 2, seed % 3 == 0, seed), {8});
        }
      } else {
        run(make_case("random", 128), {1, 64});
        run(make_case("random", 1024), {1, 64});
        run(make_case("fabric", 288), {1, 64});
        run(make_case("grid", 1024), {1, 64});
        run(make_case("chain", 256), {1});
        run(make_case("random", 1024, 8, true, true, true), {8});
        if (suite == "full") {
          run(make_case("random", 4096), {1, 64, 256});
          run(make_case("random", 16384), {1, 64});
          run(make_case("random", 65536), {1, 64});
          run(make_case("chain", 2048), {1});
          run(make_case("random", 4096, 32), {1, 64});
          run(make_case("clos", 200), {1, 64});
          run(make_case("clos", 800), {1, 64});
          run(make_case("grid", 10000), {1, 8});
          run(make_case("grid", 40000), {1});
        } else if (suite != "smoke") throw std::runtime_error("invalid suite");
      }
    } catch (const std::exception& e) {
      std::cerr << "FAIL: " << e.what() << '\n'; return 1;
    }
  }
}
