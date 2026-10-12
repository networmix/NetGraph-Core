// The allocation-failure hook of the fault-injection executable (netgraph_core_fault_tests):
// the replacement global operator new and delete (throwing and nothrow forms), the countdown and
// the leak counter declared in fault_injection.hpp, and fails_at / fails_in_workers_at /
// reaches_in_workers_at, which every fault test file sweeps its operations through.
#include <gtest/gtest.h>

#include <atomic>
#include <cstddef>
#include <cstdlib>
#include <functional>
#include <new>

#include "fault_injection.hpp"

// Allocations left before the injected failure; negative: no injection (fault_injection.hpp).
// Shared by every thread; with g_workers_only the calling test thread's allocations do not count.
std::atomic<long> g_countdown {-1};
// Live allocations of the replacement operator new (a successful allocation counts up, the
// delete of a non-null pointer counts down): the sweeps compare it before and after an operation
// to count leaks, which ASan's leak detection (off on macOS) would not.
std::atomic<long> g_live {0};
static std::atomic<bool> g_workers_only {false};
static thread_local bool t_test_thread = false;

void* operator new(std::size_t n) {
  if (g_countdown.load(std::memory_order_relaxed) >= 0 &&
      !(t_test_thread && g_workers_only.load(std::memory_order_relaxed)) &&
      g_countdown.fetch_sub(1, std::memory_order_relaxed) == 0) {
    throw std::bad_alloc();
  }
  if (void* p = std::malloc(n == 0 ? 1 : n)) {
    g_live.fetch_add(1, std::memory_order_relaxed);
    return p;
  }
  throw std::bad_alloc();
}
void* operator new[](std::size_t n) { return ::operator new(n); }
// The nothrow forms never inject: a caller that asks for memory it can do without (the temporary
// buffer of std::stable_partition) absorbs a failure by design, so an injected one would be
// reported as swallowed. Their allocations count as live like the others. (Without these
// replacements, the standard library's nothrow form would call the throwing one above and catch
// the injected std::bad_alloc.)
void* operator new(std::size_t n, const std::nothrow_t&) noexcept {
  void* p = std::malloc(n == 0 ? 1 : n);
  if (p != nullptr) g_live.fetch_add(1, std::memory_order_relaxed);
  return p;
}
void* operator new[](std::size_t n, const std::nothrow_t&) noexcept { return ::operator new(n, std::nothrow); }
void operator delete(void* p) noexcept {
  if (p != nullptr) g_live.fetch_sub(1, std::memory_order_relaxed);
  std::free(p);
}
void operator delete[](void* p) noexcept { ::operator delete(p); }
void operator delete(void* p, std::size_t) noexcept { ::operator delete(p); }
void operator delete[](void* p, std::size_t) noexcept { ::operator delete(p); }
void operator delete(void* p, const std::nothrow_t&) noexcept { ::operator delete(p); }
void operator delete[](void* p, const std::nothrow_t&) noexcept { ::operator delete(p); }

namespace {

// Runs fn with the k-th allocation failing; returns whether it was reached, and in `failed`
// whether fn raised std::bad_alloc. Any other exception of fn propagates, with the hook disarmed
// first: left armed, it would fail an allocation of whatever runs next (the test framework's
// report of that exception, a later test).
bool run(long k, bool workers_only, const std::function<void()>& fn, bool& failed) {
  const auto disarm = [] {
    g_countdown = -1;
    g_workers_only = false;
  };
  t_test_thread = true;
  g_workers_only = workers_only;
  g_countdown = k;
  failed = false;
  try {
    fn();
  } catch (const std::bad_alloc&) {
    failed = true;
  } catch (...) {
    disarm();
    throw;
  }
  const bool reached = g_countdown.load() < 0;
  disarm();
  return reached;
}

bool inject(long k, bool workers_only, const std::function<void()>& fn) {
  bool failed = false;
  const bool reached = run(k, workers_only, fn, failed);
  EXPECT_EQ(failed, reached) << "a std::bad_alloc was swallowed or replaced at allocation " << k
                             << (workers_only ? " (worker threads)" : "");
  return failed;
}

} // namespace

// Runs fn with the k-th allocation (from 0) failing; returns whether that failure happened.
bool fails_at(long k, const std::function<void()>& fn) { return inject(k, false, fn); }

bool fails_in_workers_at(long k, const std::function<void()>& fn) { return inject(k, true, fn); }

bool reaches_in_workers_at(long k, const std::function<void()>& fn, bool& bad_alloc) {
  const bool reached = run(k, true, fn, bad_alloc);
  EXPECT_TRUE(reached || !bad_alloc) << "a std::bad_alloc that was not injected, at allocation " << k
                                     << " (worker threads)";
  return reached;
}
