// The private hooks (NETGRAPH_CORE_PRIVATE_HOOKS builds only): the test controls of Python and
// native allocation failures, and the apportionment kernels.
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include "registration.hpp"

#if NETGRAPH_CORE_PRIVATE_HOOKS
#include <atomic>
#include <cstdint>
#include <cstdlib>
#include <new>
#include <thread>
#include <vector>

#include "netgraph/core/detail/kernels.hpp"
#include "netgraph/core/detail/private_hooks.hpp"
#include "netgraph/core/flow_policy.hpp"
#include "netgraph/core/types.hpp"
#include "integers.hpp"

namespace py = pybind11;

namespace netgraph::core::bindings {

namespace {

// Take the normalized exception, including its traceback, without replacing an existing
// exception instance. PyErr_GetRaisedException was added in Python 3.12.
py::object take_raised_exception() {
#if PY_VERSION_HEX >= 0x030C0000
  return py::reinterpret_steal<py::object>(PyErr_GetRaisedException());
#else
  PyObject* type = nullptr;
  PyObject* value = nullptr;
  PyObject* traceback = nullptr;
  PyErr_Fetch(&type, &value, &traceback);
  PyErr_NormalizeException(&type, &value, &traceback);
  const auto kind = py::reinterpret_steal<py::object>(type);
  auto error = py::reinterpret_steal<py::object>(value);
  const auto trace = py::reinterpret_steal<py::object>(traceback);
  if (trace && PyException_SetTraceback(error.ptr(), trace.ptr()) < 0) throw py::error_already_set();
  return error;
#endif
}

// The allocation-failure hook: the k-th object-domain allocation fails (Python's allocator
// wrapped for one call, then restored). The allocator is process-wide while the call runs, so
// the hook serves one caller at a time (the tests); its counters are atomics because the module
// does not hold the GIL for every thread (mod_gil_not_used).
PyMemAllocatorEx g_saved_allocator {};
std::atomic<long> g_fail_countdown {-1};
std::atomic<long> g_fail_seen {0};
std::atomic<bool> g_fail_hit {false};

bool fail_now() {
  g_fail_seen.fetch_add(1);
  if (g_fail_countdown.load() < 0) return false;
  if (g_fail_countdown.fetch_sub(1) == 0) {
    g_fail_hit.store(true);
    return true;
  }
  return false;
}
void* failing_malloc(void*, std::size_t n) { return fail_now() ? nullptr : g_saved_allocator.malloc(g_saved_allocator.ctx, n); }
void* failing_calloc(void*, std::size_t n, std::size_t size) {
  return fail_now() ? nullptr : g_saved_allocator.calloc(g_saved_allocator.ctx, n, size);
}
void* failing_realloc(void*, void* ptr, std::size_t n) {
  return fail_now() ? nullptr : g_saved_allocator.realloc(g_saved_allocator.ctx, ptr, n);
}
void forwarding_free(void*, void* ptr) { g_saved_allocator.free(g_saved_allocator.ctx, ptr); }

// The native allocation-failure hook: while armed for one call, the k-th allocation of
// the calling thread through the extension's global operator new (the bindings and the library
// linked into the extension) throws std::bad_alloc, and the thread's allocations and deletions
// are counted, so a test also sees what a failed call left allocated. Other threads and other
// images (libc++, Python) are not counted; the aligned forms are not replaced. The counters are
// atomics for the same reason as above.
std::atomic<bool> g_native_armed {false};
std::atomic<std::thread::id> g_native_owner {};
std::atomic<long> g_native_countdown {-1};
std::atomic<long> g_native_seen {0};
std::atomic<long> g_native_live {0};
std::atomic<bool> g_native_hit {false};

bool native_counted() noexcept {
  return g_native_armed.load(std::memory_order_relaxed) && g_native_owner.load() == std::this_thread::get_id();
}

void* native_allocate(std::size_t n) {
  if (native_counted()) {
    g_native_seen.fetch_add(1);
    if (g_native_countdown.load() >= 0 && g_native_countdown.fetch_sub(1) == 0) {
      g_native_hit.store(true);
      throw std::bad_alloc();
    }
  }
  void* p = std::malloc(n == 0 ? 1 : n);
  if (p == nullptr) throw std::bad_alloc();
  if (native_counted()) g_native_live.fetch_add(1);
  return p;
}

void* native_allocate_nothrow(std::size_t n) noexcept {
  try {
    return native_allocate(n);
  } catch (...) {
    return nullptr;
  }
}

void native_free(void* p) noexcept {
  if (p != nullptr && native_counted()) g_native_live.fetch_sub(1);
  std::free(p);
}

} // namespace

void register_private_hooks(py::module_& m) {
  m.def(
      "_fail_python_allocation",
      [](const py::object& fn, const py::tuple& args, long k) {
        // Calls fn(*args) with the k-th Python object allocation (from 0) failing, then
        // restores the allocator. Returns (hit, allocations seen, result or None, exception or
        // None). One caller at a time: the allocator is process-wide while the call runs.
        PyMem_GetAllocator(PYMEM_DOMAIN_OBJ, &g_saved_allocator);
        PyMemAllocatorEx hook {nullptr, failing_malloc, failing_calloc, failing_realloc, forwarding_free};
        g_fail_countdown.store(k);
        g_fail_seen.store(0);
        g_fail_hit.store(false);
        PyMem_SetAllocator(PYMEM_DOMAIN_OBJ, &hook);
        PyObject* r = PyObject_CallObject(fn.ptr(), args.ptr());
        PyMem_SetAllocator(PYMEM_DOMAIN_OBJ, &g_saved_allocator);
        g_fail_countdown.store(-1);
        py::object result = r != nullptr ? py::reinterpret_steal<py::object>(r) : py::object(py::none());
        py::object error = py::none();
        if (r == nullptr) error = take_raised_exception();
        return py::make_tuple(g_fail_hit.load(), g_fail_seen.load(), result, error);
      },
      py::arg("fn"), py::arg("args"), py::arg("k"),
      "private test hook: fn(*args) with the k-th Python object allocation failing -> (hit, seen, result, error)");

  m.def(
      "_fail_native_allocation",
      [](const py::object& fn, const py::tuple& args, long k) {
        // Calls fn(*args) with the calling thread's k-th native allocation (from 0) failing.
        // Returns (hit, allocations seen, result or None, exception or None, allocations of the
        // call still live). One caller at a time.
        g_native_countdown.store(k);
        g_native_seen.store(0);
        g_native_live.store(0);
        g_native_hit.store(false);
        g_native_owner.store(std::this_thread::get_id());
        g_native_armed.store(true);
        PyObject* r = PyObject_CallObject(fn.ptr(), args.ptr());
        g_native_armed.store(false);
        g_native_countdown.store(-1);
        py::object result = r != nullptr ? py::reinterpret_steal<py::object>(r) : py::object(py::none());
        py::object error = py::none();
        if (r == nullptr) error = take_raised_exception();
        return py::make_tuple(g_native_hit.load(), g_native_seen.load(), result, error, g_native_live.load());
      },
      py::arg("fn"), py::arg("args"), py::arg("k"),
      "private test hook: fn(*args) with the k-th native allocation failing -> (hit, seen, result, error, live)");

  m.def(
      "_round_robin",
      [](const py::object& total, const py::object& k) {
        const auto t = integer_value(total, 0, kDomainLimit - 1, "total");
        const auto n = integer_value(k, 1, 1 << 20, "k");
        std::vector<Amount> out(static_cast<std::size_t>(n));
        detail::round_robin(t, out);
        return out;
      },
      py::arg("total"), py::arg("k"),
      "private: round robin of total units over k members (the first total % k get one unit more)");

  m.def(
      "_highest_averages",
      [](const py::object& total, const py::object& weights) {
        const auto t = integer_value(total, 0, kDomainLimit - 1, "total");
        const auto w = wide_integers(weights, "weights");
        if (w.size() == 0) throw py::value_error("highest averages over no member");
        for (std::size_t i = 0; i < w.size(); ++i) {
          if (!w.in(i, 1, kDomainLimit - 1)) {
            throw py::value_error("weights[" + std::to_string(i) + "] = " + w.text(i) + " is outside [1, 2^62)");
          }
        }
        const auto wv = w.narrow<Amount>();
        std::vector<Amount> out(wv.size());
        detail::highest_averages(t, wv, out);
        return out;
      },
      py::arg("total"), py::arg("weights"),
      "private: highest-averages (D'Hondt) apportionment of total units over positive weights, ties by\n"
      "the lower index");
}

} // namespace netgraph::core::bindings

// The extension's global allocation functions in a hooks build: the native allocation-failure
// hook above. Unarmed they allocate as the defaults do (malloc and free, which every image's
// operator delete accepts).
using netgraph::core::bindings::native_allocate;
using netgraph::core::bindings::native_allocate_nothrow;
using netgraph::core::bindings::native_free;

void* operator new(std::size_t n) { return native_allocate(n); }
void* operator new[](std::size_t n) { return native_allocate(n); }
void* operator new(std::size_t n, const std::nothrow_t&) noexcept { return native_allocate_nothrow(n); }
void* operator new[](std::size_t n, const std::nothrow_t&) noexcept { return native_allocate_nothrow(n); }
void operator delete(void* p) noexcept { native_free(p); }
void operator delete[](void* p) noexcept { native_free(p); }
void operator delete(void* p, std::size_t) noexcept { native_free(p); }
void operator delete[](void* p, std::size_t) noexcept { native_free(p); }
void operator delete(void* p, const std::nothrow_t&) noexcept { native_free(p); }
void operator delete[](void* p, const std::nothrow_t&) noexcept { native_free(p); }

#endif
