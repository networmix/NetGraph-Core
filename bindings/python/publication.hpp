/*
  The checked construction of result objects: the Python object of a bound class is allocated
  through a checked tp_alloc, so an allocation failure raises instead of crashing, and a policy
  call makes its result object before it changes anything.
*/
#pragma once

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include <memory>
#include <type_traits>
#include <utility>

#include "netgraph/core/evaluators.hpp"
#include "netgraph/core/flow_policy.hpp"

namespace py = pybind11;

namespace netgraph::core::bindings {

// A new Python object of the bound class T owning a T made from `args`, made like pybind11's
// generic move-cast (type_caster_generic::cast) but with the object allocation checked:
// pybind11's make_new_instance uses the result of tp_alloc unchecked, so a failed Python
// allocation there is a crash rather than a MemoryError. Holder is the class's holder type
// (std::unique_ptr<T> by default, std::shared_ptr<T> for PredDAG and Path); another type is
// refused before anything is made.
// The steps, in order (new_object and object_from both use them): (1) allocate the object and
// check for null; (2) hold it in `obj` at once, so any later failure releases it; (3) allocate
// its layout; (4) make the value inside a local holder; (5) register the value and let the
// object take the holder; (6) mark the object as owner. Every failure raises and frees what was
// made exactly once. Until step 5 succeeds the object owns nothing: the local holder frees the
// value (std::shared_ptr also deletes it when its control block cannot be allocated), and the
// object is released without touching it. Taking the holder never allocates (pybind11 copies or
// moves it).
template <class T, class Holder, class... Args>
std::pair<py::object, T*> emplace_object(Args&&... args) {
  static_assert(std::is_same_v<Holder, std::unique_ptr<T>> || std::is_same_v<Holder, std::shared_ptr<T>>);
  const py::detail::type_info* tinfo = py::detail::get_type_info(typeid(T), /*throw_if_missing=*/true);
  constexpr auto holder_kind = std::is_same_v<Holder, std::shared_ptr<T>> ? py::detail::holder_enum_t::std_shared_ptr
                                                                          : py::detail::holder_enum_t::std_unique_ptr;
  if (tinfo->holder_enum_v != holder_kind) {
    py::pybind11_fail("emplace_object: Holder is not the holder of the bound class");
  }
  PyTypeObject* type = tinfo->type;
  PyObject* raw = type->tp_alloc(type, 0);
  if (raw == nullptr) throw py::error_already_set();
  auto obj = py::reinterpret_steal<py::object>(raw);
  auto* inst = reinterpret_cast<py::detail::instance*>(raw);
  inst->allocate_layout();
  inst->owned = false;
  Holder holder(new T(std::forward<Args>(args)...));
  T* value = holder.get();
  py::detail::values_and_holders(inst).begin()->value_ptr() = value;
  tinfo->init_instance(inst, &holder);
  inst->owned = true;
  return {std::move(obj), value};
}

// A new Python object of the bound class T (default holder) owning a default-constructed T
// (emplace_object). The policy calls make their result object with this before they change
// anything.
template <class T>
std::pair<py::object, T*> new_object() {
  return emplace_object<T, std::unique_ptr<T>>();
}

// The value moved into its result object after the commit: nothing may allocate there.
static_assert(std::is_nothrow_move_assignable_v<PlacementResult>);
static_assert(std::is_nothrow_move_assignable_v<SensitivityResult>);
static_assert(std::is_nothrow_move_assignable_v<py::object>);

template <class T>
py::object object_of(T value) {
  auto [obj, slot] = new_object<T>();
  *slot = std::move(value);
  return std::move(obj);
}

// A new Python object of the bound class T (held by std::shared_ptr) owning `value`, through the
// same checked allocation as new_object: for the values without a default constructor (PredDAG,
// Path) that a search returns. A failure raises and releases everything made so far.
template <class T>
py::object object_from(T value) {
  return emplace_object<T, std::shared_ptr<T>>(std::move(value)).first;
}

} // namespace netgraph::core::bindings
