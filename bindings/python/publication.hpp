/*
  The checked construction of result objects: the Python object of a bound class is allocated
  through a checked tp_alloc, so an allocation failure raises instead of crashing, and a policy
  call makes its result object before it changes anything.
*/
#pragma once

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <pybind11/typing.h>

#include <memory>
#include <optional>
#include <type_traits>
#include <utility>
#include <vector>

#include "netgraph/core/evaluators.hpp"
#include "netgraph/core/flow_policy.hpp"

namespace py = pybind11;

namespace netgraph::core::bindings {

// A new Python object of the bound class T owning a T made from `args`, made like pybind11's
// generic move-cast (type_caster_generic::cast) but with the object allocation checked:
// pybind11's make_new_instance uses the result of tp_alloc unchecked, so a failed Python
// allocation there is a crash rather than a MemoryError. Holder is the class's holder type
// (std::unique_ptr<T> by default, std::shared_ptr<T> for PredecessorDAG and Path); another type is
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
  return emplace_object<T, std::unique_ptr<T>>(std::move(value)).first;
}

// A tuple of new Python objects of the bound class T, one per row (object_of each): a failed
// allocation raises and releases what was made (pybind11's cast of a vector of bound values
// allocates each object unchecked).
template <class T>
py::tuple tuple_of(const std::vector<T>& rows) {
  PyObject* raw = PyTuple_New(static_cast<Py_ssize_t>(rows.size()));
  if (raw == nullptr) throw py::error_already_set();
  auto out = py::reinterpret_steal<py::tuple>(raw);
  for (std::size_t i = 0; i < rows.size(); ++i) {
    PyTuple_SET_ITEM(out.ptr(), static_cast<Py_ssize_t>(i), object_of<T>(rows[i]).release().ptr());
  }
  return out;
}

// A new Python object of the bound class T (held by std::shared_ptr) owning `value`, through the
// same checked allocation as new_object: for the values without a default constructor (PredecessorDAG,
// Path) that a search returns. A failure raises and releases everything made so far.
template <class T>
py::object object_from(T value) {
  return emplace_object<T, std::shared_ptr<T>>(std::move(value)).first;
}

// The Python object of a shared value of the bound class T (held by std::shared_ptr): None for a
// null pointer; the object already registered for it, if any (as pybind11's cast); else a new
// object made by emplace_object's steps with `holder` in place of a new value. pybind11's own
// cast of an existing holder marks the object as owner before it registers the value, so a
// failed registration deletes the value while the holder still owns it (a double free); here
// the object owns nothing until it holds its copy of the holder (a copy never allocates).
template <class T>
py::object shared_object(std::shared_ptr<T> holder) {
  if (!holder) return py::none();
  const py::detail::type_info* tinfo = py::detail::get_type_info(typeid(T), /*throw_if_missing=*/true);
  if (tinfo->holder_enum_v != py::detail::holder_enum_t::std_shared_ptr) {
    py::pybind11_fail("shared_object: the bound class is not held by std::shared_ptr");
  }
  if (py::handle existing = py::detail::find_registered_python_instance(holder.get(), tinfo)) {
    return py::reinterpret_steal<py::object>(existing);
  }
  PyTypeObject* type = tinfo->type;
  PyObject* raw = type->tp_alloc(type, 0);
  if (raw == nullptr) throw py::error_already_set();
  auto obj = py::reinterpret_steal<py::object>(raw);
  auto* inst = reinterpret_cast<py::detail::instance*>(raw);
  inst->allocate_layout();
  inst->owned = false;
  py::detail::values_and_holders(inst).begin()->value_ptr() = holder.get();
  tinfo->init_instance(inst, &holder);
  inst->owned = true;
  return obj;
}

template <class T>
py::object shared_object(std::shared_ptr<const T> holder) {
  return shared_object(std::const_pointer_cast<T>(std::move(holder)));
}

// A py::object known to be a T, for signatures: pybind11 names it as T (the specialization of
// handle_type_name below), so a binding that returns shared_object's result keeps its type.
template <class T>
class TypedPythonObject : public py::object {
public:
  TypedPythonObject() = default;
  TypedPythonObject(py::handle h, py::object::stolen_t) : py::object(h, py::object::stolen_t {}) {}
  TypedPythonObject(py::handle h, py::object::borrowed_t) : py::object(h, py::object::borrowed_t {}) {}
  static bool check_(py::handle h) { return py::isinstance<T>(h); }
};

// object_of(value) as a T; a list of object_of each row (as tuple_of); object_of of a present
// value or None. pybind11's cast of a bound value, a vector or an optional of them allocates the
// object unchecked.
template <class T>
TypedPythonObject<T> typed_object(T value) {
  return py::reinterpret_steal<TypedPythonObject<T>>(object_of<T>(std::move(value)).release());
}
template <class T>
py::typing::List<T> list_of(const std::vector<T>& rows) {
  PyObject* raw = PyList_New(static_cast<Py_ssize_t>(rows.size()));
  if (raw == nullptr) throw py::error_already_set();
  auto out = py::reinterpret_steal<py::typing::List<T>>(raw);
  for (std::size_t i = 0; i < rows.size(); ++i) {
    PyList_SET_ITEM(out.ptr(), static_cast<Py_ssize_t>(i), object_of<T>(rows[i]).release().ptr());
  }
  return out;
}
template <class T>
py::typing::Optional<T> optional_of(const std::optional<T>& value) {
  return py::reinterpret_steal<py::typing::Optional<T>>(value.has_value() ? object_of<T>(*value).release()
                                                                          : py::none().release());
}

// shared_object(holder) as a T (holder not null) or as T | None.
template <class T>
TypedPythonObject<std::remove_const_t<T>> typed_shared(std::shared_ptr<T> holder) {
  return py::reinterpret_steal<TypedPythonObject<std::remove_const_t<T>>>(shared_object(std::move(holder)).release());
}
template <class T>
py::typing::Optional<std::remove_const_t<T>> optional_shared(std::shared_ptr<T> holder) {
  return py::reinterpret_steal<py::typing::Optional<std::remove_const_t<T>>>(shared_object(std::move(holder)).release());
}

// The getters of read-only fields whose values are Python objects of bound classes or
// enumerations: def_readonly casts them through pybind11's unchecked allocation (a failed Python
// allocation is a crash); these make them through the checked one above. A field of a bound value
// or enumeration, an optional one, a shared one (its registered object, or None) and a vector of
// bound values (a list).
template <class C, class D>
auto field_object(D C::*member) {
  return [member](const C& c) { return typed_object(c.*member); };
}
template <class C, class D>
auto optional_field(std::optional<D> C::*member) {
  return [member](const C& c) { return optional_of(c.*member); };
}
template <class C, class D>
auto shared_field(std::shared_ptr<D> C::*member) {
  return [member](const C& c) { return optional_shared(c.*member); };
}
template <class C, class D>
auto list_field(std::vector<D> C::*member) {
  return [member](const C& c) { return list_of(c.*member); };
}

} // namespace netgraph::core::bindings

namespace pybind11::detail {
template <class T>
struct handle_type_name<netgraph::core::bindings::TypedPythonObject<T>> {
  static constexpr auto name = make_caster<T>::name;
};
} // namespace pybind11::detail
