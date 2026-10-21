// Reject an unconstructed Python policy before pybind11's ordinary value loader
// can allocate raw storage for it. Explicit __new__/failed __init__ is legal Python.
#pragma once
#include <pybind11/pybind11.h>
#include "netgraph/core/flow_policy.hpp"
#include "netgraph/core/detail/policy_engine.hpp"

namespace netgraph::core::bindings {
inline FlowPolicy& policy_value(const pybind11::handle& self) {
  const auto* type = pybind11::detail::get_type_info(typeid(FlowPolicy), true);
  // isinstance can consult a forged __class__. Only the real extension type has
  // the instance layout read below; never reinterpret an unrelated Python object.
  if (PyObject_TypeCheck(self.ptr(), type->type) == 0) throw pybind11::type_error("policy must be a FlowPolicy");
  auto* instance = reinterpret_cast<pybind11::detail::instance*>(self.ptr());
  auto value = instance->get_value_and_holder(type);
  if (!value.holder_constructed() || !value.instance_registered() || value.value_ptr() == nullptr) {
    throw pybind11::type_error("FlowPolicy is not constructed");
  }
  auto& policy = *static_cast<FlowPolicy*>(value.value_ptr());
  if (!detail::PolicyAccess::constructed(policy)) throw pybind11::type_error("FlowPolicy is not constructed");
  return policy;
}
} // namespace netgraph::core::bindings
