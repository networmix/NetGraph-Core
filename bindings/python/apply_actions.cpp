// Several demands at once: BatchMode, the actions Place, Replace and Forward, and apply_actions.
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <cstddef>
#include <memory>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include "netgraph/core/eligibility.hpp"
#include "netgraph/core/flow_ledger.hpp"
#include "netgraph/core/flow_policy.hpp"
#include "integers.hpp"
#include "publication.hpp"
#include "values.hpp"
#include "registration.hpp"

namespace py = pybind11;

namespace netgraph::core::bindings {

namespace {

// The Python actions keep their arguments as given (the objects themselves, so an action keeps its
// policy and black holes alive); the call converts them.
struct PyPlace {
  py::object policy;
  py::object offer;
};
struct PyReplace {
  py::object policy;
  py::object offer;
  py::object design_demand;
};
struct PyForward {
  py::object policy;
  py::object offer;
  py::object black_holes;
};

std::string type_name(const py::handle& obj) { return std::string(py::str(py::type::of(obj).attr("__name__"))); }

py::object policy_arg(const py::object& policy) {
  if (!py::isinstance<FlowPolicy>(policy)) throw py::type_error("policy must be a FlowPolicy (got " + type_name(policy) + ")");
  return policy;
}

// An offer: an int (an __index__ object is kept as its int) or a float, None refused.
py::object offer_arg(const py::object& offer) {
  if (offer.is_none()) throw py::type_error("offer must be an int or a float (got NoneType)");
  return units_as_given(offer, "offer");
}

py::object black_holes_arg(const py::object& black_holes) {
  (void)eligibility_or_none(black_holes, "black_holes");
  return black_holes;
}

// The offer of a converted action: +inf nullopt; any value outside the domain -1, so the call's
// own check raises the offer text in its place.
std::optional<Amount> offer_at_use(const py::object& offer) { return *units_at_use(offer); }

FlowPolicy& policy_at_use(const py::object& policy) { return policy.cast<FlowPolicy&>(); }

// The native action of entry i, or the type error of an entry that is not an action.
Action action_of(const py::handle& item, std::size_t i) {
  if (py::isinstance<PyPlace>(item)) {
    const auto& a = item.cast<const PyPlace&>();
    return Place {policy_at_use(a.policy), offer_at_use(a.offer)};
  }
  if (py::isinstance<PyReplace>(item)) {
    const auto& a = item.cast<const PyReplace&>();
    return Replace {policy_at_use(a.policy), offer_at_use(a.offer), units_at_use(a.design_demand)};
  }
  if (py::isinstance<PyForward>(item)) {
    const auto& a = item.cast<const PyForward&>();
    const Eligibility* black_holes = a.black_holes.is_none() ? nullptr : a.black_holes.cast<const Eligibility*>();
    return Forward {policy_at_use(a.policy), offer_at_use(a.offer), black_holes};
  }
  throw py::type_error("actions[" + std::to_string(i) + "] must be a Place, Replace or Forward (got " + type_name(item) + ")");
}

} // namespace

void register_apply_actions(py::module_& m) {
  py::enum_<BatchMode>(m, "BatchMode",
                       "How apply_actions runs its list. ORDERED: each action as its scalar call, in list order, on the\n"
                       "state the earlier actions left (earlier actions are served first). SIMULTANEOUS: every action forwarded at\n"
                       "the same time on one view, losses shared on links across the demands.")
      .value("ORDERED", BatchMode::Ordered)
      .value("SIMULTANEOUS", BatchMode::Simultaneous);

  py::class_<PyPlace>(m, "Place",
                      "An incremental place(offer) of policy in an apply_actions list (ORDERED only). An immutable value:\n"
                      "policy (a FlowPolicy, kept alive) and offer (int units or a float, kept as given; another type\n"
                      "raises TypeError here, a value outside the domain raises the call's ValueError at apply_actions).")
      .def(py::init([](const py::object& policy, const py::object& offer) {
             return PyPlace {policy_arg(policy), offer_arg(offer)};
           }),
           py::arg("policy"), py::arg("offer"))
      .def_readonly("policy", &PyPlace::policy, "the FlowPolicy")
      .def_readonly("offer", &PyPlace::offer, "the offer as given (int units or a float)");

  py::class_<PyReplace>(m, "Replace",
                        "A replace(offer, design_demand) of policy in an apply_actions list; SIMULTANEOUS takes it only for\n"
                        "ShortestPaths EQUAL or WEIGHTED and for Pinned. An immutable value: policy (kept alive), offer and\n"
                        "design_demand (int units or a float, kept as given; None: absent; another type raises TypeError\n"
                        "here, a value outside the domain raises the call's ValueError at apply_actions).")
      .def(py::init([](const py::object& policy, const py::object& offer, const py::object& design_demand) {
             return PyReplace {policy_arg(policy), offer_arg(offer), units_as_given(design_demand, "design_demand")};
           }),
           py::arg("policy"), py::arg("offer"), py::kw_only(), py::arg("design_demand") = py::none())
      .def_readonly("policy", &PyReplace::policy, "the FlowPolicy")
      .def_readonly("offer", &PyReplace::offer, "the offer as given (int units or a float)")
      .def_readonly("design_demand", &PyReplace::design_demand, "the design demand as given (int units or a float); None: absent");

  py::class_<PyForward>(m, "Forward",
                        "A forward(offer, black_holes) of policy over its installed program in an apply_actions list. An\n"
                        "immutable value: policy and black_holes (an Eligibility of the policy's graph or None) kept alive,\n"
                        "offer (int units or a float, kept as given; another type raises TypeError here, a value outside\n"
                        "the domain raises the call's ValueError at apply_actions).")
      .def(py::init([](const py::object& policy, const py::object& offer, const py::object& black_holes) {
             return PyForward {policy_arg(policy), offer_arg(offer), black_holes_arg(black_holes)};
           }),
           py::arg("policy"), py::arg("offer"), py::kw_only(), py::arg("black_holes") = py::none())
      .def_readonly("policy", &PyForward::policy, "the FlowPolicy")
      .def_readonly("offer", &PyForward::offer, "the offer as given (int units or a float)")
      .def_readonly("black_holes", &PyForward::black_holes, "the black holes (an Eligibility); None: none");

  m.def(
      "apply_actions",
      [](const py::object& ledger, const py::object& actions, const py::object& mode, const py::object& outputs) {
        // The arguments' types in argument order, every entry included, before any native check.
        if (!py::isinstance<FlowLedger>(ledger)) throw py::type_error("ledger must be a FlowLedger");
        if (is_invalid_sequence(actions)) {
          throw py::type_error("actions must be a sequence of Place, Replace or Forward (got " + type_name(actions) + ")");
        }
        const auto items = py::reinterpret_borrow<py::sequence>(actions);
        const Py_ssize_t count = PySequence_Size(items.ptr());
        if (count < 0) throw py::error_already_set();
        const auto n = static_cast<std::size_t>(count);
        // The entries are held for the call: the actions borrow their policies and black holes.
        std::vector<py::object> held;
        held.reserve(n);
        std::vector<Action> native;
        native.reserve(n);
        for (std::size_t i = 0; i < n; ++i) {
          PyObject* raw = PySequence_GetItem(items.ptr(), static_cast<Py_ssize_t>(i));
          if (raw == nullptr) throw py::error_already_set();
          held.push_back(py::reinterpret_steal<py::object>(raw));
          native.push_back(action_of(held.back(), i));
        }
        if (!py::isinstance<BatchMode>(mode)) throw py::type_error("mode must be a BatchMode (got " + type_name(mode) + ")");
        const auto batch_mode = mode.cast<BatchMode>();
        const auto o = outputs_of(outputs);
        // All or nothing across the boundary: the list and its n result objects exist before the
        // call, and after the commit nothing allocates (each value moves into its object). The GIL
        // stays held: the call mutates the ledger (one writer per ledger).
        auto out = py::reinterpret_steal<py::typing::List<PlacementResult>>(checked(PyList_New(count)).release());
        std::vector<PlacementResult*> slots(n);
        for (std::size_t i = 0; i < n; ++i) {
          auto [obj, slot] = new_object<PlacementResult>();
          slots[i] = slot;
          PyList_SET_ITEM(out.ptr(), static_cast<Py_ssize_t>(i), obj.release().ptr());
        }
        auto results = apply_actions(ledger.cast<FlowLedger&>(), native, batch_mode, o);
        for (std::size_t i = 0; i < n; ++i) *slots[i] = std::move(results[i]);
        return out;
      },
      py::arg("ledger"), py::arg("actions"), py::kw_only(), py::arg_v("mode", BatchMode::Ordered, "BatchMode.ORDERED"),
      py::arg("outputs") = py::none(),
      "Run the actions of several policies of ledger (a sequence of Place, Replace and Forward) as one\n"
      "transaction -> list[PlacementResult], one per action in input order, each the scalar call's fields\n"
      "for its policy. The whole list is validated first, then computed, then published once (one exchange\n"
      "of every policy's records, then their epochs): a failed call raises with the ledger and every\n"
      "policy unchanged (with the two limits in the module docstring). A policy appears at most once, and\n"
      "every policy belongs to ledger.\n"
      "mode ORDERED: each action runs in list order on a private working state exactly as its scalar call\n"
      "would on the state the earlier actions left (their carried load and holds); new records take\n"
      "FlowIds in list order. Earlier actions are served first: a reservation gives no protection. mode SIMULTANEOUS:\n"
      "Replace of ShortestPaths EQUAL or WEIGHTED (routes prepared on the view) and of Pinned, Forward of\n"
      "any policy with an installed program, every one LOSSY with a finite offer, forwarded at the same\n"
      "time on one view: the residual plus every participant's old carried load (outside owners and orphans\n"
      "stay; holds are never added). At a link that cannot carry its arrivals the loss is shared by highest\n"
      "averages over the demands (ties by policy creation order, never by list position), then over each\n"
      "demand's sources, then its tunnels; each participant's black holes are its own. Permuting the list\n"
      "changes nothing. outputs (OutputOptions; None: none) as for the scalar calls: residual and\n"
      "saturated_boundary_edges read the state after that action (ORDERED) or after the whole list\n"
      "(SIMULTANEOUS, one residual for every entry, each entry's boundary from its own sources and mask).\n"
      "Checked in this order: the arguments' types (TypeError: `ledger must be a FlowLedger`, `actions must\n"
      "be a sequence of Place, Replace or Forward (got <type>)`, `actions[<i>] must be a Place, Replace or\n"
      "Forward (got <type>)`, `mode must be a BatchMode (got <type>)`, `outputs must be an OutputOptions\n"
      "(got <type>)`); mode's value (`mode: <v> is not a value of BatchMode`); an empty list returns [];\n"
      "then each action in list order, all of its checks before the next action's: `list policy belongs\n"
      "to another ledger`, `duplicate policy in list`, SIMULTANEOUS's own rules (`SIMULTANEOUS does not\n"
      "accept Place`, `SIMULTANEOUS Replace requires ShortestPaths EQUAL or WEIGHTED, or Pinned`,\n"
      "`SIMULTANEOUS requires admission=LOSSY`, `SIMULTANEOUS requires finite offers`), then the scalar\n"
      "call's checks in its order. Only then is anything computed; an error of the computation (a\n"
      "loss-iteration error, a solver or forwarding size limit, the record identity limit counted for the\n"
      "whole list) publishes nothing.");
}

} // namespace netgraph::core::bindings
