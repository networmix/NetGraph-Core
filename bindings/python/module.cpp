/*
  Pybind11 module exposing NetGraph-Core C++ APIs to Python.

  Notes:
    - Accepts NumPy arrays (C-contiguous) and converts to spans for zero-copy
      views where possible.
    - Distances returned as float64 arrays with inf for unreachable.
    - Edge/Node masks are validated for dtype=bool and length.
*/
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <cstring>

#include "netgraph/core/k_shortest_paths.hpp"
#include "netgraph/core/max_flow.hpp"
#include "netgraph/core/flow_state.hpp"
#include "netgraph/core/shortest_paths.hpp"
#include "netgraph/core/strict_multidigraph.hpp"
#include "netgraph/core/types.hpp"
#include "netgraph/core/flow_graph.hpp"
#include "netgraph/core/flow.hpp"
#include "netgraph/core/flow_policy.hpp"

namespace py = pybind11;
using namespace netgraph::core;

// Helpers to check NumPy arrays
template <typename T>
static std::span<const T> as_span(const py::array& arr, const char* name) {
  if (!py::isinstance<py::array_t<T>>(arr)) {
    throw py::type_error(std::string(name) + ": expected numpy array of correct dtype");
  }
  if (!(arr.flags() & py::array::c_style)) {
    throw py::type_error(std::string(name) + ": array must be C-contiguous (use np.ascontiguousarray)");
  }
  auto buf = arr.request();
  return std::span<const T>(static_cast<const T*>(buf.ptr), static_cast<std::size_t>(buf.size));
}

PYBIND11_MODULE(_netgraph_core, m) {
  m.doc() = "NetGraph-Core C++ bindings";

  py::enum_<EdgeTieBreak>(m, "EdgeTieBreak")
      .value("DETERMINISTIC", EdgeTieBreak::Deterministic)
      .value("PREFER_HIGHER_RESIDUAL", EdgeTieBreak::PreferHigherResidual);

  py::class_<EdgeSelection>(m, "EdgeSelection")
      .def(py::init<>())
      .def(py::init([](bool multipath, bool require_capacity, EdgeTieBreak tie_break){
        EdgeSelection s; s.multipath = multipath; s.require_capacity = require_capacity; s.tie_break = tie_break; return s;
      }),
        py::kw_only(),
        py::arg("multipath") = true,
        py::arg("require_capacity") = false,
        py::arg("tie_break") = EdgeTieBreak::Deterministic)
      .def_readwrite("multipath", &EdgeSelection::multipath)
      .def_readwrite("require_capacity", &EdgeSelection::require_capacity)
      .def_readwrite("tie_break", &EdgeSelection::tie_break);

  py::enum_<FlowPlacement>(m, "FlowPlacement")
      .value("PROPORTIONAL", FlowPlacement::Proportional)
      .value("EQUAL_BALANCED", FlowPlacement::EqualBalanced);

  py::class_<StrictMultiDiGraph>(m, "StrictMultiDiGraph")
      .def_static(
          "from_arrays",
          [](std::int32_t num_nodes,
             py::array src, py::array dst,
             py::array capacity, py::array cost,
             bool add_reverse) {
            // public API: src/dst are int32; pass through as int32 to core
            auto src_s = as_span<std::int32_t>(src, "src");
            auto dst_s = as_span<std::int32_t>(dst, "dst");
            if (src_s.size() != dst_s.size()) throw py::type_error("src and dst must have the same length");
            auto cap_s = as_span<double>(capacity, "capacity");
            // Accept any numeric dtype for cost; if float, require finite values, then force-cast to int64
            if (py::isinstance<py::array_t<double>>(cost)) {
              auto carr = py::cast<py::array_t<double>>(cost);
              if (!(carr.flags() & py::array::c_style)) {
                throw py::type_error("cost: array must be C-contiguous (use np.ascontiguousarray)");
              }
              auto cbuf_d = carr.request();
              if (cbuf_d.ndim != 1) throw py::type_error("cost must be a 1-D array");
              const double* cd = static_cast<const double*>(cbuf_d.ptr);
              for (ssize_t i = 0; i < cbuf_d.shape[0]; ++i) {
                if (!std::isfinite(cd[i])) {
                  throw py::value_error("cost values must be finite");
                }
              }
            }
            py::array_t<std::int64_t, py::array::c_style | py::array::forcecast> cost_i64(cost);
            auto cbuf = cost_i64.request();
            std::span<const Cost> cost_s(static_cast<const Cost*>(cbuf.ptr), static_cast<std::size_t>(cbuf.size));
            return StrictMultiDiGraph::from_arrays(num_nodes,
                                                  src_s,
                                                  dst_s,
                                                  cap_s, cost_s, add_reverse);
          },
          py::arg("num_nodes"), py::arg("src"), py::arg("dst"), py::arg("capacity"), py::arg("cost"),
          py::kw_only(), py::arg("add_reverse") = false)
      .def("num_nodes", &StrictMultiDiGraph::num_nodes)
      .def("num_edges", &StrictMultiDiGraph::num_edges)
      // external link ids removed; EdgeId is the canonical id
      .def("capacity_view", [](py::object self_obj, const StrictMultiDiGraph& g){
        auto s = g.capacity_view();
        return py::array(
            py::buffer_info(
                const_cast<double*>(s.data()),
                sizeof(double),
                py::format_descriptor<double>::format(),
                1,
                { s.size() },
                { sizeof(double) }
            ),
            self_obj
        );
      })
      .def("edge_src_view", [](py::object self_obj, const StrictMultiDiGraph& g){
        auto s = g.edge_src_view();
        return py::array(
            py::buffer_info(
                const_cast<std::int32_t*>(s.data()),
                sizeof(std::int32_t),
                py::format_descriptor<std::int32_t>::format(),
                1,
                { s.size() },
                { sizeof(std::int32_t) }
            ),
            self_obj
        );
      })
      .def("edge_dst_view", [](py::object self_obj, const StrictMultiDiGraph& g){
        auto s = g.edge_dst_view();
        return py::array(
            py::buffer_info(
                const_cast<std::int32_t*>(s.data()),
                sizeof(std::int32_t),
                py::format_descriptor<std::int32_t>::format(),
                1,
                { s.size() },
                { sizeof(std::int32_t) }
            ),
            self_obj
        );
      })
      .def("cost_view", [](py::object self_obj, const StrictMultiDiGraph& g){
        auto s = g.cost_view();
        return py::array(
            py::buffer_info(
                const_cast<Cost*>(s.data()),
                sizeof(Cost),
                py::format_descriptor<Cost>::format(),
                1,
                { s.size() },
                { sizeof(Cost) }
            ),
            self_obj
        );
      })
      .def("row_offsets_view", [](const StrictMultiDiGraph& g){
        auto s = g.row_offsets_view();
        py::array_t<std::int32_t> arr(s.size());
        std::memcpy(arr.mutable_data(), s.data(), s.size()*sizeof(std::int32_t));
        return arr;
      })
      .def("col_indices_view", [](const StrictMultiDiGraph& g){
        auto s = g.col_indices_view();
        py::array_t<std::int32_t> arr(s.size());
        std::memcpy(arr.mutable_data(), s.data(), s.size()*sizeof(std::int32_t));
        return arr;
      })
      .def("adj_edge_index_view", [](const StrictMultiDiGraph& g){
        auto s = g.adj_edge_index_view();
        py::array_t<std::int32_t> arr(s.size());
        std::memcpy(arr.mutable_data(), s.data(), s.size()*sizeof(std::int32_t));
        return arr;
      })
      .def("in_row_offsets_view", [](const StrictMultiDiGraph& g){
        auto s = g.in_row_offsets_view();
        py::array_t<std::int32_t> arr(s.size());
        std::memcpy(arr.mutable_data(), s.data(), s.size()*sizeof(std::int32_t));
        return arr;
      })
      .def("in_col_indices_view", [](const StrictMultiDiGraph& g){
        auto s = g.in_col_indices_view();
        py::array_t<std::int32_t> arr(s.size());
        std::memcpy(arr.mutable_data(), s.data(), s.size()*sizeof(std::int32_t));
        return arr;
      });

  py::class_<PredDAG>(m, "PredDAG")
      .def_property_readonly("parent_offsets", [](const PredDAG& d){
        py::array_t<std::int32_t> arr(d.parent_offsets.size());
        std::memcpy(arr.mutable_data(), d.parent_offsets.data(), d.parent_offsets.size()*sizeof(std::int32_t));
        return arr;
      })
      .def_property_readonly("parents", [](const PredDAG& d){
        py::array_t<std::int32_t> arr(d.parents.size());
        std::memcpy(arr.mutable_data(), d.parents.data(), d.parents.size()*sizeof(std::int32_t));
        return arr;
      })
      .def_property_readonly("via_edges", [](const PredDAG& d){
        py::array_t<std::int32_t> arr(d.via_edges.size());
        std::memcpy(arr.mutable_data(), d.via_edges.data(), d.via_edges.size()*sizeof(std::int32_t));
        return arr;
      });

  m.def("spf",
        [](const StrictMultiDiGraph& g, std::int32_t src, py::object dst,
           py::object selection_obj, py::object residual_obj, py::object node_mask, py::object edge_mask) {
          // Basic index validation
          if (src < 0 || src >= g.num_nodes()) {
            throw py::value_error("src out of range");
          }
          std::optional<NodeId> dst_opt;
          if (!dst.is_none()) {
            auto dval = static_cast<NodeId>(py::cast<std::int32_t>(dst));
            if (dval < 0 || dval >= g.num_nodes()) {
              throw py::value_error("dst out of range");
            }
            dst_opt = dval;
          }
          // Parse selection
          EdgeSelection selection;
          if (!selection_obj.is_none()) {
            selection = py::cast<EdgeSelection>(selection_obj);
          }
          // Parse optional residual
          std::span<const double> residual_span;
          std::vector<double> residual_vec;
          if (!residual_obj.is_none()) {
            auto arr = py::cast<py::array>(residual_obj);
            if (!(arr.flags() & py::array::c_style)) throw py::type_error("residual must be C-contiguous (np.ascontiguousarray)");
            auto buf = arr.request();
            if (buf.ndim != 1) throw py::type_error("residual must be a 1-D array of float64");
            if (buf.format != py::format_descriptor<double>::format()) throw py::type_error("residual must be dtype=float64");
            if (static_cast<std::size_t>(buf.shape[0]) != static_cast<std::size_t>(g.num_edges())) throw py::type_error("residual length must equal num_edges");
            residual_vec.resize(static_cast<std::size_t>(buf.shape[0]));
            std::memcpy(residual_vec.data(), buf.ptr, residual_vec.size()*sizeof(double));
            residual_span = std::span<const double>(residual_vec.data(), residual_vec.size());
          }
          // Validate masks and obtain raw pointers (optional)
          const bool* node_ptr = nullptr;
          const bool* edge_ptr = nullptr;
          py::array node_arr;
          py::array edge_arr;
          if (!node_mask.is_none()) {
            node_arr = py::cast<py::array>(node_mask);
            if (!(node_arr.flags() & py::array::c_style)) throw py::type_error("node_mask must be C-contiguous (np.ascontiguousarray)");
            auto buf = node_arr.request();
            if (buf.ndim != 1) throw py::type_error("node_mask must be a 1-D array of bool");
            if (buf.format != py::format_descriptor<bool>::format()) throw py::type_error("node_mask must be dtype=bool");
            if (static_cast<std::size_t>(buf.shape[0]) != static_cast<std::size_t>(g.num_nodes())) throw py::type_error("node_mask length must equal num_nodes");
            node_ptr = static_cast<const bool*>(buf.ptr);
          }
          if (!edge_mask.is_none()) {
            edge_arr = py::cast<py::array>(edge_mask);
            if (!(edge_arr.flags() & py::array::c_style)) throw py::type_error("edge_mask must be C-contiguous (np.ascontiguousarray)");
            auto buf = edge_arr.request();
            if (buf.ndim != 1) throw py::type_error("edge_mask must be a 1-D array of bool");
            if (buf.format != py::format_descriptor<bool>::format()) throw py::type_error("edge_mask must be dtype=bool");
            if (static_cast<std::size_t>(buf.shape[0]) != static_cast<std::size_t>(g.num_edges())) throw py::type_error("edge_mask length must equal num_edges");
            edge_ptr = static_cast<const bool*>(buf.ptr);
          }
          py::gil_scoped_release release;
          auto res = shortest_paths(g, src, dst_opt, selection, residual_span, node_ptr, edge_ptr);
          py::gil_scoped_acquire acquire;
          py::array_t<double> dist_arr(res.first.size());
          auto* out = dist_arr.mutable_data();
          auto maxc = std::numeric_limits<Cost>::max();
          for (std::size_t i=0;i<res.first.size();++i) out[i] = (res.first[i]==maxc) ? std::numeric_limits<double>::infinity() : static_cast<double>(res.first[i]);
          return py::make_tuple(std::move(dist_arr), res.second);
        }, py::arg("g"), py::arg("src"), py::arg("dst") = py::none(), py::kw_only(), py::arg("selection") = py::none(), py::arg("residual") = py::none(), py::arg("node_mask") = py::none(), py::arg("edge_mask") = py::none());

  m.def("ksp",
        [](const StrictMultiDiGraph& g, std::int32_t src, std::int32_t dst,
           int k, py::object max_cost_factor, bool unique, py::object node_mask, py::object edge_mask) {
          if (src < 0 || src >= g.num_nodes()) throw py::value_error("src out of range");
          if (dst < 0 || dst >= g.num_nodes()) throw py::value_error("dst out of range");
          if (k <= 0) throw py::value_error("k must be >= 1");
          std::optional<double> mcf;
          if (!max_cost_factor.is_none()) mcf = py::cast<double>(max_cost_factor);
          const bool* node_ptr = nullptr;
          const bool* edge_ptr = nullptr;
          py::array node_arr;
          py::array edge_arr;
          if (!node_mask.is_none()) {
            auto arr = py::cast<py::array>(node_mask);
            if (!(arr.flags() & py::array::c_style)) throw py::type_error("node_mask must be C-contiguous (np.ascontiguousarray)");
            auto buf = arr.request();
            if (buf.ndim != 1) throw py::type_error("node_mask must be a 1-D array of bool");
            if (buf.format != py::format_descriptor<bool>::format()) throw py::type_error("node_mask must be dtype=bool");
            if (static_cast<std::size_t>(buf.shape[0]) != static_cast<std::size_t>(g.num_nodes())) throw py::type_error("node_mask length must equal num_nodes");
            node_arr = arr; node_ptr = static_cast<const bool*>(buf.ptr);
          }
          if (!edge_mask.is_none()) {
            auto arr = py::cast<py::array>(edge_mask);
            if (!(arr.flags() & py::array::c_style)) throw py::type_error("edge_mask must be C-contiguous (np.ascontiguousarray)");
            auto buf = arr.request();
            if (buf.ndim != 1) throw py::type_error("edge_mask must be a 1-D array of bool");
            if (buf.format != py::format_descriptor<bool>::format()) throw py::type_error("edge_mask must be dtype=bool");
            if (static_cast<std::size_t>(buf.shape[0]) != static_cast<std::size_t>(g.num_edges())) throw py::type_error("edge_mask length must equal num_edges");
            edge_arr = arr; edge_ptr = static_cast<const bool*>(buf.ptr);
          }
          py::gil_scoped_release release;
          auto items = k_shortest_paths(g, src, dst, k, mcf, unique, node_ptr, edge_ptr);
          py::gil_scoped_acquire acquire;
          // Convert to list of (dist ndarray, PredDAG)
          py::list out;
          for (auto& pr : items) {
            auto& dist = pr.first;
            py::array_t<double> dist_arr(dist.size());
            auto* outp = dist_arr.mutable_data();
            auto maxc = std::numeric_limits<Cost>::max();
            for (std::size_t i=0;i<dist.size();++i) outp[i] = (dist[i]==maxc) ? std::numeric_limits<double>::infinity() : static_cast<double>(dist[i]);
            out.append(py::make_tuple(std::move(dist_arr), pr.second));
          }
          return out;
        }, py::arg("g"), py::arg("src"), py::arg("dst"), py::kw_only(), py::arg("k"), py::arg("max_cost_factor") = py::none(), py::arg("unique") = true, py::arg("node_mask") = py::none(), py::arg("edge_mask") = py::none());

  // resolve_to_paths: return list of paths; each path is list of (node, tuple(edge_ids...)) ending with (dst, ())
  m.def("resolve_to_paths",
        [](const PredDAG& dag, std::int32_t src, std::int32_t dst, bool split_parallel_edges, py::object max_paths){
          std::optional<std::int64_t> mp;
          if (!max_paths.is_none()) mp = py::cast<std::int64_t>(max_paths);
          auto out = netgraph::core::resolve_to_paths(dag, src, dst, split_parallel_edges, mp);
          py::list paths;
          for (auto& path : out) {
            py::list py_elems;
            for (auto& pr : path) {
              py::tuple edge_tuple(pr.second.size());
              for (std::size_t i=0;i<pr.second.size();++i) edge_tuple[i] = py::int_(pr.second[i]);
              py_elems.append(py::make_tuple(py::int_(pr.first), edge_tuple));
            }
            paths.append(py::tuple(py_elems));
          }
          return paths;
        }, py::arg("dag"), py::arg("src"), py::arg("dst"), py::kw_only(), py::arg("split_parallel_edges") = false, py::arg("max_paths") = py::none());

  py::class_<MinCut>(m, "MinCut")
      .def_property_readonly("edges", [](const MinCut& mc){
        py::array_t<std::int32_t> arr(mc.edges.size());
        if (!mc.edges.empty()) {
          std::memcpy(arr.mutable_data(), mc.edges.data(), mc.edges.size()*sizeof(std::int32_t));
        }
        return arr;
      });
  // CostBucket/CostDistribution removed in favor of parallel arrays on FlowSummary
  py::class_<FlowSummary>(m, "FlowSummary")
      .def_readonly("total_flow", &FlowSummary::total_flow)
      .def_readonly("min_cut", &FlowSummary::min_cut)
      .def_property_readonly("costs", [](const FlowSummary& s){
        py::array_t<std::int64_t> arr(s.costs.size());
        if (!s.costs.empty()) std::memcpy(arr.mutable_data(), s.costs.data(), s.costs.size()*sizeof(std::int64_t));
        return arr;
      })
      .def_property_readonly("flows", [](const FlowSummary& s){
        py::array_t<double> arr(s.flows.size());
        if (!s.flows.empty()) std::memcpy(arr.mutable_data(), s.flows.data(), s.flows.size()*sizeof(double));
        return arr;
      })
      .def_property_readonly("edge_flows", [](const FlowSummary& s){
        py::array_t<double> arr(s.edge_flows.size());
        if (!s.edge_flows.empty()) std::memcpy(arr.mutable_data(), s.edge_flows.data(), s.edge_flows.size()*sizeof(double));
        return arr;
      })
      .def_property_readonly("residual_capacity", [](const FlowSummary& s){
        py::array_t<double> arr(s.residual_capacity.size());
        if (!s.residual_capacity.empty()) std::memcpy(arr.mutable_data(), s.residual_capacity.data(), s.residual_capacity.size()*sizeof(double));
        return arr;
      })
      .def_property_readonly("reachable_nodes", [](const FlowSummary& s){
        // Store as bool array for Python; cast bytes to bool
        py::array_t<bool> arr(s.reachable_nodes.size());
        auto* out = arr.mutable_data();
        for (std::size_t i=0;i<s.reachable_nodes.size();++i) out[i] = static_cast<bool>(s.reachable_nodes[i]);
        return arr;
      });

  // spf_residual removed; unified spf accepts optional residual and EdgeSelection

  m.def("max_flow",
        [](const StrictMultiDiGraph& g, std::int32_t src, std::int32_t dst,
           FlowPlacement placement, bool shortest_path, bool with_edge_flows,
           bool with_reachable, bool with_residuals,
           py::object node_mask, py::object edge_mask) {
          const bool* node_ptr = nullptr;
          const bool* edge_ptr = nullptr;
          py::array node_arr, edge_arr;
          if (!node_mask.is_none()) {
            node_arr = py::cast<py::array>(node_mask);
            if (!(node_arr.flags() & py::array::c_style)) throw py::type_error("node_mask must be C-contiguous (np.ascontiguousarray)");
            auto buf = node_arr.request();
            if (buf.ndim != 1) throw py::type_error("node_mask must be a 1-D array of bool");
            if (buf.format != py::format_descriptor<bool>::format()) throw py::type_error("node_mask must be dtype=bool");
            if (static_cast<std::size_t>(buf.shape[0]) != static_cast<std::size_t>(g.num_nodes())) throw py::type_error("node_mask length must equal num_nodes");
            node_ptr = static_cast<const bool*>(buf.ptr);
          }
          if (!edge_mask.is_none()) {
            edge_arr = py::cast<py::array>(edge_mask);
            if (!(edge_arr.flags() & py::array::c_style)) throw py::type_error("edge_mask must be C-contiguous (np.ascontiguousarray)");
            auto buf = edge_arr.request();
            if (buf.ndim != 1) throw py::type_error("edge_mask must be a 1-D array of bool");
            if (buf.format != py::format_descriptor<bool>::format()) throw py::type_error("edge_mask must be dtype=bool");
            if (static_cast<std::size_t>(buf.shape[0]) != static_cast<std::size_t>(g.num_edges())) throw py::type_error("edge_mask length must equal num_edges");
            edge_ptr = static_cast<const bool*>(buf.ptr);
          }
          py::gil_scoped_release release;
          auto res = calc_max_flow(g, src, dst, placement, shortest_path, with_edge_flows, with_reachable, with_residuals, node_ptr, edge_ptr);
          py::gil_scoped_acquire acquire;
          return res;
        }, py::arg("g"), py::arg("src"), py::arg("dst"), py::kw_only(), py::arg("flow_placement") = FlowPlacement::Proportional, py::arg("shortest_path") = false, py::arg("with_edge_flows") = false, py::arg("with_reachable") = false, py::arg("with_residuals") = false, py::arg("node_mask") = py::none(), py::arg("edge_mask") = py::none());

  m.def("batch_max_flow",
        [](const StrictMultiDiGraph& g, py::array pairs,
           py::object node_masks, py::object edge_masks,
           FlowPlacement placement, bool shortest_path, bool with_edge_flows,
           bool with_reachable, bool with_residuals) {
          auto buf = pairs.request();
          if (buf.ndim != 2 || buf.shape[1] != 2) throw py::type_error("pairs must be shape [B,2]");
          // accept int32 for pairs
          if (buf.format != py::format_descriptor<std::int32_t>::format()) throw py::type_error("pairs dtype must be int32");
          const std::size_t batch_size = static_cast<std::size_t>(buf.shape[0]);
          std::vector<std::pair<NodeId,NodeId>> pp;
          auto* ptr = static_cast<const std::int32_t*>(buf.ptr);
          pp.reserve(batch_size);
          for (ssize_t i=0;i<buf.shape[0];++i) {
            if (ptr[2*i] < 0 || ptr[2*i+1] < 0) throw py::type_error("pairs must be non-negative indices");
            pp.emplace_back(static_cast<std::int32_t>(ptr[2*i]), static_cast<std::int32_t>(ptr[2*i+1]));
          }
          // Parse mask lists if provided, collect raw pointers and keep arrays alive
          auto parse_mask_list = [&](py::object list_obj, std::size_t expected_len, const char* what, std::vector<const bool*>& out_ptrs, std::vector<py::array>& keep){
            if (list_obj.is_none()) return;
            auto seq = py::cast<py::sequence>(list_obj);
            if (static_cast<std::size_t>(py::len(seq)) != batch_size) throw py::type_error(std::string(what) + " length must equal number of pairs");
            for (auto item: seq) {
              auto arr = py::cast<py::array>(item);
              if (!(arr.flags() & py::array::c_style)) throw py::type_error(std::string(what) + " arrays must be C-contiguous (np.ascontiguousarray)");
              auto b = arr.request();
              if (b.ndim != 1) throw py::type_error(std::string(what) + " must be a list of 1-D arrays");
              if (b.format != py::format_descriptor<bool>::format()) throw py::type_error(std::string(what) + " arrays must be dtype=bool");
              if (static_cast<std::size_t>(b.shape[0]) != expected_len) throw py::type_error(std::string(what) + " array has wrong length");
              keep.push_back(arr);
              out_ptrs.push_back(static_cast<const bool*>(b.ptr));
            }
          };
          std::vector<const bool*> node_ptrs;
          std::vector<const bool*> edge_ptrs;
          std::vector<py::array> node_keep;
          std::vector<py::array> edge_keep;
          parse_mask_list(node_masks, static_cast<std::size_t>(g.num_nodes()), "node_masks", node_ptrs, node_keep);
          parse_mask_list(edge_masks, static_cast<std::size_t>(g.num_edges()), "edge_masks", edge_ptrs, edge_keep);
          py::gil_scoped_release release;
          auto out = netgraph::core::batch_max_flow(g, pp, placement, shortest_path, with_edge_flows, with_reachable, with_residuals, node_ptrs, edge_ptrs);
          py::gil_scoped_acquire acquire;
          return out;
        }, py::arg("g"), py::arg("pairs"), py::kw_only(), py::arg("node_masks") = py::none(), py::arg("edge_masks") = py::none(), py::arg("flow_placement") = FlowPlacement::Proportional, py::arg("shortest_path") = false, py::arg("with_edge_flows") = false, py::arg("with_reachable") = false, py::arg("with_residuals") = false);

  // FlowState bindings
  py::class_<FlowState>(m, "FlowState")
      .def(py::init<const StrictMultiDiGraph&>())
      .def(py::init([](const StrictMultiDiGraph& g, py::array residual){
        if (!(py::isinstance<py::array_t<double>>(residual))) throw py::type_error("residual must be a numpy float64 array");
        if (!(residual.flags() & py::array::c_style)) throw py::type_error("residual must be C-contiguous");
        auto buf = residual.request();
        if (buf.ndim != 1 || static_cast<std::size_t>(buf.shape[0]) != static_cast<std::size_t>(g.num_edges())) throw py::type_error("residual length must equal num_edges");
        std::span<const double> rspan(static_cast<const double*>(buf.ptr), static_cast<std::size_t>(buf.shape[0]));
        return FlowState(g, rspan);
      }))
      .def("reset", [](FlowState& fs){ fs.reset(); })
      .def("reset", [](FlowState& fs, py::array residual){
        if (!(py::isinstance<py::array_t<double>>(residual))) throw py::type_error("residual must be a numpy float64 array");
        if (!(residual.flags() & py::array::c_style)) throw py::type_error("residual must be C-contiguous");
        auto buf = residual.request();
        if (buf.ndim != 1 || static_cast<std::size_t>(buf.shape[0]) != static_cast<std::size_t>(fs.capacity_view().size())) throw py::type_error("residual length must equal num_edges");
        std::span<const double> rspan(static_cast<const double*>(buf.ptr), static_cast<std::size_t>(buf.shape[0]));
        fs.reset(rspan);
      })
      .def("capacity_view", [](py::object self_obj, const FlowState& fs){
        auto s = fs.capacity_view();
        return py::array(
            py::buffer_info(
                const_cast<double*>(s.data()), sizeof(double), py::format_descriptor<double>::format(), 1, { s.size() }, { sizeof(double) }
            ), self_obj);
      })
      .def("residual_view", [](py::object self_obj, const FlowState& fs){
        auto s = fs.residual_view();
        return py::array(
            py::buffer_info(
                const_cast<double*>(s.data()), sizeof(double), py::format_descriptor<double>::format(), 1, { s.size() }, { sizeof(double) }
            ), self_obj);
      })
      .def("edge_flow_view", [](py::object self_obj, const FlowState& fs){
        auto s = fs.edge_flow_view();
        return py::array(
            py::buffer_info(
                const_cast<double*>(s.data()), sizeof(double), py::format_descriptor<double>::format(), 1, { s.size() }, { sizeof(double) }
            ), self_obj);
      })
      .def("place_on_dag", [](FlowState& fs, std::int32_t src, std::int32_t dst, const PredDAG& dag, double requested_flow, FlowPlacement placement, bool shortest_path){
        py::gil_scoped_release rel; auto placed = fs.place_on_dag(src, dst, dag, requested_flow, placement, shortest_path); py::gil_scoped_acquire acq; return placed;
      }, py::arg("src"), py::arg("dst"), py::arg("dag"), py::arg("requested_flow") = std::numeric_limits<double>::infinity(), py::arg("flow_placement") = FlowPlacement::Proportional, py::arg("shortest_path") = false)
      .def("place_max_flow", [](FlowState& fs, std::int32_t src, std::int32_t dst, FlowPlacement placement, bool shortest_path){
        py::gil_scoped_release rel; auto total = fs.place_max_flow(src, dst, placement, shortest_path); py::gil_scoped_acquire acq; return total;
      }, py::arg("src"), py::arg("dst"), py::arg("flow_placement") = FlowPlacement::Proportional, py::arg("shortest_path") = false)
      .def("compute_min_cut", [](const FlowState& fs, std::int32_t src, py::object node_mask, py::object edge_mask){
        const bool* node_ptr = nullptr; const bool* edge_ptr = nullptr; py::array node_arr, edge_arr;
        if (!node_mask.is_none()) {
          node_arr = py::cast<py::array>(node_mask);
          if (!(node_arr.flags() & py::array::c_style)) throw py::type_error("node_mask must be C-contiguous (np.ascontiguousarray)");
          auto b = node_arr.request();
          if (b.ndim != 1 || b.format != py::format_descriptor<bool>::format()) throw py::type_error("node_mask must be 1-D bool");
          node_ptr = static_cast<const bool*>(b.ptr);
        }
        if (!edge_mask.is_none()) {
          edge_arr = py::cast<py::array>(edge_mask);
          if (!(edge_arr.flags() & py::array::c_style)) throw py::type_error("edge_mask must be C-contiguous (np.ascontiguousarray)");
          auto b = edge_arr.request();
          if (b.ndim != 1 || b.format != py::format_descriptor<bool>::format()) throw py::type_error("edge_mask must be 1-D bool");
          edge_ptr = static_cast<const bool*>(b.ptr);
        }
        py::gil_scoped_release rel; auto mc = fs.compute_min_cut(src, node_ptr, edge_ptr); py::gil_scoped_acquire acq; return mc;
      }, py::arg("src"), py::kw_only(), py::arg("node_mask") = py::none(), py::arg("edge_mask") = py::none());

  // FlowIndex and FlowGraph bindings
  py::class_<FlowIndex>(m, "FlowIndex")
      .def(py::init<NodeId, NodeId, std::int32_t, std::int64_t>())
      .def_readonly("src", &FlowIndex::src)
      .def_readonly("dst", &FlowIndex::dst)
      .def_readonly("flowClass", &FlowIndex::flowClass)
      .def_readonly("flowId", &FlowIndex::flowId);

  py::class_<FlowGraph>(m, "FlowGraph")
      .def(py::init<const StrictMultiDiGraph&>())
      .def("capacity_view", [](py::object self_obj, const FlowGraph& fg){ auto s = fg.capacity_view(); return py::array(py::buffer_info(const_cast<double*>(s.data()), sizeof(double), py::format_descriptor<double>::format(), 1, { s.size() }, { sizeof(double) }), self_obj); })
      .def("residual_view", [](py::object self_obj, const FlowGraph& fg){ auto s = fg.residual_view(); return py::array(py::buffer_info(const_cast<double*>(s.data()), sizeof(double), py::format_descriptor<double>::format(), 1, { s.size() }, { sizeof(double) }), self_obj); })
      .def("edge_flow_view", [](py::object self_obj, const FlowGraph& fg){ auto s = fg.edge_flow_view(); return py::array(py::buffer_info(const_cast<double*>(s.data()), sizeof(double), py::format_descriptor<double>::format(), 1, { s.size() }, { sizeof(double) }), self_obj); })
      .def_property_readonly("graph", [](const FlowGraph& fg){ return &fg.graph(); }, py::return_value_policy::reference)
      .def("place", [](FlowGraph& fg, const FlowIndex& idx, std::int32_t src, std::int32_t dst, const PredDAG& dag, double amount, FlowPlacement placement, bool shortest_path){ py::gil_scoped_release rel; auto placed = fg.place(idx, src, dst, dag, amount, placement, shortest_path); py::gil_scoped_acquire acq; return placed; }, py::arg("index"), py::arg("src"), py::arg("dst"), py::arg("dag"), py::arg("amount"), py::arg("flow_placement") = FlowPlacement::Proportional, py::arg("shortest_path") = false)
      .def("remove", [](FlowGraph& fg, const FlowIndex& idx){ py::gil_scoped_release rel; fg.remove(idx); py::gil_scoped_acquire acq; })
      .def("remove_by_class", [](FlowGraph& fg, std::int32_t cls){ py::gil_scoped_release rel; fg.remove_by_class(cls); py::gil_scoped_acquire acq; })
      .def("reset", [](FlowGraph& fg){ py::gil_scoped_release rel; fg.reset(); py::gil_scoped_acquire acq; })
      .def("get_flow_edges", [](const FlowGraph& fg, const FlowIndex& idx){ auto v = fg.get_flow_edges(idx); py::list out; for (auto const& pr : v) { out.append(py::make_tuple(pr.first, pr.second)); } return out; })
      .def("get_flow_path", [](const FlowGraph& fg, const FlowIndex& idx){ auto v = fg.get_flow_path(idx); py::list out; for (auto e : v) out.append(e); return out; });

  py::class_<FlowRecord>(m, "Flow")
      .def_property_readonly("index", [](const FlowRecord& f){ return f.index; })
      .def_readonly("src", &FlowRecord::src)
      .def_readonly("dst", &FlowRecord::dst)
      .def_readonly("cost", &FlowRecord::cost)
      .def_readonly("placed_flow", &FlowRecord::placed_flow);

  py::enum_<PathAlg>(m, "PathAlg")
      .value("SPF", PathAlg::SPF);

  py::class_<FlowPolicy>(m, "FlowPolicy")
      .def(py::init<PathAlg, FlowPlacement, EdgeSelection, int, std::optional<int>, std::optional<Cost>, std::optional<double>, bool, bool, int, int, bool, int, double>(),
           py::arg("path_alg") = PathAlg::SPF,
           py::arg("flow_placement") = FlowPlacement::Proportional,
           py::arg("selection") = EdgeSelection{},
           py::arg("min_flow_count") = 1,
           py::arg("max_flow_count") = py::none(),
           py::arg("max_path_cost") = py::none(),
           py::arg("max_path_cost_factor") = py::none(),
           py::arg("shortest_path") = false,
           py::arg("reoptimize_flows_on_each_placement") = false,
           py::arg("max_no_progress_iterations") = 100,
           py::arg("max_total_iterations") = 10000,
           py::arg("diminishing_returns_enabled") = true,
           py::arg("diminishing_returns_window") = 8,
           py::arg("diminishing_returns_epsilon_frac") = 1e-3)
      .def("flow_count", &FlowPolicy::flow_count)
      .def("placed_demand", &FlowPolicy::placed_demand)
      .def("place_demand", [](FlowPolicy& p, FlowGraph& fg, std::int32_t src, std::int32_t dst, std::int32_t flowClass, double volume, py::object target_per_flow, py::object min_flow){ std::optional<double> tpf; if (!target_per_flow.is_none()) tpf = py::cast<double>(target_per_flow); std::optional<double> mfl; if (!min_flow.is_none()) mfl = py::cast<double>(min_flow); py::gil_scoped_release rel; auto pr = p.place_demand(fg, src, dst, flowClass, volume, tpf, mfl); py::gil_scoped_acquire acq; return py::make_tuple(pr.first, pr.second); }, py::arg("flow_graph"), py::arg("src"), py::arg("dst"), py::arg("flowClass"), py::arg("volume"), py::arg("target_per_flow") = py::none(), py::arg("min_flow") = py::none())
      .def("rebalance_demand", [](FlowPolicy& p, FlowGraph& fg, std::int32_t src, std::int32_t dst, std::int32_t flowClass, double target){ py::gil_scoped_release rel; auto pr = p.rebalance_demand(fg, src, dst, flowClass, target); py::gil_scoped_acquire acq; return py::make_tuple(pr.first, pr.second); },
           py::arg("flow_graph"), py::arg("src"), py::arg("dst"), py::arg("flowClass"), py::arg("target"))
      .def("remove_demand", [](FlowPolicy& p, FlowGraph& fg){ py::gil_scoped_release rel; p.remove_demand(fg); py::gil_scoped_acquire acq; })
      .def_property_readonly("flows", [](const FlowPolicy& p){ py::dict out; for (auto const& kv : p.flows()) { const auto& idx = kv.first; const auto& f = kv.second; out[py::make_tuple(idx.src, idx.dst, idx.flowClass, idx.flowId)] = py::make_tuple(f.src, f.dst, f.cost, f.placed_flow); } return out; });
}
