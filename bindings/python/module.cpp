/*
  Pybind11 module exposing NetGraph-Core C++ APIs to Python.
  - Uses NumPy arrays with zero-copy spans where possible for FlowState/FlowGraph views (read-only)
  - Returns distances as float64 arrays by default (inf for unreachable); dtype='int64' available
  - Validates edge/node masks for bool dtype and length
  - StrictMultiDiGraph view methods return copies (by design) to avoid exposing lifetime
*/
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <cstring>
#include <memory>

#include "netgraph/core/k_shortest_paths.hpp"
#include "netgraph/core/max_flow.hpp"
#include "netgraph/core/flow_state.hpp"
#include "netgraph/core/shortest_paths.hpp"
#include "netgraph/core/strict_multidigraph.hpp"
#include "netgraph/core/types.hpp"
#include "netgraph/core/flow_graph.hpp"
#include "netgraph/core/flow.hpp"
#include "netgraph/core/flow_policy.hpp"
#include "netgraph/core/backend.hpp"
#include "netgraph/core/algorithms.hpp"
#include "netgraph/core/options.hpp"
#include "netgraph/core/profiling.hpp"

namespace py = pybind11;
using namespace netgraph::core;

// Helper opaque wrappers (defined at file scope to ensure stable type identity)
struct PyBackend { BackendPtr impl; };
struct PyGraph { GraphHandle handle; py::object keep_alive; std::int32_t num_nodes; std::int32_t num_edges; };

// Helpers to check NumPy arrays
namespace {
// Small helper: copy a span into a new NumPy array of the same dtype.
template <typename T>
py::array_t<T> copy_to_numpy(std::span<const T> s) {
  py::array_t<T> arr(s.size());
  if (!s.empty()) std::memcpy(arr.mutable_data(), s.data(), s.size()*sizeof(T));
  return arr;
}

// Owns a bool[] buffer and exposes a const span view onto it.
struct BoolSpan {
  std::unique_ptr<bool[]> buf;
  std::span<const bool> view;
};

// Validate 1-D C-contig NumPy bool array and deep-copy into bool[].
static BoolSpan to_bool_span_from_numpy(py::object obj, std::size_t expected_len, const char* what) {
  BoolSpan out;
  if (obj.is_none()) return out;  // empty span
  py::array arr = obj.cast<py::array>();
  if (!(arr.flags() & py::array::c_style)) {
    throw py::type_error(std::string(what) + " must be C-contiguous");
  }
  auto b = arr.request();
  if (b.ndim != 1 || b.format != py::format_descriptor<bool>::format()) {
    throw py::type_error(std::string(what) + " must be 1-D bool");
  }
  if (static_cast<std::size_t>(b.shape[0]) != expected_len) {
    throw py::type_error(std::string(what) + " length must equal " + std::to_string(expected_len));
  }
  out.buf.reset(new bool[expected_len]);
  // NumPy bool_ uses 1 byte; C++ bool is 1 byte on supported platforms.
  std::memcpy(out.buf.get(), b.ptr, expected_len * sizeof(bool));
  out.view = std::span<const bool>(out.buf.get(), expected_len);
  return out;
}
} // namespace

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
  if (buf.ndim != 1) {
    throw py::type_error(std::string(name) + ": expected 1-D array");
  }
  return std::span<const T>(static_cast<const T*>(buf.ptr), static_cast<std::size_t>(buf.size));
}

PYBIND11_MODULE(_netgraph_core, m) {
  m.doc() = "NetGraph-Core C++ bindings";

  py::enum_<EdgeTieBreak>(m, "EdgeTieBreak")
      .value("DETERMINISTIC", EdgeTieBreak::Deterministic)
      .value("PREFER_HIGHER_RESIDUAL", EdgeTieBreak::PreferHigherResidual);

  py::class_<EdgeSelection>(m, "EdgeSelection")
      .def(py::init<>())
      .def(py::init([](bool multi_edge, bool require_capacity, EdgeTieBreak tie_break){
        EdgeSelection s; s.multi_edge = multi_edge; s.require_capacity = require_capacity; s.tie_break = tie_break; return s;
      }),
        py::kw_only(),
        py::arg("multi_edge") = true,
        py::arg("require_capacity") = false,
        py::arg("tie_break") = EdgeTieBreak::Deterministic)
      .def_readwrite("multi_edge", &EdgeSelection::multi_edge)
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
             py::array ext_edge_ids) {
            // public API: src/dst are int32; pass through as int32 to core
            auto src_s = as_span<std::int32_t>(src, "src");
            auto dst_s = as_span<std::int32_t>(dst, "dst");
            if (src_s.size() != dst_s.size()) throw py::type_error("src and dst must have the same length");
            auto cap_s = as_span<double>(capacity, "capacity");
            // Cost dtype must be int64 to match internal Cost type
            auto cost_s = as_span<std::int64_t>(cost, "cost");
            // ext_edge_ids required: maps internal edge indices to caller's stable IDs
            auto ext_s = as_span<std::int64_t>(ext_edge_ids, "ext_edge_ids");
            return StrictMultiDiGraph::from_arrays(num_nodes,
                                                  src_s,
                                                  dst_s,
                                                  cap_s, cost_s, ext_s);
          },
          py::arg("num_nodes"), py::arg("src"), py::arg("dst"), py::arg("capacity"), py::arg("cost"),
          py::arg("ext_edge_ids"))
      .def("num_nodes", &StrictMultiDiGraph::num_nodes)
      .def("num_edges", &StrictMultiDiGraph::num_edges)
      .def("capacity_view", [](const StrictMultiDiGraph& g){ return copy_to_numpy<double>(g.capacity_view()); })
      .def("edge_src_view", [](const StrictMultiDiGraph& g){ return copy_to_numpy<std::int32_t>(g.edge_src_view()); })
      .def("edge_dst_view", [](const StrictMultiDiGraph& g){ return copy_to_numpy<std::int32_t>(g.edge_dst_view()); })
      .def("ext_edge_ids_view", [](const StrictMultiDiGraph& g){ return copy_to_numpy<std::int64_t>(g.ext_edge_ids_view()); })
      .def("cost_view", [](const StrictMultiDiGraph& g){ return copy_to_numpy<std::int64_t>(g.cost_view()); })
      .def("row_offsets_view", [](const StrictMultiDiGraph& g){ return copy_to_numpy<std::int32_t>(g.row_offsets_view()); })
      .def("col_indices_view", [](const StrictMultiDiGraph& g){ return copy_to_numpy<std::int32_t>(g.col_indices_view()); })
      .def("adj_edge_index_view", [](const StrictMultiDiGraph& g){ return copy_to_numpy<std::int32_t>(g.adj_edge_index_view()); })
      .def("in_row_offsets_view", [](const StrictMultiDiGraph& g){ return copy_to_numpy<std::int32_t>(g.in_row_offsets_view()); })
      .def("in_col_indices_view", [](const StrictMultiDiGraph& g){ return copy_to_numpy<std::int32_t>(g.in_col_indices_view()); })
      .def("in_adj_edge_index_view", [](const StrictMultiDiGraph& g){ return copy_to_numpy<std::int32_t>(g.in_adj_edge_index_view()); });

  // Backend and Algorithms

  // Opaque wrapper types
  py::class_<PyBackend>(m, "Backend")
      .def_static("cpu", [](){ return PyBackend{ make_cpu_backend() }; });

  py::class_<PyGraph>(m, "Graph")
      // Opaque holder; constructed via Algorithms.build_graph or build_graph_from_arrays
      .def_property_readonly("num_nodes", [](const PyGraph& pg){ return pg.num_nodes; })
      .def_property_readonly("num_edges", [](const PyGraph& pg){ return pg.num_edges; });

  py::class_<Algorithms, std::shared_ptr<Algorithms>>(m, "Algorithms")
      .def(py::init([](const PyBackend& b){ return std::make_shared<Algorithms>(b.impl); }))
      .def("build_graph", [](const Algorithms& algs, py::object graph_obj){
        const StrictMultiDiGraph& g = graph_obj.cast<const StrictMultiDiGraph&>();
        auto gh = algs.build_graph(g);
        // Keep the Python graph object alive alongside the handle when
        // referencing a non-owned graph instance.
        return PyGraph{ gh, graph_obj, g.num_nodes(), g.num_edges() };
      }, py::arg("graph"))
      .def("build_graph_from_arrays", [](const Algorithms& algs,
                                          std::int32_t num_nodes,
                                          py::array src, py::array dst,
                                          py::array capacity, py::array cost,
                                          py::array ext_edge_ids){
        // Build graph and construct shared ownership directly
        auto ext_s = as_span<std::int64_t>(ext_edge_ids, "ext_edge_ids");
        auto sp = std::make_shared<StrictMultiDiGraph>(
            StrictMultiDiGraph::from_arrays(
                num_nodes,
                as_span<std::int32_t>(src, "src"),
                as_span<std::int32_t>(dst, "dst"),
                as_span<double>(capacity, "capacity"),
                as_span<std::int64_t>(cost, "cost"),
                ext_s));
        auto gh = algs.build_graph(std::static_pointer_cast<const StrictMultiDiGraph>(sp));
        // GraphHandle holds shared ownership; no additional Python-side owner needed
        return PyGraph{ gh, py::none(), sp->num_nodes(), sp->num_edges() };
      }, py::arg("num_nodes"), py::arg("src"), py::arg("dst"), py::arg("capacity"), py::arg("cost"), py::arg("ext_edge_ids"))
      .def("spf", [](const Algorithms& algs, const PyGraph& pg, std::int32_t src,
                       py::object dst, py::object selection_obj, py::object residual_obj,
                       py::object node_mask, py::object edge_mask, bool multipath, std::string dtype){
        if (src < 0 || src >= pg.num_nodes) throw py::value_error("src out of range");
        SpfOptions opts; if (!selection_obj.is_none()) opts.selection = py::cast<EdgeSelection>(selection_obj);
        if (!dst.is_none()) { auto d = py::cast<std::int32_t>(dst); if (d < 0 || d >= pg.num_nodes) throw py::value_error("dst out of range"); opts.dst = static_cast<NodeId>(d); }
        std::vector<double> residual_vec;
        if (!residual_obj.is_none()) {
          auto arr = py::cast<py::array>(residual_obj);
          if (!(arr.flags() & py::array::c_style)) throw py::type_error("residual must be C-contiguous (np.ascontiguousarray)");
          auto buf = arr.request();
          if (buf.ndim != 1 || buf.format != py::format_descriptor<double>::format()) throw py::type_error("residual must be 1-D float64");
          residual_vec.resize(static_cast<std::size_t>(buf.shape[0]));
          std::memcpy(residual_vec.data(), buf.ptr, residual_vec.size()*sizeof(double));
          opts.residual = std::span<const double>(residual_vec.data(), residual_vec.size());
        }
        auto node_bs = to_bool_span_from_numpy(node_mask, static_cast<std::size_t>(pg.num_nodes), "node_mask");
        auto edge_bs = to_bool_span_from_numpy(edge_mask, static_cast<std::size_t>(pg.num_edges), "edge_mask");
        opts.node_mask = node_bs.view;
        opts.edge_mask = edge_bs.view;
        opts.multipath = multipath;
        py::gil_scoped_release rel; auto res = algs.spf(pg.handle, src, opts); py::gil_scoped_acquire acq;
        auto maxc = std::numeric_limits<Cost>::max();
        if (dtype == "int64") {
          py::array_t<std::int64_t> dist_arr(res.first.size());
          auto* out = dist_arr.mutable_data();
          for (std::size_t i=0;i<res.first.size();++i) out[i] = (res.first[i]==maxc) ? static_cast<std::int64_t>(maxc) : static_cast<std::int64_t>(res.first[i]);
          return py::make_tuple(std::move(dist_arr), res.second);
        } else if (dtype == "float64") {
          py::array_t<double> dist_arr(res.first.size());
          auto* out = dist_arr.mutable_data();
          for (std::size_t i=0;i<res.first.size();++i) out[i] = (res.first[i]==maxc) ? std::numeric_limits<double>::infinity() : static_cast<double>(res.first[i]);
          return py::make_tuple(std::move(dist_arr), res.second);
        } else {
          throw py::value_error("dtype must be 'float64' or 'int64'");
        }
      }, py::arg("graph"), py::arg("src"), py::arg("dst") = py::none(), py::kw_only(), py::arg("selection") = py::none(), py::arg("residual") = py::none(), py::arg("node_mask") = py::none(), py::arg("edge_mask") = py::none(), py::arg("multipath") = true, py::arg("dtype") = "float64")
      .def("ksp", [](const Algorithms& algs, const PyGraph& pg, std::int32_t src, std::int32_t dst,
                       int k, py::object max_cost_factor, bool unique, py::object node_mask, py::object edge_mask, std::string dtype){
        if (src < 0 || src >= pg.num_nodes || dst < 0 || dst >= pg.num_nodes) throw py::value_error("src/dst out of range");
        if (k <= 0) throw py::value_error("k must be >= 1");
        KspOptions opts; opts.k = k; opts.unique = unique; if (!max_cost_factor.is_none()) opts.max_cost_factor = py::cast<double>(max_cost_factor);
        auto node_bs = to_bool_span_from_numpy(node_mask, static_cast<std::size_t>(pg.num_nodes), "node_mask");
        auto edge_bs = to_bool_span_from_numpy(edge_mask, static_cast<std::size_t>(pg.num_edges), "edge_mask");
        opts.node_mask = node_bs.view; opts.edge_mask = edge_bs.view;
        py::gil_scoped_release rel; auto items = algs.ksp(pg.handle, src, dst, opts); py::gil_scoped_acquire acq;
        py::list out;
        auto maxc = std::numeric_limits<Cost>::max();
        for (auto& pr : items) {
          auto& dist = pr.first;
          if (dtype == "int64") {
            py::array_t<std::int64_t> dist_arr(dist.size());
            auto* outp = dist_arr.mutable_data();
            for (std::size_t i=0;i<dist.size();++i) outp[i] = (dist[i]==maxc) ? static_cast<std::int64_t>(maxc) : static_cast<std::int64_t>(dist[i]);
            out.append(py::make_tuple(std::move(dist_arr), pr.second));
          } else if (dtype == "float64") {
            py::array_t<double> dist_arr(dist.size());
            auto* outp = dist_arr.mutable_data();
            for (std::size_t i=0;i<dist.size();++i) outp[i] = (dist[i]==maxc) ? std::numeric_limits<double>::infinity() : static_cast<double>(dist[i]);
            out.append(py::make_tuple(std::move(dist_arr), pr.second));
          } else {
            throw py::value_error("dtype must be 'float64' or 'int64'");
          }
        }
        return out;
      }, py::arg("graph"), py::arg("src"), py::arg("dst"), py::kw_only(), py::arg("k"), py::arg("max_cost_factor") = py::none(), py::arg("unique") = true, py::arg("node_mask") = py::none(), py::arg("edge_mask") = py::none(), py::arg("dtype") = "float64")
      .def("max_flow", [](const Algorithms& algs, const PyGraph& pg, std::int32_t src, std::int32_t dst,
                            FlowPlacement placement, bool shortest_path, bool require_capacity, bool with_edge_flows, bool with_reachable, bool with_residuals,
                            py::object node_mask, py::object edge_mask){
        if (src < 0 || src >= pg.num_nodes || dst < 0 || dst >= pg.num_nodes) throw py::value_error("src/dst out of range");
        MaxFlowOptions o; o.placement = placement; o.shortest_path = shortest_path; o.require_capacity = require_capacity; o.with_edge_flows = with_edge_flows; o.with_reachable = with_reachable; o.with_residuals = with_residuals;
        auto node_bs = to_bool_span_from_numpy(node_mask, static_cast<std::size_t>(pg.num_nodes), "node_mask");
        auto edge_bs = to_bool_span_from_numpy(edge_mask, static_cast<std::size_t>(pg.num_edges), "edge_mask");
        o.node_mask = node_bs.view; o.edge_mask = edge_bs.view;
        py::gil_scoped_release rel; auto res = algs.max_flow(pg.handle, src, dst, o); py::gil_scoped_acquire acq; return res;
      }, py::arg("graph"), py::arg("src"), py::arg("dst"), py::kw_only(), py::arg("flow_placement") = FlowPlacement::Proportional, py::arg("shortest_path") = false, py::arg("require_capacity") = true, py::arg("with_edge_flows") = false, py::arg("with_reachable") = false, py::arg("with_residuals") = false, py::arg("node_mask") = py::none(), py::arg("edge_mask") = py::none())
      .def("batch_max_flow", [](const Algorithms& algs, const PyGraph& pg, py::array pairs,
                                 py::object node_masks, py::object edge_masks,
                                 FlowPlacement placement, bool shortest_path, bool require_capacity,
                                 bool with_edge_flows, bool with_reachable, bool with_residuals){
        auto buf = pairs.request();
        if (buf.ndim != 2 || buf.shape[1] != 2) throw py::type_error("pairs must be shape [B,2]");
        if (buf.format != py::format_descriptor<std::int32_t>::format()) throw py::type_error("pairs dtype must be int32");
        const std::size_t B = static_cast<std::size_t>(buf.shape[0]);
        std::vector<std::pair<NodeId,NodeId>> pp; pp.reserve(B);
        auto* p = static_cast<const std::int32_t*>(buf.ptr);
        for (std::size_t i=0;i<B;++i) { pp.emplace_back(p[2*i], p[2*i+1]); }
        std::vector<BoolSpan> node_bufs, edge_bufs;
        std::vector<std::span<const bool>> node_spans; std::vector<std::span<const bool>> edge_spans;
        auto parse_mask_list = [&](py::object list_obj, std::size_t expected_len, const char* what,
                                   std::vector<BoolSpan>& bufs, std::vector<std::span<const bool>>& spans) {
          if (list_obj.is_none()) return;
          auto seq = py::cast<py::sequence>(list_obj);
          if (static_cast<std::size_t>(py::len(seq)) != B)
            throw py::type_error(std::string(what) + " length must equal number of pairs");
          bufs.reserve(B); spans.reserve(B);
          for (auto item : seq) {
            auto bs = to_bool_span_from_numpy(item, expected_len, what);
            spans.push_back(bs.view);
            bufs.push_back(std::move(bs)); // keep storage alive across GIL release
          }
        };
        parse_mask_list(node_masks, static_cast<std::size_t>(pg.num_nodes), "node_masks", node_bufs, node_spans);
        parse_mask_list(edge_masks, static_cast<std::size_t>(pg.num_edges), "edge_masks", edge_bufs, edge_spans);
        MaxFlowOptions o; o.placement = placement; o.shortest_path = shortest_path; o.require_capacity = require_capacity; o.with_edge_flows = with_edge_flows; o.with_reachable = with_reachable; o.with_residuals = with_residuals;
        py::gil_scoped_release rel; auto out = algs.batch_max_flow(pg.handle, pp, o, node_spans, edge_spans); py::gil_scoped_acquire acq; return out;
      }, py::arg("graph"), py::arg("pairs"), py::kw_only(), py::arg("node_masks") = py::none(), py::arg("edge_masks") = py::none(), py::arg("flow_placement") = FlowPlacement::Proportional, py::arg("shortest_path") = false, py::arg("require_capacity") = true, py::arg("with_edge_flows") = false, py::arg("with_reachable") = false, py::arg("with_residuals") = false)
      .def("sensitivity_analysis", [](const Algorithms& algs, const PyGraph& pg, std::int32_t src, std::int32_t dst,
                                      FlowPlacement placement, bool shortest_path, bool require_capacity,
                                      py::object node_mask, py::object edge_mask){
        if (src < 0 || src >= pg.num_nodes || dst < 0 || dst >= pg.num_nodes) throw py::value_error("src/dst out of range");
        MaxFlowOptions o; o.placement = placement; o.shortest_path = shortest_path; o.require_capacity = require_capacity;
        auto node_bs = to_bool_span_from_numpy(node_mask, static_cast<std::size_t>(pg.num_nodes), "node_mask");
        auto edge_bs = to_bool_span_from_numpy(edge_mask, static_cast<std::size_t>(pg.num_edges), "edge_mask");
        o.node_mask = node_bs.view; o.edge_mask = edge_bs.view;
        py::gil_scoped_release rel; auto res = algs.sensitivity_analysis(pg.handle, src, dst, o); py::gil_scoped_acquire acq;
        // Return as list of tuples (EdgeId, FlowDelta)
        py::list out;
        for (auto const& pr : res) {
          out.append(py::make_tuple(pr.first, pr.second));
        }
        return out;
      }, py::arg("graph"), py::arg("src"), py::arg("dst"), py::kw_only(), py::arg("flow_placement") = FlowPlacement::Proportional, py::arg("shortest_path") = false, py::arg("require_capacity") = true, py::arg("node_mask") = py::none(), py::arg("edge_mask") = py::none());

  py::class_<PredDAG>(m, "PredDAG")
      .def_property_readonly("parent_offsets", [](const PredDAG& d){
        return copy_to_numpy<std::int32_t>(std::span<const std::int32_t>(d.parent_offsets.data(), d.parent_offsets.size()));
      })
      .def_property_readonly("parents", [](const PredDAG& d){
        return copy_to_numpy<std::int32_t>(std::span<const std::int32_t>(d.parents.data(), d.parents.size()));
      })
      .def_property_readonly("via_edges", [](const PredDAG& d){
        return copy_to_numpy<std::int32_t>(std::span<const std::int32_t>(d.via_edges.data(), d.via_edges.size()));
      })
      .def("resolve_to_paths", [](const PredDAG& dag, std::int32_t src, std::int32_t dst, bool split_parallel_edges, py::object max_paths){
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
      }, py::arg("src"), py::arg("dst"), py::kw_only(), py::arg("split_parallel_edges") = false, py::arg("max_paths") = py::none());

  py::class_<MinCut>(m, "MinCut")
      .def_property_readonly("edges", [](const MinCut& mc){
        return copy_to_numpy<std::int32_t>(std::span<const std::int32_t>(mc.edges.data(), mc.edges.size()));
      });

  py::class_<FlowSummary>(m, "FlowSummary")
      .def_readonly("total_flow", &FlowSummary::total_flow)
      .def_readonly("min_cut", &FlowSummary::min_cut)
      .def_property_readonly("costs", [](const FlowSummary& s){
        return copy_to_numpy<std::int64_t>(std::span<const std::int64_t>(s.costs.data(), s.costs.size()));
      })
      .def_property_readonly("flows", [](const FlowSummary& s){
        return copy_to_numpy<double>(std::span<const double>(s.flows.data(), s.flows.size()));
      })
      .def_property_readonly("edge_flows", [](const FlowSummary& s){
        return copy_to_numpy<double>(std::span<const double>(s.edge_flows.data(), s.edge_flows.size()));
      })
      .def_property_readonly("residual_capacity", [](const FlowSummary& s){
        return copy_to_numpy<double>(std::span<const double>(s.residual_capacity.data(), s.residual_capacity.size()));
      })
      .def_property_readonly("reachable_nodes", [](const FlowSummary& s){
        // Store as bool array for Python; cast bytes to bool
        py::array_t<bool> arr(s.reachable_nodes.size());
        auto* out = arr.mutable_data();
        for (std::size_t i=0;i<s.reachable_nodes.size();++i) out[i] = static_cast<bool>(s.reachable_nodes[i]);
        return arr;
      });

  // FlowState bindings
  py::class_<FlowState>(m, "FlowState", py::dynamic_attr())
      .def("__init__", [](py::object self, py::object graph_obj){
        const StrictMultiDiGraph& g = py::cast<const StrictMultiDiGraph&>(graph_obj);
        FlowState* fs = self.cast<FlowState*>();
        new (fs) FlowState(g);
        self.attr("_graph_ref") = graph_obj;
      }, py::arg("graph"))
      .def("__init__", [](py::object self, py::object graph_obj, py::array residual){
        const StrictMultiDiGraph& g = py::cast<const StrictMultiDiGraph&>(graph_obj);
        if (!(py::isinstance<py::array_t<double>>(residual))) throw py::type_error("residual must be a numpy float64 array");
        if (!(residual.flags() & py::array::c_style)) throw py::type_error("residual must be C-contiguous");
        auto buf = residual.request();
        if (buf.ndim != 1 || static_cast<std::size_t>(buf.shape[0]) != static_cast<std::size_t>(g.num_edges())) throw py::type_error("residual length must equal num_edges");
        std::span<const double> rspan(static_cast<const double*>(buf.ptr), static_cast<std::size_t>(buf.shape[0]));
        FlowState* fs = self.cast<FlowState*>();
        new (fs) FlowState(g, rspan);
        self.attr("_graph_ref") = graph_obj;
      }, py::arg("graph"), py::arg("residual"))
      .def("reset", [](FlowState& fs){ fs.reset(); })
      .def("reset", [](FlowState& fs, py::array residual){
        if (!(py::isinstance<py::array_t<double>>(residual))) throw py::type_error("residual must be a numpy float64 array");
        if (!(residual.flags() & py::array::c_style)) throw py::type_error("residual must be C-contiguous");
        auto buf = residual.request();
        if (buf.ndim != 1 || static_cast<std::size_t>(buf.shape[0]) != static_cast<std::size_t>(fs.capacity_view().size())) throw py::type_error("residual length must equal num_edges");
        std::span<const double> rspan(static_cast<const double*>(buf.ptr), static_cast<std::size_t>(buf.shape[0]));
        fs.reset(rspan);
      })
      .def("capacity_view", [](py::object self_obj){
        const FlowState& fs = py::cast<const FlowState&>(self_obj);
        auto s = fs.capacity_view();
        py::array view(py::buffer_info(const_cast<double*>(s.data()), sizeof(double), py::format_descriptor<double>::format(), 1, { s.size() }, { sizeof(double) }), self_obj);
        view.attr("setflags")(py::arg("write") = false);
        return view;
      })
      .def("residual_view", [](py::object self_obj){
        const FlowState& fs = py::cast<const FlowState&>(self_obj);
        auto s = fs.residual_view();
        py::array view(py::buffer_info(const_cast<double*>(s.data()), sizeof(double), py::format_descriptor<double>::format(), 1, { s.size() }, { sizeof(double) }), self_obj);
        view.attr("setflags")(py::arg("write") = false);
        return view;
      })
      .def("edge_flow_view", [](py::object self_obj){
        const FlowState& fs = py::cast<const FlowState&>(self_obj);
        auto s = fs.edge_flow_view();
        py::array view(py::buffer_info(const_cast<double*>(s.data()), sizeof(double), py::format_descriptor<double>::format(), 1, { s.size() }, { sizeof(double) }), self_obj);
        view.attr("setflags")(py::arg("write") = false);
        return view;
      })
      .def("place_on_dag", [](FlowState& fs, std::int32_t src, std::int32_t dst, const PredDAG& dag, double requested_flow, FlowPlacement placement){
        py::gil_scoped_release rel; auto placed = fs.place_on_dag(src, dst, dag, requested_flow, placement); py::gil_scoped_acquire acq; return placed;
      }, py::arg("src"), py::arg("dst"), py::arg("dag"), py::arg("requested_flow") = std::numeric_limits<double>::infinity(), py::arg("flow_placement") = FlowPlacement::Proportional)
      .def("place_max_flow", [](py::object self_obj, std::int32_t src, std::int32_t dst,
                                 FlowPlacement placement, bool shortest_path, bool require_capacity,
                                 py::object node_mask, py::object edge_mask){
        FlowState& fs = py::cast<FlowState&>(self_obj);
        // Get graph reference to validate mask lengths
        const StrictMultiDiGraph& g = py::cast<const StrictMultiDiGraph&>(self_obj.attr("_graph_ref"));
        auto node_bs = to_bool_span_from_numpy(node_mask, static_cast<std::size_t>(g.num_nodes()), "node_mask");
        auto edge_bs = to_bool_span_from_numpy(edge_mask, static_cast<std::size_t>(g.num_edges()), "edge_mask");
        py::gil_scoped_release rel;
        auto total = fs.place_max_flow(src, dst, placement, shortest_path, require_capacity,
                                       node_bs.view, edge_bs.view);
        py::gil_scoped_acquire acq;
        return total;
      }, py::arg("src"), py::arg("dst"),
         py::arg("flow_placement") = FlowPlacement::Proportional,
         py::arg("shortest_path") = false, py::arg("require_capacity") = true,
         py::kw_only(), py::arg("node_mask") = py::none(), py::arg("edge_mask") = py::none())
      .def("compute_min_cut", [](py::object self_obj, std::int32_t src, py::object node_mask, py::object edge_mask){
        const FlowState& fs = py::cast<const FlowState&>(self_obj);
        // Get graph reference to validate mask lengths
        const StrictMultiDiGraph& g = py::cast<const StrictMultiDiGraph&>(self_obj.attr("_graph_ref"));
        auto node_bs = to_bool_span_from_numpy(node_mask, static_cast<std::size_t>(g.num_nodes()), "node_mask");
        auto edge_bs = to_bool_span_from_numpy(edge_mask, static_cast<std::size_t>(g.num_edges()), "edge_mask");
        py::gil_scoped_release rel;
        auto mc = fs.compute_min_cut(src, node_bs.view, edge_bs.view);
        py::gil_scoped_acquire acq;
        return mc;
      }, py::arg("src"), py::kw_only(), py::arg("node_mask") = py::none(), py::arg("edge_mask") = py::none());

  // FlowIndex and FlowGraph bindings
  py::class_<FlowIndex>(m, "FlowIndex")
      .def(py::init<NodeId, NodeId, std::int32_t, std::int64_t>())
      .def_readonly("src", &FlowIndex::src)
      .def_readonly("dst", &FlowIndex::dst)
      .def_readonly("flowClass", &FlowIndex::flowClass)
      .def_readonly("flowId", &FlowIndex::flowId);

  py::class_<FlowGraph>(m, "FlowGraph", py::dynamic_attr())
      .def("__init__", [](py::object self, py::object graph_obj){
        const StrictMultiDiGraph& g = py::cast<const StrictMultiDiGraph&>(graph_obj);
        FlowGraph* fg = self.cast<FlowGraph*>();
        new (fg) FlowGraph(g);
        self.attr("_graph_ref") = graph_obj;
      }, py::arg("graph"))
      .def("capacity_view", [](py::object self_obj){
        const FlowGraph& fg = py::cast<const FlowGraph&>(self_obj);
        auto s = fg.capacity_view();
        py::array view(py::buffer_info(const_cast<double*>(s.data()), sizeof(double), py::format_descriptor<double>::format(), 1, { s.size() }, { sizeof(double) }), self_obj);
        view.attr("setflags")(py::arg("write") = false);
        return view;
      })
      .def("residual_view", [](py::object self_obj){
        const FlowGraph& fg = py::cast<const FlowGraph&>(self_obj);
        auto s = fg.residual_view();
        py::array view(py::buffer_info(const_cast<double*>(s.data()), sizeof(double), py::format_descriptor<double>::format(), 1, { s.size() }, { sizeof(double) }), self_obj);
        view.attr("setflags")(py::arg("write") = false);
        return view;
      })
      .def("edge_flow_view", [](py::object self_obj){
        const FlowGraph& fg = py::cast<const FlowGraph&>(self_obj);
        auto s = fg.edge_flow_view();
        py::array view(py::buffer_info(const_cast<double*>(s.data()), sizeof(double), py::format_descriptor<double>::format(), 1, { s.size() }, { sizeof(double) }), self_obj);
        view.attr("setflags")(py::arg("write") = false);
        return view;
      })
      .def_property_readonly("graph", [](const FlowGraph& fg){ return &fg.graph(); }, py::return_value_policy::reference_internal)
      .def("place", [](FlowGraph& fg, const FlowIndex& idx, std::int32_t src, std::int32_t dst, const PredDAG& dag, double amount, FlowPlacement placement){ py::gil_scoped_release rel; auto placed = fg.place(idx, src, dst, dag, amount, placement); py::gil_scoped_acquire acq; return placed; }, py::arg("index"), py::arg("src"), py::arg("dst"), py::arg("dag"), py::arg("amount"), py::arg("flow_placement") = FlowPlacement::Proportional)
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

  py::class_<FlowPolicyConfig>(m, "FlowPolicyConfig")
      .def(py::init<>())
      .def(py::init([](PathAlg path_alg,
                       FlowPlacement flow_placement,
                       EdgeSelection selection,
                       bool require_capacity,
                       bool multipath,
                       int min_flow_count,
                       py::object max_flow_count,
                       py::object max_path_cost,
                       py::object max_path_cost_factor,
                       bool shortest_path,
                       bool reoptimize_flows_on_each_placement,
                       int max_no_progress_iterations,
                       int max_total_iterations,
                       bool diminishing_returns_enabled,
                       int diminishing_returns_window,
                       double diminishing_returns_epsilon_frac){
            FlowPolicyConfig cfg;
            cfg.path_alg = path_alg;
            cfg.flow_placement = flow_placement;
            cfg.selection = selection;
            cfg.require_capacity = require_capacity;
            cfg.multipath = multipath;
            cfg.min_flow_count = min_flow_count;
            if (!max_flow_count.is_none()) cfg.max_flow_count = py::cast<int>(max_flow_count);
            if (!max_path_cost.is_none()) cfg.max_path_cost = py::cast<Cost>(max_path_cost);
            if (!max_path_cost_factor.is_none()) cfg.max_path_cost_factor = py::cast<double>(max_path_cost_factor);
            cfg.shortest_path = shortest_path;
            cfg.reoptimize_flows_on_each_placement = reoptimize_flows_on_each_placement;
            cfg.max_no_progress_iterations = max_no_progress_iterations;
            cfg.max_total_iterations = max_total_iterations;
            cfg.diminishing_returns_enabled = diminishing_returns_enabled;
            cfg.diminishing_returns_window = diminishing_returns_window;
            cfg.diminishing_returns_epsilon_frac = diminishing_returns_epsilon_frac;
            return cfg;
          }),
          py::kw_only(),
          py::arg("path_alg") = PathAlg::SPF,
          py::arg("flow_placement") = FlowPlacement::Proportional,
          py::arg("selection") = EdgeSelection{},
          py::arg("require_capacity") = true,
          py::arg("multipath") = true,
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
      .def_readwrite("path_alg", &FlowPolicyConfig::path_alg)
      .def_readwrite("flow_placement", &FlowPolicyConfig::flow_placement)
      .def_readwrite("selection", &FlowPolicyConfig::selection)
      .def_readwrite("require_capacity", &FlowPolicyConfig::require_capacity)
      .def_readwrite("multipath", &FlowPolicyConfig::multipath)
      .def_readwrite("min_flow_count", &FlowPolicyConfig::min_flow_count)
      .def_readwrite("max_flow_count", &FlowPolicyConfig::max_flow_count)
      .def_readwrite("max_path_cost", &FlowPolicyConfig::max_path_cost)
      .def_readwrite("max_path_cost_factor", &FlowPolicyConfig::max_path_cost_factor)
      .def_readwrite("shortest_path", &FlowPolicyConfig::shortest_path)
      .def_readwrite("reoptimize_flows_on_each_placement", &FlowPolicyConfig::reoptimize_flows_on_each_placement)
      .def_readwrite("max_no_progress_iterations", &FlowPolicyConfig::max_no_progress_iterations)
      .def_readwrite("max_total_iterations", &FlowPolicyConfig::max_total_iterations)
      .def_readwrite("diminishing_returns_enabled", &FlowPolicyConfig::diminishing_returns_enabled)
      .def_readwrite("diminishing_returns_window", &FlowPolicyConfig::diminishing_returns_window)
      .def_readwrite("diminishing_returns_epsilon_frac", &FlowPolicyConfig::diminishing_returns_epsilon_frac);

  py::class_<FlowPolicy>(m, "FlowPolicy", py::dynamic_attr())
      .def("__init__", [](py::object self, py::object algs_obj, py::object graph_obj, FlowPolicyConfig cfg,
                          py::object node_mask, py::object edge_mask){
            std::shared_ptr<Algorithms> algs = py::cast<std::shared_ptr<Algorithms>>(algs_obj);
            const PyGraph& pg = py::cast<const PyGraph&>(graph_obj);

            // Convert masks to spans (FlowPolicy will copy the data)
            auto node_bs = to_bool_span_from_numpy(node_mask, static_cast<std::size_t>(pg.num_nodes), "node_mask");
            auto edge_bs = to_bool_span_from_numpy(edge_mask, static_cast<std::size_t>(pg.num_edges), "edge_mask");

            // Update config with mask spans (will be copied by FlowPolicy constructor)
            cfg.node_mask = node_bs.view;
            cfg.edge_mask = edge_bs.view;

            ExecutionContext ctx(algs, pg.handle);
            FlowPolicy* fp = self.cast<FlowPolicy*>();
            new (fp) FlowPolicy(ctx, cfg);
            self.attr("_algorithms_ref") = algs_obj;
            self.attr("_graph_ref") = graph_obj;
           },
           py::arg("algorithms"), py::arg("graph"), py::arg("config"),
           py::kw_only(),
           py::arg("node_mask") = py::none(),
           py::arg("edge_mask") = py::none(),
           py::keep_alive<1, 2>()  // self keeps algorithms alive
      )
      .def("flow_count", &FlowPolicy::flow_count)
      .def("placed_demand", &FlowPolicy::placed_demand)
      .def("place_demand", [](FlowPolicy& p, FlowGraph& fg, std::int32_t src, std::int32_t dst, FlowClass flowClass, double volume, py::object target_per_flow, py::object min_flow){ std::optional<double> tpf; if (!target_per_flow.is_none()) tpf = py::cast<double>(target_per_flow); std::optional<double> mfl; if (!min_flow.is_none()) mfl = py::cast<double>(min_flow); py::gil_scoped_release rel; auto pr = p.place_demand(fg, src, dst, flowClass, volume, tpf, mfl); py::gil_scoped_acquire acq; return py::make_tuple(pr.first, pr.second); }, py::arg("flow_graph"), py::arg("src"), py::arg("dst"), py::arg("flowClass"), py::arg("volume"), py::arg("target_per_flow") = py::none(), py::arg("min_flow") = py::none())
      .def("rebalance_demand", [](FlowPolicy& p, FlowGraph& fg, std::int32_t src, std::int32_t dst, FlowClass flowClass, double target){ py::gil_scoped_release rel; auto pr = p.rebalance_demand(fg, src, dst, flowClass, target); py::gil_scoped_acquire acq; return py::make_tuple(pr.first, pr.second); },
           py::arg("flow_graph"), py::arg("src"), py::arg("dst"), py::arg("flowClass"), py::arg("target"))
      .def("remove_demand", [](FlowPolicy& p, FlowGraph& fg){ py::gil_scoped_release rel; p.remove_demand(fg); py::gil_scoped_acquire acq; })
      .def_property_readonly("flows", [](const FlowPolicy& p){ py::dict out; for (auto const& kv : p.flows()) { const auto& idx = kv.first; const auto& f = kv.second; out[py::make_tuple(idx.src, idx.dst, idx.flowClass, idx.flowId)] = py::make_tuple(f.src, f.dst, f.cost, f.placed_flow); } return out; });

  // Profiling functions (enabled via NGRAPH_CORE_PROFILE=1 environment variable)
  m.def("profiling_enabled", &profiling_enabled,
        "Return True if profiling is enabled (NGRAPH_CORE_PROFILE=1).");
  m.def("profiling_dump", []{ ProfilingStats::instance().dump(); },
        "Print profiling statistics to stderr.");
  m.def("profiling_reset", []{ ProfilingStats::instance().reset(); },
        "Clear all profiling statistics.");
}
