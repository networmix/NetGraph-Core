#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <cstring>

#include "netgraph/core/k_shortest_paths.hpp"
#include "netgraph/core/max_flow.hpp"
#include "netgraph/core/shortest_paths.hpp"
#include "netgraph/core/strict_multidigraph.hpp"
#include "netgraph/core/types.hpp"

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

  py::enum_<EdgeSelect>(m, "EdgeSelect")
      .value("ALL_MIN_COST", EdgeSelect::AllMinCost)
      .value("SINGLE_MIN_COST", EdgeSelect::SingleMinCost)
      .value("ALL_MIN_COST_WITH_CAP_REMAINING", EdgeSelect::AllMinCostWithCapRemaining);

  py::enum_<FlowPlacement>(m, "FlowPlacement")
      .value("PROPORTIONAL", FlowPlacement::Proportional)
      .value("EQUAL_BALANCED", FlowPlacement::EqualBalanced);

  py::class_<StrictMultiDiGraph>(m, "StrictMultiDiGraph")
      .def_static(
          "from_arrays",
          [](std::int32_t num_nodes,
             py::array src, py::array dst,
             py::array capacity, py::array cost,
             py::object link_ids,
             bool add_reverse) {
            // public API: src/dst are int32; pass through as int32 to core
            auto src_s = as_span<std::int32_t>(src, "src");
            auto dst_s = as_span<std::int32_t>(dst, "dst");
            if (src_s.size() != dst_s.size()) throw py::type_error("src and dst must have the same length");
            auto cap_s = as_span<double>(capacity, "capacity");
            auto cost_s = as_span<double>(cost, "cost");
            std::span<const std::int64_t> link_s;
            std::vector<std::int64_t> link_tmp;
            if (!link_ids.is_none()) {
              auto arr = py::cast<py::array>(link_ids);
              // Accept any integer dtype; force-cast to int64 and require C-contiguous
              py::array_t<std::int64_t, py::array::c_style | py::array::forcecast> arr64(arr);
              auto buf = arr64.request();
              auto* ptr = static_cast<const std::int64_t*>(buf.ptr);
              link_tmp.assign(ptr, ptr + buf.size);
              link_s = std::span<const std::int64_t>(link_tmp.data(), link_tmp.size());
            }
            return StrictMultiDiGraph::from_arrays(num_nodes,
                                                  src_s,
                                                  dst_s,
                                                  cap_s, cost_s, link_s, add_reverse);
          },
          py::arg("num_nodes"), py::arg("src"), py::arg("dst"), py::arg("capacity"), py::arg("cost"), py::arg("link_ids") = py::none(),
          py::kw_only(), py::arg("add_reverse") = false)
      .def("num_nodes", &StrictMultiDiGraph::num_nodes)
      .def("num_edges", &StrictMultiDiGraph::num_edges)
      .def("link_id_of", &StrictMultiDiGraph::link_id_of)
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
      .def("cost_view", [](py::object self_obj, const StrictMultiDiGraph& g){
        auto s = g.cost_view();
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
           EdgeSelect edge_select, bool multipath, py::object node_mask, py::object edge_mask, double eps) {
          std::optional<NodeId> dst_opt;
          if (!dst.is_none()) dst_opt = static_cast<NodeId>(py::cast<std::int32_t>(dst));
          // Validate masks and obtain raw pointers (optional)
          const bool* node_ptr = nullptr;
          const bool* edge_ptr = nullptr;
          py::array node_arr;
          py::array edge_arr;
          if (!node_mask.is_none()) {
            node_arr = py::cast<py::array>(node_mask);
            auto buf = node_arr.request();
            if (buf.ndim != 1) throw py::type_error("node_mask must be a 1-D array of bool");
            if (buf.format != py::format_descriptor<bool>::format()) throw py::type_error("node_mask must be dtype=bool");
            if (static_cast<std::size_t>(buf.shape[0]) != static_cast<std::size_t>(g.num_nodes())) throw py::type_error("node_mask length must equal num_nodes");
            node_ptr = static_cast<const bool*>(buf.ptr);
          }
          if (!edge_mask.is_none()) {
            edge_arr = py::cast<py::array>(edge_mask);
            auto buf = edge_arr.request();
            if (buf.ndim != 1) throw py::type_error("edge_mask must be a 1-D array of bool");
            if (buf.format != py::format_descriptor<bool>::format()) throw py::type_error("edge_mask must be dtype=bool");
            if (static_cast<std::size_t>(buf.shape[0]) != static_cast<std::size_t>(g.num_edges())) throw py::type_error("edge_mask length must equal num_edges");
            edge_ptr = static_cast<const bool*>(buf.ptr);
          }
          py::gil_scoped_release release;
          auto res = shortest_paths(g, src, dst_opt, edge_select, multipath, eps, node_ptr, edge_ptr);
          py::gil_scoped_acquire acquire;
          py::array_t<double> dist_arr(res.first.size());
          std::memcpy(dist_arr.mutable_data(), res.first.data(), res.first.size() * sizeof(double));
          return py::make_tuple(std::move(dist_arr), res.second);
        }, py::arg("g"), py::arg("src"), py::arg("dst") = py::none(), py::kw_only(), py::arg("edge_select") = EdgeSelect::AllMinCost, py::arg("multipath") = true, py::arg("node_mask") = py::none(), py::arg("edge_mask") = py::none(), py::arg("eps") = 1e-10);

  py::class_<Path>(m, "Path")
      .def_property_readonly("nodes", [](const Path& p){
        py::array_t<std::int32_t> arr(p.nodes.size());
        std::memcpy(arr.mutable_data(), p.nodes.data(), p.nodes.size()*sizeof(std::int32_t));
        return arr;
      })
      .def_property_readonly("edges", [](const Path& p){
        py::array_t<std::int32_t> arr(p.edges.size());
        std::memcpy(arr.mutable_data(), p.edges.data(), p.edges.size()*sizeof(std::int32_t));
        return arr;
      })
      .def_readonly("cost", &Path::cost);

  m.def("ksp",
        [](const StrictMultiDiGraph& g, std::int32_t src, std::int32_t dst,
           int k, py::object max_cost_factor, bool unique, py::object node_mask, py::object edge_mask, double eps) {
          std::optional<double> mcf;
          if (!max_cost_factor.is_none()) mcf = py::cast<double>(max_cost_factor);
          if (!node_mask.is_none()) {
            auto arr = py::cast<py::array>(node_mask);
            auto buf = arr.request();
            if (buf.ndim != 1) throw py::type_error("node_mask must be a 1-D array of bool");
            if (buf.format != py::format_descriptor<bool>::format()) throw py::type_error("node_mask must be dtype=bool");
            if (static_cast<std::size_t>(buf.shape[0]) != static_cast<std::size_t>(g.num_nodes())) throw py::type_error("node_mask length must equal num_nodes");
          }
          if (!edge_mask.is_none()) {
            auto arr = py::cast<py::array>(edge_mask);
            auto buf = arr.request();
            if (buf.ndim != 1) throw py::type_error("edge_mask must be a 1-D array of bool");
            if (buf.format != py::format_descriptor<bool>::format()) throw py::type_error("edge_mask must be dtype=bool");
            if (static_cast<std::size_t>(buf.shape[0]) != static_cast<std::size_t>(g.num_edges())) throw py::type_error("edge_mask length must equal num_edges");
          }
          py::gil_scoped_release release;
          auto paths = k_shortest_paths(g, src, dst, k, mcf, unique, eps);
          py::gil_scoped_acquire acquire;
          return paths;
        }, py::arg("g"), py::arg("src"), py::arg("dst"), py::kw_only(), py::arg("k"), py::arg("max_cost_factor") = py::none(), py::arg("unique") = true, py::arg("node_mask") = py::none(), py::arg("edge_mask") = py::none(), py::arg("eps") = 1e-10);

  py::class_<MinCut>(m, "MinCut")
      .def_readonly("edges", &MinCut::edges);
  py::class_<CostBucket>(m, "CostBucket")
      .def_readonly("cost", &CostBucket::cost)
      .def_readonly("share", &CostBucket::share);
  py::class_<CostDistribution>(m, "CostDistribution")
      .def_readonly("buckets", &CostDistribution::buckets);
  py::class_<FlowSummary>(m, "FlowSummary")
      .def_readonly("total_flow", &FlowSummary::total_flow)
      .def_readonly("min_cut", &FlowSummary::min_cut)
      .def_readonly("cost_distribution", &FlowSummary::cost_distribution)
      .def_readonly("edge_flows", &FlowSummary::edge_flows);

  m.def("calc_max_flow",
        [](const StrictMultiDiGraph& g, std::int32_t src, std::int32_t dst,
           FlowPlacement placement, bool shortest_path, double eps, bool with_edge_flows,
           py::object node_mask, py::object edge_mask) {
          const bool* node_ptr = nullptr;
          const bool* edge_ptr = nullptr;
          py::array node_arr, edge_arr;
          if (!node_mask.is_none()) {
            node_arr = py::cast<py::array>(node_mask);
            auto buf = node_arr.request();
            if (buf.ndim != 1) throw py::type_error("node_mask must be a 1-D array of bool");
            if (buf.format != py::format_descriptor<bool>::format()) throw py::type_error("node_mask must be dtype=bool");
            if (static_cast<std::size_t>(buf.shape[0]) != static_cast<std::size_t>(g.num_nodes())) throw py::type_error("node_mask length must equal num_nodes");
            node_ptr = static_cast<const bool*>(buf.ptr);
          }
          if (!edge_mask.is_none()) {
            edge_arr = py::cast<py::array>(edge_mask);
            auto buf = edge_arr.request();
            if (buf.ndim != 1) throw py::type_error("edge_mask must be a 1-D array of bool");
            if (buf.format != py::format_descriptor<bool>::format()) throw py::type_error("edge_mask must be dtype=bool");
            if (static_cast<std::size_t>(buf.shape[0]) != static_cast<std::size_t>(g.num_edges())) throw py::type_error("edge_mask length must equal num_edges");
            edge_ptr = static_cast<const bool*>(buf.ptr);
          }
          py::gil_scoped_release release;
          auto res = calc_max_flow(g, src, dst, placement, shortest_path, eps, with_edge_flows, node_ptr, edge_ptr);
          py::gil_scoped_acquire acquire;
          return res;
        }, py::arg("g"), py::arg("src"), py::arg("dst"), py::kw_only(), py::arg("flow_placement") = FlowPlacement::Proportional, py::arg("shortest_path") = false, py::arg("eps") = 1e-10, py::arg("with_edge_flows") = false, py::arg("node_mask") = py::none(), py::arg("edge_mask") = py::none());

  m.def("max_flow",
        [](const StrictMultiDiGraph& g, std::int32_t src, std::int32_t dst,
           FlowPlacement placement, bool shortest_path, double eps, bool with_edge_flows,
           py::object node_mask, py::object edge_mask) {
          const bool* node_ptr = nullptr;
          const bool* edge_ptr = nullptr;
          py::array node_arr, edge_arr;
          if (!node_mask.is_none()) {
            node_arr = py::cast<py::array>(node_mask);
            auto buf = node_arr.request();
            if (buf.ndim != 1) throw py::type_error("node_mask must be a 1-D array of bool");
            if (buf.format != py::format_descriptor<bool>::format()) throw py::type_error("node_mask must be dtype=bool");
            if (static_cast<std::size_t>(buf.shape[0]) != static_cast<std::size_t>(g.num_nodes())) throw py::type_error("node_mask length must equal num_nodes");
            node_ptr = static_cast<const bool*>(buf.ptr);
          }
          if (!edge_mask.is_none()) {
            edge_arr = py::cast<py::array>(edge_mask);
            auto buf = edge_arr.request();
            if (buf.ndim != 1) throw py::type_error("edge_mask must be a 1-D array of bool");
            if (buf.format != py::format_descriptor<bool>::format()) throw py::type_error("edge_mask must be dtype=bool");
            if (static_cast<std::size_t>(buf.shape[0]) != static_cast<std::size_t>(g.num_edges())) throw py::type_error("edge_mask length must equal num_edges");
            edge_ptr = static_cast<const bool*>(buf.ptr);
          }
          py::gil_scoped_release release;
          auto res = calc_max_flow(g, src, dst, placement, shortest_path, eps, with_edge_flows, node_ptr, edge_ptr);
          py::gil_scoped_acquire acquire;
          return res;
        }, py::arg("g"), py::arg("src"), py::arg("dst"), py::kw_only(), py::arg("flow_placement") = FlowPlacement::Proportional, py::arg("shortest_path") = false, py::arg("eps") = 1e-10, py::arg("with_edge_flows") = false, py::arg("node_mask") = py::none(), py::arg("edge_mask") = py::none());

  m.def("batch_max_flow",
        [](const StrictMultiDiGraph& g, py::array pairs,
           py::object node_masks, py::object edge_masks,
           FlowPlacement placement, bool shortest_path, double eps, bool with_edge_flows, py::object threads, py::object seed) {
          auto buf = pairs.request();
          if (buf.ndim != 2 || buf.shape[1] != 2) throw py::type_error("pairs must be shape [B,2]");
          // accept int32 for pairs
          if (buf.format != py::format_descriptor<std::int32_t>::format()) throw py::type_error("pairs dtype must be int32");
          std::vector<std::pair<NodeId,NodeId>> pp;
          auto* ptr = static_cast<const std::int32_t*>(buf.ptr);
          pp.reserve(static_cast<std::size_t>(buf.shape[0]));
          for (ssize_t i=0;i<buf.shape[0];++i) {
            if (ptr[2*i] < 0 || ptr[2*i+1] < 0) throw py::type_error("pairs must be non-negative indices");
            pp.emplace_back(static_cast<std::int32_t>(ptr[2*i]), static_cast<std::int32_t>(ptr[2*i+1]));
          }
          // Parse mask lists if provided, collect raw pointers and keep arrays alive
          auto parse_mask_list = [&](py::object list_obj, std::size_t expected_len, const char* what, std::vector<const bool*>& out_ptrs, std::vector<py::array>& keep){
            if (list_obj.is_none()) return;
            auto seq = py::cast<py::sequence>(list_obj);
            for (auto item: seq) {
              auto arr = py::cast<py::array>(item);
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
          int threads_val = threads.is_none() ? 0 : py::cast<int>(threads);
          std::optional<std::uint64_t> seed_val;
          if (!seed.is_none()) seed_val = py::cast<std::uint64_t>(seed);
          py::gil_scoped_release release;
          auto out = netgraph::core::batch_max_flow(g, pp, placement, shortest_path, eps, with_edge_flows, threads_val, seed_val, node_ptrs, edge_ptrs);
          py::gil_scoped_acquire acquire;
          return out;
        }, py::arg("g"), py::arg("pairs"), py::kw_only(), py::arg("node_masks") = py::none(), py::arg("edge_masks") = py::none(), py::arg("flow_placement") = FlowPlacement::Proportional, py::arg("shortest_path") = false, py::arg("eps") = 1e-10, py::arg("with_edge_flows") = false, py::arg("threads") = py::none(), py::arg("seed") = py::none());
}
