#include "netgraph/core/shortest_paths_reuse.hpp"

#include <algorithm>
#include <optional>
#include <stdexcept>
#include <string>

#include "netgraph/core/detail/ledger_access.hpp"
#include "netgraph/core/detail/policy_engine.hpp"
#include "netgraph/core/detail/program_data.hpp"
#include "netgraph/core/detail/route_builders.hpp"

namespace netgraph::core {
namespace detail {
namespace {
bool eligible(const FlowPolicyConfig& config) {
  const auto* sp = std::get_if<ShortestPaths>(&config.routing);
  return sp && config.admission == Admission::Lossless && sp->allocation == RouteAllocation::Equal &&
      (!sp->source_allocation || sp->source_allocation == SourceAllocation::Equal) && !sp->weight_basis && !sp->edge_weights &&
      sp->tie_break == RouteTieBreak::FewerHops && sp->multipath && sp->parallel_edges;
}
bool node(const Eligibility* mask, NodeId v) { return mask == nullptr || mask->allows_node(v); }
bool edge(const Eligibility* mask, EdgeId e) { return mask == nullptr || mask->allows_edge(e); }
std::uint64_t uid(const Eligibility* mask) { return mask ? mask->uid() : 0; }
bool same_mask(const Graph& graph, const Eligibility* a, const Eligibility* b) {
  if (uid(a) == uid(b)) return true;
  for (NodeId v = 0; v < graph.num_nodes(); ++v) if (node(a, v) != node(b, v)) return false;
  for (EdgeId e = 0; e < graph.num_edges(); ++e) if (edge(a, e) != edge(b, e)) return false;
  return true;
}
}

struct ReusePreparation {
  std::shared_ptr<const Graph> graph;
  std::optional<Eligibility> mask;
  std::vector<NodeId> sources, destinations, roster;
  std::shared_ptr<const ForwardingProgram> program;
};

struct ReuseDemand {
  std::vector<NodeId> sources, destinations, roster;
  std::shared_ptr<const ForwardingProgram> program;
  std::shared_ptr<const PredecessorDAG> members;
  std::vector<std::uint8_t> member, target;
  std::vector<std::int32_t> degree, offsets, incoming;
};

struct ReuseBaseline {
  std::shared_ptr<const Graph> graph;
  std::optional<Eligibility> mask;
  std::vector<std::optional<ReuseDemand>> demands;
};

std::uint64_t ShortestPathsReuseAccess::graph_uid(const ShortestPathsPreparation& preparation) noexcept {
  return preparation.data_->graph->uid();
}

void ShortestPathsReuseAccess::check(const ShortestPathsPreparation& preparation, const PolicyEngine& engine) {
  const auto& p = *preparation.data_;
  if (!same_mask(*p.graph, p.mask ? &*p.mask : nullptr, engine.eligibility)) {
    throw std::invalid_argument("preparation Eligibility does not match policy Eligibility");
  }
  if (p.sources != engine.sources || p.destinations != engine.destinations) {
    throw std::invalid_argument("preparation endpoints do not match policy endpoints");
  }
  if (!eligible(engine.config)) throw std::invalid_argument("preparation configuration does not match policy configuration");
}

void ShortestPathsReuseAccess::install(const ShortestPathsPreparation& preparation, PolicyEngine& engine) {
  const auto& p = *preparation.data_;
  engine.reachable_sources = p.roster;
  engine.prepared_program = p.program;
  engine.preparation_consumed = true;
}
} // namespace detail

ShortestPathsReuse::ShortestPathsReuse(std::span<FlowPolicy* const> policies) {
  using namespace detail;
  if (policies.empty()) throw std::invalid_argument("baseline policies must not be empty");
  for (const auto* p : policies) {
    if (p == nullptr || !PolicyAccess::constructed(*p)) throw std::invalid_argument("FlowPolicy is not constructed");
  }
  auto& first = PolicyAccess::engine(*policies.front());
  for (const auto* p : policies) {
    if (&PolicyAccess::engine(*p).ledger != &first.ledger) throw std::invalid_argument("baseline policies must share one ledger");
  }
  for (const auto* p : policies) {
    if (!same_mask(*first.graph, first.eligibility, PolicyAccess::engine(*p).eligibility)) {
      throw std::invalid_argument("baseline policies must share one Eligibility");
    }
  }
  auto out = std::make_shared<ReuseBaseline>();
  out->graph = first.graph;
  if (first.eligibility) out->mask = *first.eligibility;
  const auto gate = LedgerAccess::lock(first.ledger);
  for (std::size_t i = 0; i < policies.size(); ++i) {
    const auto& pe = PolicyAccess::engine(*policies[i]);
    if (!eligible(pe.config)) { out->demands.emplace_back(); continue; }
    if (!pe.epoch.program) throw std::invalid_argument("baseline policy " + std::to_string(i) + " has no prepared program");
    ReuseDemand d;
    d.sources = pe.sources;
    d.destinations = pe.destinations;
    d.roster = pe.reachable_sources;
    d.program = pe.epoch.program;
    const auto& data = ProgramAccess::data(*d.program);
    const auto members = data.common_replay_support ? std::span<const EdgeId>(*data.common_replay_support) : std::span<const EdgeId> {};
    d.members = std::make_shared<const PredecessorDAG>(PredecessorDAGBuilder::of_edges(*out->graph, members));
    const auto n = static_cast<std::size_t>(out->graph->num_nodes());
    d.member.assign(static_cast<std::size_t>(out->graph->num_edges()), 0);
    d.target.assign(n, 0);
    d.degree.assign(n, 0);
    d.offsets.assign(n + 1, 0);
    for (NodeId t : pe.allowed_destinations) d.target[static_cast<std::size_t>(t)] = 1;
    const auto tails = out->graph->edge_sources(), heads = out->graph->edge_destinations();
    for (EdgeId e : members) {
      d.member[static_cast<std::size_t>(e)] = 1;
      ++d.degree[static_cast<std::size_t>(tails[static_cast<std::size_t>(e)])];
      ++d.offsets[static_cast<std::size_t>(heads[static_cast<std::size_t>(e)]) + 1];
    }
    for (std::size_t v = 1; v < d.offsets.size(); ++v) d.offsets[v] += d.offsets[v - 1];
    auto fill = d.offsets;
    d.incoming.resize(members.size());
    for (EdgeId e : members) d.incoming[static_cast<std::size_t>(fill[static_cast<std::size_t>(heads[static_cast<std::size_t>(e)])]++)] = e;
    out->demands.emplace_back(std::move(d));
  }
  data_ = std::move(out);
}

std::vector<std::shared_ptr<ShortestPathsPreparation>> ShortestPathsReuse::prepare(const Eligibility* mask) const {
  using namespace detail;
  const auto& g = *data_->graph;
  if (mask && mask->graph_uid() != g.uid()) throw std::invalid_argument("foreign Eligibility: built for another graph");
  std::vector<std::shared_ptr<ShortestPathsPreparation>> out(data_->demands.size());
  const auto* baseline = data_->mask ? &*data_->mask : nullptr;
  std::vector<NodeId> failed_nodes;
  std::vector<EdgeId> failed_edges;
  for (NodeId v = 0; v < g.num_nodes(); ++v) {
    if (node(mask, v) && !node(baseline, v)) return out;
    if (node(baseline, v) && !node(mask, v)) failed_nodes.push_back(v);
  }
  for (EdgeId e = 0; e < g.num_edges(); ++e) {
    if (edge(mask, e) && !edge(baseline, e)) return out;
    if (edge(baseline, e) && !edge(mask, e)) failed_edges.push_back(e);
  }
  const auto tails = g.edge_sources();
  for (std::size_t i = 0; i < data_->demands.size(); ++i) {
    if (!data_->demands[i]) continue;
    const auto& d = *data_->demands[i];
    auto degree = d.degree;
    std::vector<std::uint8_t> affected(degree.size(), 0), removed(d.member.size(), 0);
    std::vector<NodeId> queue;
    auto mark = [&](NodeId v) {
      if (!affected[static_cast<std::size_t>(v)]) { affected[static_cast<std::size_t>(v)] = 1; queue.push_back(v); }
    };
    auto kill = [&](EdgeId e) {
      const auto ei = static_cast<std::size_t>(e);
      if (removed[ei]) return;
      removed[ei] = 1;
      const auto u = static_cast<std::size_t>(tails[ei]);
      if (affected[u]) return;
      if (--degree[u] == 0 && !(d.target[u] && node(mask, static_cast<NodeId>(u)))) mark(static_cast<NodeId>(u));
    };
    for (NodeId v : failed_nodes) mark(v);
    for (EdgeId e : failed_edges) if (d.member[static_cast<std::size_t>(e)]) kill(e);
    for (std::size_t q = 0; q < queue.size(); ++q) {
      const auto v = static_cast<std::size_t>(queue[q]);
      for (auto k = d.offsets[v]; k < d.offsets[v + 1]; ++k) kill(d.incoming[static_cast<std::size_t>(k)]);
    }
    bool accepted = true;
    std::vector<NodeId> roots;
    for (NodeId r : d.roster) if (node(mask, r)) {
      if (affected[static_cast<std::size_t>(r)]) accepted = false;
      roots.push_back(r);
    }
    if (!accepted) continue;
    auto hint = std::make_shared<ReusePreparation>();
    hint->graph = data_->graph;
    if (mask) hint->mask = *mask;
    hint->sources = d.sources;
    hint->destinations = d.destinations;
    hint->roster = roots;
    if (roots.empty()) {
      hint->program = empty_program(g.uid(), d.sources, d.destinations, false);
    } else {
      auto plan = restrict_shortest_path_route(g, *d.members, roots, d.destinations, mask);
      const auto& old = ProgramAccess::data(*d.program);
      std::vector<ProgramEntry> entries(roots.size());
      for (std::size_t r = 0; r < roots.size(); ++r) {
        const auto k = static_cast<std::size_t>(std::lower_bound(d.roster.begin(), d.roster.end(), roots[r]) - d.roster.begin());
        entries[r].min_route_cost = old.entries[k].min_route_cost;
      }
      hint->program = shortest_path_program(g, d.sources, d.destinations, roots,
          std::make_shared<const PreparedRoute>(std::move(plan)), std::move(entries));
    }
    out[i] = std::shared_ptr<ShortestPathsPreparation>(new ShortestPathsPreparation(std::move(hint)));
  }
  return out;
}
} // namespace netgraph::core
