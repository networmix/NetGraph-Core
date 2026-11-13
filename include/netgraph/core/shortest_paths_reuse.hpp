/* Immutable baseline route preparation for deletion-only ShortestPaths cases. */
#pragma once

#include <cstdint>
#include <memory>
#include <span>
#include <vector>

#include "netgraph/core/flow_policy.hpp"

namespace netgraph::core {
namespace detail {
struct ReuseBaseline;
struct ReusePreparation;
struct ShortestPathsReuseAccess;
}

// One case's topology-derived preparation. It owns no ledger or policy, and never
// contains placed amounts. Its immutable route arrays may be shared by policies.
class ShortestPathsPreparation {
private:
  friend class ShortestPathsReuse;
  friend struct detail::ShortestPathsReuseAccess;
  explicit ShortestPathsPreparation(std::shared_ptr<const detail::ReusePreparation> data) : data_(std::move(data)) {}
  std::shared_ptr<const detail::ReusePreparation> data_;
};

// Capture prepared default Equal/Lossless ShortestPaths programs of one ledger.
// Unsupported configurations keep their list position but prepare a null hint.
// All input policies must have the same effective mask; eligible ones must have
// installed a program. Capturing retains the graph and route values, not policies.
class ShortestPathsReuse {
public:
  explicit ShortestPathsReuse(std::span<FlowPolicy* const> policies);
  // One dense mask comparison shared by all demands; no hint for a noncontained
  // mask or for a demand with a surviving source whose distance changes.
  [[nodiscard]] std::vector<std::shared_ptr<ShortestPathsPreparation>> prepare(const Eligibility* eligibility = nullptr) const;

private:
  std::shared_ptr<const detail::ReuseBaseline> data_;
};

namespace detail {
struct ShortestPathsReuseAccess {
  [[nodiscard]] static std::uint64_t graph_uid(const ShortestPathsPreparation& preparation) noexcept;
  static void check(const ShortestPathsPreparation& preparation, const PolicyEngine& engine);
  // Install an invisible template before binding. Every precommit computation
  // creates a fresh public program; only successful publication retains one.
  static void install(const ShortestPathsPreparation& preparation, PolicyEngine& engine);
};
}
} // namespace netgraph::core
