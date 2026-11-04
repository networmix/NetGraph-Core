// One policy epoch's identities, load and stored record data. The ledger supplies its uid
// when materializing a public snapshot, so copies never retain an old ledger identity.
#pragma once

#include <cstdint>
#include <memory>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "netgraph/core/flow_record.hpp"
#include "netgraph/core/graph.hpp"
#include "netgraph/core/forwarding_program.hpp"

namespace netgraph::core::detail {

struct PreparedRoute;
struct ForwardingView;
struct ProgramFlow;

struct RecordInventory {
  NodeId source {0};
  std::optional<std::int32_t> entry_id {};
};

struct DerivedInputs {
  std::shared_ptr<const ForwardingProgram> program;
  std::shared_ptr<const PreparedRoute> plan;
  std::vector<std::uint64_t> roster, basis;
  std::vector<NodeId> destinations;
  std::vector<RecordInventory> inventory;
};
struct ForwardingHistory {
  std::shared_ptr<const ForwardingHistory> previous;
  std::optional<Amount> offer;
  Counter128 offset;
  std::vector<std::uint64_t> dark_heads;
  bool raw_head_offers {false};
  std::vector<Amount> carried, drops;
  ~ForwardingHistory() {
    // Retire an unshared suffix iteratively; a long epoch must not consume the call stack.
    while (previous && previous.use_count() == 1) {
      auto next = std::move(const_cast<ForwardingHistory&>(*previous).previous);
      previous = std::move(next);
    }
  }
};
struct DerivedRow {
  Counter128 offered, source_drops;
  Amount delivered {0};
};
struct DerivedRecords {
  std::shared_ptr<const DerivedInputs> inputs;
  std::shared_ptr<const ForwardingHistory> history;
  std::vector<DerivedRow> rows;
};

struct StoredRecord {
  std::shared_ptr<const PredecessorDAG> route {};
  Cost min_route_cost {kUnreachable};
  std::optional<Counter128> offered {};
  Amount delivered {0};
  std::vector<std::pair<EdgeId, Amount>> edge_flows {};
  std::vector<std::pair<EdgeId, Counter128>> edge_drops {};
  std::vector<std::pair<EdgeId, Amount>> edge_reservations {};
};

// One policy epoch. Use the accessors for identity and aggregate data: their storage
// depends on the record kind and count.
struct LedgerEntry {
  std::int64_t owner {0};
  Admission admission {Admission::Lossless};
  // Stored-record identities; derived entries keep these in derived->inputs instead.
  std::vector<NodeId> destinations {};
  std::int64_t first_counter {0};
  std::vector<RecordInventory> inventory {};
  // Carried, held and dropped vectors stay empty when one stored record owns
  // these aggregates in its edge_* fields.
  std::vector<std::pair<EdgeId, Amount>> carried {};
  std::vector<std::pair<EdgeId, Amount>> holds {};
  Amount delivered {0};
  std::vector<std::pair<EdgeId, Counter128>> drops {};
  std::vector<StoredRecord> stored {};
  // Source records are materialized from derived inputs. With one root, the
  // aggregate is its only record's traffic, so no forwarding history is needed.
  std::optional<DerivedRecords> derived {};

  [[nodiscard]] const auto& aggregate_carried() const noexcept {
    return stored.size() == 1 ? stored.front().edge_flows : carried;
  }
  [[nodiscard]] const auto& aggregate_holds() const noexcept {
    return stored.size() == 1 ? stored.front().edge_reservations : holds;
  }
  [[nodiscard]] const auto& aggregate_drops() const noexcept {
    return stored.size() == 1 ? stored.front().edge_drops : drops;
  }
  [[nodiscard]] const std::vector<NodeId>& destination_nodes() const noexcept {
    return derived ? derived->inputs->destinations : destinations;
  }
  [[nodiscard]] const std::vector<RecordInventory>& record_inventory() const noexcept {
    return derived ? derived->inputs->inventory : inventory;
  }

  [[nodiscard]] std::vector<FlowRecord> records(const Graph&, std::uint64_t uid, std::size_t first, std::size_t count) const;
  [[nodiscard]] FlowRecord materialize(const Graph& graph, std::uint64_t uid, std::size_t position) const {
    return std::move(records(graph, uid, position, 1).front());
  }
};

inline void accumulate_counter(Counter128& value, const Counter128& increment, const char* label) {
  const std::uint64_t lo = value.lo + increment.lo;
  const std::uint64_t hi = value.hi + increment.hi + (lo < value.lo ? 1U : 0U);
  if (hi >= (std::uint64_t {1} << 63)) throw std::invalid_argument(std::string(label) + " limit exceeded");
  value = Counter128 {hi, lo};
}

// Consumes a complete set of records, keeping their identity fields only in the common entry.
// Every allocation precedes publication; the caller owns the pending value until then.
[[nodiscard]] LedgerEntry stored_entry(const Graph& graph, std::int64_t owner, std::span<FlowRecord> records);
[[nodiscard]] LedgerEntry derived_entry(const Graph&, std::int64_t owner, Admission, std::span<const NodeId> destinations,
    const ForwardingView&, bool entry_records, const LedgerEntry* previous, std::optional<Amount> offer,
    Counter128 offset, ProgramFlow& flow);

inline bool bit_at(std::span<const std::uint64_t> bits, std::size_t i) { return (bits[i / 64] >> (i % 64)) & 1U; }
inline void set_bit(std::vector<std::uint64_t>& bits, std::size_t i) { bits[i / 64] |= std::uint64_t {1} << (i % 64); }

#if NETGRAPH_CORE_PRIVATE_HOOKS
struct EntryWork {
  std::size_t forwards {0}, attribution {0}, joint {0};
  std::size_t common_rosters {0}, derivations {0}, routes {0};
};
EntryWork& thread_entry_work() noexcept;
#endif

} // namespace netgraph::core::detail
