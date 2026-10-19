/* The process-wide uid counter of graphs, Eligibility values and ledgers (internal). */
#pragma once

#include <atomic>
#include <cstdint>

namespace netgraph::core::detail {

// A fresh uid from one process-wide 64-bit counter: never 0, never reused, at most the largest
// int64 (2^63 - 1 identities); after it the counter is spent: std::runtime_error `object identity
// limit exceeded`, nothing changed.
[[nodiscard]] std::uint64_t next_uid();
// The rule of next_uid on a given counter (the counter holds the last uid taken, 0 for none).
[[nodiscard]] std::uint64_t take_uid(std::atomic<std::uint64_t>& counter);

} // namespace netgraph::core::detail
