/* The process-wide uid counter of graphs, Eligibility values and ledgers (internal). */
#pragma once

#include <cstdint>

namespace netgraph::core::detail {

// A fresh uid from one process-wide 64-bit counter: never 0, never reused.
[[nodiscard]] std::uint64_t next_uid() noexcept;

} // namespace netgraph::core::detail
