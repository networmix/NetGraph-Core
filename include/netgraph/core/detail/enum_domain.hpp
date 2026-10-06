/* Membership of enumeration values (internal; rule 13). An enum-typed field can hold a value
   outside its enumeration (a C++ cast, a Python enum built from an integer); such a value is
   outside the domain and is rejected (std::invalid_argument naming the field) before it
   reaches any rule or record. */
#pragma once

#include <stdexcept>
#include <string>

#include "netgraph/core/options.hpp"

namespace netgraph::core::detail {

[[noreturn]] inline void not_a_member(const char* field, int value, const char* enumeration) {
  throw std::invalid_argument(std::string(field) + ": " + std::to_string(value) + " is not a value of " + enumeration);
}

inline void check_member(FlowPlacement v, const char* field) {
  switch (v) {
    case FlowPlacement::MaxFlow:
    case FlowPlacement::Weighted:
    case FlowPlacement::EqualBalanced: return;
  }
  not_a_member(field, static_cast<int>(v), "FlowPlacement");
}

inline void check_member(Admission v, const char* field) {
  switch (v) {
    case Admission::Lossless:
    case Admission::Lossy: return;
  }
  not_a_member(field, static_cast<int>(v), "Admission");
}

inline void check_member(EdgeTieBreak v, const char* field) {
  switch (v) {
    case EdgeTieBreak::Deterministic:
    case EdgeTieBreak::PreferHigherPriority: return;
  }
  not_a_member(field, static_cast<int>(v), "EdgeTieBreak");
}

inline void check_member(EdgePriority v, const char* field) {
  switch (v) {
    case EdgePriority::Residual:
    case EdgePriority::Capacity:
    case EdgePriority::Headroom: return;
  }
  not_a_member(field, static_cast<int>(v), "EdgePriority");
}

inline void check_member(Direction v, const char* field) {
  switch (v) {
    case Direction::Forward:
    case Direction::Reverse: return;
  }
  not_a_member(field, static_cast<int>(v), "Direction");
}

// The two enum-typed fields of a selection.
inline void check_selection(const PathSelection& selection) {
  check_member(selection.tie_break, "selection.tie_break");
  check_member(selection.priority, "selection.priority");
}

} // namespace netgraph::core::detail
