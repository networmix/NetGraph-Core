/* Membership of enumeration values (internal). An enum-typed field can hold a value
   outside its enumeration (a C++ cast, a Python enum built from an integer); such a value is
   outside the domain and is rejected (std::invalid_argument naming the field) before it
   reaches any rule or record. */
#pragma once

#include <stdexcept>
#include <string>

#include "netgraph/core/options.hpp"
#include "netgraph/core/search.hpp"

namespace netgraph::core::detail {

[[noreturn]] inline void not_a_member(const char* field, int value, const char* enumeration) {
  throw std::invalid_argument(std::string(field) + ": " + std::to_string(value) + " is not a value of " + enumeration);
}

inline void check_member(RouteAllocation v, const char* field) {
  switch (v) {
    case RouteAllocation::Equal:
    case RouteAllocation::Weighted:
    case RouteAllocation::Optimized: return;
  }
  not_a_member(field, static_cast<int>(v), "RouteAllocation");
}

inline void check_member(Admission v, const char* field) {
  switch (v) {
    case Admission::Lossless:
    case Admission::Lossy: return;
  }
  not_a_member(field, static_cast<int>(v), "Admission");
}

inline void check_member(SourceAllocation v, const char* field) {
  switch (v) {
    case SourceAllocation::Equal:
    case SourceAllocation::Optimized: return;
  }
  not_a_member(field, static_cast<int>(v), "SourceAllocation");
}

inline void check_member(WeightBasis v, const char* field) {
  switch (v) {
    case WeightBasis::Capacity:
    case WeightBasis::Explicit:
    case WeightBasis::ResidualSnapshot: return;
  }
  not_a_member(field, static_cast<int>(v), "WeightBasis");
}

inline void check_member(RouteTieBreak v, const char* field) {
  switch (v) {
    case RouteTieBreak::FewerHops:
    case RouteTieBreak::HigherCapacity:
    case RouteTieBreak::HigherAvailableCapacity:
    case RouteTieBreak::HigherHeadroom: return;
  }
  not_a_member(field, static_cast<int>(v), "RouteTieBreak");
}

inline void check_member(Direction v, const char* field) {
  switch (v) {
    case Direction::Forward:
    case Direction::Reverse: return;
  }
  not_a_member(field, static_cast<int>(v), "Direction");
}

inline void check_member(SpfQueue v, const char* field) {
  switch (v) {
    case SpfQueue::Auto:
    case SpfQueue::Heap: return;
  }
  not_a_member(field, static_cast<int>(v), "SpfQueue");
}

} // namespace netgraph::core::detail
