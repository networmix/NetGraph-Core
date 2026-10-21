// Finite-state test controls. These models never disable synchronization in the
// library or create an unsynchronized native data race.
#pragma once

#include <algorithm>
#include <array>
#include <cstddef>
#include <set>
#include <stdexcept>
#include <string_view>
#include <vector>

namespace netgraph::core::detail {

struct PublicationModelResult {
  std::size_t states {0}, mixed {0};
};
inline PublicationModelResult publication_model(bool locked) {
  // pc(writer), pc(reader), owner, load, reservation, program.
  using State = std::array<int, 6>;
  const State initial {0, 0, -1, 0, 0, 0};
  std::set<State> visited;
  std::vector<State> pending {initial};
  PublicationModelResult out;
  for (std::size_t n = 0; n < pending.size(); ++n) {
    const auto current = pending[n];
    for (std::size_t thread = 0; thread < 2; ++thread) {
      State next = current;
      const int pc = current[thread];
      const int end = thread == 0 ? 5 : 3;
      if (pc == end) continue;
      if (pc == 0) {
        if (locked && current[2] != -1) continue;
        if (locked) next[2] = static_cast<int>(thread);
      } else if (pc == end - 1) {
        if (locked) next[2] = -1;
      } else if (thread == 0) {
        next[static_cast<std::size_t>(pc + 2)] = 6;
      } else if (!(current[3] == current[4] && current[4] == current[5])) {
        ++out.mixed;
      }
      ++next[thread];
      if (visited.insert(next).second) pending.push_back(next);
    }
  }
  out.states = visited.size();
  return out;
}

struct PhaseModelResult {
  std::size_t states {0}, partial {0}, python_under_lock {0}, late_hits {0};
};
inline PhaseModelResult phase_model(std::string_view protocol) {
  enum Op { Py, Gate, Ungate, Compute, Publish, Snapshot, Shared, Unshare, Lookup, Hits, Exclusive, Thaw, Unexclusive, Read };
  std::vector<Op> a;
  if (protocol == "publish-result-first") a = {Py, Gate, Compute, Publish, Ungate};
  else if (protocol == "publish-result-after") a = {Gate, Compute, Publish, Ungate, Py};
  else if (protocol == "read-result-failure") a = {Gate, Snapshot, Ungate, Py};
  else if (protocol == "read-result-under-gate") a = {Gate, Snapshot, Py, Ungate};
  else if (protocol == "frozen-hits") a = {Py, Shared, Lookup, Lookup, Hits, Unshare};
  else if (protocol == "frozen-hits-after-lease") a = {Py, Shared, Lookup, Lookup, Unshare, Hits};
  else throw std::invalid_argument("unknown phase protocol");
  const std::vector<Op> b = std::find(a.begin(), a.end(), Shared) != a.end()
      ? std::vector<Op>{Exclusive, Thaw, Unexclusive} : std::vector<Op>{Gate, Read, Ungate};
  // pc A/B, gate owner, exclusive owner, shared count, mode, ledger,
  // published hits, private hits, failed, mode at lease acquisition.
  using State = std::array<int, 11>;
  const State initial {0, 0, -1, -1, 0, 0, 0, 0, 0, 0, -1};
  std::set<State> seen {initial};
  std::vector<State> pending {initial};
  PhaseModelResult out;
  while (!pending.empty()) {
    const State current = pending.back();
    pending.pop_back();
    bool advances = false;
    auto enqueue = [&](const State& next) {
      advances = true;
      if (seen.insert(next).second) pending.push_back(next);
    };
    for (std::size_t t = 0; t < 2; ++t) {
      const auto& code = t == 0 ? a : b;
      if (static_cast<std::size_t>(current[t]) == code.size()) continue;
      const Op op = code[static_cast<std::size_t>(current[t])];
      State next = current;
      ++next[t];
      switch (op) {
        case Py: {
          if (current[2] == 0 || current[4] > 0) ++out.python_under_lock;
          enqueue(next);
          next = current;
          next[0] = static_cast<int>(a.size());
          if (next[2] == 0) next[2] = -1;
          next[4] = 0;
          next[9] = 1;
          enqueue(next);
          break;
        }
        case Gate:
          if (current[2] == -1) { next[2] = static_cast<int>(t); enqueue(next); }
          break;
        case Ungate: next[2] = -1; enqueue(next); break;
        case Shared:
          if (current[3] == -1) { ++next[4]; next[10] = current[5]; enqueue(next); }
          break;
        case Unshare: --next[4]; enqueue(next); break;
        case Publish: next[6] = 1; enqueue(next); break;
        case Lookup: ++next[8]; enqueue(next); break;
        case Hits:
          if (current[5] != current[10]) ++out.late_hits;
          next[7] += current[8]; next[8] = 0; enqueue(next); break;
        case Exclusive:
          if (current[3] == -1 && current[4] == 0) { next[3] = 1; enqueue(next); }
          break;
        case Unexclusive: next[3] = -1; enqueue(next); break;
        case Thaw: next[5] = 1; enqueue(next); break;
        case Compute: case Snapshot: case Read: enqueue(next); break;
      }
    }
    if (!advances) {
      if (static_cast<std::size_t>(current[0]) != a.size() || static_cast<std::size_t>(current[1]) != b.size()) {
        throw std::logic_error("phase model deadlocked");
      }
      if (current[9] != 0 && (current[6] != 0 || current[7] != 0)) ++out.partial;
    }
  }
  out.states = seen.size();
  return out;
}
} // namespace netgraph::core::detail
