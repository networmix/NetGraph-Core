#pragma once

#include <stdexcept>
#include <string>

namespace netgraph::core {

struct TypeError : public std::runtime_error {
  using std::runtime_error::runtime_error;
};

struct ValueError : public std::runtime_error {
  using std::runtime_error::runtime_error;
};

struct RuntimeError : public std::runtime_error {
  using std::runtime_error::runtime_error;
};

} // namespace netgraph::core
