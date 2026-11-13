"""Prepare an isolated pinned source export and experiment-only wheel environment."""

import difflib
import json
from pathlib import Path
import shutil
import subprocess
import tarfile

ROOT = Path(__file__).resolve().parents[2]
TASK = Path(__file__).resolve().parent
SCRATCH = ROOT / "scratch/split"


def replace(text, old, new):
    assert text.count(old) == 1, old
    return text.replace(old, new)


def main():
    source, venv = SCRATCH / "core", SCRATCH / "venv"
    assert not source.exists() and not venv.exists()
    archive = SCRATCH / "core-ref.tar"
    subprocess.run(
        [
            "git",
            "-C",
            str(ROOT / "worktrees/core-ref"),
            "archive",
            "9cd7642",
            "-o",
            str(archive),
        ],
        check=True,
    )
    source.mkdir()
    with tarfile.open(archive) as stream:
        stream.extractall(source, filter="data")
    patch = ""
    for relative in ("bindings/python/native_call.hpp", "bindings/python/module.cpp"):
        path = source / relative
        old = path.read_text()
        if path.name == "native_call.hpp":
            new = replace(
                old,
                "#include <chrono>",
                "#include <chrono>\n#include <array>\n#include <cstdint>\n#include <ctime>",
            )
            new = replace(
                new,
                "namespace netgraph::core::bindings {",
                """namespace netgraph::core::bindings {
// Experiment-only cumulative counters; no allocation or cross-thread writes.
inline thread_local std::array<std::uint64_t, 7> split_totals{};
inline std::uint64_t split_wall() noexcept {
  return std::chrono::duration_cast<std::chrono::nanoseconds>(
      std::chrono::steady_clock::now().time_since_epoch()).count();
}
inline std::uint64_t split_cpu() noexcept {
  timespec value{};
  if (clock_gettime(CLOCK_THREAD_CPUTIME_ID, &value) != 0) split_totals[6] = 1;
  return std::uint64_t(value.tv_sec) * 1000000000ULL + value.tv_nsec;
}""",
            )
            new = replace(
                new,
                "  Detached() : initialized_(Py_IsInitialized() != 0), state_(PyEval_SaveThread()) {}",
                """  explicit Detached(bool deletion = true)
      : initialized_(Py_IsInitialized() != 0), state_(PyEval_SaveThread()),
        wall_(split_wall()), cpu_(split_cpu()), deletion_(deletion) {}""",
            )
            new = replace(
                new,
                "    attach(state_);",
                """    const auto c1 = split_cpu();
    const auto w1 = split_wall();
    attach(state_);
    const auto w2 = split_wall();
    const auto c2 = split_cpu();
    split_totals[0] += w1 - wall_;
    split_totals[1] += c1 - cpu_;
    split_totals[2] += w2 - w1;
    split_totals[3] += c2 - c1;
    ++split_totals[4];
    split_totals[5] += deletion_;""",
            )
            new = replace(
                new,
                "  PyThreadState* state_;",
                "  PyThreadState* state_;\n  std::uint64_t wall_, cpu_;\n  bool deletion_;",
            )
            new = replace(
                new, "  const Detached detached;", "  const Detached detached(false);"
            )
        else:
            new = replace(
                old,
                '#include "registration.hpp"',
                '#include "registration.hpp"\n#include "native_call.hpp"\n#include <pybind11/stl.h>',
            )
            extra = """
  m.def("_split_snapshot", [] { return bindings::split_totals; });
  m.def("_split_probe", [](int count, bool fail) {
    for (int i = 0; i < count; ++i) bindings::native_call([&] {
      if (fail) throw std::invalid_argument("split exception control");
    });
  });
  m.def("_split_clock_reads", [](int count) {
    std::array<std::uint64_t, 6> out{};
    const auto begin_w = bindings::split_wall(), begin_c = bindings::split_cpu();
    for (int i = 0; i < count; ++i) {
      const auto w0 = bindings::split_wall(), c0 = bindings::split_cpu();
      const auto c1 = bindings::split_cpu(), w1 = bindings::split_wall();
      const auto w2 = bindings::split_wall(), c2 = bindings::split_cpu();
      out[0] += w1-w0; out[1] += c1-c0;
      out[2] += w2-w1; out[3] += c2-c1;
    }
    out[5] = bindings::split_cpu()-begin_c;
    out[4] = bindings::split_wall()-begin_w;
    return out;
  });
"""
            position = new.rfind("\n}")
            assert position > 0
            new = new[:position] + extra + new[position:]
        path.write_text(new)
        patch += "".join(
            difflib.unified_diff(
                old.splitlines(True),
                new.splitlines(True),
                fromfile="a/" + relative,
                tofile="b/" + relative,
            )
        )
    (TASK / "split.patch").write_text(patch)
    shutil.copytree(ROOT / "scratch/consumer/ref/venv", venv, symlinks=True)
    print(
        json.dumps(
            {
                "source": str(source),
                "venv": str(venv),
                "summary": "PASS: isolated pinned export patched and reference venv copied",
            }
        )
    )


if __name__ == "__main__":
    main()
