"""Export pinned source, apply the experiment-only binding patch, copy the venv."""

import hashlib
import json
import shutil
import subprocess
import tarfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SCRATCH = ROOT / "scratch/clocks"
SOURCE = SCRATCH / "core-delay"
VENV = SCRATCH / "venv-delay"


def main():
    assert not SOURCE.exists() and not VENV.exists(), "owned build already exists"
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
    SOURCE.mkdir()
    with tarfile.open(archive) as stream:
        stream.extractall(SOURCE, filter="data")
    module = SOURCE / "bindings/python/module.cpp"
    old = module.read_text()
    extra = """
  m.def("_clocks_delay", [](std::uint64_t ns) {
    std::uint64_t elapsed;
    {
      const py::gil_scoped_release release;
      const auto start = std::chrono::steady_clock::now();
      if (ns != 0) std::this_thread::sleep_for(std::chrono::nanoseconds(ns));
      elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(
          std::chrono::steady_clock::now() - start).count();
    }
    return elapsed;
  });
  m.def("_clocks_native_cpu", [](std::uint64_t ns) {
    std::uint64_t consumed = 0;
    {
      const py::gil_scoped_release release;
      std::thread worker([&] {
        auto cpu = [] {
          timespec value{};
          if (clock_gettime(CLOCK_THREAD_CPUTIME_ID, &value) != 0)
            throw std::runtime_error("thread CPU clock failed");
          return std::uint64_t(value.tv_sec) * 1000000000ULL + value.tv_nsec;
        };
        const auto start = cpu();
        while ((consumed = cpu() - start) < ns) {}
      });
      worker.join();
    }
    return consumed;
  });
"""
    new = (
        "#include <chrono>\n#include <thread>\n#include <ctime>\n#include <cstdint>\n"
        + old
    )
    position = new.rfind("\n}")
    assert position > 0
    new = new[:position] + extra + new[position:]
    module.write_text(new)
    import difflib

    patch = "".join(
        difflib.unified_diff(
            old.splitlines(True),
            new.splitlines(True),
            fromfile="a/bindings/python/module.cpp",
            tofile="b/bindings/python/module.cpp",
        )
    )
    (Path(__file__).parent / "delay.patch").write_text(patch)
    shutil.copytree(ROOT / "scratch/consumer/ref/venv", VENV, symlinks=True)
    site = VENV / "lib/python3.14/site-packages"
    for path in site.glob("*netgraph*"):
        assert path.resolve().is_relative_to(VENV), path
    metadata = {
        "source": str(SOURCE),
        "venv": str(VENV),
        "patch_sha256": hashlib.sha256(patch.encode()).hexdigest(),
        "summary": "clean pinned export patched; reference venv copied without modifying the original",
    }
    (SCRATCH / "build-preparation.json").write_text(
        json.dumps(metadata, indent=2) + "\n"
    )
    print(json.dumps(metadata), flush=True)


if __name__ == "__main__":
    main()
