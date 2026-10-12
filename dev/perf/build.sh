#!/usr/bin/env bash
# Build one A/B/A timing harness and its build record (README.md, "Building the two sides").
#
#   bash dev/perf/build.sh REV OUT_DIR
#
# REV is required: `worktree` compiles this checkout's include/ and src/ (uncommitted changes
# included); anything else is a git revision whose include/ and src/ are extracted with
# `git archive` into OUT_DIR/src. The driver (harness.cpp, aba_rows.hpp, aba_api.hpp) always
# comes from this checkout, so both sides of a comparison run the same rows, each compiled
# against its own Core. OUT_DIR must be absent or empty; it receives `harness` and the build
# record `harness.json`.
#
# The compiler is clang++ from PATH (Xcode's shim adds the SDK itself) with exactly the flags of
# the first flags=(...) line below. CXX overrides the compiler; on macOS a set CXX also gets
# -isysroot "$(xcrun --show-sdk-path)" (a toolchain compiler such as "$(xcrun --find clang++)"
# finds no system headers without it), and a set CXX without an SDK is refused.
set -euo pipefail
flags=(-std=c++20 -O3 -DNDEBUG -funroll-loops -fno-math-errno -fno-trapping-math -Wall -Wextra)
[ $# -eq 2 ] && [ -n "$1" ] && [ -n "$2" ] || { echo "usage: build.sh REV|worktree OUT_DIR" >&2; exit 2; }
rev=$1
out=$2
root=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd -P)
case "$out" in /*) ;; *) out=$PWD/$out ;; esac
if [ -e "$out" ] && { [ ! -d "$out" ] || [ -n "$(ls -A "$out")" ]; }; then
  echo "build.sh: $out must be absent or an empty directory" >&2
  exit 2
fi
if [ "$rev" = worktree ]; then
  commit=$(git -C "$root" rev-parse HEAD)
else
  commit=$(git -C "$root" rev-parse --verify --quiet "$rev^{commit}") || { echo "build.sh: no revision $rev" >&2; exit 2; }
fi
# The compiler and its SDK are settled before OUT_DIR is touched: a refusal leaves it as it was.
cxx=${CXX:-clang++}
sdk_path=
sdk_version=
if command -v xcrun > /dev/null 2>&1; then
  sdk_path=$(xcrun --show-sdk-path 2> /dev/null || true)
  sdk_version=$(xcrun --show-sdk-version 2> /dev/null || true)
fi
isysroot=
if [ -n "${CXX:-}" ] && [ "$(uname -s)" = Darwin ]; then
  if [ -z "$sdk_path" ] || [ ! -d "$sdk_path" ]; then
    echo "build.sh: CXX=$CXX is set but xcrun names no macOS SDK, so this compiler would not find the system headers; install Xcode or the Command Line Tools, or unset CXX to build with the clang++ shim" >&2
    exit 2
  fi
  isysroot=$sdk_path
  flags+=(-isysroot "$isysroot")
fi
mkdir -p "$out"
out=$(cd "$out" && pwd -P)
cd "$root"
if [ "$rev" = worktree ]; then
  src=$root
  dirty=$(git status --short -- include src | wc -l | tr -d ' ')
else
  src=$out/src
  mkdir -p "$src"
  git archive "$commit" include src | tar -x -C "$src"
  dirty=0
fi
# A harness without its record would be timed unchecked (aba.py checks a record's driver), so
# a failure from the compile to the written record removes both.
trap 'rm -f "$out/harness" "$out/harness.json"' EXIT
"$cxx" "${flags[@]}" -I"$src/include" dev/perf/harness.cpp "$src"/src/*.cpp -o "$out/harness"
buildv=
if command -v vtool > /dev/null 2>&1; then buildv=$(vtool -show-build "$out/harness" 2> /dev/null || true); fi
OUT=$out/harness SRC=$src REV=$rev COMMIT=$commit DIRTY=$dirty CXXV="$("$cxx" --version | head -1)" \
  FLAGS="${flags[*]}" CXXPATH="$(command -v "$cxx" || true)" CXXSET="${CXX-}" CXXISSET="${CXX+1}" \
  ISYSROOT="$isysroot" SDKPATH="$sdk_path" SDKVERSION="$sdk_version" SDKROOTENV="${SDKROOT-}" \
  SDKROOTSET="${SDKROOT+1}" DTSET="${MACOSX_DEPLOYMENT_TARGET+1}" DT="${MACOSX_DEPLOYMENT_TARGET-}" \
  BUILDV="$buildv" python3 - << 'EOF'
import glob
import hashlib
import json
import os


def sha(path):
    with open(path, "rb") as f:
        return hashlib.sha256(f.read()).hexdigest()


def env_or_null(value, is_set):
    return os.environ[value] if os.environ[is_set] else None


src = os.environ["SRC"]
lines = [
    (os.path.relpath(p, src), sha(p))
    for sub in ("include", "src")
    for p in glob.glob(os.path.join(src, sub, "**", "*"), recursive=True)
    if os.path.isfile(p)
]
driver = sorted(glob.glob("dev/perf/*.hpp") + ["dev/perf/harness.cpp"])
# LC_BUILD_VERSION as vtool prints it: "platform MACOS", "minos 15.0", "sdk 27.0", ...
build = {}
raw = []
for ln in os.environ["BUILDV"].splitlines()[1:]:
    parts = ln.split(None, 1)
    if len(parts) == 2:
        raw.append(f"{parts[0]} {parts[1].strip()}")
        build.setdefault(parts[0], parts[1].strip())
info = {
    "binary": os.environ["OUT"],
    "binary_sha256": sha(os.environ["OUT"]),
    "rev": os.environ["REV"],
    "commit": os.environ["COMMIT"],
    "dirty_core_paths": int(os.environ["DIRTY"]),
    # sha256 of "<sha256>  <path>\n" lines over include/ and src/, sorted by path.
    "core_sources_sha256": hashlib.sha256(
        "".join(f"{h}  {p}\n" for p, h in sorted(lines)).encode()
    ).hexdigest(),
    "driver_sha256": {p: sha(p) for p in driver},
    "compiler": os.environ["CXXV"],
    "flags": os.environ["FLAGS"],
    # The compiler binary the name resolved to, and CXX as given (null: unset, the PATH clang++).
    "compiler_path": os.environ["CXXPATH"] or None,
    "cxx_env": env_or_null("CXXSET", "CXXISSET"),
    # The SDK passed with -isysroot (a set CXX on macOS), else null (the shim's own choice).
    "isysroot": os.environ["ISYSROOT"] or None,
    # xcrun's SDK: the one the clang++ shim and -isysroot use unless SDKROOT names another.
    "sdk": {
        "path": os.environ["SDKPATH"] or None,
        "version": os.environ["SDKVERSION"] or None,
        "sdkroot_env": env_or_null("SDKROOTENV", "SDKROOTSET"),
        "source": "xcrun --show-sdk-path / --show-sdk-version"
        if os.environ["SDKPATH"]
        else "no xcrun SDK on this host: the compiler's own system headers",
    },
    # Requested: MACOSX_DEPLOYMENT_TARGET (null when unset: the compiler's default, the SDK's
    # OS version). Effective: the minimum OS and SDK the binary carries (LC_BUILD_VERSION).
    "deployment_target": {
        "requested": env_or_null("DT", "DTSET"),
        "effective_minos": build.get("minos"),
        "binary_sdk": build.get("sdk"),
        "platform": build.get("platform"),
        "raw": "; ".join(raw) or None,
        "source": "vtool -show-build" if raw else "no LC_BUILD_VERSION read (no vtool, or it printed none)",
    },
    "built_by": "dev/perf/build.sh",
}
with open(os.environ["OUT"] + ".json", "w") as f:
    json.dump(info, f, indent=2)
print(
    f"{os.environ['OUT']}: {info['rev']} {info['commit'][:12]} core "
    f"{info['core_sources_sha256'][:16]} sha256 {info['binary_sha256'][:16]}"
)
EOF
trap - EXIT
