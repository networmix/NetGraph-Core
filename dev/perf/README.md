# Timing tool (A/B/A)

Developer-only tool for timing two builds of NetGraph-Core against each other, as `AGENTS.md`
asks for every performance claim. Nothing here ships in the wheel or the sdist, and no timing
runs in CI; only the tool's own tests do (below), and the repository gates (`make check-ci`,
`make sanitize-test`, `make tsan-test`) stay the acceptance tests. Verify correctness first
(bit-identical results where the library is exact, certificates for flows); time second.

**A/B/A** means: the reference build A, the candidate build B, then A again, each block in a
fresh process, repeated seven times; the candidate's time is compared with the mean of the two
reference blocks around it, so slow drift of the machine cancels out. A bootstrap interval over
the seven triples says whether a difference is real.

| File | Role |
| --- | --- |
| `build.sh` | builds one harness binary against the Core of a revision or of the working tree, with its build record |
| `harness.cpp`, `aba_rows.hpp`, `aba_api.hpp` | the harness: the rows (workloads, semantic work per iteration, digests) and the adapter that maps them onto the current public Core API (plus two private functions of the queue rule, below) |
| `aba.py` | runs one A/B/A comparison: blocks, statistics, verdicts, `report.json`, `summary.md`, `blocks.csv`, `meta.txt` |
| `queries/` | the frozen row inputs (`*.json`, `storage/inputs.json`) with their `SHA256SUMS` |
| `test_tool.py` | the tool's own tests (`make perf-tool-test`) |

## Quick start

Compare the working tree (candidate) with a revision (reference); both binaries carry the same
rows, each compiled against its own Core. Run on a quiet machine, nothing else running.

```sh
bash dev/perf/build.sh <ref-rev> build/perf/aba/ref       # e.g. HEAD, or main
bash dev/perf/build.sh worktree build/perf/aba/cand
venv/bin/python dev/perf/aba.py \
    --old-bin build/perf/aba/ref/harness --new-bin build/perf/aba/cand/harness \
    --old-label ref --new-label cand --rows "kernel.*,maxflow.*" --out build/perf/aba/run
```

`aba.py` needs Python 3.10 or newer (the venv's; macOS's system `python3` 3.9 is refused with a
message). The run prints a plan (the expected duration of every row) before the first timed
block, then one line per row, and writes `report.json`, `summary.md` (the table to read),
`blocks.csv` (every block; in an aborted run, those of the finished warm-up pairs and triples)
and `meta.txt` (every hash) to `--out` (default `dev/perf/results/<timestamp>/`, ignored by
git). The run's verdict is the worst row verdict, and a row whose results differ between the two
builds fails: the tool compares two builds of the current Core, so a change it times (an
optimisation) must leave every result identical. A change that is meant to change results is
checked for correctness by the test suites, not timed against the old results. The adapter
`aba_api.hpp` is written against the current public API, apart from the queue check of the
storage search rows, which reads the library's queue rule through the private
`detail::bucket_width` and `detail::uses_buckets`: a reference revision with another API needs
the adapter adapted first.

**Across an API change.** To compare with a reference revision whose API this adapter does not
compile against, build the reference harness in a checkout of that revision instead, with its
own driver adapted to its API and the same rows (`bash dev/perf/build.sh worktree OUT_DIR` in that
checkout), move its `harness.json` aside so that the run does not refuse it for carrying another
driver, and pass the binary as `--old-bin`. The run still requires both harnesses to list the same
rows, record for record, and `report.json` records that side's build record as empty.

**Exit statuses and acceptance.** `aba.py` exits 0 when the run completed, whatever its verdict
and conformance; 1 when it aborted after starting blocks; 2 when it refused an input or failed
before the first block; 130 for Ctrl-C before the first block. Every failure writes its reason
to `report.json` and `summary.md`, with verdict `aborted` and conformance false. A reason names
the side, mode and row when it comes from a block. A refusal before row selection writes a
minimal report with no rows. A failed harness, malformed response, changed or unreadable binary,
non-finite derived value, operating-system error, internal exception or Ctrl-C ends the run.
The first Ctrl-C stops the active harness; later interrupts are held while the outputs are
written. An interrupt during final writing records an abort and rewrites both outputs.

Both files are created before the first block. Each final output is attempted even if another
fails; text that cannot be formed is replaced by a minimal abort record. Complete output text
replaces its file in one operation; if directory permissions no longer allow replacement, an
existing writable file is retried directly with the abort reason. An invalid, ambiguous,
occupied or unwritable `--out` cannot receive a new report; existing unrelated files are never
overwritten. A persistent write failure (for example a full disk), or a process killed by an
uncatchable signal, can prevent output. Write failures are also reported on stderr. `--help`
prints help and writes no run files. Both outputs list every selected row with the step it reached: `measured`, and in an
aborted run `not run`, `warm-up` (its warm-up step did not finish: a warm-up block failed, or a
count or the expected duration derived from them is not finite), `planned` (warmed up, not
timed) or `timing` (its timed step did not finish: a timed block failed, or a row statistic
derived from them is not finite). A run killed by another signal writes nothing more. Its exit
status is no acceptance: a run is accepted when it is conforming and its verdict is `pass`, the
first two items of `summary.md` (and `report.json`'s `conforming` and `verdict`).

## Building the two sides

`build.sh` requires `bash`, `git`, `tar`, `python3` and the compiler on `PATH` (on macOS,
`xcrun` for an explicitly selected compiler). `vtool` supplies optional build diagnostics.

`bash dev/perf/build.sh REV OUT_DIR` compiles this checkout's `harness.cpp` (with `aba_rows.hpp`
and `aba_api.hpp`) against the Core sources of `REV`: `worktree` takes this checkout's `include/`
and `src/` with their uncommitted changes; a git revision is extracted with `git archive` into
`OUT_DIR/src`. `REV` is required (there is no default reference), and a revision git does not
know is refused. `OUT_DIR` must be absent or empty; it receives `harness` and the build record
`harness.json`. The compiler is `clang++` from `PATH` with the first `flags=(...)` line of
`build.sh` (`-std=c++20 -O3 -DNDEBUG -funroll-loops -fno-math-errno -fno-trapping-math
-falign-functions=64 -Wall -Wextra`). `CXX` overrides the compiler; on macOS a set `CXX` also gets `-isysroot "$(xcrun
--show-sdk-path)"` (a toolchain compiler such as `$(xcrun --find clang++)` finds no system headers
without it), and a set `CXX` without an SDK is refused. A refusal exits 2. A failed compile exits
with the compiler's status, and a failure while the build record is written exits non-zero;
either leaves neither `harness` nor `harness.json` (a revision's extracted sources stay in
`OUT_DIR/src`, so a second attempt needs another or an emptied `OUT_DIR`).

Function alignment reduces timing variation due to code placement. Build both sides with
the same flags. The script enables neither `-march=native` nor PGO. These harness flags do not
change the package wheels' Release flags.

### Build records

`harness.json` names everything that changes the binary's bytes:

| Field | Meaning |
| --- | --- |
| `binary`, `binary_sha256` | the output path and the hash of the file written there; the file name is part of the bytes (below) |
| `rev`, `commit`, `dirty_core_paths` | the revision asked for, its commit (`HEAD` for `worktree`), and for `worktree` the number of paths `git status --short` lists under `include/` and `src/` (0 for a revision) |
| `core_sources_sha256` | sha256 of `"<sha256>  <path>\n"` lines over every file under the `include/` and `src/` compiled against (the checkout's or the extracted ones), sorted by path; a path with a part that starts with a dot is left out |
| `driver_sha256` | the sha256 of `dev/perf/harness.cpp` and of every `dev/perf/*.hpp`, by path |
| `compiler`, `compiler_path`, `cxx_env` | the compiler's `--version` line, the program the name resolved to on `PATH`, and `CXX` as given (null when unset) |
| `flags`, `isysroot` | the flags passed (`-isysroot <SDK>` included when added), and that SDK path (null when none was passed) |
| `sdk.path`, `sdk.version`, `sdk.sdkroot_env`, `sdk.source` | the SDK `xcrun` names (the one the `clang++` shim and `-isysroot` use) and `SDKROOT`, which overrides the shim's choice; path and version null when `xcrun` names no SDK |
| `deployment_target.requested` | `MACOSX_DEPLOYMENT_TARGET` (null when unset: the compiler then targets the SDK's own OS version) |
| `deployment_target.effective_minos`, `binary_sdk`, `platform`, `raw`, `source` | the minimum OS and SDK version the binary carries, read back from its `LC_BUILD_VERSION` load command with `vtool -show-build` (the first value of each when the binary carries several; `raw` has every line) |
| `built_by` | `dev/perf/build.sh` |

Besides the sources and the flags, the bytes depend on the **deployment target** (the same
sources built with `MACOSX_DEPLOYMENT_TARGET=15.0` and without it differ), the **output file
name** (the linker signs every arm64 binary ad hoc with the file name as the signature's
identifier, and the binary's UUID follows; copying a binary keeps its bytes, the name at link
time is what counts), the **SDK** and the **compiler**. A pinned binary hash is reproducible
only with all of them equal. `aba.py` refuses a binary whose hash differs from its record. Every field in the table is
required and no other key is accepted, including inside `sdk` and `deployment_target`.
All values are strings except `dirty_core_paths` (a finite nonnegative integer), the three
objects `sdk`, `deployment_target`, `driver_sha256`, and these nullable strings: `compiler_path`, `cxx_env`, `isysroot`, `sdk.path`, `sdk.version`,
`sdk.sdkroot_env`, and `deployment_target.requested`, `effective_minos`, `binary_sdk`, `platform`,
`raw`. A null means the value was unavailable or the environment variable was unset. `binary_sha256` and `core_sources_sha256` are 64 lowercase hexadecimal digits;
`driver_sha256` has exactly this checkout's driver paths and matching hashes. Descriptive
strings are recorded as text; the old build path, revision and compiler are not re-executed or
used to locate inputs.

## Protocol

**Blocks.** Each block is a fresh process of one side's harness running one row: setup
(untimed), one untimed iteration that fixes the reference digest, then a fixed iteration count
with every iteration timed; the harness reports the median iteration time and `block_s`, its own
wall-clock measurement of the timed loop. Every iteration must return the reference digest (else
the harness exits 3 and the run aborts), which the block line carries as 16 lowercase hexadecimal
digits. A block's whole output is its response, and the run
aborts unless it is exactly one block line for the requested row ended by a newline, with
nothing before or after it (a usable line is never picked out of a broken response). A block
line is malformed, and aborts the run, when it is not one JSON object, names a key twice or a
key outside the listed fields, or a required number is missing, of another type, not finite (an integer beyond the largest
float counts as not finite) or outside its domain: the iteration count is an integer in [1, 2^63 - 1]; the median iteration
time and the harness's own peak RSS are positive; `block_s`, `setup_s` and the children's peak
RSS may be 0 (the harness prints the two times rounded to microseconds, and starts no child
process); the statistics no rule reads (`mean_s`, `min_s`, `p10_s`, `p90_s`, `max_s`) are
nonnegative when present. A timed block must report the iteration count it was given. The
numeric options are checked the same way (`--tolerance`, `--rss-floor-mib` and
`--calibration-margin` may be 0, `--seed` is an integer of either sign within the finite float range, the others are positive), and the
values derived from them must be finite (`--target-seconds` times one plus the margin, the floor
in bytes). "Inputs" below lists every input's rule.

**Derived values.** Every number of a block line is finite, and every one but the iteration
count is read as a float, but arithmetic on huge ones need not be finite: the mean of two
medians near 1e308 is infinite, and a ratio against it would read 0, a false pass. So a value
derived from the blocks that is not finite (a mean of A1 and A2, a ratio, the statistic, its
interval, the drift, a median, a peak sum, the memory allowance, a calibrated count, an expected
duration or a finish time beyond any date) aborts the run before the row's statistics or verdict
are recorded or printed.

**Per row.** One warm-up block per side (discarded; it calibrates the iteration counts), then
seven triples A1, B, A2 in that order (fewer only with `--allow-short`, which labels the run
"not gate evidence"). The counts are **per side**: `ceil(target * (1 + margin) / median_side)`
from that side's own warm-up block (a calibrated count above 2^63 - 1 aborts with the side and row named) (`--target-seconds`, default 2, `--calibration-margin`,
default 0.10), so every block on either side lasts about the target even when the two sides'
speeds differ by orders of magnitude. A1 and A2 share side A's count. `--iters` pins the counts
instead (below). The block value is the median iteration time, so the ratio does not depend on
the counts. The warm-up blocks of every selected row run first; the expected duration of every
row and of the run is printed and stored (`report.json` `plan`, each row's `expected`) before
the first timed block.

**Statistics.** `r_t = B_t / mean(A1_t, A2_t)`; the statistic is `median_t(r_t - 1)`; its 95
percent percentile interval comes from `--resamples` (default 10,000) bootstrap resamples of
whole triples, seeded with `--seed`; drift is `median_t |A1_t / A2_t - 1|`.

**Verdicts per row.**

* timing: `inconclusive` when the drift exceeds the tolerance (`--tolerance`, default 0.03),
  else `pass` when the interval's upper bound is within the tolerance, else `fail`;
* blocks: `fail` when any timed block's `block_s` is below `--target-seconds`;
* same-side clusters: per side, the timed blocks' per-iteration times sorted; a gap between two
  consecutive blocks above the tolerance with at least two blocks on each side marks the side
  `bimodal` and fails the row even when the median passes;
* memory (below);
* digests: every timed block of both sides must report one and the same digest.

The row is `fail` if any check fails, else `inconclusive` if timing is, else `pass`. The run's
verdict is `aborted` when the run aborted, else the worst row verdict. A run is **conforming**
(gate evidence) only when it did not abort, `--allow-short` was not given (without it fewer than
seven triples or a target below 2 s are refused; a measured block shorter than the target fails
its row), no timed block was short and the
allocator environment is normal; the summary's first line says which.

**Memory.** Each block's peak RSS (resident set size: the memory the process holds in RAM) of the
harness process tree is the larger of the harness's own `getrusage` (self + children)
`ru_maxrss` and `ru_maxrss` from the orchestrator's `os.wait4` (macOS reports bytes, Linux KiB,
converted). The harness's own record of the iteration times is the same in every block, whatever
its iteration count: allocated and written before the row's setup, it keeps a count per duration
below 65,536 clock ticks (65.5 µs with the nanosecond clocks of macOS and Linux) and a list of
the longer ones with 65,536 slots, which grows only past 4.3 s of such iterations in one block
(1 MiB of bookkeeping in all; its statistics are those of the sorted list of every iteration's
time). A faster side runs more iterations, and a record that grew with them would count in its
peak. Reference = the largest A1 block peak, candidate = the largest B block peak (over the
timed triples); the row fails when the candidate exceeds `max(1.10 * reference, reference + 8
MiB)` bytes (`--budget`, `--rss-floor-mib`). A1, the first reference block of every triple, is
fixed in advance: a larger A2 peak or a larger warm-up peak never rescues a failure. Why the
maxima: macOS's allocator keeps freed large blocks in a cache whose pages stay resident, so on
several rows the block peaks of any build land on one of two levels about 9 to 10 MiB apart at
random; medians can fall on different levels by chance and fail identical binaries, while the
largest peaks of both sides reach the upper level. A row records `rss_a1_max_bytes`,
`rss_b_max_bytes`, `rss_allowance_bytes`, `rss_balanced_ratio`, `rss_a2_max_bytes` (a
diagnostic) and `rss_verdict`. The check compares the largest observed peaks; it is no absolute
memory bound.

**Allocator environment.** `aba.py` neither sets nor clears an allocator variable
(`MallocLargeCache=0` would hide the cache levels and with them a candidate's genuine peak). It
records the harness's environment in `report.json` `allocator_env` and on the summary's
"allocator environment" line: whether `MallocLargeCache` is set, and every variable whose name
starts with `Malloc` or `MALLOC_` or is `DYLD_INSERT_LIBRARIES`, `LD_PRELOAD` or
`GLIBC_TUNABLES`. A run with any such variable is not conforming.

**Provenance.** `meta.txt` lists the sha256 of both harness binaries (each checked against its
build record when it has one) and of the orchestrator and driver sources (`aba.py`,
`harness.cpp`, `aba_rows.hpp`, `aba_api.hpp`, hashed at the start). A side's binary is
re-hashed before and after every block and at the end of the run; a difference, or a binary
that can no longer be read, aborts the run. `--expect-meta meta.txt` refuses to start unless
every hash equals a previous run's (the refusal names the hashes that differ). Both sides'
harnesses run in the repository root with the environment `aba.py` was started in, unchanged:
`aba.py` sets and clears no variable.

**Machine state.** Each block records the 1-minute load before and after; the summary has a
machine-state line (load averages, top CPU consumers) at the start and the end of the run.

### The sides

Every selected row runs as A1 = `--old-bin`, B = `--new-bin`, A2 = `--old-bin` with every rule
above: two builds of the current Core. Both must run this checkout's driver (`build.sh` builds
them so): the two harnesses must list the same rows, record for record, and a side's build
record, when it has one, must carry this checkout's driver as its `driver_sha256`: exactly the
driver files `meta.txt` names (`harness.cpp`, `aba_rows.hpp`, `aba_api.hpp`), each with the
sha256 it has at the start of the run; else the run is refused, so a harness built before the
rows or the adapter changed is not timed. A side without a build record is not checked this way
(`report.json` records its build record as empty).

### Selecting rows

`--rows` takes comma-separated globs over row names and groups (`all`, `construct.*`,
`kernel`); spaces around a glob are ignored, and no text means `all`. **A glob that matches no
row or group refuses the run** before any block, with the list of the groups that exist, even
beside globs that match. The report records, per glob, the rows it selected and the rows of its
groups it did not, and per group the rows the whole selection left out.

Other options: `--iters row=N,...` (both sides) or `row=NA:NB,...` (per side) pins the timed
blocks' iteration counts of a row in place of the calibrated ones (every item names a selected
row, once; a count is at most 2^63 - 1, the largest the harness reads); `--triples` (default 7,
fewer only with `--allow-short`, "Per row" above); `--seed` (default 20261003) and
`--resamples` of the bootstrap; `--old-label` and `--new-label` (defaults `old` and `new`; two
different labels, below); `--out` (absent or empty).

**Names.** A row's name and its group are 1 to 64 of the characters `a`-`z`, `0`-`9`, `.`,
`_` and `-`; a side's label is 1 to 32 of `a`-`z`, `0`-`9`, `_` and `-` (no dot). They make the
stderr logs' names, `blocks/<row>.<label>.<mode>.stderr`, and the fields of `meta.txt`'s lines:
lower case only, so no two names share a log on a file system that ignores case, and no dot in a
label, so a log's name splits back into one row, label and mode. A log collects the stderr of
every block of its row, side and mode, in run order; a log that is still empty after a block is
removed, and an abort for a block names the log only when that block wrote to it ("nothing on
stderr" otherwise). A path the run records (the binaries', `--out`, `--expect-meta`, this
directory's) has printable characters only: no tab, line break or other control or separator
character but the space.

### Inputs

Everything `aba.py` reads that it did not produce in the same process has a rule, checked where
it is read; anything else is refused (exit 2 before the first block) or aborts the run (once
blocks have run). Only the keys listed below are accepted, with plain scalar values except the build record's
three explicitly listed objects. JSON is read with no key twice in an object (Python's reader would otherwise
keep the last value of a repeated key).

| Input | Rule |
| --- | --- |
| a harness response (every mode) | exit 0; UTF-8 lines, each ended by a newline and each one JSON object; nothing else |
| `aba-list` | at least one row; each a string `row` and a string `group`, each a name ("Names" above) (`desc`, `queries` strings when present); no other keys; no row twice; the same list, record for record, from both harnesses |
| `aba-run`, `aba-calibrate` | exactly one record; required keys `row`, `iterations`, `median_s`, `block_s`, `setup_s`, `digest`, `rss_self_bytes`, `rss_children_bytes`; optional `mean_s`, `min_s`, `p10_s`, `p90_s`, `max_s`; no others. The row is the requested one, its numbers as in "Blocks" above, its digest 16 lowercase hexadecimal digits; `aba-run`'s count the requested one; the values derived from them finite ("Derived values" above) |
| `--rows` | the characters of names, the glob characters `*`, `?`, `[`, `]` and `!`, commas and spaces; every glob matches a row or a group on its own; at least one glob |
| command line | full option names only, each at most once; argparse checks required options and value types; raw arguments are strings, recorded only after all checks pass |
| numeric options | of their type, finite, in their domain ("Blocks" above); the values derived from them finite |
| `--iters` | every item `row=N` or `row=NA:NB` (positive integers up to 2^63 - 1) naming a selected row, once |
| `--expect-meta` | a printable path; UTF-8, newline-ended lines: `#` comments; `label`, side A/B, valid label (tab-separated, side at most once); or exactly five tab-separated hash fields: A/B, `binary`, `harness`, sha256, printable path; or `run`, `source`, one of the four source files, sha256, printable path. No hash named twice, at least one; exactly this run's hash names and values |
| `--old-bin`, `--new-bin` | a readable file whose resolved path is printable; its build record (`<resolved path>.json`; a `.json` beside a symbolic link to the binary that is not that file is refused), when there is one, the complete typed object of "Build records", whose `binary_sha256` is the binary's and whose `driver_sha256` is this checkout's driver ("The sides" above) |
| `--old-label`, `--new-label` | two different labels ("Names" above; they name the sides in every block, stderr log and `meta.txt`) |
| `--out` | a printable path, absent or an empty directory (the stderr logs append; an earlier run's files would mix in) |
| operating-system errors | before the first block a refusal, from the first block on an abort ("Exit statuses" above) |
| the environment | OS-provided string keys and values, passed to the harness unchanged; allocator variables recorded by the prefix/name rule above (a run with one is not conforming); `PATH` resolves diagnostic commands; no other environment value is interpreted by the driver |
| system readings | block loads, elapsed times and wait4 peak: finite nonnegative numbers, checked before storage; CPU count: positive integer or null; machine load: three finite nonnegative numbers, or null when unavailable. CPU/platform names and timestamps are text |
| diagnostic commands | `ps`: decoded text with replacement for bad UTF-8; among its first five lines only a finite nonnegative CPU share and nonempty command name are kept; unavailable output means no consumers. `sysctl`: decoded text, or empty when unavailable |
| other file reads | binary and source contents are arbitrary readable bytes, hashed by sha256; log sizes are OS file sizes used only to name a nonempty stderr log; `--out` directory entries only test emptiness. The query JSON and hash files are not read by `aba.py`: their frozen values are compiled into the harness, and the tool tests check them against those files |
| the Python running `aba.py` | 3.10 or newer |

## Rows

`harness aba-list` prints every row as JSON (`row`, `group`, `desc`, and `queries` for a row
with a frozen query list: the file and its sha256). The harness's other modes are `aba-run ROW N`
(`N` a positive integer) and `aba-calibrate ROW SECONDS` (`SECONDS` a positive number). Exit
status 1 is a failure inside an iteration (the untimed first one included), 2 a usage error or
an unknown row, 3 a digest that changed between iterations, 4 a row whose input or setup check
failed.

| Group | Rows | Work per iteration |
| --- | --- | --- |
| `spf` | `spf.full.{random16k,clos,grid10k,elig_residual16k,wide_fallback16k}` | 8 (grid: 4) full forward searches; `elig_residual16k` with static exclusions and a capacity-gated residual; `wide_fallback16k` with costs up to 2^20 (the heap) |
| `spf` | `spf.dag.{near65k,far16k}` | 64 / 8 searches to frozen targets (`queries/near65k.json`, `far16k.json`) |
| `spf` | `spf.to.random16k`, `spf.to.targetset_clos` | 8 reverse searches; 16 reverse searches rooted at a destination set (`queries/targetset_clos.json`) |
| `spf` | `spf.single.canonical16k` | 8 canonical single paths to frozen targets (`queries/canonical16k.json`) on the residual |
| `construct` | `construct.random65k`, `construct.clos` | graph construction from arrays, unit conversion included |
| `kernel` | `kernel.{proportional,eb,lossy}.clos`, `kernel.{proportional,eb}.grid`, `kernel.weighted.clos` | placements on fixed cost-only DAGs into one ledger |
| `capacity` | `capacity.costonly_{eb,prop}.failures` | ShortestPaths capacity over 32 failure states x 8 pairs |
| `maxflow` | `maxflow.exact.{totals,detailed}` | exact maximum flow over 8 failure states x 4 pairs |
| `sensitivity` | `sens.exact`, `sens.policy_all_edges`, `sens.equal_global` | exact-model sensitivity; policy-model loss for every eligible edge; equal-source-allocation CentralizedTE, every edge solved afresh |
| `sensitivity` | `sens.exact.{cut,partial,branch,reroute,capped}` | full exact-model sensitivity calls on a 64-link saturated chain, 32 partial-repair branches, 32 branches with a second upstream route and a shared downstream cut, 32 branches with a parallel upstream member and a downstream cut, and the chain at binding offer 1; nonuniform costs, baseline/classification/copies/results included |
| `pooled` | `pooled.{ecmp,lossy}.clos`, `pooled.from_dag.clos`, `pooled.from_dag_build.clos` | 8 demands of 8 pooled sources on one ledger, Lossless / Lossy above capacity; the same demands by Pinned over prebuilt `from_dag` programs, and building those programs |
| `global` | `global.equal.random` | equal-source-allocation CentralizedTE capacity, 2 demands of 4 sources |
| `program` | `program.forward.global`, `program.replay.lossy`, `program.stale.lossy`, `program.replace.costonly`, `program.ring.lossy` | forwarding programs: a CentralizedTE program forwarded, replayed under failures (Lossy), stale under black holes, ShortestPaths re-prepared by replace; a 64-head program on a ring of 64 links whose 64 barriers form one cyclic group (6 rounds) |
| `batch` | `batch.{scalar,ordered,simultaneous}.clos`, `batch.simultaneous.ring` | several demands at once: the 8 pooled Lossy demands as 8 scalar `replace` calls, one ORDERED `apply_actions` list (the same work and digest) and one SIMULTANEOUS list (the joint evaluator; the demands share no link, so the same digest); the ring's 64 tunnels as 64 policies in one SIMULTANEOUS list (the same traffic amounts and digest as `program.ring.lossy`; 64 publishing policies against one, so a 2x operating budget that includes their per-policy bookkeeping) |
| `t1` | `t1.setup.fabric`, `t1.forward.fabric` | T1 with 16 tunnels per head, 12 demands on one ledger: the setup with its holds; the installed programs forwarded under 32 offers |
| `t2` | `t2.exact16.headroom`, `t2.exact16.cost`, `t2.upto256.headroom`, `t2.size.mesh` | T2 on the same fabric and demands: 16 tunnels per head with and without `prefer_headroom`, up to 256 with it; 16 tunnels for 10 heads on a 500-node mesh (2000 links, 5/8 to 8/8 free), sizing-dominated |
| `msd` | `msd.loop` | the maximum-supported-demand bracket and bisection with a per-series DAG cache |
| `ksp` | `ksp.fabric3_k8` | k = 8 shortest paths for 8 pairs |
| `ksp` | `ksp.longprefix_k64` | 96-edge common prefix and 64 branches, k = 64; all prefix and result work timed |
| `ksp` | `ksp.weighted_k32` | random 128 x6, seed64, four pairs, k = 32 |
| `ksp` | `ksp.collapse` | 12-hop parallel prefix, 32 branches, HIGHER_CAPACITY Collapse |
| `ksp` | `ksp.ties` | 64 equal-cost branches, each of the four tie rules |
| `ksp` | `ksp.k1` | the long-prefix population with k = 1 |
| `cache` | `cache.cold`, `cache.warm`, `cache.frozen.reads` | 64 direct lookups from a cleared, prepopulated building, or prepopulated frozen cache |
| `concurrency` | `concurrency.place_zero`, `delivered`, `program`, `config` | 10,000 calls of one operation on a prepared one-link policy (capacity 16, cost 3, q=1, replace(4)); config is the immutable control |
| `concurrency` | `concurrency.two_ledgers_serial`, `concurrency.two_ledgers_parallel` | 32 evaluate(32) calls on each independent ledger; 384-node/2816-edge layered graph, q=1, capacity 64, two thread launches and joins in both modes |
| `eligibility` | `eligibility.build_intersect` | 64 failure states built and intersected |
| `construct` | the fourteen storage construction rows (`construct.wgrid256`, `fanout_clos`, `random1m_nat`, `random1m_scat`, `wide65k`, `k65536`, `k65537`, `star1m`, `sparse1m`, `empty1m`, `random65k_bprime`, `random65k_reversed`, `star1m_smallk`, `instar1m_smallk`) | the constructor on prepared input: many short rows, the sort fallback above 65,536 distinct costs and both sides of that limit, very long rows, mostly isolated nodes, an empty graph, input already in or reversed from edge-id order |
| `spf` | the six storage search rows (`spf.to.{wgrid256,random1m,wide_fallback16k,width65536,width65537}`, `spf.full.random1m`) | 8 full searches each: the bucket queue, a large scattered graph, the heap, and both sides of the bucket-width limit 65,536 |

Every row's `desc` (in `aba-list`) describes its input and its work; the generators, seeds and
sources are in `aba_rows.hpp`. The row names are fixed (the frozen inputs and earlier reports
name them): `proportional` and `prop` in a name mean the OPTIMIZED route allocation, and `lossy`
means the EQUAL route allocation with LOSSY admission.

The tiny native operation rows use added time per call against the same row on the reference:
at most 0.1 microseconds for a gated call. The two direct-hit cache rows have the same
0.1-microsecond allowance per lookup; `cache.cold` keeps the ordinary 3% tolerance. The immutable
configuration control keeps the ordinary tolerance. For two ledgers compare candidate serial
against candidate parallel on identical work, with a speedup target of at least 1.5; both modes
include two thread launches/joins. These budgets are interpreted separately from the default
percentage verdict of `aba.py`. Correctness checks assert full delivery (one unit per head,
below every link's capacity) before accepting a time.

**The storage rows** exercise the graph's storage and construction. Every input comes from a
generator in `aba_rows.hpp` (`storage_inputs()`) with fixed seeds, one random draw per
statement, explicit costs and positive integral capacities. Their timed digests read only
API-independent values (node-indexed counts and distances, costs, physical capacities, edges
named by their input position) on 16 sampled nodes. A search row also checks its queue in setup
(exit 4 otherwise): the library's choice (the bucket width W and buckets or heap) and the width
read from the input (max cost / cost GCD + 1) must both equal the frozen expectation, and its
roots must be the rule's `(i * 997 + 1) mod N`. A sampled digest is not correctness evidence.

### Frozen inputs

The destination rows' targets, the endpoint set and the storage inputs are frozen, so every
build works on exactly the same queries: `queries/*.json` with `queries/SHA256SUMS` (each query
list with its row, its input generator and digest, and its rule), `queries/storage/inputs.json`
with `queries/storage/SHA256SUMS` (each storage input's shape, full-content digest and
`input_sha256`, and each storage search row's roots and queue), and the two frozen blocks and the
storage search table of `aba_rows.hpp` that compile their values in. A row refuses an input its
generator no longer reproduces (exit 4): the query lists by the input digest, the storage inputs
by their shape and full-content digest. The byte hashes inside the JSON files (`arrays_sha256`,
`input_sha256`) record the inputs as they were frozen; no tool in the checkout recomputes them or
regenerates the frozen files. A deliberate change of a generator needs new values in the JSON
files, their `SHA256SUMS` and `aba_rows.hpp` together, and the tool's tests check that they agree
(and that every frozen width up to 65,536 takes the buckets and every wider one the heap).

## The tool's own tests

```sh
make perf-tool-test          # or: venv/bin/python -B dev/perf/test_tool.py -v
```

Run them whenever a file under `dev/perf/` changes. CI runs them on one job (macOS, Python
3.14): the adapter compiles against the library's API and two private functions, so a library
change that breaks the tool fails there. They need Python 3.10 or newer, `git`, `tar`, `python3` and `clang++` on
`PATH` for one unoptimised build of the harness's driver and Core (a few seconds); they time
nothing. They test the rules of this file:

* `aba.py` end to end with fake harnesses (short Python scripts that answer like the harness):
  the protocol's block order and counts, every verdict rule, the exit statuses and the summary's
  first lines, every input of "Inputs" (each refusal and each abort, operating-system errors
  before and after the first block included), every abort again with several rows selected
  (each row's status in both outputs), an interrupt during the blocks, a second one while the
  first is handled, and one while the final outputs are written (a fake `ps` sends it), the
  child environment, the allocator record and the stderr logs; and its parts directly: the
  memory rule's boundary, the bootstrap interval, the cluster rule, the row selection, the
  block-line parser, a block's stderr log, derived values that are not finite, the SIGINT handler
  and the hashes of the tool's sources;
* `build.sh` with a fake compiler (no real build) in a scratch repository of two revisions:
  its refusals, a failing compile, a failing build-record step, the sources of the revision asked
  for, and every field of its build record;
* the harness's driver compiled with a few extra test rows: `aba-list`, one iteration of every
  row, the block line (its own peak RSS in bytes, the same for 1,000 and 2,000,000 iterations;
  its statistics against the sorted list of the same times), every exit status, and the
  frozen-input, root and queue checks;
* the frozen files: both `SHA256SUMS` files, and the values compiled into `aba_rows.hpp` against
  the JSON files (each storage row on its input).

The `sens.exact`, cut/repair population and `sens.policy_all_edges` rows also have
`.w2` counterparts: the identical full query on two workers. Unsuffixed rows use one worker.
