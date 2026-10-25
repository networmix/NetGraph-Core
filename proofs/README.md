# NetGraph-Core proofs

This is a Lean 4 project. Lean is a proof assistant: a program that checks a mathematical proof
step by step. The project proves properties of the library's round-robin rule, which splits a
number of whole units among members in order (10 units over 3 members give 4, 3 and 3), and of
the increments it gives when an offer grows over several calls. Each theorem and its proof are
explained in plain words, with worked examples, in
[docs/formal/theorems/round-robin.md](../docs/formal/theorems/round-robin.md).

The proofs are about the mathematical rule. Separate tests compare it exactly with the native
package on chosen inputs; nothing here proves that the C++ code implements the rule on every
input. Installing, building or running the package never needs Lean, and the wheel contains
none of this directory.

## What you need

- Lean's toolchain manager, [Elan](https://lean-lang.org/install/). It installs the Lean version
  named in `lean-toolchain`; keep that file as it is.
- Git, which Lake (Lean's build tool) uses to fetch mathlib, the mathematics library the proofs
  use. `lakefile.toml` pins mathlib to one Git revision and `lake-manifest.json` pins every
  dependency of it, so every check uses the same sources.
- Python 3.11 or newer, for the inventory checker `Checks/audit.py` (standard library only).

No native build and no Git history of this repository are needed.

## Checking the proofs from the source archive

From the root of the extracted source archive (sdist), with an absolute output directory of
your choice:

```bash
cd proofs
lake exe cache get
lake --wfail build
lake exe proof-audit --output /absolute/output/check.json
python3 -I -B Checks/audit.py --verify-report /absolute/output/check.json
lake env leanchecker --fresh Checks
```

- `lake exe cache get` downloads mathlib's prebuilt files for the pinned revision, so the build
  does not recompile mathlib; it only saves time.
- `lake --wfail build` checks every proof, treating warnings as errors.
- `lake exe proof-audit` compares each theorem in the built project with its record in
  `theorems.json`: the exact statement, the definitions it depends on, the axioms it rests on
  (only `propext`, `Classical.choice` and `Quot.sound` are accepted), its explanation, the
  sections of `docs/semantics.md` and the native source files it describes, and the pinned
  toolchain and dependencies. It also runs the examples against the expected results in
  `cases.json`. Any difference fails the check. `Checks/audit.py --verify-report` checks that the
  report is complete and matches these inputs.
- `lake env leanchecker --fresh Checks` checks every declaration again in a fresh environment
  with Lean's kernel, the small core that accepts proofs.

The archive holds every local file these commands read; Lean and mathlib are downloaded.

## Checking in a Git checkout

The repository wraps the same commands in `dev/proofs/run.sh` (the `dev/` directory is not in the
source archive). From the repository root, after `cd proofs && lake exe cache get`:

```bash
bash dev/proofs/run.sh check     # build, audit and examples; writes build/proofs/check.json
PROOFS_DEPENDENCY_ORIGIN=cache+local bash dev/proofs/run.sh replay   # check, then the fresh replay
bash dev/proofs/run.sh test      # the gate refuses planted errors, in a disposable copy
```

`PROOFS_DEPENDENCY_ORIGIN` records where mathlib's built files came from: `cache+local` for the
downloaded cache plus local builds, `source` for a build from source, `mixed` for both. Set
`PROOFS_PYTHON` to an absolute interpreter path if `python3` is not the one to use.

To compare the rule with a native wheel, install that wheel in a clean environment first:

```bash
PROOFS_NATIVE_PYTHON=/absolute/path/to/venv/bin/python \
PROOFS_WHEEL=/absolute/path/to/netgraph_core.whl \
bash dev/proofs/run.sh bridge
```

The bridge checks that the installed package is that wheel, then compares results exactly with
`cases.json`: sources passed out of node order come back in canonical order, and repeated calls
through `ShortestPaths` and `Pinned` add the right increments, including after the total offer
passes 2^64. Cases near 2^127 are checked on the model only; the public calls cannot reach that
total. The report is `build/proofs/bridge.json`.

## Changing a theorem

`theorems.json` is a reviewed record, not a list of whatever compiles: a changed statement,
definition, explanation or mapped input fails the audit until the change is reviewed and its
digests are recorded again. The rules for adding or changing a theorem are in the Git
repository's `dev/proofs/README.md`.
