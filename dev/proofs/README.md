# Lean proofs: rules for contributors

`proofs/` is a Lean 4 project that states parts of the library's mathematical model and proves
properties of it; Lean checks every step of each proof. This page gives the rules a theorem
meets before it is accepted, and how the checks and the release use the proofs. The commands
are in [proofs/README.md](../../proofs/README.md); the explained theorems are under
[docs/formal/theorems/](../../docs/formal/theorems/round-robin.md). Ordinary native and Python
development does not need Lean, and `make check` does not run the proofs.

## What a proof shows

[docs/semantics.md](../../docs/semantics.md) is the statement of what the library computes. A
Lean definition expresses a selected part of it, and each theorem states its assumptions
exactly. When the contract, the model and the code disagree, review the three together; a proof
never changes a default, an order, a supported domain or a result on its own.

Three different claims need three different kinds of evidence:

| Claim | Evidence it needs |
| --- | --- |
| A property holds for the mathematical model | A reviewed statement and a Lean proof checked by Lean's kernel, with audited dependencies |
| A solver's result for one input meets its objective | A certificate for that exact input, accepted by a checker that is itself proved sound; the compiled checker and its input decoding stay trusted assumptions |
| The C++ code implements the model on its whole supported domain | A refinement proof connecting the code and its machine semantics to the model |

The project gives the first kind only. Comparisons of the model with the native package (the
bridge, below) check the tested inputs exactly; they do not give the third claim. A proof about
the model says nothing about native exception safety, memory safety, locking or caching. Never
describe the package as "formally verified": name the model property that is proved.

## Files

```text
proofs/README.md                     commands, for a checkout and for the source archive
proofs/lean-toolchain                exact Lean release (mathlib's own pin, copied byte for byte)
proofs/lakefile.toml                 Lake package, mathlib at an exact Git revision
proofs/lake-manifest.json            every transitive dependency revision
proofs/NetGraphCore.lean             imports every model and proof module
proofs/NetGraphCore/Apportionment/   the round-robin model and its theorems
proofs/Checks.lean                   imports every check module
proofs/Checks/                       axiom audit, examples, Python source-safety checker (audit.py)
proofs/cases.json                    example inputs and expected outputs, derived independently
docs/formal/theorems/                plain explanation of each theorem and its proof
dev/proofs/run.sh                    entry point: check, replay, test, bridge
dev/proofs/test_tool.py              negative controls of the gate
dev/proofs/bridge.py                 exact comparisons with an installed native wheel
dev/proofs/replay_report.py          evidence file of the kernel replay
dev/proofs/package.py                source archive and wheel contents, published-file identity
.github/workflows/proofs.yml         CI jobs: Lean proofs, Lean distribution, Lean bridge
```

A new module directory under `proofs/NetGraphCore/` is created with its first definitions. Scratch
arguments, review notes, timing logs and abandoned models stay out of the project.

## What ships where

| Distribution | Proof content |
| --- | --- |
| Git checkout and version tag | Everything above |
| Source archive (sdist) | `proofs/` without `proofs/.lake/`, and `docs/`, including the explanations; no `dev/` |
| Wheel | No proof sources, Lean, mathlib or proof build products |
| `pip install`, CMake, run time | Lean is never installed or run; building the library does not build the proofs |
| CI evidence | Reports are workflow artifacts of one commit, never shipped |

The source archive's file rule is `sdist.include` and `sdist.exclude` in `pyproject.toml`;
`dev/stubs/check_contents.py` repeats it independently, and `make stubs-check` and
`dev/proofs/package.py` compare the built archive with it. Because the archive's file list is
explicit, the `/proofs/.lake/` line in `.gitignore` alone does not keep Lake's dependencies and
build products out: the sdist exclusion does. Reports go to `build/proofs/`, which neither git
nor the archive takes.

The archive holds every local file its proof checks read: the Lean sources, the checker and
`proofs/cases.json`. It needs nothing from `dev/`, `.git` or a sibling checkout. Lean and the
dependencies selected by Lake are downloaded.

## Adding or changing a theorem

1. Find the rule in `docs/semantics.md`, who calls it and the intended claim. List the exact
   assumptions, the degenerate cases (zero, empty, all-zero weights) and known counterexamples.
   Settle an ambiguity in the contract before choosing a convenient formal statement.
2. Write the plain explanation of the statement and of the proof first (next section), and have
   the definitions and statement reviewed before the Lean proof. Check the expanded types and
   assumptions, that the intended domain is not empty, and member identity and exact integer
   rules. Use concrete cases that satisfy the assumptions and cases that break a necessary one.
3. Write the simple model, then the supporting lemmas, then the theorem. Keep the definitions
   that state what is computed separate from any faster evaluator and from the proof itself.
   No unfinished proof enters the project.
4. Import the theorem's module from its aggregation root and add executable examples with
   expected values derived by hand in `proofs/cases.json`, never generated from Lean or from the
   native code. A comparison with the native code gets a plausible wrong variant that the
   comparison rejects. The audit discovers project theorems from the built Lean modules.
5. Run `bash dev/proofs/run.sh check`, `replay` and `test`, and `bridge` for a connected theorem.
   A change to native code still runs its own checks under `AGENTS.md`; a proof replaces none.
6. Get an independent review of the statement, definitions, domain and evidence, by a reviewer
   able to check the mathematical argument, and the readability review described below. The
   review keeps what is proved about the model apart from evidence about the implementation.
7. Land the definitions, proofs, explanations, cases and any contract change together.

A theorem is accepted when it has all of the following:

- a statement reviewed against `docs/semantics.md`, with its domain stated in the explanation;
- a Lean proof that builds with warnings as errors and depends only on the standard axioms;
- a plain explanation of the statement and the proof, with a worked example and a boundary
  example, and a failing case for every assumption that can fail;
- doc comments on every project declaration, private helpers included;
- the gate passing (`check`, `replay`, `test`, and `bridge` for native comparisons).

A proof checks the statement actually written in Lean. Review changes to the contract, statement,
assumptions and native implementation together: a true theorem can describe the wrong rule, and
finite native comparisons do not establish equivalence on every supported input.

## Plain explanations

Every theorem and its proof, supporting lemmas included, is explained so that an engineer who
knows the library, has no mathematical training and does not read Lean can follow it. This is a
condition of acceptance, not documentation added afterwards. Each theorem has its own section in
a page under `docs/formal/theorems/`. [round-robin.md](../../docs/formal/theorems/round-robin.md) shows
the expected form. A section explains:

- the practical question and the result in ordinary words: what is split, carried, bounded or
  kept, and why it matters to a caller;
- every symbol and term at first use, with units, and each assumption as a concrete operating
  condition, zero and empty cases included;
- a small worked example with its actual inputs, intermediate steps and result, checked by an
  executable case;
- the proof as a sequence of steps: what is known, why the next conclusion follows and which
  earlier fact it uses; for an induction, a bound or a transformation, the reasoning itself, not
  only the name of the technique;
- each supporting lemma's meaning and why it is needed, and each mathlib fact used, in plain
  words where it is used;
- the formal declaration and the library behaviour it describes, the conditions under which it
  applies, and a concrete failing case for each necessary restriction.

Examples answer a concrete question. Each theorem has a representative example and a boundary
example; where an assumption can fail, a contrasting example shows the failure. Choose small
inputs that expose the rule (unequal weights, an explicit member order where ties matter, exact
units), explain every input, step and output so a reader can redo the calculation by hand, derive
the expected values before running any code, and compare outputs exactly. A formula, a code block
or "it is easy to see" without a walkthrough does not count.

In Lean, every definition and theorem has a doc comment saying what it means in ordinary words,
and each substantive proof block has a comment naming its reasoning step. Explanations and
comments follow the actual proof; they never state a different rule or hide an assumption behind
an analogy.

Executable examples are checked against `proofs/cases.json`. A human reviewer checks the text,
doc comments and correspondence with the model: an independent reader restates the result, its
conditions and why it holds without
reading Lean, recalculates the examples and says what each one shows. Missing, misleading or
unexplained examples block acceptance even when Lean accepts the proof.

## Modelling rules

- Prove properties over unbounded integers and, separately, that the machine values stay in
  range: amounts, cumulative counters, products, intermediate sums, sentinels and refusal limits.
  A proof over natural numbers does not show that a C++ expression cannot overflow. Rounding of
  floating-point inputs and outputs is not covered unless a theorem states it.
- A graph model keeps edge identity: parallel edges, self-loops, the EdgeId order, and the
  difference between an edge of capacity 0 and an excluded edge.
- Do not add theorems that hold only because of the choice of a type.

## Proof checks

`bash dev/proofs/run.sh check` builds every proof with warnings as errors, then runs the axiom
and example audit. Local modules are discovered from the Lean source files, and each must be
imported by its aggregation root. The audit loads those built modules and checks the actual
transitive axioms of every mathematical declaration, including helpers and generated declarations.
It refuses project axioms and dependencies on any axiom except `propext`, `Classical.choice` and
`Quot.sound`. No hand-maintained list selects which theorems receive these checks.

The source preflight refuses `sorry`, `admit`, project axioms, `native_decide`, unsafe code,
kernel-skip options and commands that add declarations unchecked in the mathematical modules
and executable examples. Audit metaprogramming lives only in `proofs/Checks/`. The source checks
supplement Lean's build and replay; they do not establish kernel validity on their own.

The executable model cases must equal the independently derived inputs and exact outputs in
`proofs/cases.json`. The generated report names the discovered theorems and their transitive
axioms and contains the model's actual examples. The command refuses a missing or incomplete
report, so disabling the audit cannot make the check pass.

`bash dev/proofs/run.sh replay` also replays every declaration in a fresh environment with
`leanchecker --fresh`, which reuses Lean's own kernel (it is not an independent kernel); the
release requires it. Set `PROOFS_DEPENDENCY_ORIGIN` truthfully (`cache+local`, `source` or `mixed`).

`bash dev/proofs/run.sh test` runs negative controls in a disposable project: forbidden source
commands, an imported axiom, an imported unfinished proof and a disabled audit must fail. A false
theorem inserted with kernel checking skipped inside a checker module must fail the fresh replay.
The controls also exercise exact bridge comparisons and the publication artifact identity check.
`bash dev/proofs/run.sh test --preflight-only` runs the source, comparison and artifact controls
without Lean; it does not run the build, transitive-axiom or replay controls.

Review establishes that the proved statements express the intended rules. The bridge compares
the model with the selected native wheel on the stated finite cases; these comparisons establish
agreement on those inputs only.

## Pins

The project pins one mathlib release by Git revision, Lean at the version that revision names,
and every transitive revision in `lake-manifest.json`. Routine checks never run `lake update` or
follow a branch; a mathlib cache download is a build aid, must match the pins, and is not proof
evidence. Upgrading a pin is a reviewed change that rebuilds every proof and reruns the audit,
the replay and the controls.

## CI and release

`.github/workflows/proofs.yml` ("Lean verification") runs on every pull request and every branch
push, with no path filters, and can be started by hand. Its jobs:

- **Lean proofs**: the check, the replay and the negative controls, on Ubuntu;
- **Lean distribution**: builds the sdist and wheel with planted files in `proofs/.lake/` and
  `build/proofs/`, checks their exact contents, extracts the archive into a directory without
  `.git` or `dev/`, builds and audits the proofs there and replays them, and records the
  archive's and wheels' SHA-256;
- **Lean bridge**: installs that wheel and compares it exactly with the model cases.

A pull request merges only when these three pass. For a version tag, `release.yml` calls the same
workflow on the release's own sdist and wheels, and publication waits for it; the publish job
checks again that the files it uploads are the ones whose proofs were checked. Push the release
branch first and let the workflow pass there before pushing the tag
([CONTRIBUTING.md](../../CONTRIBUTING.md), "Release Process").

## References

- Lean axioms and incomplete proofs: <https://lean-lang.org/doc/reference/latest/Axioms/>
- Lake configuration and manifest: <https://lean-lang.org/doc/reference/latest/Build-Tools-and-Distribution/Lake/>
- mathlib projects and cache: <https://leanprover-community.github.io/install/project>
- Kernel replay checker: <https://github.com/leanprover/lean4checker>
- Reusable GitHub workflows: <https://docs.github.com/en/actions/how-tos/reuse-automations/reuse-workflows>
