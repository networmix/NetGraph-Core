"""All-root exclusive sample weights; unmatched inline phases remain unresolved."""

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import re
import sys

from run_workload import ROOT, equal
from capture import REF

sys.path.insert(0, str(ROOT / "scratch/orient-add/scripts"))
from original_analyze_native_sample import parse  # noqa: E402
from native_attribution import raw_pc  # noqa: E402 -- parser only; no old ranges

RANGES = [
    (0xC4EF0, 0xC51CC, "solver_search"),
    (0xC51CC, 0xC52E4, "solver_potentials"),
    (0xC5338, 0xC5C44, "solver_blocking_flow"),
    (0xC5CC8, 0xC6140, "solver_search"),
]


def verify_ranges():
    path = ROOT / "scratch/clocks/reference-arm64.txt"
    instructions = {
        int(m[1], 16): m[2]
        for line in path.read_text().splitlines()
        if (m := re.match(r"([0-9a-f]{16})\s+(.*)", line))
    }
    witnesses = {
        0xC2D68: "bl\t0xc4e24",
        0xC4F08: "add\tx8, x8, #0x1",
        0xC4F28: "str\tx22, [x8, x10, lsl #3]",
        0xC51F4: "add\tx12, x14, x12",
        0xC5200: "strb\tw15, [x12, x11]",
        0xC5368: "str\tw8, [x19, #0x1528]",
        0xC5494: "cmp\tx10, x21",
        0xC5588: "add\tx8, x8, #0x1",
        0xC5870: "ldrsw\tx10, [x8, #-0x4]!",
        0xC5BDC: "sub\tx14, x14, x9",
        0xC5BF0: "str\tx14, [x11, x12, lsl #3]",
        0xC5DF0: "str\tx10, [x13]",
    }
    for pc, instruction in witnesses.items():
        assert instruction in instructions[pc], (hex(pc), instructions[pc])
    return {hex(pc): instructions[pc] for pc in witnesses}


WAIT = (
    "__psynch_",
    "__ulock_wait",
    "semaphore_wait",
    "mach_msg",
    "pthread_join",
    "__workq_kernreturn",
    "__sigsuspend",
    "__select",
    "kevent",
    "nanosleep",
    "__semwait_signal",
)


def classify(frames):
    native = any("(in _netgraph_core." in f for f in frames)
    leaf = frames[-1].split("  (in ", 1)[0]
    if any(token in leaf for token in WAIT):
        return "blocked_wait", native
    if not native:
        return "python_runtime_or_unresolved_outside_core", False
    for frame in reversed(frames):
        pc = raw_pc(frame)
        if pc is not None:
            for left, right, category in RANGES:
                if left <= pc < right:
                    return category, True
    mapping = [
        ("solver_search", ("::dijkstra(", "::dijkstra<")),
        ("solver_blocking_flow", ("::build_levels<", "::find_path<", "::drain<")),
        ("solver_closure_cut", ("::close_sources(", "::cut_edges(")),
        (
            "result_record_program_construction",
            (
                "::derived_entry(",
                "::stored_entry(",
                "LedgerAccess::commit_batch(",
                "::decompose(",
                "::optimized_flow_program(",
                "ProgramData::initialize_storage(",
                "PredecessorDAGBuilder::of_edges(",
                "PlacementResult::",
                "FlowRecord::",
            ),
        ),
        (
            "roster_support_search",
            (
                "::static_costs(",
                "::static_support(",
                "::tunnel_roster(",
                "::run_search(",
            ),
        ),
        (
            "ledger_policy_fixed_cost",
            (
                "FlowLedger::FlowLedger(",
                "FlowPolicy::FlowPolicy(",
                "PolicyAccess::prepare(",
                "PolicyAccess::construct(",
                "FlowPolicy::~FlowPolicy(",
                "FlowLedger::~FlowLedger(",
                "LedgerAccess::detach(",
                "FlowPolicy::remove(",
            ),
        ),
    ]
    for category, tokens in mapping:
        if any(token in f for f in frames for token in tokens):
            return category, True
    if any("::solve_min_cost_flow(" in f for f in frames):
        return "solver_inline_phase_unresolved", True
    return "unresolved_native_bindings", True


def analyze(source, expected, out):
    index = json.loads((source / "index.json").read_text())
    assert (
        index["native_sha256"] == REF
        and index["sample_exit"] == index["worker_exit"] == 0
    )
    witnesses = verify_ranges()
    equal(
        json.loads((source / "results.json").read_text()),
        json.loads((expected / "results.json").read_text()),
    )
    text = (source / "sample.txt").read_text()
    nodes = parse(text)
    roots = [n for n in nodes if n.parent is None]
    totals, byroot, leaves, functions = Counter(), {}, [], Counter()
    for i, node in enumerate(nodes):
        own = node.count - sum(nodes[c].count for c in node.children)
        if not own:
            continue
        chain, cursor = [], i
        while cursor is not None:
            chain.append(nodes[cursor])
            cursor = nodes[cursor].parent
        chain.reverse()
        frames = [n.frame for n in chain]
        category, native = classify(frames)
        root = str(chain[0].line)
        byroot.setdefault(root, Counter())[category] += own
        totals[category] += own
        functions[node.frame.split("  (in ", 1)[0]] += own
        leaves.append(
            {
                "line": node.line,
                "root": root,
                "self_samples": own,
                "category": category,
                "native_context": native,
                "frames": frames,
            }
        )
    for root in roots:
        assert sum(byroot[str(root.line)].values()) == root.count
    assert sum(totals.values()) == sum(r.count for r in roots)
    top = text.split("Sort by top of stack, same collapsed", 1)[1].split(
        "Binary Images:", 1
    )[0]
    crosschecks = {}
    for line in top.splitlines():
        match = re.fullmatch(r"\s*(netgraph::.+?)  \(in .+?\)\s+(\d+)\s*", line)
        if match:
            assert functions[match[1]] == int(match[2])
            crosschecks[match[1]] = int(match[2])
    report = {
        "exploratory": True,
        "native_sha256": REF,
        "sample_sha256": hashlib.sha256(text.encode()).hexdigest(),
        "roots": [
            {
                "line": r.line,
                "frame": r.frame,
                "samples": r.count,
                "categories": byroot[str(r.line)],
            }
            for r in roots
        ],
        "categories": totals,
        "reference_instruction_witnesses": witnesses,
        "phase_ranges": RANGES,
        "top_checks": crosschecks,
        "unresolved_note": "Four partial solver ranges have same-reference instruction/source witnesses. Other inline solver, binding and stripped worker frames remain unresolved. Counts are sampled stacks, not CPU time.",
        "summary": "PASS: all-root exclusive weights conserved, named top counts and full export equal",
    }
    out.mkdir()
    (out / "summary.json").write_text(json.dumps(report, indent=2) + "\n")
    with (out / "self-weights.jsonl").open("x") as stream:
        for row in leaves:
            stream.write(json.dumps(row) + "\n")
        stream.write(
            json.dumps(
                {
                    "summary": "each self weight attributed once",
                    "samples": sum(totals.values()),
                }
            )
            + "\n"
        )
    print(report["summary"])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    parser.add_argument("expected", type=Path)
    parser.add_argument("out", type=Path)
    args = parser.parse_args()
    analyze(args.source, args.expected, args.out)
