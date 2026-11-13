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

sys.path.insert(0, str(ROOT / "tasks/orient-consumers/probes"))
from analyze_native_sample import classify as orient_classify  # noqa: E402

RANGES = []  # No absolute PCs from a different wheel are used.


def verify_reference(expected_sha=REF):
    assert expected_sha == REF, "wrong binary for static_support offsets"
    extension = next(
        (ROOT / "scratch/consumer/ref/venv/lib/python3.14/site-packages").glob(
            "_netgraph_core*.so"
        )
    )
    assert hashlib.sha256(extension.read_bytes()).hexdigest() == REF
    import subprocess

    text = subprocess.run(
        ["/usr/bin/otool", "-arch", "arm64", "-tvV", str(extension)],
        check=True,
        capture_output=True,
        text=True,
        timeout=30,
    ).stdout
    rows = text.splitlines()
    start = next(
        i
        for i, row in enumerate(rows)
        if row.startswith("__ZN8netgraph4core6detail14static_support")
    )
    excerpt = rows[start : start + 49]
    assert "00000000000d38fc\tbl\t0x142b14" in excerpt
    assert "00000000000d3924\tbl\t0x14414c" in excerpt
    source = (ROOT / "worktrees/core-ref/src/search.cpp").read_text()
    body = source[
        source.index("PredecessorDAG static_support(") : source.index(
            "std::vector<Cost> static_costs("
        )
    ]
    assert body.index("run_search(ws, in)") < body.index("emit_full(ws, em)")
    return {
        "extension_sha256": REF,
        "disassembly_sha256": hashlib.sha256(text.encode()).hexdigest(),
        "static_support": excerpt,
        "source_body": body,
        "return_offsets": [140, 180],
        "qualification": "Targets stripped: branch names inferred from source order and call-site argument/result context, as in ORIENT; no old absolute PC range reused.",
    }


verify_ranges = verify_reference


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
    from types import SimpleNamespace

    native = any("(in _netgraph_core." in f for f in frames)
    leaf = frames[-1].split("  (in ", 1)[0]
    if any(token in leaf for token in WAIT):
        return "blocked_wait", native
    return orient_classify([SimpleNamespace(frame=f) for f in frames])


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
        "core_samples_including_wait": sum(
            r["self_samples"] for r in leaves if r["native_context"]
        ),
        "core_samples_excluding_wait": sum(
            r["self_samples"]
            for r in leaves
            if r["native_context"] and r["category"] != "blocked_wait"
        ),
        "blocked_in_core": sum(
            r["self_samples"]
            for r in leaves
            if r["native_context"] and r["category"] == "blocked_wait"
        ),
        "blocked_outside_core": sum(
            r["self_samples"]
            for r in leaves
            if not r["native_context"] and r["category"] == "blocked_wait"
        ),
        "top_checks": crosschecks,
        "unresolved_note": "P1 ORIENT categories; static_support return offsets revalidated on the reference wheel. Unknown branches stay unresolved; explicit blocked leaves excluded from phase denominator. Counts are sampled stacks, not CPU time.",
        "summary": "PASS: all-root exclusive weights conserved, named top counts and full export equal",
    }
    denominator = report["core_samples_excluding_wait"]
    assert denominator > 0
    report["native_phase_share_percent"] = {
        category: 100
        * sum(
            r["self_samples"]
            for r in leaves
            if r["native_context"] and r["category"] == category
        )
        / denominator
        for category in (
            "support_and_route_program_construction",
            "placement_admission_forwarding",
            "native_ledger_result_materialization",
            "shortest_path_search_branch",
            "unresolved_native_and_bindings",
            "native_setup_masks_cleanup",
        )
    }
    assert abs(sum(report["native_phase_share_percent"].values()) - 100) < 1e-8
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
