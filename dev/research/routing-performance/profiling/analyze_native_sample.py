"""Partition one macOS sample call tree without adding inclusive counts.

At each node, self weight is inclusive weight minus immediate children.
Every self weight is assigned once by its enclosing measured native branch.
Unknown symbols outside a named branch remain unresolved native/binding work.
"""

import argparse
import hashlib
import json
import re
import subprocess
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path

EFFORT = Path(__file__).resolve().parents[3]
SCRATCH = EFFORT / "scratch/orient-consumers"
CORE = EFFORT / "worktrees/core-orient"


@dataclass
class Node:
    line: int
    depth: int
    count: int
    frame: str
    parent: int | None
    children: list[int] = field(default_factory=list)


def parse(text):
    nodes, stack = [], []
    active = False
    for line_number, line in enumerate(text.splitlines(), 1):
        if line == "Call graph:":
            active = True
            continue
        if line.startswith("Total number in stack"):
            break
        if not active or not line.strip():
            continue
        match = re.fullmatch(r"([ +!:|]*)(\d+) (.+)", line)
        assert match is not None, (line_number, line)
        depth = len(match[1])
        while stack and nodes[stack[-1]].depth >= depth:
            stack.pop()
        parent = stack[-1] if stack else None
        node = Node(line_number, depth, int(match[2]), match[3], parent)
        if parent is not None:
            assert depth == nodes[parent].depth + 2, (line_number, depth)
            nodes[parent].children.append(len(nodes))
        nodes.append(node)
        stack.append(len(nodes) - 1)
    assert nodes, "empty call graph"
    for node in nodes:
        assert sum(nodes[c].count for c in node.children) <= node.count, node
    return nodes


def name(frame):
    return frame.split("  (in ", 1)[0]


def classify(chain):
    frames = [node.frame for node in chain]
    native = any("(in _netgraph_core." in frame for frame in frames)
    if not native:
        if name(frames[-1]) == "read" and any(
            "textiowrapper_readline" in frame.lower() for frame in frames
        ):
            return "readiness_wait", False
        return "python_and_runtime_outside_core", False
    # Pinned-wheel disassembly ties return PCs +140 and +180 to the two
    # calls in static_support: run_search and emit_full, respectively.
    for frame in frames:
        if "netgraph::core::detail::static_support(" in frame:
            if re.search(r"\) \+ 140\s", frame):
                return "shortest_path_search_branch", True
            if re.search(r"\) \+ 180\s", frame):
                return "support_and_route_program_construction", True
    if any(
        token in frame
        for frame in frames
        for token in (
            "netgraph::core::detail::prepare_shortest_path_routes(",
            "netgraph::core::detail::shortest_path_program(",
            "netgraph::core::detail::ProgramData::initialize_storage(",
            "netgraph::core::detail::make_forwarding_view(",
        )
    ):
        return "support_and_route_program_construction", True
    if any(
        token in frame
        for frame in frames
        for token in (
            "netgraph::core::detail::lossless_step(",
            "netgraph::core::detail::allocate(",
        )
    ):
        return "placement_admission_forwarding", True
    if any(
        token in frame
        for frame in frames
        for token in (
            "netgraph::core::detail::derived_entry(",
            "netgraph::core::detail::LedgerAccess::commit_batch(",
        )
    ):
        return "native_ledger_result_materialization", True
    if any(
        token in frame
        for frame in frames
        for token in (
            "netgraph::core::FlowLedger::FlowLedger(",
            "netgraph::core::Eligibility::intersect(",
            "netgraph::core::Graph::eligibility(",
            "netgraph::core::FlowPolicy::remove(",
            "netgraph::core::detail::LedgerAccess::detach(",
            "netgraph::core::FlowPolicy::~FlowPolicy(",
            "netgraph::core::FlowLedger::~FlowLedger(",
        )
    ):
        return "native_setup_masks_cleanup", True
    return "unresolved_native_and_bindings", True


def normalize(value, parent=None):
    if isinstance(value, dict):
        return {
            key: None
            if key == "workers" and parent == "metadata"
            else normalize(item, key)
            for key, item in value.items()
        }
    if isinstance(value, list):
        return [normalize(item) for item in value]
    return value


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--capture", type=Path, default=SCRATCH / "native-outside-01/1")
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    sample_path = args.capture / "sample.txt"
    text = sample_path.read_text()
    nodes = parse(text)
    roots = [node for node in nodes if node.parent is None]
    assert len(roots) == 1 and "com.apple.main-thread" in roots[0].frame
    total = roots[0].count
    grouped, native_total, leaf_rows, function_self = Counter(), 0, [], Counter()
    for i, node in enumerate(nodes):
        own = node.count - sum(nodes[c].count for c in node.children)
        if not own:
            continue
        chain = []
        current = i
        while current is not None:
            chain.append(nodes[current])
            current = nodes[current].parent
        chain.reverse()
        category, native = classify(chain)
        grouped[category] += own
        native_total += own if native else 0
        function_self[name(node.frame)] += own
        leaf_rows.append(
            {
                "sample_line": node.line,
                "self_samples": own,
                "category": category,
                "native_context": native,
                "stack_lines": [n.line for n in chain],
                "leaf": node.frame,
            }
        )
    assert sum(grouped.values()) == total
    # Independent supplied 'top of stack' table verifies tree depth/self math.
    top_section = text.split("Sort by top of stack, same collapsed", 1)[1].split(
        "Binary Images:", 1
    )[0]
    top_checks = {}
    for prefix in (
        "netgraph::core::detail::allocate(",
        "netgraph::core::detail::prepare_shortest_path_routes(",
        "netgraph::core::detail::LedgerAccess::commit_batch(",
        "netgraph::core::detail::derived_entry(",
    ):
        row = next(
            line for line in top_section.splitlines() if line.strip().startswith(prefix)
        )
        expected = int(re.search(r"(\d+)\s*$", row)[1])
        found = sum(
            count
            for function, count in function_self.items()
            if function.startswith(prefix)
        )
        assert found == expected, (prefix, expected, found)
        top_checks[prefix] = found
    identity = json.loads((SCRATCH / "source-inventory.json").read_text())[
        "identities"
    ]["runtime"]
    extension = Path(identity["extension_path"])
    extension_hash = hashlib.sha256(extension.read_bytes()).hexdigest()
    assert extension_hash == identity["extension_sha256"]
    disassembly = subprocess.run(
        ["/usr/bin/otool", "-arch", "arm64", "-tvV", str(extension)],
        capture_output=True,
        text=True,
        check=True,
        timeout=30,
    ).stdout
    instructions = {}
    for line in disassembly.splitlines():
        match = re.match(r"^([0-9a-f]{16})\s", line)
        if match:
            address = int(match[1], 16)
            if 0x1CF1C <= address < 0x1D074:
                instructions[address] = line
    assert "bl\t0x19494" in instructions[0x1CFA4]
    assert "bl\t0x1d2ac" in instructions[0x1CFCC]
    source = (CORE / "src/search.cpp").read_text().splitlines()
    static_source = "\n".join(source[1173:1190])
    assert static_source.index("run_search(ws, in)") < static_source.index(
        "emit_full(ws, em)"
    )
    expected = normalize(json.loads((SCRATCH / "dcbb96-1-02/results.json").read_text()))
    actual_raw = json.loads((args.capture / "results.json").read_text())
    actual = normalize(actual_raw)
    assert json.dumps(actual) == json.dumps(expected), "sampled export differs"
    assert actual_raw["metadata"]["core_extension_sha256"] == extension_hash
    capture_index = json.loads((args.capture / "index.json").read_text())
    assert capture_index["sample_exit"] == capture_index["worker_exit"] == 0
    rows = [
        {
            "category": category,
            "samples": count,
            "share_all_percent": count / total * 100,
            "share_native_context_percent": count / native_total * 100
            if category not in ("readiness_wait", "python_and_runtime_outside_core")
            else None,
        }
        for category, count in grouped.items()
    ]
    summary = {
        "measured_exploratory": True,
        "captures": 1,
        "requested_duration_s": 15,
        "workers": 1,
        "thread_roots": len(roots),
        "all_main_thread_samples": total,
        "native_context_samples": native_total,
        "categories": rows,
        "source_sample": str(sample_path),
        "sample_sha256": hashlib.sha256(sample_path.read_bytes()).hexdigest(),
        "extension_sha256": extension_hash,
        "top_of_stack_crosschecks": top_checks,
        "worker": json.loads((args.capture / "worker.json").read_text()),
        "interpretation": "exclusive self weights attributed by enclosing branches; percentages are stack samples, not end-to-end seconds; static_support+140/+180 names inferred from matching ARM64 call sites and source",
        "export_comparison": "all values/types/order equal except metadata.workers",
        "summary": "sample tree conserved, top-of-stack checks passed, pinned export matched, one exploratory capture partitioned",
    }
    args.out.mkdir()
    with (args.out / "summary.json").open("x") as stream:
        json.dump(summary, stream, indent=2)
    with (args.out / "self-weights.jsonl").open("x") as stream:
        for row in leaf_rows:
            stream.write(json.dumps(row) + "\n")
        stream.write(
            json.dumps({"summary": "every self weight assigned once", "samples": total})
            + "\n"
        )
    with (args.out / "static-support-arm64.txt").open("x") as stream:
        stream.write(
            "\n".join(instructions.values())
            + "\nSUMMARY matching ARM64 static_support call sites verified\n"
        )
    print(
        json.dumps(
            {
                "summary": summary["summary"],
                "all_samples": total,
                "native_samples": native_total,
                "categories": rows,
            }
        )
    )


if __name__ == "__main__":
    main()
