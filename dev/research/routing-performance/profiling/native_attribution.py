"""Wheel-specific conservative native attribution; verified ARM64 instruction witnesses.

Ranges are image-relative PCs of the pinned wheel, not symbols of a different
harness. The source and compiler control flow identify inlined solver phases.
Unknown work stays unknown. Only complete evaluate_pairs binding subtrees are
used to scale the separately timed per-evaluation native API duration.
"""

import hashlib
import re
from pathlib import Path

BASE = Path(__file__).resolve().parents[1]
EFFORT = BASE.parents[1]
WHEEL_SHA256 = "6db37900d83800b0b5c227104889f47259c5761dd640c1d0d8684fac61009ff1"
# Half-open PC ranges. Calls are attributed using the enclosing return PC as in
# the predecessor; detailed descendants take priority over enclosing phases.
PHASE_RANGES = [
    (0xA3D8, 0xA734, "solver_residual_dijkstra"),
    (0xA734, 0xA7A8, "solver_potentials"),
    (0xA7A8, 0xAC0C, "solver_residual_dijkstra"),
    (0xAC0C, 0xACAC, "solver_potentials"),
    (0xACE4, 0xB66C, "solver_blocking_flow_dinic"),
    (0xB830, 0xB854, "solver_residual_closure_cut"),
    (0xB880, 0xBAD0, "solver_potentials"),
    (0xE1EC, 0xE548, "solver_residual_closure_cut"),
    (0xE548, 0xE6F8, "solver_residual_closure_cut"),
]


def raw_pc(frame):
    match = re.search(
        r"\(in _netgraph_core\.[^)]*\).*?load address 0x[0-9a-f]+ \+ (0x[0-9a-f]+)",
        frame,
    )
    return int(match[1], 16) if match else None


def batch_context(chain):
    for node in chain:
        if "detail::evaluate_pairs_in_call(" in node.frame:
            return True
        pc = raw_pc(node.frame)
        if pc is not None and (0x1331E4 <= pc < 0x133344 or 0x134850 <= pc < 0x135D98):
            return True
    return False


def apply_actions_context(chain):
    return any(
        "netgraph::core::apply_actions(" in node.frame
        or (
            (pc := raw_pc(node.frame)) is not None
            and (0x12CAAC <= pc < 0x12CBE4 or 0x12CCF4 <= pc < 0x12E00C)
        )
        for node in chain
    )


def classify(chain):
    frames = [node.frame for node in chain]
    native = any("(in _netgraph_core." in f for f in frames)
    if not native:
        if (
            any("textiowrapper_readline" in f.lower() for f in frames)
            and frames[-1].split("  (in ", 1)[0] == "read"
        ):
            return "readiness_wait", False
        return "python_runtime_outside_core", False
    solver = any("detail::solve_min_cost_flow(" in f for f in frames)
    if solver:
        for node in reversed(chain):
            pc = raw_pc(node.frame)
            if pc is not None:
                for lo, hi, category in PHASE_RANGES:
                    if lo <= pc < hi:
                        return category, True
        return "solver_setup_output_other", True
    for category, tokens in [
        (
            "roster_support_search",
            (
                "detail::static_costs(",
                "detail::static_support(",
                "detail::tunnel_roster(",
                "detail::run_search(",
            ),
        ),
        (
            "result_record_program_construction",
            (
                "detail::derived_entry(",
                "detail::stored_entry(",
                "LedgerAccess::commit_batch(",
                "detail::decompose(",
                "detail::optimized_flow_program(",
                "detail::ProgramData::initialize_storage(",
                "detail::PredecessorDAGBuilder::of_edges(",
                "detail::min_route_cost(",
                "detail::compare_programs(",
                "PlacementResult::",
                "FlowRecord::",
            ),
        ),
        (
            "ledger_policy_construction_cleanup",
            (
                "FlowLedger::FlowLedger(",
                "FlowPolicy::FlowPolicy(",
                "PolicyAccess::prepare(",
                "PolicyAccess::construct(",
                "FlowPolicy::~FlowPolicy(",
                "FlowLedger::~FlowLedger(",
                "LedgerAccess::detach(",
                "FlowPolicy::remove(",
                "Eligibility::intersect(",
                "Graph::eligibility(",
            ),
        ),
    ]:
        if any(token in frame for frame in frames for token in tokens):
            return category, True
    # The known batch wrapper includes argument conversion/validation, Python
    # result slots, publication and cleanup. Do not sweep unknown C++ evaluation
    # work into bindings merely because its enclosing wrapper is on the stack.
    if (
        batch_context(chain)
        and not any("detail::evaluate_pairs_in_call(" in f for f in frames)
    ) or (
        apply_actions_context(chain)
        and not any("netgraph::core::apply_actions(" in f for f in frames)
    ):
        return "binding_conversion_publication", True
    return "unresolved_native_bindings", True


def verify_witnesses():
    import _netgraph_core as extension

    assert (
        hashlib.sha256(Path(extension.__file__).read_bytes()).hexdigest()
        == WHEEL_SHA256
    )
    lines = (BASE / "wheel-arm64.txt").read_text().splitlines()
    instructions = {
        int(m[1], 16): (number, line)
        for number, line in enumerate(lines, 1)
        if (m := re.match(r"([0-9a-f]{16})\s", line))
    }
    # The call chain enters build(), then run(). The latter inlines Dijkstra,
    # settle(), drain(), edge-flow output and published potentials.
    witnesses = {
        0x12CB04: "bl\t0x12ccf4",
        0x12D4B8: "bl\t__ZN8netgraph4core13apply_actions",
        0x7B18: "bl\t0x9a38",
        0x7B54: "bl\t0xa31c",
        0xA3F8: "add\tx8, x8, #0x1",  # dijkstra_runs
        0xA41C: "str\tx11, [x8, x10, lsl #3]",  # infinite distances
        0xA648: "str\txzr, [x8, x19, lsl #3]",  # source distance=0
        0xAB34: "sub\tx19, x8, x9",  # reduced distance: du+pi[u]+cost-pi[v]
        0xAB48: "str\tx19, [x8, w27, sxtw #3]",  # relax distance
        0xA75C: "add\tx12, x14, x12",  # potentials[v]+=dist[v]
        0xA768: "strb\tw15, [x12, x11]",  # ever_settled[v]=1
        0xAD30: "str\tw8, [x27, #0x1528]",  # Dinic reach generation
        0xAED8: "cmp\tx10, x21",  # zero reduced cost in reverse BFS
        0xAFA0: "bl\t0x142f44",  # level_.assign(-1), memset
        0xAFAC: "str\tx8, [x22, #0x90]",  # dinic_phases++
        0xB18C: "ldrsw\tx9, [x24, #-0x4]!",  # DFS current stack head
        0xB83C: "bl\t0xe1ec",
        0xB850: "bl\t0xe548",  # close_sources / cut_edges
        0xB8D0: "ldr\tx9, [x27, #0x1d8]",  # ever_settled mask -> output potentials
        0x13528C: "bl\t__ZN8netgraph4core6detail22evaluate_pairs_in_call",
        0x133248: "add\tx2, x2, #0x850",  # specialized thunk chooses lambda at 0x134850
    }
    for address, expected in witnesses.items():
        assert expected in instructions[address][1], (
            hex(address),
            expected,
            instructions[address],
        )
    source = (
        (EFFORT / "worktrees/core-orient/src/min_cost_flow.cpp")
        .read_text()
        .splitlines()
    )
    source_witnesses = {
        "dijkstra": [
            (i + 1, line)
            for i, line in enumerate(source)
            if "bool dijkstra(" in line or "const P nd = base +" in line
        ],
        "potentials": [
            (i + 1, line)
            for i, line in enumerate(source)
            if "potentials_[vi] +=" in line or "ever_settled_[vi] = 1" in line
        ],
        "dinic": [
            (i + 1, line)
            for i, line in enumerate(source)
            if "Amount drain(" in line
            or "bool build_levels(" in line
            or "bool find_path(" in line
        ],
        "run": [
            (i + 1, line)
            for i, line in enumerate(source)
            if "Output run(" in line
            or "dijkstra<Cost>(c)" in line
            or "total = drain<Cost>" in line
        ],
    }
    assert all(source_witnesses.values())
    return {
        "extension_sha256": WHEEL_SHA256,
        "instruction_witnesses": [
            {
                "pc": hex(pc),
                "disassembly_line": instructions[pc][0],
                "instruction": instructions[pc][1],
            }
            for pc in witnesses
        ],
        "source_witnesses": source_witnesses,
        "phase_ranges": [
            {"start": hex(lo), "end": hex(hi), "category": cat}
            for lo, hi, cat in PHASE_RANGES
        ],
        "method": "closest identified native ancestor; solver ranges inferred from the same wheel control flow and source; binding range verified by evaluate_pairs_in_call call site; unknowns remain unresolved",
    }
