"""Small separating controls for the new counters and sample/export attribution."""

import copy
import hashlib
import json
from pathlib import Path
import sys

import _netgraph_core as native
import netgraph_core as core

from analyze import delta
from analyze_sample import classify, parse, verify_reference
from reuse import ROOT, load
from run_workload import equal


def refused(operation):
    try:
        operation()
    except AssertionError:
        return
    raise AssertionError("poison accepted")


def main(out):
    checks = []
    extension = Path(native.__file__).resolve()
    assert extension.is_relative_to(Path(sys.prefix).resolve())
    assert (
        hashlib.sha256(extension.read_bytes()).hexdigest()
        == json.loads((ROOT / "scratch/split/build-identity.json").read_text())[
            "scratch_sha256"
        ]
    )
    before = native._split_snapshot()
    assert core.to_units(2.5, 1) == 3
    assert delta(before, native._split_snapshot()) == [0] * 6
    assert core.to_units([0.0, 1.2, 2.5], 1).tolist() == [0, 1, 3]
    vector = delta(before, native._split_snapshot())
    assert vector[4:] == [1, 0] and vector[0] > 0
    checks.append({"scalar_attached_vector_detached_cross_translation_unit": vector})
    before = native._split_snapshot()
    native._split_probe(2, False)
    multiple = delta(before, native._split_snapshot())
    assert multiple[4:] == [2, 0] and multiple[0] > 0
    before = native._split_snapshot()
    try:
        native._split_probe(1, True)
    except ValueError as error:
        assert str(error) == "split exception control"
    else:
        raise AssertionError("exception smoke did not raise")
    assert delta(before, native._split_snapshot())[4:] == [1, 0]
    checks.append(
        {"two_detaches": multiple, "exception": "PASS: RAII retake counted once"}
    )
    assert 12 - (3 + 4) - (1 + 1) == 3
    assert (90 - 80 - 2) == 8  # CPU retake is CPU, not wall latency.
    refused(lambda: delta([0] * 7, [8, 8, 1, 1, 0, 0, 0]))
    refused(lambda: delta([0] * 7, [0, 0, 0, 0, 0, 0, 1]))
    t4 = load("analyze_clocks")
    assert t4.union([(2, 7), (10, 14)]) == 9
    assert t4.complement([(2, 7), (10, 14)], 0, 20) == 11
    checks.append(
        "arithmetic, independent complement, missing-count and invalid-clock poisons rejected"
    )
    expected = {"metadata": {"core_extension_sha256": "ref"}, "rows": [1, 2, False]}
    equal(expected, copy.deepcopy(expected))
    for mutate in (
        lambda d: d["rows"].__setitem__(0, 3),
        lambda d: d["rows"].__setitem__(2, 0),
        lambda d: d.__delitem__("rows"),
        lambda d: d.__setitem__("extra", 1),
        lambda d: d["rows"].reverse(),
    ):
        poisoned = copy.deepcopy(expected)
        mutate(poisoned)
        refused(lambda: equal(expected, poisoned, identity=True))
    checks.append("export scalar/type/missing/extra/order poisons rejected")
    nodes = parse(
        "Call graph:\n    10 root1\n    + 6 child\n    7 root2\nTotal number in stack\n"
    )
    assert sum(n.count - sum(nodes[c].count for c in n.children) for n in nodes) == 17
    assert (
        classify(
            [
                "netgraph::core::detail::allocate( (in _netgraph_core.test)",
                "__psynch_cvwait",
            ]
        )[0]
        == "blocked_wait"
    )
    refused(lambda: verify_reference("wrong binary"))
    witness = verify_reference()
    for offset, category in (
        (140, "shortest_path_search_branch"),
        (180, "support_and_route_program_construction"),
    ):
        assert (
            classify(
                [
                    f"netgraph::core::detail::static_support(...) + {offset}  (in _netgraph_core.test)"
                ]
            )[0]
            == category
        )
    checks.append(
        {"sample_roots_wait_and_wrong_binary": "PASS", "reference_witness": witness}
    )
    record = {
        "checks": checks,
        "summary": "PASS: live zero/multiple/exception detaches, cross-TU counters, arithmetic/complement, export poisons, all-root wait classification and reference-only offsets",
    }
    with out.open("x") as stream:
        json.dump(record, stream, indent=2)
        stream.write("\n")
    print(record["summary"], flush=True)


if __name__ == "__main__":
    main(Path(sys.argv[1]))
