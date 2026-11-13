"""Short, unqualified controls for compilation, lifetimes, statistics and launchers."""

from array import array
import argparse
import copy
import json
from pathlib import Path
from types import SimpleNamespace
import sys
import time
from unittest.mock import patch

from kernel import Plan, interval_summary
import experiment as experiment
from events import BASE, Blobs, digest, read_stream
from replay import Replay
from values import first_difference


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("out", type=Path)
    args = parser.parse_args()
    out = args.out.resolve()
    out.mkdir()
    rows = []
    intervals = interval_summary(
        [
            array("q", [10**9, 4 * 10**9, 8 * 10**9, 9 * 10**9]),
            array("q", [3 * 10**9, 5 * 10**9]),
        ]
    )
    assert intervals == {"union_s": 5.0, "elapsed_s": 8.0}
    try:
        interval_summary([array("q", [2, 1])])
    except ValueError:
        pass
    else:
        raise AssertionError("backward clock accepted")
    rows.append("overlap union and backward-clock refusal")
    triple = [
        {"union_s": 10.0, "elapsed_s": 20.0, "repeats": 1, "rss_bytes": 20 * 1024**2}
        for _ in range(3)
    ]
    stat = experiment.statistics_for(
        [copy.deepcopy(triple) for _ in range(7)], "union_s"
    )
    assert (
        stat["median_ratio"] == 1 and stat["ci95"] == (1, 1) and stat["qualified_noise"]
    )
    shifted = copy.deepcopy(triple)
    shifted[1]["union_s"] = 9.0
    assert experiment.statistics_for([shifted] * 7, "union_s")["ci95"] == (0.9, 0.9)
    assert not experiment.statistics_for([triple], "union_s")["qualified_noise"]
    drift = copy.deepcopy(triple)
    drift[0]["union_s"] = 10.4
    assert not experiment.statistics_for([drift] * 7, "union_s")["qualified_noise"]
    try:
        experiment.equal_export({"x": True}, {"x": 1})
    except AssertionError:
        pass
    else:
        raise AssertionError("export bool/int poison accepted")
    rows.append(
        "hand-derived A/A, shifted ratio, short-sample/drift refusals, export type poison"
    )
    try:
        experiment.execute(
            [sys.executable, "-c", "import time; time.sleep(1)"], out / "timeout", 0.03
        )
    except ValueError:
        assert json.loads((out / "timeout/process.json").read_text())["timeout"]
    else:
        raise AssertionError("timeout not caught")
    rows.append("owned child group bounded and reaped")

    import netgraph_core as core

    recording = BASE / "scratch/replay/control-record-r12-01"
    reference = out / "reference"
    Replay(core, recording, reference, 2, "a", "0/1").run()
    plan = Plan(core, recording, 2)
    report = plan.block(2, dry=True, observation_root=out / "compiled")
    for repeat in range(2):
        candidate = out / "compiled" / str(repeat)
        for job in json.loads((reference / "manifest.json").read_text())["jobs"]:
            left, right = (
                list(read_stream(reference, job)),
                list(read_stream(candidate, job)),
            )
            assert len(left) == len(right)
            for expected, actual in zip(left, right, strict=True):
                difference = first_difference(
                    expected,
                    actual,
                    Blobs(reference / "payloads"),
                    Blobs(candidate / "payloads"),
                )
                assert difference is None, (repeat, job, difference)
    assert report["dry"] and report["repeats"] == 2
    rows.append(
        "two fresh compiled native sequences equal accepted check-mode fields/state/lifetimes"
    )

    import yaml

    source = BASE / "worktrees/netgraph-orient/scenarios/readme_example.yml"
    document = yaml.safe_load(source.read_text())
    for step in document["workflow"]:
        if step["type"] == "TrafficMatrixPlacement":
            step["iterations"], step["workers"] = 4, 1
    toy = out / "toy.yml"
    toy.write_text(yaml.safe_dump(document, sort_keys=False))
    document["workflow"] = []
    empty = out / "empty.yml"
    empty.write_text(yaml.safe_dump(document, sort_keys=False))
    worker = Path(__file__).with_name("worker.py")
    experiment.execute(
        [
            sys.executable,
            "-B",
            worker,
            "record",
            toy,
            out / "toy-record",
            "--workers",
            "1",
        ],
        out / "record-process",
    )
    grouped = out / "toy-record/recording"
    grouped_reference = out / "grouped-reference"
    Replay(core, grouped, grouped_reference, 2, "a", "0/1").run()
    grouped_report = Plan(core, grouped, 2).block(
        2, dry=True, observation_root=out / "grouped-compiled"
    )
    assert grouped_report["python_workers"] == 2
    for repeat in range(2):
        candidate = out / "grouped-compiled" / str(repeat)
        for job in json.loads((grouped_reference / "manifest.json").read_text())[
            "jobs"
        ]:
            left, right = (
                list(read_stream(grouped_reference, job)),
                list(read_stream(candidate, job)),
            )
            assert len(left) == len(right)
            for expected, actual in zip(left, right, strict=True):
                assert (
                    first_difference(
                        expected,
                        actual,
                        Blobs(grouped_reference / "payloads"),
                        Blobs(candidate / "payloads"),
                    )
                    is None
                )
    rows.append(
        "two parallel compiled sequences preserve every grouped result/state/lifetime"
    )
    manifest = out / "inputs.json"
    experiment.save(
        manifest,
        {
            "rows": [
                {
                    "workload": "control",
                    "workers": 1,
                    "input": str(toy),
                    "input_sha256": digest(toy),
                    "empty_input": str(empty),
                    "recording": str(out / "toy-record/recording"),
                }
            ]
        },
    )
    script = Path(__file__).with_name("experiment.py")
    experiment.execute(
        [sys.executable, "-B", script, "pilot", manifest, out / "pilot", "--smoke"],
        out / "pilot-process",
    )
    experiment.execute(
        [
            sys.executable,
            "-B",
            script,
            "run",
            out / "pilot/plan.json",
            out / "window",
            "--smoke",
        ],
        out / "window-process",
    )
    experiment.execute(
        [
            sys.executable,
            "-B",
            script,
            "analyse",
            out / "window/window.json",
            out / "analysis",
            "--smoke",
        ],
        out / "analysis-process",
    )
    analysis = json.loads((out / "analysis/analysis.json").read_text())
    assert len(analysis["rows"]) == 3 and all(
        not row["qualified_noise"] for row in analysis["rows"]
    )
    assert all(
        not kind["qualified"]
        for row in analysis["rows"]
        for kind in row["kinds"].values()
    )
    assert json.loads((out / "pilot/export-poison-result.json").read_text())["caught"]
    rows.append(
        "real small NetGraph CLI/consumer/replay window, export poison, lossless archives; smoke never qualifies"
    )

    bad_plan = json.loads((out / "pilot/plan.json").read_text())
    bad_plan.update(
        {
            "smoke": False,
            "projected_s": 1000,
            "started_monotonic": time.monotonic(),
            "holder": "owned-test",
        }
    )
    bad_file = out / "over-budget.json"
    experiment.save(bad_file, bad_plan)
    with (
        patch.dict(
            experiment.os.environ,
            {"HEAVY_SLOT_HELD": "1", "HEAVY_OWNER": "replay-timing"},
        ),
        patch.object(
            experiment, "HOLDER", SimpleNamespace(read_text=lambda: "owned-test")
        ),
    ):
        try:
            experiment.run(
                SimpleNamespace(
                    input=bad_file,
                    out=out / "never-started",
                    smoke=False,
                    minutes=0.001,
                )
            )
        except ValueError as error:
            assert "budget" in str(error)
        else:
            raise AssertionError("over-budget window accepted")
    rows.append("over-budget plan refuses before launching a block")
    fake = {
        "smoke": False,
        "plan": {
            "rows": [
                {
                    "workload": "statistics-control",
                    "workers": 1,
                    "null": {"rows": [{"in_span_ns": 1}]},
                }
            ]
        },
        "blocks": [],
    }
    for kind in ("replay", "cli", "matched"):
        for triple in range(7):
            for label in ("A1", "B", "A2"):
                block = {
                    "row": 0,
                    "experiment": kind,
                    "triple": triple,
                    "slot": label,
                    "union_s": 9.0 if label == "B" else 10.0,
                    "elapsed_s": 22.0 if label == "B" else 20.0,
                    "repeats": 1,
                    "rss_bytes": 20 * 1024**2,
                    "released_wrappers": {},
                    "before": None,
                    "after": None,
                    "kinds": {
                        "run:heavy": {
                            "calls": 1,
                            "wall_s": 9.0 if label == "B" else 10.0,
                            "thread_cpu_s": 1.0,
                        }
                    },
                }
                if label == "B":
                    block["kinds"]["run:missing"] = {
                        "calls": 1,
                        "wall_s": 1.0,
                        "thread_cpu_s": 1.0,
                    }
                fake["blocks"].append(block)
    fake_file = out / "statistics-fixture.json"
    experiment.save(fake_file, fake)
    fixture_out = out / "statistics-analysis"
    fixture_out.mkdir()
    experiment.analyse(SimpleNamespace(input=fake_file, out=fixture_out))
    result = json.loads((fixture_out / "analysis.json").read_text())
    assert all(
        row["union_elapsed_disagree"] and not row["qualified_noise"]
        for row in result["rows"]
    )
    assert not result["rows"][2]["kinds"]["run:missing"]["qualified"]
    rows.append(
        "opposite union/elapsed directions and missing E7 phase/kind stay inconclusive"
    )
    experiment.save(
        out / "summary.json",
        {
            "rows": rows,
            "qualification": "none: short untimed controls only",
            "summary": "PASS stage-2 controls",
        },
    )
    print(json.dumps({"rows": rows, "summary": "PASS stage-2 controls"}), flush=True)


if __name__ == "__main__":
    main()
