"""Prove that the actual stream comparator rejects a wrong cut and ledger state."""

import argparse
import gzip
import json
import shutil
import subprocess
from pathlib import Path

from events import REFERENCE, Blobs, read_stream
from values import dumps


def result_cut(node):
    if isinstance(node, dict):
        if (
            node.get("$type") == "PlacementResult"
            and isinstance(node.get("min_cut"), dict)
            and node["min_cut"].get("$list")
        ):
            node["min_cut"]["$list"][0] += 1
            return True
        return any(result_cut(v) for v in node.values())
    if isinstance(node, list):
        return any(result_cut(v) for v in node)
    return False


def ledger_unit(record, blobs):
    for state in record.get("state", {}).values():
        if state and "full" in state and state["full"].get("$type") == "FlowLedger":
            spec = state["full"]["carried"]["$ndarray"]
            data = blobs.get(spec)
            data.reshape(-1)[0] += 1
            state["full"]["carried"]["$ndarray"] = blobs.put(data)
            return True
    return False


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("reference", type=Path)
    parser.add_argument("out", type=Path)
    args = parser.parse_args()
    args.out.mkdir()
    rows = []
    for poison in ("min_cut", "ledger_carried"):
        candidate = args.out / poison
        shutil.copytree(args.reference, candidate)
        manifest = json.loads((candidate / "manifest.json").read_text())
        caught, witness = False, None
        for job in sorted(manifest["jobs"]):
            records = list(read_stream(candidate, job))
            blobs = Blobs(candidate / "payloads")
            for record in records:
                changed = (
                    result_cut(record.get("result"))
                    if poison == "min_cut"
                    else ledger_unit(record, blobs)
                )
                if changed:
                    witness = {"job": job, "call": record["call"]}
                    break
            if witness:
                with gzip.open(
                    candidate / (job + ".jsonl.gz"), "wt", encoding="utf-8"
                ) as stream:
                    for record in records:
                        stream.write(dumps(record) + "\n")
                    stream.write(dumps({"footer": len(records)}) + "\n")
                break
        if witness is None:
            raise AssertionError(f"preflight never reached {poison}")
        command = [
            str(REFERENCE / "bin/python"),
            "-B",
            str(Path(__file__).with_name("replay.py")),
            "compare",
            str(args.reference),
            str(candidate),
        ]
        result = subprocess.run(command, capture_output=True, text=True, timeout=30)
        caught = result.returncode != 0 and all(
            token in result.stderr
            for token in ("job=", "call=", "field=", "reference=", "candidate=")
        )
        with (args.out / (poison + ".log")).open("x") as log:
            log.write(
                result.stdout
                + result.stderr
                + f"\nSUMMARY: poison={poison} caught={caught} witness={witness}\n"
            )
        if not caught:
            raise AssertionError(
                f"comparator did not reject {poison} with a field witness"
            )
        rows.append(
            {"poison": poison, **witness, "exit": result.returncode, "caught": caught}
        )
    report = {
        "rows": rows,
        "summary": "PASS both real-result and real-ledger poisonings rejected",
    }
    (args.out / "summary.json").write_text(json.dumps(report, indent=2) + "\n")
    print(dumps(report))


if __name__ == "__main__":
    main()
