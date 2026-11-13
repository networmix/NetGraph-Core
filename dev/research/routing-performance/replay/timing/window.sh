#!/bin/bash
# Run calibration and all rows inside the coordinator's single heavy-slot hold.
set -euo pipefail
script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
effort_dir=$(cd -- "$script_dir/../../.." && pwd)
task_python="$effort_dir/scratch/consumer/ref/venv/bin/python"
inputs="$1"
output="$2"
mkdir "$output"
"$task_python" -B "$script_dir/experiment.py" pilot "$inputs" "$output/pilot"
exec "$task_python" -B "$script_dir/experiment.py" run "$output/pilot/plan.json" "$output/run" --minutes 120
