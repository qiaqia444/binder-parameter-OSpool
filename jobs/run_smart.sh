#!/usr/bin/env bash
set -euo pipefail
mkdir -p output logs

# Check if we should run manual correlators based on output prefix
if [[ "${8:-}" == manual_* ]]; then
    echo "Running with manual correlators..."
    exec julia --project=. run_manual.jl "$@"
else
    echo "Running with standard ITensorCorrelators..."
    exec julia --project=. run_standard.jl "$@"
fi
