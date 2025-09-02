#!/usr/bin/env bash
set -euo pipefail

# Create output and logs directories
mkdir -p output logs

# Detect which script to run based on output prefix
if [[ "${8:-}" == manual_* ]]; then
    SCRIPT="run_manual.jl"
    echo "Running manual correlator script: $SCRIPT"
else
    SCRIPT="run_standard.jl"
    echo "Running standard correlator script: $SCRIPT"
fi

# In container environment, julia should be available directly
echo "Executing: julia --project=. $SCRIPT $*"
exec julia --project=. "$SCRIPT" "$@"
