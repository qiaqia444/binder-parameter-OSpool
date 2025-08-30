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

# Try different Julia paths commonly available on OSG
JULIA_PATHS=(
    "/cvmfs/oasis.opensciencegrid.org/mis/apptainer/images/julia/bin/julia"
    "/cvmfs/singularity.opensciencegrid.org/opensciencegrid/julia/bin/julia"
    "/usr/local/bin/julia"
    "/usr/bin/julia"
    "julia"
)

JULIA_CMD=""
for path in "${JULIA_PATHS[@]}"; do
    if command -v "$path" &> /dev/null; then
        JULIA_CMD="$path"
        echo "Found Julia at: $JULIA_CMD"
        break
    fi
done

if [[ -z "$JULIA_CMD" ]]; then
    echo "ERROR: Julia not found in any of the expected locations"
    echo "Tried paths: ${JULIA_PATHS[*]}"
    exit 1
fi

# Execute the appropriate script
echo "Executing: $JULIA_CMD --project=. $SCRIPT $*"
exec "$JULIA_CMD" --project=. "$SCRIPT" "$@"
