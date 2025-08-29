#!/usr/bin/env bash
set -euo pipefail

# Find Singularity executable
if command -v singularity >/dev/null 2>&1; then
    SINGULARITY_CMD="singularity"
elif command -v apptainer >/dev/null 2>&1; then
    SINGULARITY_CMD="apptainer"
elif [ -f "/usr/bin/singularity" ]; then
    SINGULARITY_CMD="/usr/bin/singularity"
elif [ -f "/usr/local/bin/singularity" ]; then
    SINGULARITY_CMD="/usr/local/bin/singularity"
else
    echo "ERROR: Neither singularity nor apptainer found!"
    echo "Available commands:"
    which -a singularity apptainer || echo "None found"
    echo "PATH: $PATH"
    exit 1
fi

echo "Using container runtime: $SINGULARITY_CMD"

# Use transferred container file (transferred as part of job)
CONTAINER_PATH="./container.sif"

# Check if container file exists
if [ ! -f "$CONTAINER_PATH" ]; then
    echo "ERROR: Container file not found at $CONTAINER_PATH"
    echo "Current directory contents:"
    ls -la
    exit 1
fi

echo "Using container: $CONTAINER_PATH"

# Create project directory structure if needed
mkdir -p output
mkdir -p .julia_local

# Set up environment for writable Julia depot
export JULIA_DEPOT_PATH="$(pwd)/.julia_local:/tmp/julia"

# All arguments are passed directly to Julia
# HTCondor will pass them as individual arguments
$SINGULARITY_CMD exec --bind $(pwd):/work --pwd /work $CONTAINER_PATH julia --project=. run.jl "$@"
