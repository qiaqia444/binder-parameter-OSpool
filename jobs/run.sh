#!/usr/bin/env bash
set -euo pipefail

# Use Singularity container with Julia
CONTAINER_PATH="/ospool/ap40/data/qia.wang/container.sif"

# Create project directory structure if needed
mkdir -p output

# All arguments are passed directly to Julia
# HTCondor will pass them as individual arguments
singularity exec --bind $(pwd):/work --pwd /work $CONTAINER_PATH julia --project=. run.jl "$@"
