#!/usr/bin/env bash
set -euo pipefail

# Use Singularity container with Julia
CONTAINER_PATH="/ospool/ap40/data/qia.wang/container.sif"

# Create project directory structure if needed
mkdir -p output

# Run the main simulation script using Singularity
singularity exec --bind $(pwd):/work --pwd /work $CONTAINER_PATH julia --project=. run.jl ${PARAMS_JSON:+--params "$PARAMS_JSON"} --outdir output
