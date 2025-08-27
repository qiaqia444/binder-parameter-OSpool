#!/usr/bin/env bash
set -euo pipefail

# Use the Julia available on the execution node
JULIA_CMD="julia"

# Set up Julia environment
export JULIA_DEPOT_PATH="$(pwd)/.julia"
mkdir -p "$JULIA_DEPOT_PATH"

# Create project directory structure if needed
mkdir -p output

# Run the main simulation script
$JULIA_CMD --project=. run.jl ${PARAMS_JSON:+--params "$PARAMS_JSON"} --outdir output
