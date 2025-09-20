#!/bin/bash

# HTCondor job script for forced +1 measurement Binder parameter calculation
# This script runs a single forced measurement simulation

# Get arguments from HTCondor
ARGS="$@"

echo "=== Forced +1 Measurement Job Start ==="
echo "Job started at: $(date)"
echo "Running on: $(hostname)"
echo "Arguments: $ARGS"
echo "Working directory: $(pwd)"

# List available files
echo "Available files:"
ls -la

# Check Julia installation
echo "Julia version:"
julia --version

# Activate Julia environment and run the simulation
echo "Starting forced +1 measurement simulation..."
echo "Command: julia --project=. run_forced.jl $ARGS"

julia --project=. run_forced.jl $ARGS
EXIT_CODE=$?

echo "Simulation completed with exit code: $EXIT_CODE"

# List output files
echo "Output files created:"
ls -la output/ 2>/dev/null || echo "No output directory found"

echo "Job finished at: $(date)"
echo "=== Forced +1 Measurement Job End ==="

exit $EXIT_CODE