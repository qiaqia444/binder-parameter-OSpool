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

# Create output directory
mkdir -p output

# Check if Apptainer/Singularity is available and use container
if command -v apptainer >/dev/null 2>&1; then
    echo "Using Apptainer with Julia container..."
    # Use pre-pulled Julia container or pull if needed
    if [ ! -f julia.sif ]; then
        echo "Pulling Julia container..."
        apptainer pull julia.sif docker://julia:1.11
    fi
    # Install packages in container
    echo "Setting up Julia environment in container..."
    apptainer exec julia.sif julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'
    # Run simulation
    echo "Running forced simulation in container..."
    echo "Command: apptainer exec julia.sif julia --project=. run_forced.jl $ARGS"
    apptainer exec julia.sif julia --project=. run_forced.jl $ARGS
    EXIT_CODE=$?
elif command -v singularity >/dev/null 2>&1; then
    echo "Using Singularity with Julia container..."
    # Use Singularity (older version)
    if [ ! -f julia.sif ]; then
        echo "Pulling Julia container..."
        singularity pull julia.sif docker://julia:1.11
    fi
    # Install packages in container
    echo "Setting up Julia environment in container..."
    singularity exec julia.sif julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'
    # Run simulation
    echo "Running forced simulation in container..."
    echo "Command: singularity exec julia.sif julia --project=. run_forced.jl $ARGS"
    singularity exec julia.sif julia --project=. run_forced.jl $ARGS
    EXIT_CODE=$?
else
    echo "No container runtime found, this job requires Apptainer/Singularity"
    exit 1
fi

echo "Simulation completed with exit code: $EXIT_CODE"

# List output files
echo "Output files created:"
ls -la output/ 2>/dev/null || echo "No output directory found"

echo "Job finished at: $(date)"
echo "=== Forced +1 Measurement Job End ==="

exit $EXIT_CODE