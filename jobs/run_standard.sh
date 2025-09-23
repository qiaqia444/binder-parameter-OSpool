#!/bin/bash

# Standard BinderSim job runner
# This script runs the standard quantum trajectory simulation

echo "=========================================="
echo "Starting Standard BinderSim Job"
echo "Arguments: $@"
echo "Working directory: $(pwd)"
echo "=========================================="

# Load container if available
if [ -f "image.def" ]; then
    echo "Container definition found, using Apptainer..."
    # Run with container
    apptainer exec --bind $(pwd):/work --pwd /work image.sif julia --project=. run.jl "$@"
else
    echo "No container found, running directly..."
    # Run directly with Julia
    julia --project=. run.jl "$@"
fi

exit_code=$?
echo "=========================================="
echo "Standard BinderSim job completed with exit code: $exit_code"
echo "=========================================="
exit $exit_code