#!/bin/bash

# BinderSimWithDummy job runner
# This script runs the dummy site version of BinderSim

echo "=========================================="
echo "Starting BinderSim with Dummy Site Job"
echo "Arguments: $@"
echo "Working directory: $(pwd)"
echo "=========================================="

# Load container if available
if [ -f "image.def" ]; then
    echo "Container definition found, using Apptainer..."
    # Run with container
    apptainer exec --bind $(pwd):/work --pwd /work image.sif julia --project=. run_dummy.jl "$@"
else
    echo "No container found, running directly..."
    # Run directly with Julia
    julia --project=. run_dummy.jl "$@"
fi

exit_code=$?
echo "=========================================="
echo "BinderSim with Dummy job completed with exit code: $exit_code"
echo "=========================================="
exit $exit_code