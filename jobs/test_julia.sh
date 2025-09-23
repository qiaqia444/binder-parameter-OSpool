#!/bin/bash

# Test script to verify Julia container functionality
echo "=== Julia Container Test Start ==="
echo "Job started at: $(date)"
echo "Running on: $(hostname)"
echo "Working directory: $(pwd)"

# List available files
echo "Available files:"
ls -la

# Check available container runtimes
echo "Container runtime check:"
echo "Apptainer available: $(command -v apptainer >/dev/null 2>&1 && echo 'YES' || echo 'NO')"
echo "Singularity available: $(command -v singularity >/dev/null 2>&1 && echo 'YES' || echo 'NO')"

# Test container functionality
if command -v apptainer >/dev/null 2>&1; then
    echo "Testing Apptainer with Julia container..."
    # Pull Julia container
    echo "Pulling Julia container..."
    apptainer pull julia.sif docker://julia:1.11
    # Test Julia in container
    echo "Testing Julia in container..."
    apptainer exec julia.sif julia --version
    apptainer exec julia.sif julia -e 'println("Julia container test successful!")'
    # Test package installation
    echo "Testing package environment setup..."
    apptainer exec julia.sif julia --project=. -e 'using Pkg; Pkg.status(); println("Package environment test successful!")'
    EXIT_CODE=0
elif command -v singularity >/dev/null 2>&1; then
    echo "Testing Singularity with Julia container..."
    # Pull Julia container
    echo "Pulling Julia container..."
    singularity pull julia.sif docker://julia:1.11
    # Test Julia in container
    echo "Testing Julia in container..."
    singularity exec julia.sif julia --version
    singularity exec julia.sif julia -e 'println("Julia container test successful!")'
    # Test package installation
    echo "Testing package environment setup..."
    singularity exec julia.sif julia --project=. -e 'using Pkg; Pkg.status(); println("Package environment test successful!")'
    EXIT_CODE=0
else
    echo "ERROR: No container runtime found!"
    EXIT_CODE=1
fi

echo "Test completed with exit code: $EXIT_CODE"
echo "Job finished at: $(date)"
echo "=== Julia Container Test End ==="

exit $EXIT_CODE