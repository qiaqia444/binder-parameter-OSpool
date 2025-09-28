#!/bin/bash

# Test script to verify Julia container functionality
echo "=== Julia Container Test Start ==="
echo "Job started at: $(date)"
echo "Running on: $(hostname)"
echo "Working directory: $(pwd)"

# List available files
echo "Available files:"
ls -la

# Test Julia functionality
echo "Testing Julia..."
julia --version
julia -e 'println("Julia container test successful!")'

# Test package installation
echo "Testing package environment setup..."
julia --project=. -e 'using Pkg; Pkg.status(); println("Package environment test successful!")'
EXIT_CODE=$?

echo "Test completed with exit code: $EXIT_CODE"
echo "Job finished at: $(date)"
echo "=== Julia Container Test End ==="

exit $EXIT_CODE