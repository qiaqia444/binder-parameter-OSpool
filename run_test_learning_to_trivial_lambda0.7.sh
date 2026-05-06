#!/bin/bash

# Test script for Learning-to-Trivial Lambda0.7
# Quick verification before submitting 2,200 jobs

echo "=========================================="
echo "Testing Learning-to-Trivial Lambda0.7"
echo "=========================================="
echo ""
echo "This will run a quick test with:"
echo "  - L = 8"
echo "  - P_x = 0.3 (mid-transition region)"
echo "  - ntrials = 50 (fast)"
echo "  - Expected runtime: ~20-30 seconds"
echo ""

# Check Julia
if ! command -v julia &> /dev/null; then
    echo "ERROR: Julia not found in PATH"
    exit 1
fi

echo "Running test..."
julia --project=. test_learning_to_trivial_lambda0.7.jl

EXIT_CODE=$?

echo ""
if [ $EXIT_CODE -eq 0 ]; then
    echo "✓ Test completed successfully!"
    echo "Results saved to: test_learning_to_trivial_lambda0.7_results/"
    echo ""
    echo "Next step: Submit full batch"
    echo "  ./submit_learning_to_trivial_lambda0.7.sh"
else
    echo "✗ Test failed with exit code $EXIT_CODE"
    exit $EXIT_CODE
fi
