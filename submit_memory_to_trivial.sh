#!/bin/bash

# Submit Memory-to-Trivial Transition Scan jobs

echo "=========================================="
echo "Submitting Memory-to-Trivial Scan Jobs"
echo "=========================================="

# Check if params file exists
if [ ! -f "jobs/params_memory_to_trivial.txt" ]; then
    echo "Error: jobs/params_memory_to_trivial.txt not found!"
    echo "Run: julia jobs/make_params_memory_to_trivial.jl"
    exit 1
fi

# Count number of jobs
NJOBS=$(wc -l < jobs/params_memory_to_trivial.txt)
echo "Found $NJOBS parameter sets in params_memory_to_trivial.txt"

echo ""
echo "Physics setup:"
echo "  - Memory-to-trivial transition scan"
echo "  - Fixed λ_x = 0.1 (X measurements)"
echo "  - Fixed λ_zz = 0.7 (ZZ measurement strength)"
echo "  - Scanning P_x from 0 to 0.5"
echo "  - P_zz = P_x (coupled dephasing)"
echo "  - Using DENSITY MATRIX evolution (correct physics!)"
echo ""

echo "Submitting Memory-to-Trivial Scan jobs..."
cd jobs
condor_submit jobs_memory_to_trivial.submit
cd ..

echo ""
echo "Memory-to-Trivial Scan jobs submitted successfully!"
echo "   Jobs submitted: $NJOBS"
echo "   Output files: output/memory_to_trivial_*.json"
echo ""
echo "Monitor with: condor_q"
echo "Check logs in: jobs/logs/"
