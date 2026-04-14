#!/bin/bash

# Submit Memory-to-Trivial Transition Scan jobs (λ_x = 0.21, λ_zz = 0.49)

echo "=========================================="
echo "Submitting Memory-to-Trivial Scan Jobs (λ_x = 0.21, λ_zz = 0.49)"
echo "=========================================="

# Check if params file exists
if [ ! -f "jobs/params_memory_to_trivial_lx0.3.txt" ]; then
    echo "Error: jobs/params_memory_to_trivial_lx0.3.txt not found!"
    echo "Run: julia jobs/make_params_memory_to_trivial_lx0.3.jl"
    exit 1
fi

# Count number of jobs
NJOBS=$(wc -l < jobs/params_memory_to_trivial_lx0.3.txt)
echo "Found $NJOBS parameter sets in params_memory_to_trivial_lx0.3.txt"

echo ""
echo "Physics setup:"
echo "  - Memory-to-trivial transition scan"
echo "  - Fixed λ_x = 0.21 (X measurements)"
echo "  - Fixed λ_zz = 0.49 (ZZ measurement strength)"
echo "  - Scanning P_x from 0 to 0.5"
echo "  - P_zz = P_x (coupled dephasing)"
echo "  - Using DENSITY MATRIX evolution (optimized implementation)"
echo "  - 24 hours per trial for improved statistics"
echo ""

echo "Submitting Memory-to-Trivial Scan jobs..."
cd jobs
condor_submit jobs_memory_to_trivial_lx0.3.submit
cd ..

echo ""
echo "Memory-to-Trivial Scan jobs (λ_x = 0.21, λ_zz = 0.49) submitted successfully!"
echo "   Jobs submitted: $NJOBS"
echo "   Output files: output/memory_to_trivial_L*_lx0.21_lzz0.49_*.json"
echo ""
echo "Monitor with: condor_q"
echo "Check logs in: jobs/logs/"
