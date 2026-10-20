#!/bin/bash

# Submit Memory-to-Trivial EA Binder FINE-grid scan jobs
# (lambda = 0.3: lambda_x = 0.21, lambda_zz = 0.49; P = 0.025:0.005:0.075)

echo "=========================================="
echo "Submitting Memory-to-Trivial EA Binder Scan Jobs (lambda = 0.3: λ_x = 0.21, λ_zz = 0.49)"
echo "=========================================="

# Check if params file exists
if [ ! -f "jobs/params_memory_to_trivial_lambda0.3_fine.txt" ]; then
    echo "Error: jobs/params_memory_to_trivial_lambda0.3_fine.txt not found!"
    echo "Run: cd jobs && julia make_params_memory_to_trivial_lambda0.3_fine.jl"
    exit 1
fi

# Count number of jobs
NJOBS=$(wc -l < jobs/params_memory_to_trivial_lambda0.3_fine.txt)
echo "Found $NJOBS parameter sets in params_memory_to_trivial_lambda0.3_fine.txt"

echo ""
echo "Physics setup:"
echo "  - Memory-to-trivial transition scan, plain EA Binder (not dual EA)"
echo "  - Fixed λ_x = 0.21 (X measurements)"
echo "  - Fixed λ_zz = 0.49 (ZZ measurement strength)"
echo "  - Scanning P_x from 0.025 to 0.075 in steps of 0.005"
echo "  - P_zz = P_x (coupled dephasing)"
echo "  - Using DENSITY MATRIX evolution (optimized implementation)"
echo "  - 24 hours per job"
echo ""

echo "Submitting Memory-to-Trivial EA Binder Scan jobs..."
mkdir -p jobs/logs
cd jobs
condor_submit jobs_memory_to_trivial_lambda0.3_fine.submit
cd ..

echo ""
echo "Memory-to-Trivial EA Binder Scan jobs (lambda = 0.3) submitted successfully!"
echo "   Jobs submitted: $NJOBS"
echo "   Output files: jobs/output/memory_to_trivial_lambda0.3_fine_L*_lx0.21_lzz0.49_*.json"
echo ""
echo "Monitor with: condor_q"
echo "Check logs in: jobs/logs/"
