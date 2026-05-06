#!/bin/bash

# Submit Learning-to-Trivial Transition Scan jobs (λ_x = 0.49, λ_zz = 0.21)

echo "=========================================="
echo "Submitting Learning-to-Trivial Scan Jobs (λ_x = 0.49, λ_zz = 0.21)"
echo "=========================================="

# Check if params file exists
if [ ! -f "jobs/params_learning_to_trivial_lambda0.7.txt" ]; then
    echo "Error: jobs/params_learning_to_trivial_lambda0.7.txt not found!"
    echo "Run: julia jobs/make_params_learning_to_trivial_lambda0.7.jl"
    exit 1
fi

# Count number of jobs
NJOBS=$(wc -l < jobs/params_learning_to_trivial_lambda0.7.txt)
echo "Found $NJOBS parameter sets in params_learning_to_trivial_lambda0.7.txt"

echo ""
echo "Physics setup:"
echo "  - Learning-to-trivial transition scan (lambda=0.7 family)"
echo "  - Fixed λ_x = 0.49 (X measurements)"
echo "  - Fixed λ_zz = 0.21 (ZZ measurement strength)"
echo "  - Scanning P_x from 0 to 0.5"
echo "  - P_zz = P_x (coupled dephasing)"
echo "  - Using DENSITY MATRIX evolution with Rényi-2 Binder"
echo "  - 24 hours per trial for improved statistics"
echo ""

echo "Submitting Learning-to-Trivial Scan jobs..."
cd jobs
condor_submit jobs_learning_to_trivial_lambda0.7.submit
cd ..

echo ""
echo "Learning-to-Trivial Scan jobs (λ_x = 0.49, λ_zz = 0.21) submitted successfully!"
echo "   Jobs submitted: $NJOBS"
echo "   Output files: output/learning_to_trivial_lambda0.7_L*_lx0.49_lzz0.21_*.json"
echo ""
echo "Monitor with: condor_q"
echo "Check logs in: jobs/logs/"
