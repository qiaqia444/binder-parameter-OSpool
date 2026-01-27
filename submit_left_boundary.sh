#!/bin/bash

# Submit Left Boundary Scan jobs

echo "=========================================="
echo "Submitting Left Boundary Scan Jobs"
echo "=========================================="

# Check if params file exists
if [ ! -f "jobs/params_left_boundary.txt" ]; then
    echo "Error: jobs/params_left_boundary.txt not found!"
    echo "Run: julia jobs/make_params_left_boundary.jl"
    exit 1
fi

# Count number of jobs
NJOBS=$(wc -l < jobs/params_left_boundary.txt)
echo "Found $NJOBS parameter sets in params_left_boundary.txt"

echo ""
echo "Physics setup:"
echo "  - Left boundary scan (dephasing-induced transition)"
echo "  - Fixed λ_x = 0.3 (X measurement strength)"
echo "  - Fixed λ_zz = 0.0 (no ZZ measurements)"
echo "  - Scanning P_x from 0 to 0.5"
echo "  - Fixed P_zz = 0.0 (no ZZ dephasing)"
echo "  - Using DENSITY MATRIX evolution (correct physics!)"
echo ""

echo "Submitting Left Boundary Scan jobs..."
cd jobs
condor_submit jobs_left_boundary.submit
cd ..

echo ""
echo "Left Boundary Scan jobs submitted successfully!"
echo "   Jobs submitted: $NJOBS"
echo "   Output files: output/left_boundary_*.json"
echo ""
echo "Monitor with: condor_q"
echo "Check logs in: jobs/logs/"
