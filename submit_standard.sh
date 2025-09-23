#!/bin/bash

# Submit Standard BinderSim jobs

echo "=========================================="
echo "Submitting Standard BinderSim Jobs"
echo "=========================================="

# Check if params file exists
if [ ! -f "jobs/params_adaptive.txt" ]; then
    echo "Error: jobs/params_adaptive.txt not found!"
    exit 1
fi

# Count number of jobs
NJOBS=$(wc -l < jobs/params_adaptive.txt)
echo "Found $NJOBS parameter sets in params_adaptive.txt"

echo ""
echo "Submitting Standard BinderSim jobs..."
cd jobs
condor_submit jobs_standard.submit
cd ..

echo ""
echo "Standard BinderSim jobs submitted successfully!"
echo "   Jobs submitted: $NJOBS"
echo "   Output files: output/standard_*.json"
echo ""
echo "Monitor with: condor_q"