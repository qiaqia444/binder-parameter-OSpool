#!/bin/bash

# Submit Forced BinderSim jobs

echo "=========================================="
echo "Submitting Forced BinderSim Jobs"
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
echo "Submitting Forced BinderSim jobs..."
cd jobs
condor_submit jobs_forced.submit
cd ..

echo ""
echo "Forced BinderSim jobs submitted successfully!"
echo "   Jobs submitted: $NJOBS"
echo "   Output files: output/forced_*.json"
echo ""
echo "Monitor with: condor_q"