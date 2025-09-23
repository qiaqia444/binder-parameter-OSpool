#!/bin/bash

# Submit Dummy Site BinderSim jobs

echo "=========================================="
echo "Submitting Dummy Site BinderSim Jobs"
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
echo "Submitting Dummy Site BinderSim jobs..."
cd jobs
condor_submit jobs_dummy.submit
cd ..

echo ""
echo "Dummy Site BinderSim jobs submitted successfully!"
echo "   Jobs submitted: $NJOBS"
echo "   Output files: output/dummy_*.json"
echo ""
echo "Monitor with: condor_q"