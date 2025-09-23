#!/bin/bash

# Submit all three Binder simulation methods
# Uses the same params_adaptive.txt for all three

echo "=========================================="
echo "Submitting Three Binder Simulation Methods"
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
echo "SUBMITTING JOBS:"

# 1. Standard BinderSim
echo "1. Submitting Standard BinderSim jobs..."
cd jobs
condor_submit jobs_standard.submit
echo "   Standard jobs submitted"

# 2. Forced BinderSim  
echo "2. Submitting Forced BinderSim jobs..."
condor_submit jobs_forced.submit
echo "   Forced jobs submitted"

# 3. Dummy Site BinderSim
echo "3. Submitting Dummy Site BinderSim jobs..."
condor_submit jobs_dummy.submit
echo "   Dummy jobs submitted"

cd ..

echo ""
echo "SUMMARY:"
echo "   Standard BinderSim:     $NJOBS jobs submitted"
echo "   Forced BinderSim:       $NJOBS jobs submitted" 
echo "   Dummy Site BinderSim:   $NJOBS jobs submitted"
echo "   Total jobs:             $((NJOBS * 3))"

echo ""
echo "Output will be saved to:"
echo "   Standard:  output/standard_*.json"
echo "   Forced:    output/forced_*.json"
echo "   Dummy:     output/dummy_*.json"

echo ""
echo "Monitor jobs with:"
echo "   condor_q"
echo "   condor_q -better-analyze"

echo ""
echo "All jobs submitted successfully!"