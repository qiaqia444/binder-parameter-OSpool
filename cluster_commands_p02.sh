#!/bin/bash

# Complete workflow for running P=0.2 dephasing simulations on cluster
# Copy and paste these commands one by one, or run this entire script

echo "=========================================="
echo "DEPHASING P=0.2 CLUSTER WORKFLOW"
echo "=========================================="

# Step 1: Navigate to project directory and pull latest code
echo ""
echo "Step 1: Pulling latest code from GitHub..."
cd ~/binder-parameter-OSpool
git pull origin main

# Step 2: Verify files are present
echo ""
echo "Step 2: Verifying P=0.2 files..."
ls -lh run_dephasing_p02.jl
ls -lh jobs/run_dephasing_p02.sh
ls -lh jobs/dephasing_p02_jobs.submit
ls -lh jobs/params_dephasing_p02.txt
wc -l jobs/params_dephasing_p02.txt

# Step 3: Submit jobs to HTCondor
echo ""
echo "Step 3: Submitting 153 jobs to HTCondor..."
cd jobs
condor_submit dephasing_p02_jobs.submit

# Step 4: Check job status
echo ""
echo "Step 4: Checking job status..."
echo "Run this command periodically to monitor progress:"
echo "  condor_q"
echo ""
echo "To see detailed status:"
echo "  condor_q -nobatch"
echo ""
echo "To check for held jobs:"
echo "  condor_q -hold"
echo ""

echo "=========================================="
echo "Jobs submitted! Monitor with: condor_q"
echo "=========================================="
echo ""
echo "After jobs complete, run:"
echo "  cd ~/binder-parameter-OSpool"
echo "  ./collect_dephasing_p02_results.sh"
echo "  wormhole send dephasing_p02_results_*.tar.gz"
echo "=========================================="

