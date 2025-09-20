#!/bin/bash
# Setup script for forced measurement jobs on cluster
# Run this once before submitting forced measurement jobs

echo "=== Setting up Forced Measurement Job Environment ==="

# Ensure we're in the right directory
cd /home/qia.wang/binder-parameter-OSpool

# Sync with GitHub to get latest code
echo "Syncing with GitHub..."
git pull origin main

# Verify forced measurement files exist
echo "Verifying forced measurement files..."
required_files=(
    "src/BinderSimForced.jl"
    "run_forced.jl"
    "jobs/jobs_forced.submit"
    "jobs/run_forced.sh"
    "jobs/params_forced.txt"
)

for file in "${required_files[@]}"; do
    if [ -f "$file" ]; then
        echo "✓ $file exists"
    else
        echo "✗ Missing: $file"
        exit 1
    fi
done

# Create clean output directory
echo "Setting up output directory..."
if [ -d "output" ]; then
    backup_dir="output_backup_$(date +%Y%m%d_%H%M%S)"
    echo "Backing up existing output to $backup_dir"
    mv output "$backup_dir"
fi

mkdir -p output
mkdir -p logs
echo "✓ Created clean output and logs directories"

# Verify parameter file
param_count=$(wc -l < jobs/params_forced.txt)
echo "Parameter file has $param_count jobs (expecting 153)"

if [ "$param_count" -eq 153 ]; then
    echo "✓ Parameter file is correct"
else
    echo "✗ Parameter file has wrong number of jobs"
    exit 1
fi

echo ""
echo "=== Setup Complete ==="
echo "Ready to submit forced measurement jobs!"
echo ""
echo "To submit jobs, run:"
echo "  cd jobs"
echo "  condor_submit jobs_forced.submit"
echo ""
echo "Expected output:"
echo "  153 files in output/forced_*.json format"
echo "  Files will be: forced_L{8,12,16}_lam{0.1-0.9}_s{1,2,3}.json"