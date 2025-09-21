#!/bin/bash
# Setup script for adaptive random sampling jobs on cluster
# Run this once before submitting adaptive measurement jobs

echo "=== Setting up Adaptive Random Sampling Job Environment ==="

# Ensure we're in the right directory
cd /home/qia.wang/binder-parameter-OSpool

# Sync with GitHub to get latest code with simultaneous ZZ measurements
echo "Syncing with GitHub..."
git pull origin main

# Verify adaptive measurement files exist
echo "Verifying adaptive measurement files..."
required_files=(
    "src/BinderSim.jl"
    "run_adaptive.jl"
    "jobs/jobs_adaptive.submit"
    "jobs/run_adaptive.sh"
    "jobs/params_adaptive.txt"
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
param_count=$(wc -l < jobs/params_adaptive.txt)
echo "Parameter file has $param_count jobs (expecting 153)"

if [ "$param_count" -eq 153 ]; then
    echo "✓ Parameter file is correct"
else
    echo "✗ Parameter file has wrong number of jobs"
    exit 1
fi

# Show parameter summary
echo ""
echo "=== Parameter Summary ==="
echo "System sizes: $(grep -o 'L [0-9]*' jobs/params_adaptive.txt | sort -u)"
echo "Lambda values: $(grep -o 'lambda_x [0-9.]*' jobs/params_adaptive.txt | sort -u | wc -l) values"
echo "Trials per job: $(grep -o 'ntrials [0-9]*' jobs/params_adaptive.txt | head -1)"
echo "Total jobs: $param_count"

echo ""
echo "=== Setup Complete ==="
echo "Ready to submit adaptive random sampling jobs!"
echo ""
echo "To submit jobs, run:"
echo "  cd jobs"
echo "  condor_submit jobs_adaptive.submit"
echo ""
echo "Expected output:"
echo "  153 files in output/L{8,12,16}_lam{0.1-0.9}_s{1,2,3}.json format"
echo "  Each job runs 2000 random measurement trials"
echo ""
echo "Monitor with:"
echo "  condor_q"
echo "  watch -n 30 'condor_q | tail -10'"