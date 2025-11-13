#!/bin/bash

# Script to collect results from completed P=0.5 dephasing jobs
# Run this on the submit node after jobs complete

echo "Collecting P=0.5 dephasing simulation results..."

# Create results directory
mkdir -p results_dephasing_p05

# Count total jobs
total_jobs=$(wc -l < params_dephasing_p05.txt)
echo "Total jobs submitted: $total_jobs"

# Collect all completed result files
count=0
for dir in dephasing_p05_L*_P0.5_s*/; do
    if [ -d "$dir" ]; then
        if [ -d "$dir/output" ]; then
            # Copy all JSON files from output directory
            for file in "$dir"/output/*.json; do
                if [ -f "$file" ]; then
                    cp "$file" results_dephasing_p05/
                    ((count++))
                fi
            done
        fi
    fi
done

echo "Collected $count result files"
echo "Results saved in: results_dephasing_p05/"

# List some examples
echo ""
echo "Sample of collected files:"
ls results_dephasing_p05/ | head -10

# Summary
echo ""
echo "Collection complete!"
echo "  Expected: $total_jobs"
echo "  Found:    $count"

if [ $count -eq $total_jobs ]; then
    echo "  Status:   ✓ All jobs completed successfully"
else
    echo "  Status:   ⚠ Some jobs may be missing"
    missing=$((total_jobs - count))
    echo "  Missing:  $missing jobs"
fi
