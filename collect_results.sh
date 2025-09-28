#!/bin/bash

# Collect and organize forced measurement Binder parameter results
# Run this script after all HTCondor jobs complete

echo "=== Binder Parameter Results Collection ==="
echo "Starting collection at: $(date)"

# Create organized output directories
mkdir -p results/forced_measurements/L8
mkdir -p results/forced_measurements/L12  
mkdir -p results/forced_measurements/L16
mkdir -p results/analysis

# Navigate to jobs directory
cd jobs

# Check if output directory exists
if [ ! -d "output" ]; then
    echo "ERROR: No output directory found. Jobs may not have completed yet."
    exit 1
fi

echo "Found output directory with $(ls output/*.json 2>/dev/null | wc -l) result files"

# Organize by system size
echo "Organizing results by system size..."
cp output/forced_L8_*.json ../results/forced_measurements/L8/ 2>/dev/null
cp output/forced_L12_*.json ../results/forced_measurements/L12/ 2>/dev/null  
cp output/forced_L16_*.json ../results/forced_measurements/L16/ 2>/dev/null

# Count results
L8_count=$(ls ../results/forced_measurements/L8/*.json 2>/dev/null | wc -l)
L12_count=$(ls ../results/forced_measurements/L12/*.json 2>/dev/null | wc -l)
L16_count=$(ls ../results/forced_measurements/L16/*.json 2>/dev/null | wc -l)

echo "Results organized:"
echo "  L=8:  $L8_count files"
echo "  L=12: $L12_count files" 
echo "  L=16: $L16_count files"
echo "  Total: $((L8_count + L12_count + L16_count)) files"

# Check for failures
failure_count=$(ls output/*_FAILED.json 2>/dev/null | wc -l)
if [ $failure_count -gt 0 ]; then
    echo "WARNING: Found $failure_count failed jobs"
    cp output/*_FAILED.json ../results/analysis/ 2>/dev/null
    echo "Failure files copied to results/analysis/"
fi

# Create summary report
echo "Creating summary report..."
cat > ../results/analysis/collection_summary.txt << EOF
Binder Parameter Forced Measurement Results Summary
Collection Date: $(date)
=================================================

Total Expected Jobs: 153 (51 per system size)
Total Completed Jobs: $((L8_count + L12_count + L16_count))
Total Failed Jobs: $failure_count

System Size Breakdown:
- L=8:  $L8_count/51 completed ($(( 100 * L8_count / 51 ))%)
- L=12: $L12_count/51 completed ($(( 100 * L12_count / 51 ))%)  
- L=16: $L16_count/51 completed ($(( 100 * L16_count / 51 ))%)

Lambda Values: 0.1, 0.2, 0.3, 0.4, 0.46-0.54 (17 values)
Samples per lambda: 3 (2000 trials each forced measurement)

File Format: Each JSON contains:
- Complete simulation parameters (L, lambda_x, lambda_zz, lambda)
- Binder parameter result and statistics  
- Computational metadata (maxdim, cutoff, seed)
- Success/failure status

Next Steps for Analysis:
1. Verify all expected files are present
2. Load JSON data for statistical analysis
3. Plot Binder parameter vs lambda for each L
4. Analyze finite-size scaling behavior
5. Compare with standard quantum trajectory results

Archive Command for Transfer:
tar -czf forced_binder_results_\$(date +%Y%m%d).tar.gz results/
EOF

echo "Summary report created: results/analysis/collection_summary.txt"

# Create transfer archive
cd ..
echo "Creating transfer archive..."
tar -czf "forced_binder_results_$(date +%Y%m%d_%H%M).tar.gz" results/

archive_size=$(du -h "forced_binder_results_$(date +%Y%m%d_%H%M).tar.gz" | cut -f1)
echo "Archive created: forced_binder_results_$(date +%Y%m%d_%H%M).tar.gz ($archive_size)"

echo ""
echo "=== Collection Complete ==="
echo "To transfer to your Mac using Magic Wormhole:"
echo "  1. On OSPool: wormhole send forced_binder_results_$(date +%Y%m%d_%H%M).tar.gz"
echo "  2. On Mac: wormhole receive [code]"
echo ""
echo "After transfer, extract with:"
echo "  tar -xzf forced_binder_results_*.tar.gz"
echo ""
echo "Results will be organized in results/ directory ready for analysis!"