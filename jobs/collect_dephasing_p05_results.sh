#!/bin/bash

# Script to collect results from completed P=0.5 dephasing jobs (CORRECTED VERSION)
# Run this on the submit node after jobs complete
# L=8,12,16,20 (small systems)

echo "Collecting P=0.5 dephasing simulation results (CORRECTED)..."
echo "Small systems: L=8,12,16,20"

# Create results directory with timestamp
TIMESTAMP=$(date +%Y%m%d_%H%M)
RESULTS_DIR="dephasing_p05_results_${TIMESTAMP}"
mkdir -p ${RESULTS_DIR}

# Count total jobs
total_jobs=$(wc -l < params_dephasing_p05.txt)
echo "Total jobs submitted: $total_jobs (expected 680)"

# Find and copy result files
echo "Searching for dephasing_p05_*.json files..."
find . -name "dephasing_p05_L8_*.json" -o -name "dephasing_p05_L12_*.json" -o -name "dephasing_p05_L16_*.json" -o -name "dephasing_p05_L20_*.json" | while read file; do cp "$file" ${RESULTS_DIR}/; done

# Count files
NUM_FILES=$(ls -1 ${RESULTS_DIR}/*.json 2>/dev/null | wc -l)
echo "Collected ${NUM_FILES} result files"

if [ ${NUM_FILES} -eq 0 ]; then
    echo "WARNING: No result files found!"
    echo "Expected 680 files (4 L values × 17 λ values × 10 samples)"
else
    echo "Expected: 680 files (L=8,12,16,20 × 17 λ values × 10 samples)"
    
    # Check by system size
    for L in 8 12 16 20; do
        count=$(ls -1 ${RESULTS_DIR}/dephasing_p05_L${L}_*.json 2>/dev/null | wc -l)
        echo "  L=${L}: ${count} files (expected 170)"
    done
fi

# Create tarball
TARBALL="${RESULTS_DIR}.tar.gz"
echo ""
echo "Creating tarball: ${TARBALL}"
tar -czf ${TARBALL} ${RESULTS_DIR}/

# Get tarball size
SIZE=$(ls -lh ${TARBALL} | awk '{print $5}')
echo "Tarball size: ${SIZE}"

echo ""
echo "================================================"
echo "Collection complete!"
echo "================================================"
echo "Transfer with: wormhole send ${TARBALL}"
echo ""
