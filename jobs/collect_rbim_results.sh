#!/bin/bash

# Script to collect RBIM simulation results
# Run this on the submit node after jobs complete

echo "Collecting RBIM simulation results..."
echo "Random Bond Ising Model (λ_x = 0)"
echo "Scanning strategy: P vs Binder at fixed λ_zz to find Nishimori point"

# Create results directory with timestamp
TIMESTAMP=$(date +%Y%m%d_%H%M)
RESULTS_DIR="rbim_results_${TIMESTAMP}"
mkdir -p ${RESULTS_DIR}

# Count total jobs
total_jobs=$(wc -l < params_rbim.txt)
echo "Total jobs submitted: $total_jobs (expected 3840)"

# Find and copy result files
echo "Searching for rbim_*.json files..."
find . -name "rbim_L*.json" | while read file; do cp "$file" ${RESULTS_DIR}/; done

# Count files
NUM_FILES=$(ls -1 ${RESULTS_DIR}/*.json 2>/dev/null | wc -l)
echo "Collected ${NUM_FILES} result files"

if [ ${NUM_FILES} -eq 0 ]; then
    echo "WARNING: No result files found!"
    echo "Expected 3840 files (8 L × 16 λ_zz × 3 P_x × 10 samples)"
else
    echo "Expected: 3840 files (8 L values × 16 λ_zz × 3 P_x × 10 samples)"
    
    # Check by dephasing strength
    for P in 0.2 0.5 0.8; do
        count=$(ls -1 ${RESULTS_DIR}/rbim_*_P${P}_*.json 2>/dev/null | wc -l)
        expected=$((8 * 16 * 10))
        echo "  P=$P: ${count} files (expected ${expected})"
    done
    
    echo ""
    # Check by system size
    for L in 8 12 16 20 24 28 32 36; do
        count=$(ls -1 ${RESULTS_DIR}/rbim_L${L}_*.json 2>/dev/null | wc -l)
        expected=$((16 * 3 * 10))
        echo "  L=${L}: ${count} files (expected ${expected})"
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
