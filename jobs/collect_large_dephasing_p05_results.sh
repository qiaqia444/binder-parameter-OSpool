#!/bin/bash

# Collection script for Large Dephasing P=0.5 results (L=20,24,28,32,36) CORRECTED
# Run this on the cluster after jobs complete
# CORRECTED VERSION: Proper quantum channel implementation

# Get timestamp for unique directory name
TIMESTAMP=$(date +%Y%m%d_%H%M)
RESULTS_DIR="large_dephasing_p05_results_${TIMESTAMP}"

echo "================================================"
echo "Collecting Large Dephasing P=0.5 Results (CORRECTED)"
echo "================================================"
echo "Creating directory: ${RESULTS_DIR}"

# Create results directory
mkdir -p ${RESULTS_DIR}

# Find and copy large dephasing P=0.5 result files (L=24,28,32,36 only)
echo "Searching for large_dephasing_p05_*.json files (L=24,28,32,36)..."
find . -name "large_dephasing_p05_L24_*.json" -o -name "large_dephasing_p05_L28_*.json" -o -name "large_dephasing_p05_L32_*.json" -o -name "large_dephasing_p05_L36_*.json" | while read file; do cp "$file" ${RESULTS_DIR}/; done

# Count files
NUM_FILES=$(ls -1 ${RESULTS_DIR}/*.json 2>/dev/null | wc -l)
echo "Collected ${NUM_FILES} result files"

if [ ${NUM_FILES} -eq 0 ]; then
    echo "WARNING: No result files found!"
    echo "Expected 1800 files (4 L values × 15 λ values × 30 samples)"
else
    echo "Expected: 1800 files (L=24,28,32,36 × 15 λ values × 30 samples)"
    echo "Note: L=20 is in small systems jobs"
    
    # Check by system size
    for L in 24 28 32 36; do
        count=$(ls -1 ${RESULTS_DIR}/large_dephasing_p05_L${L}_*.json 2>/dev/null | wc -l)
        echo "  L=${L}: ${count} files (expected 450)"
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
