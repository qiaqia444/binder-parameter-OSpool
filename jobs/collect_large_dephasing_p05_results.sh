#!/bin/bash

# Collection script for Large Dephasing P=0.5 results (L=20-36)
# Run this on the cluster after jobs complete

# Get timestamp for unique directory name
TIMESTAMP=$(date +%Y%m%d_%H%M)
RESULTS_DIR="large_dephasing_p05_results_${TIMESTAMP}"

echo "================================================"
echo "Collecting Large Dephasing P=0.5 Results"
echo "================================================"
echo "Creating directory: ${RESULTS_DIR}"

# Create results directory
mkdir -p ${RESULTS_DIR}

# Find and copy all large dephasing P=0.5 result files
echo "Searching for large_dephasing_p05_*.json files..."
find output -name "large_dephasing_p05_*.json" -type f -exec cp {} ${RESULTS_DIR}/ \;

# Count files
NUM_FILES=$(ls -1 ${RESULTS_DIR}/*.json 2>/dev/null | wc -l)
echo "Collected ${NUM_FILES} result files"

if [ ${NUM_FILES} -eq 0 ]; then
    echo "WARNING: No result files found!"
    echo "Expected 225 files (5 L values × 15 λ values × 3 samples)"
else
    echo "Expected: 225 files (L=20,24,28,32,36 × 15 λ values × 3 samples)"
    
    # Check by system size
    for L in 20 24 28 32 36; do
        count=$(ls -1 ${RESULTS_DIR}/large_dephasing_p05_L${L}_*.json 2>/dev/null | wc -l)
        echo "  L=${L}: ${count} files (expected 45)"
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
