#!/bin/bash

# Collect and organize right-boundary finite-time test results

echo "=== Right Boundary T_max Test Results Collection ==="
echo "Starting collection at: $(date)"

TIMESTAMP=$(date +%Y%m%d_%H%M)
RESULTS_DIR="right_boundary_tmax_results_${TIMESTAMP}"

echo "Creating results directory: $RESULTS_DIR"
mkdir -p "$RESULTS_DIR"

for L in 16 24 32; do
    for Tf in 4 8; do
        mkdir -p "$RESULTS_DIR/L${L}/Tf${Tf}"
    done
done

cd jobs

if [ ! -d "output" ]; then
    echo "ERROR: No output directory found. Jobs may not have completed yet."
    exit 1
fi

echo "Found output directory with $(ls output/*.json 2>/dev/null | wc -l) result files"

echo "Collecting right_boundary_tmax files..."
for L in 16 24 32; do
    for Tf in 4 8; do
        find output -name "right_boundary_tmax_L${L}_Tf${Tf}_lx0.70_lzz0.00_Px*.json" ! -name "*FAILED*" -exec cp {} "../${RESULTS_DIR}/L${L}/Tf${Tf}/" \; 2>/dev/null
        count=$(ls "../${RESULTS_DIR}/L${L}/Tf${Tf}/"*.json 2>/dev/null | wc -l)
        echo "  L=$L Tf=$Tf: $count files"
    done
done

cd ..
total_count=$(find "$RESULTS_DIR" -name "*.json" | wc -l)
echo "Total results collected: $total_count files"

cd jobs
failure_count=$(ls output/right_boundary_tmax_*_FAILED.json 2>/dev/null | wc -l)
if [ $failure_count -gt 0 ]; then
    echo "WARNING: Found $failure_count failed jobs"
    mkdir -p "../${RESULTS_DIR}/failed"
    cp output/*_FAILED.json "../${RESULTS_DIR}/failed/" 2>/dev/null
    echo "Failure files copied to ${RESULTS_DIR}/failed/"
else
    echo "No failed jobs detected"
fi

cd ..

echo "Creating compressed archive..."
tar -czf "${RESULTS_DIR}.tar.gz" "$RESULTS_DIR"
ARCHIVE_SIZE=$(du -h "${RESULTS_DIR}.tar.gz" | cut -f1)
echo "Archive created: ${RESULTS_DIR}.tar.gz (${ARCHIVE_SIZE})"

echo ""
echo "=== Collection Summary ==="
echo "Results directory: $RESULTS_DIR"
echo "Archive: ${RESULTS_DIR}.tar.gz"
echo "Archive size: $ARCHIVE_SIZE"
echo "Total files: $total_count"
echo "Failed jobs: $failure_count"
echo ""
echo "Then analyze:"
echo "  julia analyze_right_boundary_tmax.jl"
echo ""
echo "Collection completed at: $(date)"
