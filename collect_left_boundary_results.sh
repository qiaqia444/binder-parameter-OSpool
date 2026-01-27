#!/bin/bash

# Collect and organize left boundary scan results
# Run this script after all HTCondor jobs complete

echo "=== Left Boundary Scan Results Collection ==="
echo "Starting collection at: $(date)"

# Create timestamped results directory
TIMESTAMP=$(date +%Y%m%d_%H%M)
RESULTS_DIR="left_boundary_results_${TIMESTAMP}"

echo "Creating results directory: $RESULTS_DIR"
mkdir -p "$RESULTS_DIR"

# Create subdirectories for each system size
for L in 8 10 12 14 16; do
    mkdir -p "$RESULTS_DIR/L${L}"
done

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
for L in 8 10 12 14 16; do
    cp output/left_boundary_L${L}_*.json "../${RESULTS_DIR}/L${L}/" 2>/dev/null
    count=$(ls "../${RESULTS_DIR}/L${L}/"*.json 2>/dev/null | wc -l)
    echo "  L=$L: $count files"
done

# Count total results
cd ..
total_count=$(find "$RESULTS_DIR" -name "*.json" | wc -l)
echo "Total results collected: $total_count files"

# Check for failures
cd jobs
failure_count=$(ls output/left_boundary_*_FAILED.json 2>/dev/null | wc -l)
if [ $failure_count -gt 0 ]; then
    echo "WARNING: Found $failure_count failed jobs"
    mkdir -p "../${RESULTS_DIR}/failed"
    cp output/*_FAILED.json "../${RESULTS_DIR}/failed/" 2>/dev/null
    echo "Failure files copied to ${RESULTS_DIR}/failed/"
fi
cd ..

# Create archive
echo "Creating compressed archive..."
tar -czf "${RESULTS_DIR}.tar.gz" "$RESULTS_DIR"
ARCHIVE_SIZE=$(du -h "${RESULTS_DIR}.tar.gz" | cut -f1)
echo "Archive created: ${RESULTS_DIR}.tar.gz (${ARCHIVE_SIZE})"

# Print summary
echo ""
echo "=== Collection Summary ==="
echo "Results directory: $RESULTS_DIR"
echo "Archive: ${RESULTS_DIR}.tar.gz"
echo "Archive size: $ARCHIVE_SIZE"
echo "Total files: $total_count"
echo "Failed jobs: $failure_count"
echo ""
echo "=== Transfer to Mac with Magic Wormhole ==="
echo "On cluster, run:"
echo "  wormhole send ${RESULTS_DIR}.tar.gz"
echo ""
echo "On your Mac, run:"
echo "  wormhole receive"
echo "  # Enter the wormhole code when prompted"
echo ""
echo "Then extract and analyze:"
echo "  tar -xzf ${RESULTS_DIR}.tar.gz"
echo "  julia analyze_left_boundary.jl"
echo ""
echo "Collection completed at: $(date)"
