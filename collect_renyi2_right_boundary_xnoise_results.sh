#!/bin/bash

# Collect and organize right boundary Renyi-2 (after last X dephasing) scan results
# Run this script after all HTCondor jobs complete

echo "=== Right Boundary Renyi-2 (after last X dephasing) Results Collection ==="
echo "Starting collection at: $(date)"

# Create timestamped results directory
TIMESTAMP=$(date +%Y%m%d_%H%M)
RESULTS_DIR="renyi2_xnoise_right_boundary_results_${TIMESTAMP}"

echo "Creating results directory: $RESULTS_DIR"
mkdir -p "$RESULTS_DIR"

# Create subdirectories for each system size
for L in 8 16 24 32 40 48 56; do
    mkdir -p "$RESULTS_DIR/L${L}"
done

# Navigate to jobs directory
cd jobs

# Check if output directory exists
if [ ! -d "output" ]; then
    echo "ERROR: No output directory found. Jobs may not have completed yet."
    exit 1
fi

echo "Found output directory with $(ls output/renyi2_xnoise_right_boundary_*.json 2>/dev/null | wc -l) result files"

# Organize by system size (only files of THIS family; production Renyi-2 / dual-EA files are not matched)
echo "Collecting ALL λ_x=0.7, λ_zz=0.0 renyi2_xnoise files..."
for L in 8 16 24 32 40 48 56; do
    find output -name "renyi2_xnoise_right_boundary_L${L}_lx0.70_lzz0.00_*.json" ! -name "*FAILED*" -exec cp {} "../${RESULTS_DIR}/L${L}/" \; 2>/dev/null
    count=$(ls "../${RESULTS_DIR}/L${L}/"*.json 2>/dev/null | wc -l)
    echo "  L=$L: $count files"
done

# Count total results
cd ..
total_count=$(find "$RESULTS_DIR" -name "*.json" | wc -l)
echo "Total results collected: $total_count files"

# Check for failures
cd jobs
failure_count=$(ls output/renyi2_xnoise_right_boundary_*_FAILED.json 2>/dev/null | wc -l)
if [ $failure_count -gt 0 ]; then
    echo "WARNING: Found $failure_count failed jobs"
    mkdir -p "../${RESULTS_DIR}/failed"
    cp output/renyi2_xnoise_right_boundary_*_FAILED.json "../${RESULTS_DIR}/failed/" 2>/dev/null
    echo "Failure files copied to ${RESULTS_DIR}/failed/"
else
    echo "✓ No failed jobs detected"
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
echo "  julia --project=. analyze_renyi2_right_boundary_xnoise.jl ${RESULTS_DIR} [production_right_boundary_results_dir]"
echo ""
echo "Collection completed at: $(date)"
