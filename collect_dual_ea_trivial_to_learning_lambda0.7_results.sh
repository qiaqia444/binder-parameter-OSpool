#!/bin/bash

# Collect and organize trivial-to-learning (lambda = 0.7) dual EA Binder scan results
# Run this script after all HTCondor jobs complete

echo "=== Trivial-to-Learning (lambda=0.7) dual EA Scan Results Collection ==="
echo "Starting collection at: $(date)"

# Create timestamped results directory
TIMESTAMP=$(date +%Y%m%d_%H%M)
RESULTS_DIR="dual_ea_trivial_to_learning_lambda0.7_results_${TIMESTAMP}"

echo "Creating results directory: $RESULTS_DIR"
mkdir -p "$RESULTS_DIR"

SIZES="16 24 32 40"

# Create subdirectories for each system size
for L in $SIZES; do
    mkdir -p "$RESULTS_DIR/L${L}"
done

# Navigate to jobs directory
cd jobs

# Check if output directory exists
if [ ! -d "output" ]; then
    echo "ERROR: No output directory found. Jobs may not have completed yet."
    exit 1
fi

echo "Found output directory with $(ls output/dual_ea_trivial_to_learning_lambda0.7_*.json 2>/dev/null | wc -l) dual-EA result files"

# Organize by system size. Only dual_ea_trivial_to_learning_lambda0.7_* files with
# lx0.49 / lzz0.21 are matched, so Renyi-2 learning_to_trivial_* outputs are never mixed in.
echo "Collecting ALL lambda=0.7 (λ_x=0.49, λ_zz=0.21) dual-EA files..."
for L in $SIZES; do
    find output -name "dual_ea_trivial_to_learning_lambda0.7_L${L}_lx0.49_lzz0.21_*.json" ! -name "*FAILED*" -exec cp {} "../${RESULTS_DIR}/L${L}/" \; 2>/dev/null
    count=$(ls "../${RESULTS_DIR}/L${L}/"*.json 2>/dev/null | wc -l)
    echo "  L=$L: $count files"
done

# Count total results
cd ..
total_count=$(find "$RESULTS_DIR" -name "*.json" | wc -l)
echo "Total results collected: $total_count files"

# Check for failures
cd jobs
failure_count=$(ls output/dual_ea_trivial_to_learning_lambda0.7_*_FAILED.json 2>/dev/null | wc -l)
if [ $failure_count -gt 0 ]; then
    echo "WARNING: Found $failure_count failed jobs"
    mkdir -p "../${RESULTS_DIR}/failed"
    cp output/dual_ea_trivial_to_learning_lambda0.7_*_FAILED.json "../${RESULTS_DIR}/failed/" 2>/dev/null
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
echo "  julia analyze_dual_ea_trivial_to_learning_lambda0.7.jl"
echo ""
echo "Collection completed at: $(date)"
