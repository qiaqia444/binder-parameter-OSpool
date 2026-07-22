#!/bin/bash

# Collect and organize learning-to-trivial transition scan results
# (ITensorCorrelators.jl variant, λ_x = 0.56, λ_zz = 0.14)
# Run this script after all HTCondor jobs complete

echo "=== Learning-to-Trivial Transition Scan Results Collection (ITensorCorrelators.jl, λ_x = 0.56, λ_zz = 0.14) ==="
echo "Starting collection at: $(date)"

# Create timestamped results directory
TIMESTAMP=$(date +%Y%m%d_%H%M)
RESULTS_DIR="learning_to_trivial_lambda0.8_itensorcorrelators_results_${TIMESTAMP}"

echo "Creating results directory: $RESULTS_DIR"
mkdir -p "$RESULTS_DIR"

# Create subdirectories for each system size
for L in 8 16 24 32; do
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

# Validate that we're collecting from correct parameter set
echo "Validating parameter set..."
EXPECTED_LAMBDA_X="0.56"
EXPECTED_LAMBDA_ZZ="0.14"
SAMPLE_FILE=$(ls output/learning_to_trivial_lambda0.8_itensorcorrelators_L*_lx0.56_lzz0.14_*.json 2>/dev/null | head -1)

if [ -z "$SAMPLE_FILE" ]; then
    echo "ERROR: No files found matching pattern learning_to_trivial_lambda0.8_itensorcorrelators_L*_lx0.56_lzz0.14_*.json"
    echo "Please check that jobs have completed and output files exist."
    exit 1
fi

echo "✓ Found output files with corrected parameters (λ_x=$EXPECTED_LAMBDA_X, λ_zz=$EXPECTED_LAMBDA_ZZ)"

# Organize by system size (collect ALL lx0.56 lzz0.14 files)
echo "Collecting ALL λ_x=$EXPECTED_LAMBDA_X, λ_zz=$EXPECTED_LAMBDA_ZZ files..."
for L in 8 16 24 32; do
    # Collect all files matching pattern (no -mtime filter)
    find output -name "learning_to_trivial_lambda0.8_itensorcorrelators_L${L}_lx0.56_lzz0.14_*.json" ! -name "*FAILED*" -exec cp {} "../${RESULTS_DIR}/L${L}/" \; 2>/dev/null
    count=$(ls "../${RESULTS_DIR}/L${L}/"*.json 2>/dev/null | wc -l)
    echo "  L=$L: $count files"
done

# Count total results
cd ..
total_count=$(find "$RESULTS_DIR" -name "*.json" | wc -l)
echo "Total results collected: $total_count files"

# Check for failures
cd jobs
failure_count=$(ls output/learning_to_trivial_lambda0.8_itensorcorrelators_*_FAILED.json 2>/dev/null | wc -l)
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

echo ""
echo "✓ Results collection complete!"
echo "   Results directory: $RESULTS_DIR"
echo "   Compressed: ${RESULTS_DIR}.tar.gz (${ARCHIVE_SIZE})"
echo ""
