#!/bin/bash

# Collect and organize memory-to-trivial transition scan results (λ_x = 0.21, λ_zz = 0.49)
# Run this script after all HTCondor jobs complete

echo "=== Memory-to-Trivial Transition Scan Results Collection (λ_x = 0.21, λ_zz = 0.49) ==="
echo "Starting collection at: $(date)"

# Create timestamped results directory
TIMESTAMP=$(date +%Y%m%d_%H%M)
RESULTS_DIR="memory_to_trivial_lx0.3_results_${TIMESTAMP}"

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
EXPECTED_LAMBDA_X="0.21"
EXPECTED_LAMBDA_ZZ="0.49"
SAMPLE_FILE=$(ls output/memory_to_trivial_L*_lx0.21_lzz0.49_*.json 2>/dev/null | head -1)

if [ -z "$SAMPLE_FILE" ]; then
    echo "ERROR: No files found matching pattern memory_to_trivial_L*_lx0.21_lzz0.49_*.json"
    echo "Please check that jobs have completed and output files exist."
    exit 1
fi

echo "✓ Found output files with corrected parameters (λ_x=$EXPECTED_LAMBDA_X, λ_zz=$EXPECTED_LAMBDA_ZZ)"

# Organize by system size (collect ALL memory_to_trivial_L_P files, no time filter)
echo "Collecting ALL memory_to_trivial_L*_P*.json files..."
for L in 8 16 24 32; do
    # Collect all files matching pattern: memory_to_trivial_L{L}_P*.json (no -mtime filter)
    find output -name "memory_to_trivial_L${L}_P*.json" ! -name "*FAILED*" -exec cp {} "../${RESULTS_DIR}/L${L}/" \; 2>/dev/null
    count=$(ls "../${RESULTS_DIR}/L${L}/"*.json 2>/dev/null | wc -l)
    if [ $count -gt 0 ]; then
        echo "  L=$L: $count files collected"
    else
        echo "  L=$L: 0 files (jobs may still be running or not yet started)"
    fi
done

# Count total results
cd ..
total_count=$(find "$RESULTS_DIR" -name "*.json" | wc -l)
echo "Total results collected: $total_count files"

# Check for failures
cd jobs
failure_files=$(ls output/memory_to_trivial_*_FAILED.json 2>/dev/null)
failure_count=$(echo "$failure_files" | grep -c "FAILED" 2>/dev/null || echo "0")
if [ $failure_count -gt 0 ]; then
    echo "WARNING: Found $failure_count failed jobs"
    mkdir -p "../${RESULTS_DIR}/failed"
    cp output/*_FAILED.json "../${RESULTS_DIR}/failed/" 2>/dev/null
    echo "Failure files copied to ${RESULTS_DIR}/failed/"
fi
cd ..

# Validate collected data
echo ""
echo "=== Data Validation ==="
validation_errors=0

# Check if we have the expected parameter values in collected files
if [ $total_count -gt 0 ]; then
    # Sample check: verify first file has correct parameters
    sample=$(find "$RESULTS_DIR" -name "*.json" | head -1)
    if [ ! -z "$sample" ]; then
        # Extract job name from filename to verify parameters
        filename=$(basename "$sample")
        if [[ $filename == *"lx0.21_lzz0.49"* ]]; then
            echo "✓ Collected files contain correct parameter set (λ_x=0.21, λ_zz=0.49)"
        else
            echo "✗ WARNING: Filename pattern mismatch in $filename"
            validation_errors=$((validation_errors + 1))
        fi
    fi
fi

# Check if collection is complete
expected_total=3200  # Should have 3200 parameter sets (4 L × 20 P × 40 samples)
if [ $total_count -lt $expected_total ]; then
    missing=$((expected_total - total_count))
    echo "⏳ Incomplete collection: $missing/$expected_total files still missing"
    echo "   (Jobs are still running. 24-hour time limit per trial)"
    validation_errors=$((validation_errors + 1))
else
    echo "✓ Collection complete: All $total_count files collected"
fi

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
echo "Total files collected: $total_count / 3200 (expected)"
echo "Failed jobs: $failure_count"
echo "Parameters: λ_x = 0.21 (X measurements), λ_zz = 0.49 (ZZ measurements)"
echo "Timeout per trial: 24 hours (86400 seconds)"
echo ""

# Per-system-size breakdown
echo "=== Files by System Size ==="
for L in 8 16 24 32; do
    count=$(find "$RESULTS_DIR/L${L}" -name "*.json" 2>/dev/null | wc -l)
    printf "  L=%-2d: %3d files\n" $L $count
done
echo ""

if [ $validation_errors -eq 0 ]; then
    echo "✓ Data validation passed!"
else
    echo "⚠ Data validation warnings detected (see above)"
fi
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
echo "  # Run analysis scripts on the results"
echo ""
echo "=== Results Collection Complete ==="
