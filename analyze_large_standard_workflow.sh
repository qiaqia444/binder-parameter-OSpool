#!/bin/bash

# Complete workflow for standard Binder parameter analysis
# Run this script on your Mac after receiving data from cluster

set -e

echo "=================================================="
echo "STANDARD BINDER PARAMETER ANALYSIS WORKFLOW"
echo "=================================================="
echo "Started at: $(date)"
echo ""

# Check if we have an archive file as argument
if [ $# -eq 0 ]; then
    echo "Usage: $0 <standard_results_archive.tar.gz>"
    echo "Example: $0 standard_results_20251003_1500.tar.gz"
    echo ""
    echo "Steps to run this workflow:"
    echo "1. On cluster: ./collect_large_standard_results.sh"
    echo "2. Transfer: wormhole send standard_results_YYYYMMDD_HHMM.tar.gz"
    echo "3. On Mac: wormhole receive"
    echo "4. On Mac: ./analyze_large_standard_workflow.sh standard_results_YYYYMMDD_HHMM.tar.gz"
    exit 1
fi

ARCHIVE_FILE="$1"

# Check if archive exists
if [ ! -f "$ARCHIVE_FILE" ]; then
    echo "Error: Archive file '$ARCHIVE_FILE' not found!"
    exit 1
fi

echo "Processing archive: $ARCHIVE_FILE"

# Extract archive
echo "Extracting archive..."
tar -xzf "$ARCHIVE_FILE"

# Find the extracted directory
RESULTS_DIR=$(basename "$ARCHIVE_FILE" .tar.gz)

if [ ! -d "$RESULTS_DIR" ]; then
    echo "Error: Expected directory '$RESULTS_DIR' not found after extraction!"
    exit 1
fi

echo "Extracted to directory: $RESULTS_DIR"

# Count results with correct naming pattern
JSON_COUNT=$(find "$RESULTS_DIR" -name "standard_L*_lam*_s*.json" 2>/dev/null | wc -l)
if [ $JSON_COUNT -eq 0 ]; then
    # Try backup pattern
    JSON_COUNT=$(find "$RESULTS_DIR" -name "standard_*.json" | grep -v FAILED | wc -l)
fi

echo "Found $JSON_COUNT result files"

if [ $JSON_COUNT -eq 0 ]; then
    echo "Error: No result files found in $RESULTS_DIR"
    exit 1
fi

# Run analysis
echo ""
echo "Running standard Binder parameter analysis..."
echo "Command: julia analyze_standard.jl $RESULTS_DIR"
echo ""

# Check if Julia is available at the known location
JULIA_PATH="/Applications/Julia-1.11.app/Contents/Resources/julia/bin/julia"

if [ -f "$JULIA_PATH" ]; then
    echo "Using Julia at: $JULIA_PATH"
    $JULIA_PATH analyze_standard.jl "$RESULTS_DIR"
else
    echo "Using system Julia..."
    julia analyze_standard.jl "$RESULTS_DIR"
fi

# List generated files
echo ""
echo "=================================================="
echo "ANALYSIS COMPLETE"
echo "=================================================="
echo "Generated files:"
ls -la *.png *.txt 2>/dev/null | grep -E "(standard|binder)" || echo "No analysis files generated"

echo ""
echo "Original archive: $ARCHIVE_FILE"
echo "Extracted directory: $RESULTS_DIR"
echo "Completed at: $(date)"
echo "=================================================="