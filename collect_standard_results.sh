#!/bin/bash

# Script to collect and package standard job results for transfer
# Run this on the cluster after standard jobs complete

echo "=== Standard Job Results Collection ==="
echo "Started at: $(date)"

# Create results directory with timestamp
TIMESTAMP=$(date +"%Y%m%d_%H%M")
RESULTS_DIR="standard_results_${TIMESTAMP}"
mkdir -p $RESULTS_DIR

echo "Collecting results into: $RESULTS_DIR"

# Copy all standard job output files with correct naming pattern
echo "Copying JSON result files..."
find output -name "standard_L*_lam*_s*.json" -type f | while read file; do
    if [ -f "$file" ]; then
        cp "$file" "$RESULTS_DIR/"
        echo "  Copied: $file"
    fi
done

# Also check for any other standard result patterns as backup
find output -name "standard_*.json" -type f | while read file; do
    if [ -f "$file" ]; then
        cp "$file" "$RESULTS_DIR/"
        echo "  Copied: $file"
    fi
done

# Copy parameter files for reference
echo "Copying parameter files..."
cp jobs/params.txt "$RESULTS_DIR/" 2>/dev/null || echo "  Warning: params.txt not found"
cp jobs/jobs.submit "$RESULTS_DIR/standard_jobs.submit" 2>/dev/null || echo "  Warning: jobs.submit not found"

# Copy a few sample log files for debugging if needed
echo "Copying sample log files..."
mkdir -p "$RESULTS_DIR/sample_logs"
find logs -name "*standard*.out" -type f | head -5 | while read file; do
    if [ -f "$file" ]; then
        cp "$file" "$RESULTS_DIR/sample_logs/"
    fi
done

# Create summary of collected files
echo "Creating collection summary..."
cat > "$RESULTS_DIR/collection_summary.txt" << EOF
Standard Job Results Collection Summary
=====================================
Collection Date: $(date)
Collection Directory: $RESULTS_DIR

Files Collected:
$(find "$RESULTS_DIR" -type f | wc -l) total files

JSON Results:
$(find "$RESULTS_DIR" -name "*.json" | grep -v summary | wc -l) result files

File List:
$(find "$RESULTS_DIR" -type f | sort)

Disk Usage:
$(du -sh "$RESULTS_DIR")
EOF

# Display summary
echo ""
echo "=== Collection Summary ==="
cat "$RESULTS_DIR/collection_summary.txt"

# Create compressed archive
ARCHIVE_NAME="${RESULTS_DIR}.tar.gz"
echo ""
echo "Creating compressed archive: $ARCHIVE_NAME"
tar -czf "$ARCHIVE_NAME" "$RESULTS_DIR"

echo ""
echo "=== Ready for Transfer ==="
echo "Archive created: $ARCHIVE_NAME"
echo "Archive size: $(ls -lh $ARCHIVE_NAME | awk '{print $5}')"
echo ""
echo "To transfer to your Mac using magic-wormhole:"
echo "1. On cluster: wormhole send $ARCHIVE_NAME"
echo "2. On Mac: wormhole receive"
echo ""
echo "Collection completed at: $(date)"