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
echo "Searching for standard result files..."
echo "Looking for pattern: standard_L*_lam*_s*.json"

# Check multiple possible output locations
SEARCH_DIRS="output jobs/output . ./output ./jobs/output"

FOUND_FILES=""
for dir in $SEARCH_DIRS; do
    if [ -d "$dir" ]; then
        echo "Searching in: $dir"
        FILES_IN_DIR=$(find "$dir" -name "standard_L*_lam*_s*.json" -type f 2>/dev/null)
        if [ ! -z "$FILES_IN_DIR" ]; then
            FOUND_FILES="$FOUND_FILES $FILES_IN_DIR"
            echo "  Found $(echo $FILES_IN_DIR | wc -w) files in $dir"
        fi
    fi
done

if [ -z "$FOUND_FILES" ]; then
    echo "No files found with pattern standard_L*_lam*_s*.json"
    echo "Searching for any standard*.json files in all directories..."
    for dir in $SEARCH_DIRS; do
        if [ -d "$dir" ]; then
            FILES_IN_DIR=$(find "$dir" -name "standard*.json" -type f 2>/dev/null)
            if [ ! -z "$FILES_IN_DIR" ]; then
                FOUND_FILES="$FOUND_FILES $FILES_IN_DIR"
                echo "  Found standard*.json files in $dir"
            fi
        fi
    done
fi

if [ -z "$FOUND_FILES" ]; then
    echo "No standard result files found at all!"
    echo ""
    echo "Directory diagnostics:"
    for dir in $SEARCH_DIRS; do
        if [ -d "$dir" ]; then
            echo "Contents of $dir:"
            ls -la "$dir" | head -5
            echo ""
        else
            echo "$dir: directory not found"
        fi
    done
else
    echo ""
    echo "Found files:"
    for file in $FOUND_FILES; do
        echo "  $file"
    done
    
    echo ""
    echo "Copying files..."
    for file in $FOUND_FILES; do
        if [ -f "$file" ]; then
            cp "$file" "$RESULTS_DIR/"
            echo "  Copied: $file"
        fi
    done
fi

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
$(find "$RESULTS_DIR" -name "standard*.json" | wc -l) standard result files
$(find "$RESULTS_DIR" -name "*.json" | grep -v summary | wc -l) total JSON files

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