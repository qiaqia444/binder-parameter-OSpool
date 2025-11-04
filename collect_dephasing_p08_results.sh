#!/bin/bash

# Script to collect and package dephasing P=0.8 job results for transfer
# Run this on the cluster after dephasing P=0.8 jobs complete

echo "=== Dephasing P=0.8 Job Results Collection ==="
echo "Started at: $(date)"

# Create results directory with timestamp
TIMESTAMP=$(date +"%Y%m%d_%H%M")
RESULTS_DIR="dephasing_p08_results_${TIMESTAMP}"
mkdir -p $RESULTS_DIR

echo "Collecting results into: $RESULTS_DIR"

# Copy all dephasing P=0.8 job output files with correct naming pattern
echo "Searching for dephasing P=0.8 result files..."
echo "Looking for pattern: dephasing_p08_L*_lam*_P*_s*.json"

# Check multiple possible output locations
SEARCH_DIRS="output jobs/output . ./output ./jobs/output"

FOUND_FILES=""
for dir in $SEARCH_DIRS; do
    if [ -d "$dir" ]; then
        echo "Searching in: $dir"
        FILES_IN_DIR=$(find "$dir" -name "dephasing_p08_L*_lam*_P*_s*.json" -type f 2>/dev/null)
        if [ ! -z "$FILES_IN_DIR" ]; then
            FOUND_FILES="$FOUND_FILES $FILES_IN_DIR"
            echo "  Found $(echo $FILES_IN_DIR | wc -w) files in $dir"
        fi
    fi
done

if [ -z "$FOUND_FILES" ]; then
    echo "No files found with pattern dephasing_p08_L*_lam*_P*_s*.json"
    echo "Searching for any dephasing_p08*.json files in all directories..."
    for dir in $SEARCH_DIRS; do
        if [ -d "$dir" ]; then
            FILES_IN_DIR=$(find "$dir" -name "dephasing_p08*.json" -type f 2>/dev/null)
            if [ ! -z "$FILES_IN_DIR" ]; then
                FOUND_FILES="$FOUND_FILES $FILES_IN_DIR"
                echo "  Found dephasing_p08*.json files in $dir"
            fi
        fi
    done
fi

if [ -z "$FOUND_FILES" ]; then
    echo "No dephasing P=0.8 result files found at all!"
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
cp jobs/params_dephasing_p08.txt "$RESULTS_DIR/" 2>/dev/null || echo "  Warning: params_dephasing_p08.txt not found"
cp jobs/dephasing_p08_jobs.submit "$RESULTS_DIR/" 2>/dev/null || echo "  Warning: dephasing_p08_jobs.submit not found"

# Copy a few sample log files for debugging if needed
echo "Copying sample log files..."
mkdir -p "$RESULTS_DIR/sample_logs"
find logs -name "*dephasing_p08*.out" -type f | head -5 | while read file; do
    if [ -f "$file" ]; then
        cp "$file" "$RESULTS_DIR/sample_logs/"
    fi
done

# Copy a few error logs if they exist
find logs -name "*dephasing_p08*.err" -type f -size +0 | head -5 | while read file; do
    if [ -f "$file" ]; then
        cp "$file" "$RESULTS_DIR/sample_logs/"
    fi
done

# Create summary of collected files
echo "Creating collection summary..."
cat > "$RESULTS_DIR/collection_summary.txt" << EOF
Dephasing P=0.8 Job Results Collection Summary
=====================================
Collection Date: $(date)
Collection Directory: $RESULTS_DIR
Dephasing Strength: P_x = P_zz = 0.8

Files Collected:
$(find "$RESULTS_DIR" -type f | wc -l) total files

JSON Results:
$(find "$RESULTS_DIR" -name "dephasing_p08*.json" | wc -l) dephasing P=0.8 result files
$(find "$RESULTS_DIR" -name "*.json" | grep -v summary | wc -l) total JSON files

Parameter Combinations:
Expected: 153 jobs (L=8,12,16 × 17 lambda values × 3 samples)
Collected: $(find "$RESULTS_DIR" -name "dephasing_p08*.json" | wc -l) result files

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
echo "2. On Mac: wormhole receive <code-from-step-1>"
echo ""
echo "Alternative transfer methods:"
echo "  scp: scp $ARCHIVE_NAME your-mac:~/Desktop/"
echo "  rsync: rsync -avz $ARCHIVE_NAME your-mac:~/Desktop/"
echo ""
echo "Collection completed at: $(date)"
