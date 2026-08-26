#!/bin/bash

# Collect and organize right-boundary defect-insertion results.
# Run this script after all HTCondor jobs complete.

echo "=== Right Boundary Defect Insertion Results Collection ==="
echo "Starting collection at: $(date)"

TIMESTAMP=$(date +%Y%m%d_%H%M)
RESULTS_DIR="defect_insertion_results_${TIMESTAMP}"
echo "Creating results directory: $RESULTS_DIR"
mkdir -p "$RESULTS_DIR"

for L in 16 24 32; do
	mkdir -p "$RESULTS_DIR/L${L}"
done

cd jobs
if [ ! -d "output" ]; then
	echo "ERROR: No output directory found. Jobs may not have completed yet."
	exit 1
fi

echo "Found $(find output -maxdepth 1 -name '*.json' | wc -l) result files"
for L in 16 24 32; do
	find output -maxdepth 1 \
		-name "defect_insertion_L${L}_*.json" \
		! -name "*_FAILED.json" \
		-exec cp {} "../${RESULTS_DIR}/L${L}/" \; 2>/dev/null
	count=$(find "../${RESULTS_DIR}/L${L}" -name '*.json' | wc -l)
	echo "  L=$L: $count files"
done

cd ..
total_count=$(find "$RESULTS_DIR" -name '*.json' | wc -l)
echo "Total results collected: $total_count files"

cd jobs
failure_count=$(find output -maxdepth 1 -name 'defect_insertion_*_FAILED.json' | wc -l)
if [ "$failure_count" -gt 0 ]; then
	echo "WARNING: Found $failure_count failed jobs"
	mkdir -p "../${RESULTS_DIR}/failed"
	cp output/defect_insertion_*_FAILED.json "../${RESULTS_DIR}/failed/" 2>/dev/null || true
else
	echo "No failed jobs detected"
fi
cd ..

echo "Creating compressed archive..."
tar -czf "${RESULTS_DIR}.tar.gz" "$RESULTS_DIR"
ARCHIVE_SIZE=$(du -h "${RESULTS_DIR}.tar.gz" | cut -f1)

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
echo "  julia defect_insertion_analyze.jl ${RESULTS_DIR}/L16/*.json ${RESULTS_DIR}/L24/*.json ${RESULTS_DIR}/L32/*.json"
echo ""
echo "Collection completed at: $(date)"
