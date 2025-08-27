#!/bin/bash
# Complete workflow: Download results and create plots

echo "Binder Parameter Analysis Workflow"
echo ""

# Step 1: Download results
echo "Downloading results from cluster..."
./download_results.sh

echo ""
echo "Checking download status..."
result_count=$(ls binder-simulation-results/output/*.json 2>/dev/null | wc -l)
echo "Downloaded $result_count result files"

if [ $result_count -eq 0 ]; then
    echo "No result files found. Make sure your cluster jobs have completed."
    exit 1
fi

echo ""
echo "Creating plots and analysis..."

# Add Plots to the environment if not already there
julia -e "using Pkg; Pkg.add(\"Plots\")" 2>/dev/null || true

# Run the analysis
julia analyze_and_plot.jl

echo ""
echo "Analysis complete."
