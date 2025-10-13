#!/bin/bash

# HTCondor job script for CRITICAL scaling analysis Binder parameter calculation
# This script runs standard quantum trajectory simulations for critical lambda=0.5 jobs (L=24,28,32,36)
# Uses enhanced resources: 72h runtime, 16GB memory, 20GB disk

# Get individual arguments from HTCondor
L=$1
lambda_x=$2
lambda_zz=$3
lambda=$4
ntrials=$5
seed=$6
sample=$7
out_prefix=$8

echo "=== Critical Scaling Binder Parameter Job Start ==="
echo "Job started at: $(date)"
echo "Running on: $(hostname)"
echo "Parameters: L=$L lambda_x=$lambda_x lambda_zz=$lambda_zz lambda=$lambda ntrials=$ntrials seed=$seed sample=$sample out_prefix=$out_prefix"
echo "Working directory: $(pwd)"

# List available files
echo "Available files:"
ls -la

# Create output directory
mkdir -p output

# Check Julia version
echo "Julia version:"
julia --version

# Install packages
echo "Setting up Julia environment..."
julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'

# Run simulation for critical scaling analysis
echo "Running critical scaling simulation (L=$L, lambda=$lambda)..."
echo "Command: julia --project=. run.jl --L $L --lambda_x $lambda_x --lambda_zz $lambda_zz --lambda $lambda --ntrials $ntrials --seed $seed --sample $sample --out_prefix $out_prefix --outdir output"

# Set optimal threading and memory for large systems with enhanced resources
export JULIA_NUM_THREADS=4
export OPENBLAS_NUM_THREADS=4
export MKL_NUM_THREADS=4
export JULIA_MAX_NUM_PRECOMPILE_FILES=10  # Reduce memory usage

# Add progress reporting for large L with critical lambda=0.5
echo "Starting critical scaling simulation: L=$L at lambda=$lambda (expected runtime with enhanced resources: $([ $L -le 24 ] && echo "4-12 hours" || [ $L -le 32 ] && echo "12-48 hours" || echo "24-72 hours"))"

# Run with optimized settings for large L with reduced chunk size for L>=28
if [ $L -ge 28 ]; then
    julia --project=. -t 4 run.jl --L $L --lambda_x $lambda_x --lambda_zz $lambda_zz --lambda $lambda --ntrials $ntrials --seed $seed --sample $sample --out_prefix $out_prefix --outdir output --chunk4 10000
else
    julia --project=. -t 4 run.jl --L $L --lambda_x $lambda_x --lambda_zz $lambda_zz --lambda $lambda --ntrials $ntrials --seed $seed --sample $sample --out_prefix $out_prefix --outdir output
fi
EXIT_CODE=$?

echo "Simulation completed with exit code: $EXIT_CODE"

# List output files
echo "Output files created:"
ls -la output/ 2>/dev/null || echo "No output directory found"

echo "Job finished at: $(date)"
echo "=== Critical Scaling Binder Parameter Job End ==="

exit $EXIT_CODE