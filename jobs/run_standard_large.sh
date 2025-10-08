#!/bin/bash

# HTCondor job script for LARGE standard Binder parameter calculation
# This script runs standard quantum trajectory simulations for large system sizes (L=20,24,28,32,36)

# Get individual arguments from HTCondor
L=$1
lambda_x=$2
lambda_zz=$3
lambda=$4
ntrials=$5
seed=$6
sample=$7
out_prefix=$8

echo "=== Large Standard Binder Parameter Job Start ==="
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

# Run simulation for large system
echo "Running large standard simulation (L=$L)..."
echo "Command: julia --project=. run.jl --L $L --lambda_x $lambda_x --lambda_zz $lambda_zz --lambda $lambda --ntrials $ntrials --seed $seed --sample $sample --out_prefix $out_prefix --outdir output"

# Set optimal threading for large systems
export JULIA_NUM_THREADS=4
export OPENBLAS_NUM_THREADS=4
export MKL_NUM_THREADS=4

# Run with optimized settings for large L
julia --project=. -t 4 run.jl --L $L --lambda_x $lambda_x --lambda_zz $lambda_zz --lambda $lambda --ntrials $ntrials --seed $seed --sample $sample --out_prefix $out_prefix --outdir output
EXIT_CODE=$?

echo "Simulation completed with exit code: $EXIT_CODE"

# List output files
echo "Output files created:"
ls -la output/ 2>/dev/null || echo "No output directory found"

echo "Job finished at: $(date)"
echo "=== Large Standard Binder Parameter Job End ==="

exit $EXIT_CODE