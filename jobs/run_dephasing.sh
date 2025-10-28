#!/bin/bash

# HTCondor job script for dephasing channel Binder parameter calculation
# This script runs quantum trajectory simulations with weak measurements + dephasing

# Get individual arguments from HTCondor
L=$1
P_x=$2
lambda_x=$3
lambda_zz=$4
lambda=$5
ntrials=$6
seed=$7
sample=$8
out_prefix=$9

echo "=== Dephasing Channel Binder Parameter Job Start ==="
echo "Job started at: $(date)"
echo "Running on: $(hostname)"
echo "Parameters: L=$L P_x=$P_x lambda_x=$lambda_x lambda_zz=$lambda_zz lambda=$lambda ntrials=$ntrials seed=$seed sample=$sample out_prefix=$out_prefix"
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

# Run simulation
echo "Running dephasing channel simulation..."
echo "Command: julia --project=. run_dephasing.jl --L $L --P_x $P_x --lambda_x $lambda_x --lambda_zz $lambda_zz --lambda $lambda --ntrials $ntrials --seed $seed --sample $sample --out_prefix $out_prefix --outdir output"
julia --project=. run_dephasing.jl --L $L --P_x $P_x --lambda_x $lambda_x --lambda_zz $lambda_zz --lambda $lambda --ntrials $ntrials --seed $seed --sample $sample --out_prefix $out_prefix --outdir output
EXIT_CODE=$?

echo "Simulation completed with exit code: $EXIT_CODE"

# List output files
echo "Output files created:"
ls -la output/ 2>/dev/null || echo "No output directory found"

echo "Job finished at: $(date)"
echo "=== Dephasing Channel Binder Parameter Job End ==="

exit $EXIT_CODE
