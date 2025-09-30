#!/bin/bash

# HTCondor job script for bipartite entropy calculation
# This script runs a single bipartite entropy simulation

# Get individual arguments from HTCondor
L=$1
lambda_x=$2
lambda_zz=$3
lambda=$4
ntrials=$5
seed=$6
sample=$7
out_prefix=$8

echo "=== Bipartite Entropy Job Start ==="
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

# Run simulation
echo "Running bipartite entropy simulation..."
echo "Command: julia --project=. run_bipartite.jl --L $L --lambda_x $lambda_x --lambda_zz $lambda_zz --lambda $lambda --ntrials $ntrials --seed $seed --sample $sample --out_prefix $out_prefix --outdir output"
julia --project=. run_bipartite.jl --L $L --lambda_x $lambda_x --lambda_zz $lambda_zz --lambda $lambda --ntrials $ntrials --seed $seed --sample $sample --out_prefix $out_prefix --outdir output

exit_code=$?

echo "=== Job completed at: $(date) ==="
echo "Exit code: $exit_code"

# List output files
echo "Generated files:"
ls -la output/

exit $exit_code