#!/bin/bash

# HTCondor job script for learning-to-trivial transition scan (λ_x = 0.49, λ_zz = 0.21)
# Density matrix evolution with Rényi-2 Binder
# Time limit: 24 hours per trial

# Get individual arguments from HTCondor
L=$1
lambda_x=$2
lambda_zz=$3
P_x=$4
P_zz=$5
ntrials=$6
seed=$7
sample=$8
out_prefix=$9

echo "=== Learning-to-Trivial Transition Scan Job Start (λ_x = 0.49, λ_zz = 0.21) ==="
echo "Job started at: $(date)"
echo "Running on: $(hostname)"
echo "Parameters: L=$L lambda_x=$lambda_x lambda_zz=$lambda_zz P_x=$P_x P_zz=$P_zz ntrials=$ntrials seed=$seed sample=$sample"
echo "Working directory: $(pwd)"

# List available files
echo "Available files:"
ls -la

# Create output directory
mkdir -p output

# Check Julia version
echo "Julia version:"
julia --version

# Set threading environment variables to use all 4 CPUs
export JULIA_NUM_THREADS=4
export OPENBLAS_NUM_THREADS=4
export MKL_NUM_THREADS=4
export BLAS_NUM_THREADS=4

echo "Threading enabled: JULIA_NUM_THREADS=$JULIA_NUM_THREADS"

# Install packages
echo "Setting up Julia environment..."
julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'

# Run learning-to-trivial transition scan simulation
echo "Running learning-to-trivial transition scan with Rényi-2 Binder..."
echo "Command: julia run_learning_to_trivial_lambda0.7_scan.jl $L $lambda_x $lambda_zz $P_x $P_zz $ntrials $seed $sample $out_prefix"

julia --project=. run_learning_to_trivial_lambda0.7_scan.jl \
    $L \
    $lambda_x \
    $lambda_zz \
    $P_x \
    $P_zz \
    $ntrials \
    $seed \
    $sample \
    $out_prefix

EXIT_CODE=$?

echo "Simulation completed with exit code: $EXIT_CODE"

# List output files
echo "Output files created:"
ls -la output/ 2>/dev/null || echo "No output directory found"

echo "Job finished at: $(date)"
echo "=== Learning-to-Trivial Transition Scan Job End ==="

exit $EXIT_CODE
