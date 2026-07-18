#!/bin/bash

# HTCondor job script for lambda-path scan (Rényi-2 / Edwards-Anderson
# overlap susceptibilities, ITensorCorrelators.jl)
# Sweeps lambda_x = delta*lambda, lambda_zz = delta*(1-lambda) at fixed q

# Get individual arguments from HTCondor
L=$1
lambda=$2
delta=$3
q=$4
ntrials=$5
seed=$6
sample=$7
out_prefix=$8

echo "=== Lambda-Path Susceptibility Scan Job Start (ITensorCorrelators.jl) ==="
echo "Job started at: $(date)"
echo "Running on: $(hostname)"
echo "Parameters: L=$L lambda=$lambda delta=$delta q=$q ntrials=$ntrials seed=$seed sample=$sample"
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

# Run lambda-path scan simulation (single lambda point per job)
echo "Running lambda-path scan (doubled-MPS Born sampling, ITensorCorrelators.jl)..."
echo "Command: julia --project=. run_lambda_scan_susceptibility_scan.jl --L $L --lambda_min $lambda --lambda_max $lambda --lambda_steps 1 --delta $delta --q $q --ntrials $ntrials --seed $seed --output_dir output --output_file ${out_prefix}.json"

julia --project=. run_lambda_scan_susceptibility_scan.jl \
    --L $L \
    --lambda_min $lambda \
    --lambda_max $lambda \
    --lambda_steps 1 \
    --delta $delta \
    --q $q \
    --ntrials $ntrials \
    --seed $seed \
    --output_dir output \
    --output_file "${out_prefix}.json"

exit_code=$?

echo "Job completed with exit code: $exit_code"
echo "Job ended at: $(date)"

if [ $exit_code -ne 0 ]; then
    echo "ERROR: Job failed with exit code $exit_code"
    # Create FAILED file marker
    touch output/${out_prefix}_FAILED.json
fi

exit $exit_code
