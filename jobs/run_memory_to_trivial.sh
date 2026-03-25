#!/bin/bash

# HTCondor job script for memory-to-trivial transition scan
# Density matrix evolution with dephasing channels

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

echo "=== Memory-to-Trivial Transition Scan Job Start ==="
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

# Run memory-to-trivial transition scan simulation
echo "Running memory-to-trivial transition scan with density matrix evolution..."
echo "Command: julia --project=. run_memory_to_trivial_scan.jl --L $L --lambda_x $lambda_x --lambda_zz $lambda_zz --P_min $P_x --P_max $P_x --P_steps 1 --P_zz_mode fixed --P_zz_value $P_zz --ntrials $ntrials --seed $seed --output_dir output --output_file ${out_prefix}.json"

julia --project=. run_memory_to_trivial_scan.jl \
    --L $L \
    --lambda_x $lambda_x \
    --lambda_zz $lambda_zz \
    --P_min $P_x \
    --P_max $P_x \
    --P_steps 1 \
    --P_zz_mode fixed \
    --P_zz_value $P_zz \
    --ntrials $ntrials \
    --seed $seed \
    --maxdim 512 \
    --output_dir output \
    --output_file ${out_prefix}.json

EXIT_CODE=$?

echo "Simulation completed with exit code: $EXIT_CODE"

# List output files
echo "Output files created:"
ls -la output/ 2>/dev/null || echo "No output directory found"

echo "Job finished at: $(date)"
echo "=== Memory-to-Trivial Transition Scan Job End ==="

exit $EXIT_CODE
