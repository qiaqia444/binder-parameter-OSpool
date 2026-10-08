#!/bin/bash

# HTCondor job script for the right boundary Renyi-2 Binder scan measured AFTER THE LAST X DEPHASING
# (before the last ZZ-dephasing layer)
# L=48,56 near-critical zoom: 2 CPUs, maxdim 512, cutoff 1e-13, obs_maxdim_factor 2, obs_cutoff 1e-14, T_max_factor 4
# (identical settings to run_right_boundary_L48_L56.sh)
# Same dynamics and seeds as run_right_boundary_L48_L56.sh; only the observation time differs.

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

echo "=== Right Boundary Renyi-2 (after last X dephasing) Scan Job Start ==="
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

# Set threading environment variables to match request_cpus=2
export JULIA_NUM_THREADS=2
export OPENBLAS_NUM_THREADS=2
export MKL_NUM_THREADS=2
export BLAS_NUM_THREADS=2

echo "Threading enabled: JULIA_NUM_THREADS=$JULIA_NUM_THREADS"

# Install packages
echo "Setting up Julia environment..."
julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'

echo "Running right boundary Renyi-2 scan (observable after the last X dephasing)..."
echo "Command: julia --project=. run_renyi2_right_boundary_xnoise_scan.jl --L $L --lambda_x $lambda_x --lambda_zz $lambda_zz --P_min $P_x --P_max $P_x --P_steps 1 --ntrials $ntrials --maxdim 512 --cutoff 1e-13 --obs_maxdim_factor 2 --obs_cutoff 1e-14 --T_max_factor 4 --seed $seed --output_dir output --output_file ${out_prefix}.json"

julia --project=. run_renyi2_right_boundary_xnoise_scan.jl \
    --L $L \
    --lambda_x $lambda_x \
    --lambda_zz $lambda_zz \
    --P_min $P_x \
    --P_max $P_x \
    --P_steps 1 \
    --ntrials $ntrials \
    --maxdim 512 \
    --cutoff 1e-13 \
    --obs_maxdim_factor 2 \
    --obs_cutoff 1e-14 \
    --T_max_factor 4 \
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
