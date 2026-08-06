#!/bin/bash

# HTCondor job script for right-boundary finite-time test

L=$1
lambda_x=$2
lambda_zz=$3
P_x=$4
P_zz=$5
T_max_factor=$6
ntrials=$7
seed=$8
sample=$9
out_prefix=${10}

echo "=== Right Boundary T_max Test Job Start ==="
echo "Job started at: $(date)"
echo "Running on: $(hostname)"
echo "Parameters: L=$L lambda_x=$lambda_x lambda_zz=$lambda_zz P_x=$P_x P_zz=$P_zz T_max_factor=$T_max_factor ntrials=$ntrials seed=$seed sample=$sample"
echo "Working directory: $(pwd)"

mkdir -p output

echo "Julia version:"
julia --version

export JULIA_NUM_THREADS=4
export OPENBLAS_NUM_THREADS=4
export MKL_NUM_THREADS=4
export BLAS_NUM_THREADS=4

echo "Threading enabled: JULIA_NUM_THREADS=$JULIA_NUM_THREADS"

echo "Setting up Julia environment..."
julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'

echo "Running finite-time right-boundary scan point..."

julia --project=. run_right_boundary_tmax_scan.jl \
    --L $L \
    --lambda_x $lambda_x \
    --lambda_zz $lambda_zz \
    --P_min $P_x \
    --P_max $P_x \
    --P_steps 1 \
    --T_max_factor $T_max_factor \
    --ntrials $ntrials \
    --seed $seed \
    --output_dir output \
    --output_file "${out_prefix}.json"

exit_code=$?

echo "Job completed with exit code: $exit_code"
echo "Job ended at: $(date)"

if [ $exit_code -ne 0 ]; then
    echo "ERROR: Job failed with exit code $exit_code"
    touch output/${out_prefix}_FAILED.json
fi

exit $exit_code
