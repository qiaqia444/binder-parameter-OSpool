#!/bin/bash

# HTCondor worker for the standalone right-boundary defect diagnostic.

L=$1
lambda_x=$2
lambda_zz=$3
P_x=$4
P_zz=$5
ntrials=$6
seed=$7
sample=$8
out_prefix=$9

echo "=== Right Boundary Defect Insertion Job Start ==="
echo "Job started at: $(date)"
echo "Running on: $(hostname)"
echo "Parameters: L=$L lambda_x=$lambda_x lambda_zz=$lambda_zz P_x=$P_x P_zz=$P_zz ntrials=$ntrials seed=$seed sample=$sample"
echo "Working directory: $(pwd)"

echo "Available files:"
ls -la

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

echo "Running right-boundary defect-insertion diagnostic..."
echo "Command: julia --project=. defect_insertion_run.jl --L $L --lambda_x $lambda_x --lambda_zz $lambda_zz --P $P_x --ntrials $ntrials --seed $seed --output_dir output --output_file ${out_prefix}.csv"

julia --project=. defect_insertion_run.jl \
    --L "$L" \
    --lambda_x "$lambda_x" \
    --lambda_zz "$lambda_zz" \
    --P "$P_x" \
    --ntrials "$ntrials" \
    --seed "$seed" \
    --output_dir output \
    --output_file "${out_prefix}.csv"

exit_code=$?
echo "Job completed with exit code: $exit_code"
echo "Job ended at: $(date)"

if [ "$exit_code" -ne 0 ]; then
    echo "ERROR: Job failed with exit code $exit_code"
    touch "output/${out_prefix}_FAILED.json"
fi

exit "$exit_code"
