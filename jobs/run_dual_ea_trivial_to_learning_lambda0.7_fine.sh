#!/bin/bash

# HTCondor job script for the trivial-to-learning dual EA Binder scan, FINE grid
# (lambda = 0.7: lambda_x = 0.49, lambda_zz = 0.21; both weak-measurement layers active).
# Same as run_dual_ea_trivial_to_learning_lambda0.7.sh, except only 2 CPUs are requested/used
# (request_cpus = 2 in jobs_dual_ea_trivial_to_learning_lambda0.7_fine.submit), so the thread
# counts below are 2 instead of 4. The driver and its defaults are unchanged.
# Reuses the learning-to-trivial doubled-MPS Born-sampling dynamics (lt_evolve_one_trial),
# observable = dual Edwards-Anderson Binder (open-chain endpoints a = 2..L)

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

echo "=== Trivial-to-Learning dual EA Scan (fine grid) Job Start ==="
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

# Run trivial-to-learning dual EA scan (P_x = P_zz = q). Driver defaults:
# maxdim=256, cutoff=1e-12, T_max_factor=4.
echo "Running trivial-to-learning dual EA scan with doubled-MPS Born sampling..."
echo "Command: julia --project=. run_dual_ea_trivial_to_learning_lambda0.7_scan.jl --L $L --lambda_x $lambda_x --lambda_zz $lambda_zz --P_min $P_x --P_max $P_x --P_steps 1 --ntrials $ntrials --seed $seed --output_dir output --output_file ${out_prefix}.json"

julia --project=. run_dual_ea_trivial_to_learning_lambda0.7_scan.jl \
    --L $L \
    --lambda_x $lambda_x \
    --lambda_zz $lambda_zz \
    --P_min $P_x \
    --P_max $P_x \
    --P_steps 1 \
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
