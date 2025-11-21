#!/bin/bash

# Wrapper script for HTCondor: RBIM simulation
# Random Bond Ising Model in strong dephasing limit (λ_x = 0)

# Arguments from HTCondor
L=$1
P_x=$2
lambda_zz=$3
ntrials=$4
seed=$5
sample=$6
out_prefix=$7

echo "========================================"
echo "RBIM Simulation"
echo "Strong dephasing limit: λ_x = 0"
echo "========================================"
echo "System size:        L = $L"
echo "Dephasing prob:     P_x = $P_x"
echo "ZZ measurement:     lambda_zz = $lambda_zz"
echo "X measurement:      lambda_x = 0 (RBIM)"
echo "Trials:             ntrials = $ntrials"
echo "Seed:               seed = $seed"
echo "Sample:             sample = $sample"
echo "Output prefix:      $out_prefix"
echo "EA correlation: [⟨Z_i Z_j⟩²]_J"
echo "========================================"

# Create output directory
mkdir -p output

# Check Julia version
echo "Julia version:"
julia --version

# Install packages
echo "Setting up Julia environment..."
julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'

# Run the Julia script
julia --project=. run_rbim.jl \
    --L $L \
    --P_x $P_x \
    --lambda_zz $lambda_zz \
    --ntrials $ntrials \
    --seed $seed \
    --sample $sample \
    --out_prefix $out_prefix

exit_code=$?

# List output files
echo "Output files created:"
ls -la output/ 2>/dev/null || echo "No output directory found"

echo "========================================"
echo "Simulation completed with exit code: $exit_code"
echo "========================================"

exit $exit_code
