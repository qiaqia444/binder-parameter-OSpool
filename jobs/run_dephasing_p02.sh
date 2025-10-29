#!/bin/bash

# Wrapper script for HTCondor: Dephasing P=0.2 simulation
# This runs the Julia simulation with dephasing channels

# Arguments from HTCondor
L=$1
P_x=$2
lambda_x=$3
lambda_zz=$4
lambda=$5
ntrials=$6
seed=$7
sample=$8
out_prefix=$9

echo "========================================"
echo "Dephasing P=0.2 Simulation Starting"
echo "========================================"
echo "System size:        L = $L"
echo "Dephasing prob:     P_x = $P_x"
echo "Lambda X:           lambda_x = $lambda_x"
echo "Lambda ZZ:          lambda_zz = $lambda_zz"
echo "Lambda:             lambda = $lambda"
echo "Trials:             ntrials = $ntrials"
echo "Seed:               seed = $seed"
echo "Sample:             sample = $sample"
echo "Output prefix:      $out_prefix"
echo "========================================"

# Run the Julia script
julia --project=. run_dephasing_p02.jl \
    --L $L \
    --P_x $P_x \
    --lambda_x $lambda_x \
    --lambda_zz $lambda_zz \
    --lambda $lambda \
    --ntrials $ntrials \
    --seed $seed \
    --sample $sample \
    --out_prefix $out_prefix

exit_code=$?

echo "========================================"
echo "Simulation completed with exit code: $exit_code"
echo "========================================"

exit $exit_code
