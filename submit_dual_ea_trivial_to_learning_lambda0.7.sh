#!/bin/bash

# Run a small LOCAL trivial-to-learning (lambda = 0.7) dual EA Binder scan (smoke test / quick look).
# Cluster production runs use jobs/jobs_dual_ea_trivial_to_learning_lambda0.7.submit instead.

echo "=========================================="
echo "Trivial-to-Learning (lambda=0.7) dual EA Binder Scan (local)"
echo "=========================================="
echo ""
echo "Physics setup:"
echo "  - lambda = 0.7 (delta = 0.7): λ_x = 0.49, λ_zz = 0.21, P_x = P_zz = q"
echo "  - Scanning q from 0 to 0.5"
echo "  - Dynamics: doubled-MPS Born sampling (renyi2_learning_to_trivial_lambda0.7_core.jl)"
echo "  - Observable: dual EA Binder, open-chain endpoints a = 2..L"
echo ""

echo "Validating the observable against dense references..."
julia --project=. dual_ea_trivial_to_learning_lambda0.7_validate.jl || exit 1

echo ""
echo "Running dual EA trivial-to-learning (lambda=0.7) scan locally..."
julia --project=. run_dual_ea_trivial_to_learning_lambda0.7_scan.jl \
    --L ${1:-8} \
    --P_min 0.0 --P_max 0.5 --P_steps ${2:-11} \
    --ntrials ${3:-20} \
    --output_dir dual_ea_trivial_to_learning_lambda0.7_results

echo ""
echo "✓ Dual EA trivial-to-learning (lambda=0.7) scan completed!"
echo "   Results saved in: dual_ea_trivial_to_learning_lambda0.7_results/"
