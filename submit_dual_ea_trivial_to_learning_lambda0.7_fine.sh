#!/bin/bash

# Run a small LOCAL trivial-to-learning (lambda = 0.7) dual EA Binder FINE-grid scan
# (smoke test / quick look). Cluster production runs use
# jobs/jobs_dual_ea_trivial_to_learning_lambda0.7_fine.submit instead.

echo "=========================================="
echo "Trivial-to-Learning (lambda=0.7) dual EA Binder Scan, fine grid (local)"
echo "=========================================="
echo ""
echo "Physics setup:"
echo "  - lambda = 0.7 (delta = 0.7): λ_x = 0.49, λ_zz = 0.21, P_x = P_zz = q"
echo "  - Scanning q from 0.05 to 0.10 in 11 points (step 0.005)"
echo "  - Dynamics: doubled-MPS Born sampling (renyi2_learning_to_trivial_lambda0.7_core.jl)"
echo "  - Observable: dual EA Binder, open-chain endpoints a = 2..L"
echo ""

echo "Validating the observable against dense references..."
julia --project=. dual_ea_trivial_to_learning_lambda0.7_validate.jl || exit 1

L=${1:-8}
echo ""
echo "Running dual EA trivial-to-learning (lambda=0.7) fine-grid scan locally..."
julia --project=. run_dual_ea_trivial_to_learning_lambda0.7_scan.jl \
    --L $L \
    --P_min 0.05 --P_max 0.10 --P_steps ${2:-11} \
    --ntrials ${3:-20} \
    --output_dir dual_ea_trivial_to_learning_lambda0.7_fine_local_results/L${L} \
    --output_file dual_ea_trivial_to_learning_lambda0.7_fine_L${L}_local.json

echo ""
echo "✓ Dual EA trivial-to-learning (lambda=0.7) fine-grid scan completed!"
echo "   Results saved in: dual_ea_trivial_to_learning_lambda0.7_fine_local_results/L${L}/"
echo "   Analyze with: julia analyze_dual_ea_trivial_to_learning_lambda0.7_fine.jl dual_ea_trivial_to_learning_lambda0.7_fine_local_results 2"
