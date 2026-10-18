#!/bin/bash

# Run a small LOCAL right boundary dual EA Binder FINE-grid scan (smoke test / quick look).
# Cluster production runs use jobs/jobs_dual_ea_right_boundary_fine.submit instead.

echo "=========================================="
echo "Right Boundary dual EA Binder Scan, fine grid (local)"
echo "=========================================="
echo ""
echo "Physics setup:"
echo "  - Right boundary: λ_x = 0.7, λ_zz = 0.0, P_x = P_zz = q"
echo "  - Scanning q from 0.25 to 0.35 in 11 points (step 0.01)"
echo "  - Dynamics: existing doubled-MPS Born sampling (renyi2_right_boundary_core.jl)"
echo "  - Observable: dual EA Binder, open-chain endpoints a = 2..L"
echo ""

echo "Validating the observable against dense references..."
julia --project=. dual_ea_right_boundary_validate.jl || exit 1

L=${1:-8}
echo ""
echo "Running dual EA right boundary fine-grid scan locally..."
julia --project=. run_dual_ea_right_boundary_scan.jl \
    --L $L \
    --P_min 0.25 --P_max 0.35 --P_steps ${2:-11} \
    --ntrials ${3:-20} \
    --output_dir dual_ea_right_boundary_fine_local_results/L${L} \
    --output_file dual_ea_right_boundary_fine_L${L}_local.json

echo ""
echo "✓ Dual EA right boundary fine-grid scan completed!"
echo "   Results saved in: dual_ea_right_boundary_fine_local_results/L${L}/"
echo "   Analyze with: julia analyze_dual_ea_right_boundary_fine.jl dual_ea_right_boundary_fine_local_results 2"
