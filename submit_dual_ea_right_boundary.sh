#!/bin/bash

# Run a small LOCAL right boundary dual EA Binder scan (smoke test / quick look).
# Cluster production runs use jobs/jobs_dual_ea_right_boundary*.submit instead.

echo "=========================================="
echo "Right Boundary dual EA Binder Scan (local)"
echo "=========================================="
echo ""
echo "Physics setup:"
echo "  - Right boundary: λ_x = 0.7, λ_zz = 0.0, P_x = P_zz = q"
echo "  - Scanning q from 0 to 0.5"
echo "  - Dynamics: existing doubled-MPS Born sampling (renyi2_right_boundary_core.jl)"
echo "  - Observable: dual EA Binder, open-chain endpoints a = 2..L"
echo ""

echo "Validating the observable against dense references..."
julia --project=. dual_ea_right_boundary_validate.jl || exit 1

echo ""
echo "Running dual EA right boundary scan locally..."
julia --project=. run_dual_ea_right_boundary_scan.jl \
    --L ${1:-8} \
    --P_min 0.0 --P_max 0.5 --P_steps ${2:-11} \
    --ntrials ${3:-20} \
    --output_dir dual_ea_right_boundary_results

echo ""
echo "✓ Dual EA right boundary scan completed!"
echo "   Results saved in: dual_ea_right_boundary_results/"
