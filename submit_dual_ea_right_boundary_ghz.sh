#!/bin/bash

# Run a small LOCAL right boundary dual EA Binder scan with a GHZ initial state (smoke test).
# Cluster production runs use jobs/jobs_dual_ea_right_boundary_ghz.submit instead.

echo "=========================================="
echo "Right Boundary dual EA Binder Scan, GHZ start (local)"
echo "=========================================="
echo ""
echo "Physics setup:"
echo "  - Right boundary: λ_x = 0.7, λ_zz = 0.0, P_x = P_zz = q"
echo "  - Initial state: |GHZ> = (|0...0> + |1...1>)/sqrt(2)"
echo "  - Scanning q from 0 to 0.5"
echo "  - Observable: dual EA Binder, open-chain endpoints a = 2..L"
echo ""

echo "Validating the observable and the GHZ dynamics against dense references..."
julia --project=. dual_ea_right_boundary_ghz_validate.jl || exit 1

L=${1:-8}
echo ""
echo "Running dual EA right boundary GHZ scan locally..."
julia --project=. run_dual_ea_right_boundary_ghz_scan.jl \
    --L $L \
    --P_min 0.0 --P_max 0.5 --P_steps ${2:-11} \
    --ntrials ${3:-20} \
    --output_dir dual_ea_right_boundary_ghz_local_results/L${L} \
    --output_file dual_ea_right_boundary_ghz_L${L}_local.json

echo ""
echo "✓ Dual EA right boundary GHZ scan completed!"
echo "   Results saved in: dual_ea_right_boundary_ghz_local_results/L${L}/"
echo "   Analyze with: julia analyze_dual_ea_right_boundary_ghz.jl dual_ea_right_boundary_ghz_local_results"
