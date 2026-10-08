#!/bin/bash

# Run a small LOCAL right boundary Renyi-2 scan measured after the last X dephasing
# (smoke test / quick look). Cluster production runs use jobs/jobs_renyi2_right_boundary_xnoise*.submit.

echo "=============================================================="
echo "Right Boundary Renyi-2 Binder, after the last X dephasing (local)"
echo "=============================================================="
echo ""
echo "Physics setup:"
echo "  - Right boundary: λ_x = 0.7, λ_zz = 0.0, P_x = P_zz = q, q from 0 to 0.5"
echo "  - Dynamics: production doubled-MPS Born sampling (renyi2_right_boundary_core.jl)"
echo "  - Observable: Rényi-2 Binder on the state after the X dephasing of the LAST step"
echo "    (before the last ZZ-dephasing layer), as for the left-boundary EA code"
echo ""

echo "Validating against the production evolution and a dense reference..."
julia --project=. renyi2_right_boundary_xnoise_validate.jl || exit 1

echo ""
echo "Running the scan locally..."
julia --project=. run_renyi2_right_boundary_xnoise_scan.jl \
    --L ${1:-8} \
    --P_min 0.0 --P_max 0.5 --P_steps ${2:-11} \
    --ntrials ${3:-20} \
    --output_dir renyi2_right_boundary_xnoise_results

echo ""
echo "✓ Done. Results saved in: renyi2_right_boundary_xnoise_results/"
