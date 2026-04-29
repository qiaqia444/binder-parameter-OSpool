#!/bin/bash

# Submit Right Boundary Scan jobs

echo "=========================================="
echo "Submitting Right Boundary Scan Jobs"
echo "=========================================="
echo ""
echo "Physics setup:"
echo "  - Right boundary scan (measurement-induced transition)"
echo "  - Fixed λ_x = 0.7 (X measurement strength)"
echo "  - Fixed λ_zz = 0.0 (no ZZ measurements)"
echo "  - Scanning P_x from 0 to 0.5"
echo "  - P_zz = 0 (no ZZ dephasing)"
echo "  - Focus: Rényi-2 Binder (M2, M4 moments)"
echo "  - Using DENSITY MATRIX evolution"
echo ""

echo "Running right boundary scan locally..."
julia run_right_boundary_simple.jl

echo ""
echo "✓ Right Boundary Scan completed!"
echo "   Results saved in: right_boundary_results/"
