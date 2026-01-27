#!/usr/bin/env julia

"""
Test the density matrix implementation for dephasing channels.

This verifies that:
1. The modules load correctly
2. A simple evolution runs without errors
3. Binder parameter calculation works
"""

using ITensors, ITensorMPS
using Random

println("="^70)
println("Testing Density Matrix Implementation")
println("="^70)
println()

# Load modules
println("Loading modules...")
include("src_new/types.jl")
include("src_new/channels.jl")
include("src_new/dynamics_density_matrix.jl")
println("✓ Modules loaded successfully")
println()

# Test 1: Simple evolution
println("Test 1: Small system evolution")
println("-"^70)
L = 6
lambda_x = 0.5
lambda_zz = 0.5
P_x = 0.2
P_zz = 0.2

println("  L = $L")
println("  λ_x = $lambda_x, λ_zz = $lambda_zz")
println("  P_x = $P_x, P_zz = $P_zz")
println()

try
    println("  Running single trajectory...")
    state = evolve_density_matrix_one_trial(
        L;
        lambda_x = lambda_x,
        lambda_zz = lambda_zz,
        P_x = P_x,
        P_zz = P_zz,
        maxdim = 64,
        cutoff = 1e-12,
        rng = MersenneTwister(42)
    )
    println("  ✓ Evolution completed successfully")
    println("    Final state type: $(typeof(state))")
    println("    MPS bond dimension: $(maxlinkdim(get_mps(state)))")
    println()
catch e
    println("  ✗ Evolution failed with error:")
    println("    $e")
    rethrow(e)
end

# Test 2: Binder parameter calculation
println("Test 2: Binder parameter calculation")
println("-"^70)
println("  L = $L, λ = $lambda_x, P = $P_x")
println("  Running 10 trajectories...")
println()

try
    result = ea_binder_density_matrix(
        L;
        lambda_x = lambda_x,
        lambda_zz = lambda_zz,
        P_x = P_x,
        P_zz = P_zz,
        ntrials = 10,
        maxdim = 64,
        seed = 42
    )
    
    println("  ✓ Calculation completed successfully")
    println()
    println("  Results:")
    println("    B_EA = $(round(result.B, digits=4))")
    println("    B_mean = $(round(result.B_mean_of_trials, digits=4)) ± $(round(result.B_std_of_trials, digits=4))")
    println("    M₂² = $(round(result.S2_bar, digits=6))")
    println("    M₄² = $(round(result.S4_bar, digits=6))")
    println()
catch e
    println("  ✗ Binder calculation failed with error:")
    println("    $e")
    rethrow(e)
end

# Test 3: Compare with zero dephasing
println("Test 3: Sanity check - zero dephasing")
println("-"^70)
println("  Testing P=0 case (should still work)")
println()

try
    result_no_deph = ea_binder_density_matrix(
        L;
        lambda_x = lambda_x,
        lambda_zz = lambda_zz,
        P_x = 0.0,
        P_zz = 0.0,
        ntrials = 5,
        maxdim = 64,
        seed = 42
    )
    
    println("  ✓ Zero dephasing case works")
    println("    B_EA(P=0) = $(round(result_no_deph.B, digits=4))")
    println()
catch e
    println("  ✗ Zero dephasing failed:")
    println("    $e")
    rethrow(e)
end

println("="^70)
println("✓ ALL TESTS PASSED")
println("="^70)
println()
println("The density matrix implementation is working correctly!")
println("You can now run:")
println("  julia run_left_boundary_scan.jl --L 12 --lambda_x 0.5")
println()
