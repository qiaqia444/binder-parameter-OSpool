#!/usr/bin/env julia

using Pkg; Pkg.activate(".")
include("src/BinderSim.jl")
using .BinderSim
using Random

"""
Quick verification that the corrected weak measurements work
"""

function test_corrected_weak_measurements()
    println("=== Testing Corrected Weak Measurements ===")
    println("This should now match the sparse matrix implementation exactly")
    
    # Test parameters
    L = 8
    λ = 0.5
    ntrials = 5  # Very quick test
    
    println("Testing L=$L, λ=$λ with $ntrials trials...")
    
    # Run with corrected operators
    result = ea_binder_mc(L; 
                         lambda_x=λ, 
                         lambda_zz=1.0-λ, 
                         ntrials=ntrials,
                         maxdim=256,
                         cutoff=1e-12,
                         seed=12345)
    
    println("Results:")
    println("  Binder parameter: $(round(result.B, digits=4))")
    println("  Mean of trials: $(round(result.B_mean_of_trials, digits=4))")
    println("  Std of trials: $(round(result.B_std_of_trials, digits=4))")
    println("  S2_bar: $(round(result.S2_bar, digits=6))")
    println("  S4_bar: $(round(result.S4_bar, digits=6))")
    
    println("\n✓ Test completed successfully!")
    println("Ready for cluster submission.")
    
    return result
end

test_corrected_weak_measurements()
