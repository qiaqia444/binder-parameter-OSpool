#!/usr/bin/env julia

# Quick test of dephasing channel simulation
# Run with: julia test_dephasing.jl

using Pkg; Pkg.activate(".")
using ITensors, ITensorMPS
using Random

include("src/BinderSimDephasing.jl")
using .BinderSimDephasing

println("========================================")
println("Testing Dephasing Channel Functions")
println("========================================")

# Test parameters
L = 8
P_x = 0.1
P_zz = 0.1
ntrials = 10
maxdim = 128
cutoff = 1e-12

Random.seed!(1234)

println("\nTest 1: Single trial evolution with dephasing")
println("L=$L, P_x=$P_x, P_zz=$P_zz")

try
    ψ, sites = evolve_one_trial_dephasing(L; P_x=P_x, P_zz=P_zz, 
                                          maxdim=maxdim, cutoff=cutoff)
    println("✓ Single trial completed successfully")
    println("  Bond dimension: $(maxlinkdim(ψ))")
    println("  Norm: $(norm(ψ))")
catch e
    println("✗ Error in single trial: $e")
    rethrow(e)
end

println("\nTest 2: Full Binder parameter calculation with dephasing")
println("Running $ntrials trials...")

try
    result = ea_binder_mc_dephasing(
        L;
        P_x = P_x,
        P_zz = P_zz,
        maxdim = maxdim,
        cutoff = cutoff,
        ntrials = ntrials
    )
    
    println("✓ Simulation completed successfully!")
    println("\nResults:")
    println("  Binder parameter (B_EA): $(round(result.B, digits=6))")
    println("  Mean B over trials: $(round(result.B_mean_of_trials, digits=6))")
    println("  Std B over trials: $(round(result.B_std_of_trials, digits=6))")
    println("  S2_bar: $(round(result.S2_bar, digits=6))")
    println("  S4_bar: $(round(result.S4_bar, digits=6))")
    println("  Completed trials: $(result.ntrials)")
    
catch e
    println("✗ Error in full simulation: $e")
    rethrow(e)
end

println("\nTest 3: Comparing no dephasing (P=0) vs dephasing (P>0)")

try
    println("  Running with P_x=0.0, P_zz=0.0...")
    result_no_deph = ea_binder_mc_dephasing(L; P_x=0.0, P_zz=0.0, 
                                            maxdim=maxdim, cutoff=cutoff, ntrials=10)
    
    println("  Running with P_x=0.2, P_zz=0.2...")
    result_with_deph = ea_binder_mc_dephasing(L; P_x=0.2, P_zz=0.2, 
                                              maxdim=maxdim, cutoff=cutoff, ntrials=10)
    
    println("\n  Comparison:")
    println("    No dephasing: B = $(round(result_no_deph.B, digits=4))")
    println("    With dephasing: B = $(round(result_with_deph.B, digits=4))")
    println("    Difference: $(round(abs(result_no_deph.B - result_with_deph.B), digits=4))")
    
    println("\n✓ All tests passed!")
    
catch e
    println("✗ Error in comparison: $e")
    rethrow(e)
end

println("\n========================================")
println("Dephasing channel implementation ready!")
println("========================================")
