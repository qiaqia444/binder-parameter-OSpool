#!/usr/bin/env julia

using Pkg
Pkg.activate(".")

include("src/BinderSimManual.jl")
using .BinderSimManual

println("Testing manual tensor contraction approach...")
println("=" ^60)

# Test small system first (L=12) to verify it matches original approach
println("Testing L=12 (should work with both approaches)...")
try
    result_L12 = ea_binder_mc_manual(12; 
                                    lambda_x=0.5, lambda_zz=0.5, 
                                    ntrials=3, 
                                    maxdim=128, cutoff=1e-10,
                                    manual_maxdim=64, manual_cutoff=1e-8,
                                    chunk_size=50,
                                    seed=1234)
    
    println("L=12 SUCCESS!")
    println("  Binder parameter: $(result_L12.B)")
    println("  Trials completed: $(result_L12.ntrials_completed)/$(result_L12.ntrials)")
    println("  Mean ± std: $(result_L12.B_mean_of_trials) ± $(result_L12.B_std_of_trials)")
    
catch e
    println("L=12 FAILED: $e")
end

println()
println("=" ^60)

# Test medium system (L=16) - this should now work!
println("Testing L=16 (previously failed with ITensorCorrelators)...")
try
    result_L16 = ea_binder_mc_manual(16; 
                                    lambda_x=0.5, lambda_zz=0.5, 
                                    ntrials=2, 
                                    maxdim=64, cutoff=1e-8,
                                    manual_maxdim=32, manual_cutoff=1e-6,
                                    chunk_size=25,
                                    seed=1235)
    
    println("L=16 SUCCESS!")
    println("  Binder parameter: $(result_L16.B)")
    println("  Trials completed: $(result_L16.ntrials_completed)/$(result_L16.ntrials)")
    println("  Mean ± std: $(result_L16.B_mean_of_trials) ± $(result_L16.B_std_of_trials)")
    
catch e
    println("L=16 FAILED: $e")
    println("Error details: ", sprint(showerror, e))
end

println()
println("=" ^60)

# Test large system (L=20) - the ultimate test!
println("Testing L=20 (ultimate challenge)...")
try
    result_L20 = ea_binder_mc_manual(20; 
                                    lambda_x=0.5, lambda_zz=0.5, 
                                    ntrials=1, 
                                    maxdim=32, cutoff=1e-6,
                                    manual_maxdim=16, manual_cutoff=1e-4,
                                    chunk_size=10,
                                    seed=1236)
    
    println("L=20 SUCCESS!")
    println("  Binder parameter: $(result_L20.B)")
    println("  Trials completed: $(result_L20.ntrials_completed)/$(result_L20.ntrials)")
    
catch e
    println("L=20 FAILED: $e")
    println("Error details: ", sprint(showerror, e))
end

println()
println("=" ^60)
println("Manual tensor contraction testing complete!")
