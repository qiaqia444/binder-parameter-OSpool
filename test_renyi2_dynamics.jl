#!/usr/bin/env julia

"""
Test the new renyi2_dynamics pipeline with realistic parameters
"""

using JSON
using Random
using Statistics
using Dates
using ITensors, ITensorMPS

# Load the new dynamics module
include("renyi2_dynamics.jl")

function main()
    println("="^70)
    println("RENYI2-DYNAMICS PIPELINE TEST")
    println("="^70)
    println()
    
    L = 8
    lambda_x = 0.49
    lambda_zz = 0.21
    P_x = 0.1
    P_zz = 0.1
    ntrials = 10
    
    println("Configuration:")
    println("  L = $L")
    println("  λ_x = $lambda_x, λ_zz = $lambda_zz")
    println("  P_x = $P_x, P_zz = $P_zz")
    println("  ntrials = $ntrials")
    println()
    
    println("Running ensemble calculation...")
    t_start = time()
    
    result = renyi2_binder_density_matrix_separate(
        L;
        lambda_x = lambda_x,
        lambda_zz = lambda_zz,
        P_x = P_x,
        P_zz = P_zz,
        ntrials = ntrials,
        maxdim = 256,
        cutoff = 1e-12,
        seed = 42,
        T_max = 2*L,
        strobe = :after_full_layer,
        verbose = true
    )
    
    t_elapsed = time() - t_start
    
    println()
    println("✓ Completed in $(round(t_elapsed, digits=1)) seconds")
    println()
    
    println("Results:")
    println("  B (ensemble) = $(round(result.B, digits=4))")
    println("  B (mean of trials) = $(round(result.B_mean_of_trials, digits=4)) ± $(round(result.B_std_of_trials, digits=4))")
    println("  M₂_bar = $(round(result.M2_bar, digits=6))")
    println("  M₄_bar = $(round(result.M4_bar, digits=6))")
    println("  Purity = $(round(result.purity_bar, digits=4))")
    println("  Valid trials = $(result.n_valid)/$(result.ntrials)")
    println("  Max link dim = $(result.max_linkdim)")
    println()
    
    if result.n_valid < result.ntrials
        println("⚠ Warning: $(result.n_invalid) invalid trials (NaN Binder)")
    end
    
    # Save results
    result_dict = Dict(
        "L" => L,
        "lambda_x" => lambda_x,
        "lambda_zz" => lambda_zz,
        "P_x" => P_x,
        "P_zz" => P_zz,
        "B" => result.B,
        "B_mean_of_trials" => result.B_mean_of_trials,
        "B_std_of_trials" => result.B_std_of_trials,
        "M2_bar" => result.M2_bar,
        "M4_bar" => result.M4_bar,
        "purity_bar" => result.purity_bar,
        "ntrials" => result.ntrials,
        "n_valid" => result.n_valid,
        "n_invalid" => result.n_invalid,
        "max_linkdim" => result.max_linkdim,
        "time_seconds" => t_elapsed,
        "timestamp" => Dates.format(now(), "yyyy-mm-dd HH:MM:SS")
    )
    
    mkpath("test_renyi2_dynamics")
    outpath = "test_renyi2_dynamics/result.json"
    open(outpath, "w") do f
        JSON.print(f, result_dict, 2)
    end
    
    println("✓ Results saved to: $outpath")
end

main()
