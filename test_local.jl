#!/usr/bin/env julia
"""
Test script to run a single simulation locally before submitting to cluster
"""

using Pkg; Pkg.activate(".")
push!(LOAD_PATH, joinpath(@__DIR__, "src"))
using BinderSim
using JSON

function test_single_run()
    # Test parameters similar to your cluster jobs but smaller for quick test
    params = Dict(
        "L" => 12,
        "lambda_x" => 0.5,
        "lambda_zz" => 0.5,
        "ntrials" => 20,  # Small number for quick test (cluster will use 1000)
        "maxdim" => 256,
        "cutoff" => 1e-12,
        "chunk4" => 50_000,
        "seed" => 1234
    )
    
    println("Testing Binder simulation with parameters:")
    println(JSON.json(params, 2))
    println()
    
    # Run the simulation
    println("Running simulation...")
    println("Note: Using T_max = 2L = $(2*params["L"]) time steps")
    result = ea_binder_mc(params["L"]; 
                         lambda_x=params["lambda_x"], 
                         lambda_zz=params["lambda_zz"],
                         ntrials=params["ntrials"],
                         maxdim=params["maxdim"],
                         cutoff=params["cutoff"],
                         chunk4=params["chunk4"],
                         seed=params["seed"])
    
    println("Results:")
    println("  Binder parameter (EA): $(result.B)")
    println("  Binder parameter (mean of trials): $(result.B_mean_of_trials)")
    println("  Binder parameter (std of trials): $(result.B_std_of_trials)")
    println("  S2_bar: $(result.S2_bar)")
    println("  S4_bar: $(result.S4_bar)")
    println("  Number of trials: $(result.ntrials)")
    println()
    println("Ready for cluster deployment!")
    println("  - System sizes L: [12, 16, 20, 24, 28]")
    println("  - Lambda range: 0.1-0.9 with fine grid around 0.5")
    println("  - Trials per job: 1000")
    println("  - Total jobs: 255")
    println("  - Total Monte Carlo trials: 255,000")
    
    return result
end

if !isinteractive()
    test_single_run()
end
