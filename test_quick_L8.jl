#!/usr/bin/env julia

using JSON
using Random
using Statistics
using Dates
using ITensors, ITensorMPS

# Load modules
include("src_new/types.jl")
include("src_new/channels.jl")
include("src_new/dynamics_density_matrix.jl")
include("src_new/renyi2_dynamics_density_matrix_1.jl")

function main()
    println("="^70)
    println("QUICK TEST: L=8, λ_x=0.49, λ_zz=0.21, P_x=0.1, P_zz=0.1, ntrials=1")
    println("="^70)
    println()
    
    L = 8
    lambda_x = 0.49
    lambda_zz = 0.21
    P_x = 0.1
    P_zz = 0.1
    ntrials = 1
    maxdim = 256
    cutoff = 1e-12
    seed = 12345
    output_dir = "test_quick_L8_results"
    
    mkpath(output_dir)
    
    println("Configuration:")
    println("  L = $L")
    println("  λ_x = $lambda_x")
    println("  λ_zz = $lambda_zz")
    println("  P_x = $P_x")
    println("  P_zz = $P_zz")
    println("  ntrials = $ntrials")
    println()
    
    println("Starting calculation...")
    t_start = time()
    
    result = renyi2_binder_density_matrix(
        L;
        lambda_x = lambda_x,
        lambda_zz = lambda_zz,
        P_x = P_x,
        P_zz = P_zz,
        ntrials = ntrials,
        maxdim = maxdim,
        cutoff = cutoff,
        seed = seed,
        use_optimized = false,
        verbose = true
    )
    
    t_elapsed = time() - t_start
    
    println()
    println("✓ Completed in $(round(t_elapsed, digits=1)) seconds")
    println()
    
    println("Results:")
    println("  B = $(round(result.B, digits=4))")
    println("  M₂_bar = $(round(result.M2_bar, digits=6))")
    println("  M₄_bar = $(round(result.M4_bar, digits=6))")
    println("  Purity = $(round(result.purity_bar, digits=4))")
    println("  Valid/Total = $(result.n_valid)/$(result.ntrials)")
    println()
    
    # Save results
    result_dict = Dict(
        "L" => L,
        "lambda_x" => lambda_x,
        "lambda_zz" => lambda_zz,
        "P_x" => P_x,
        "P_zz" => P_zz,
        "B" => result.B,
        "M2_bar" => result.M2_bar,
        "M4_bar" => result.M4_bar,
        "purity_bar" => result.purity_bar,
        "ntrials" => result.ntrials,
        "n_valid" => result.n_valid,
        "n_invalid" => result.n_invalid,
        "time_seconds" => t_elapsed,
        "timestamp" => Dates.format(now(), "yyyy-mm-dd HH:MM:SS")
    )
    
    outpath = joinpath(output_dir, "result.json")
    open(outpath, "w") do f
        JSON.print(f, result_dict, 2)
    end
    
    println("✓ Results saved to: $outpath")
end

main()
