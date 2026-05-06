#!/usr/bin/env julia

"""
Learning-to-Trivial Lambda0.7 - Single Point Calculation

Computes Rényi-2 Binder for one configuration using the CORRECT dynamics pipeline.

Usage (from HTCondor job):
    julia run_learning_to_trivial_lambda0.7_scan.jl \\
        L lambda_x lambda_zz P_x P_zz ntrials seed sample out_prefix

Example:
    julia run_learning_to_trivial_lambda0.7_scan.jl \\
        8 0.49 0.21 0.3 0.3 100 12345 1 learning_to_trivial_lambda0.7_L8_Px0.30
"""

using JSON
using Random
using Statistics
using Dates
using ITensors, ITensorMPS

# Load the CORRECTED separate pipeline (with channel averaging dephasing fix)
include("src_new/types.jl")
include("src_new/channels.jl")
include("src_new/dynamics_density_matrix.jl")
include("src_new/renyi2_binder.jl")
include("src_new/renyi2_dynamics_separate.jl")

function main()
    # Parse command-line arguments
    if length(ARGS) < 9
        error("Usage: julia run_learning_to_trivial_lambda0.7_scan.jl L lambda_x lambda_zz P_x P_zz ntrials seed sample out_prefix")
    end
    
    L = parse(Int, ARGS[1])
    lambda_x = parse(Float64, ARGS[2])
    lambda_zz = parse(Float64, ARGS[3])
    P_x = parse(Float64, ARGS[4])
    P_zz = parse(Float64, ARGS[5])
    ntrials = parse(Int, ARGS[6])
    seed = parse(Int, ARGS[7])
    sample = parse(Int, ARGS[8])
    out_prefix = ARGS[9]
    
    # Create output directory
    mkpath("output")
    
    println("="^70)
    println("Learning-to-Trivial Transition: λ_x=0.49, λ_zz=0.21")
    println("="^70)
    println()
    println("Configuration:")
    println("  L = $L")
    println("  λ_x = $lambda_x (fixed)")
    println("  λ_zz = $lambda_zz (fixed)")
    println("  P_x = $P_x")
    println("  P_zz = $P_zz")
    println("  ntrials = $ntrials")
    println("  seed = $seed")
    println("  sample = $sample")
    println()
    
    println("Starting Rényi-2 Binder calculation (using separate pipeline with corrected dephasing)...")
    t_start = time()
    
    # Run calculation using CORRECTED separate pipeline with channel averaging dephasing
    result = renyi2_binder_density_matrix_separate(
        L;
        lambda_x = lambda_x,
        lambda_zz = lambda_zz,
        P_x = P_x,
        P_zz = P_zz,
        ntrials = ntrials,
        maxdim = 256,
        cutoff = 1e-12,
        seed = seed,
        verbose = false,
    )
    
    t_elapsed = time() - t_start
    
    println("✓ Completed in $(round(t_elapsed, digits=1)) seconds")
    println()
    println("Results:")
    println("  B (ensemble) = $(round(result.B, digits=6))")
    println("  B (mean ± std) = $(round(result.B_mean_of_trials, digits=6)) ± $(round(result.B_std_of_trials, digits=6))")
    println("  M₂_bar = $(round(result.M2_bar, digits=8))")
    println("  M₄_bar = $(round(result.M4_bar, digits=8))")
    println("  Purity = $(round(result.purity_bar, digits=6))")
    println("  Valid/Total = $(result.n_valid)/$(result.ntrials)")
    println()
    
    # Prepare output dictionary
    output_dict = Dict(
        "L" => L,
        "lambda_x" => lambda_x,
        "lambda_zz" => lambda_zz,
        "P_x" => P_x,
        "P_zz" => P_zz,
        "sample" => sample,
        "B" => result.B,
        "B_mean_of_trials" => result.B_mean_of_trials,
        "B_std_of_trials" => result.B_std_of_trials,
        "M2_bar" => result.M2_bar,
        "M4_bar" => result.M4_bar,
        "purity_bar" => result.purity_bar,
        "M2s" => result.M2s,
        "M4s" => result.M4s,
        "Bs" => result.Bs,
        "purities" => result.purities,
        "reliabilities" => result.reliabilities,
        "ntrials" => result.ntrials,
        "n_valid" => result.n_valid,
        "n_invalid" => result.n_invalid,
        "time_seconds" => t_elapsed,
        "timestamp" => Dates.format(now(), "yyyy-mm-dd HH:MM:SS"),
    )
    
    # Write output JSON
    outpath = "output/$(out_prefix).json"
    open(outpath, "w") do f
        JSON.print(f, output_dict, 2)
    end
    
    println("✓ Results saved to: $outpath")
    println()
    println("="^70)
end

main()
