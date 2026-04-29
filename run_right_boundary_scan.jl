#!/usr/bin/env julia

"""
Right Boundary Scan: Measurement-Induced Phase Transition with Rényi-2 Binder

Scan dephasing probability P_x from 0 to 1 at FIXED measurement strength λ_x.
This explores the right boundary of the phase diagram.

CRITICAL: Uses density matrix evolution (DiagonalStateMPS) for correct dephasing physics.
Focus: Rényi-2 Binder moments (M2, M4) instead of entanglement entropy proxy.
"""

using ArgParse
using JSON
using Random
using Statistics
using Dates
using ITensors, ITensorMPS

# Load modules in correct order
include("src_new/types.jl")
include("src_new/channels.jl")
include("src_new/dynamics_density_matrix.jl")
include("src_new/renyi2_binder.jl")

function parse_commandline()
    s = ArgParseSettings(description = "Right boundary scan: vary P_x at fixed λ_x with Rényi-2 Binder")
    
    @add_arg_table! s begin
        "--L"
            help = "System size"
            arg_type = Int
            default = 12
        
        "--lambda_x"
            help = "X measurement strength (FIXED)"
            arg_type = Float64
            default = 0.7
        
        "--lambda_zz"
            help = "ZZ measurement strength (FIXED)"
            arg_type = Float64
            default = 0.0
        
        "--P_min"
            help = "Minimum X dephasing probability"
            arg_type = Float64
            default = 0.0
        
        "--P_max"
            help = "Maximum X dephasing probability"
            arg_type = Float64
            default = 0.5
        
        "--P_steps"
            help = "Number of P_x values to scan"
            arg_type = Int
            default = 11
        
        "--ntrials"
            help = "Number of Monte Carlo trajectories"
            arg_type = Int
            default = 1000
        
        "--maxdim"
            help = "Maximum bond dimension"
            arg_type = Int
            default = 256
        
        "--cutoff"
            help = "SVD truncation cutoff"
            arg_type = Float64
            default = 1e-12
        
        "--seed"
            help = "Random seed"
            arg_type = Int
            default = 42
        
        "--output_dir"
            help = "Output directory"
            arg_type = String
            default = "right_boundary_results"
        
        "--output_file"
            help = "Output filename (optional, for cluster jobs)"
            arg_type = String
            default = ""
    end
    
    return parse_args(s)
end

function main()
    args = parse_commandline()
    
    L = args["L"]
    lambda_x = args["lambda_x"]
    lambda_zz = args["lambda_zz"]
    P_min = args["P_min"]
    P_max = args["P_max"]
    P_steps = args["P_steps"]
    ntrials = args["ntrials"]
    maxdim = args["maxdim"]
    cutoff = args["cutoff"]
    seed = args["seed"]
    output_dir = args["output_dir"]
    output_file = args["output_file"]
    
    println("="^70)
    println("RIGHT BOUNDARY SCAN: Measurement-Induced Phase Transition")
    println("="^70)
    println("Physics: Varying X dephasing P_x at fixed X measurement λ_x")
    println("Method: Density matrix evolution (DiagonalStateMPS)")
    println("Observable: Rényi-2 Binder (M2, M4 moments)")
    println()
    println("Parameters:")
    println("  System size L = $L")
    println("  X measurement λ_x = $lambda_x (FIXED)")
    println("  ZZ measurement λ_zz = $lambda_zz (FIXED - no ZZ)")
    println("  X dephasing scan: P_x ∈ [$P_min, $P_max] with $P_steps points")
    println("  ZZ dephasing: P_zz = 0 (no ZZ dephasing)")
    println("  Trajectories: $ntrials")
    println("  Max bond dimension: $maxdim")
    println("="^70)
    println()
    
    # Create P values to scan
    P_values = range(P_min, P_max, length=P_steps)
    
    # Create output directory
    mkpath(output_dir)
    
    results = []
    
    for (i, P) in enumerate(P_values)
        P_zz = 0.0  # Always zero for right boundary
        
        println("\n[$i/$P_steps] Running P_x = $(round(P, digits=3)), P_zz = $(round(P_zz, digits=3))")
        println("  (λ_x = $lambda_x, λ_zz = $lambda_zz)")
        println("-"^70)
        
        t_start = time()
        
        # Use Rényi-2 Binder with density matrix evolution
        result = renyi2_binder_density_matrix(
            L;
            lambda_x = lambda_x,
            lambda_zz = lambda_zz,
            P_x = P,
            P_zz = P_zz,
            ntrials = ntrials,
            maxdim = maxdim,
            cutoff = cutoff,
            seed = seed + i,  # Different seed for each P
            use_optimized = true,
            verbose = false
        )
        
        t_elapsed = time() - t_start
        
        # Store results
        result_dict = Dict(
            "L" => L,
            "lambda_x" => lambda_x,
            "lambda_zz" => lambda_zz,
            "P_x" => P,
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
            "maxdim" => maxdim,
            "cutoff" => cutoff,
            "time_seconds" => t_elapsed
        )
        
        push!(results, result_dict)
        
        # Print results
        println("  B (ensemble) = $(round(result.B, digits=4))")
        println("  B (mean of trials) = $(round(result.B_mean_of_trials, digits=4)) ± $(round(result.B_std_of_trials, digits=4))")
        println("  M₂_bar = $(round(result.M2_bar, digits=6))")
        println("  M₄_bar = $(round(result.M4_bar, digits=6))")
        println("  Purity = $(round(result.purity_bar, digits=4))")
        println("  Valid/Total = $(result.n_valid)/$(result.ntrials)")
        println("  Time: $(round(t_elapsed, digits=1)) seconds")
        
        # Save intermediate results
        if !isempty(output_file)
            # Use specified filename for cluster jobs
            outpath = joinpath(output_dir, output_file)
        else
            # Use timestamped filename for local runs
            timestamp = Dates.format(now(), "yyyymmdd_HHMM")
            outpath = joinpath(output_dir, "right_boundary_L$(L)_lambda$(lambda_x)_$(timestamp).json")
        end
        
        open(outpath, "w") do f
            JSON.print(f, results, 4)
        end
    end
    
    println("\n" * "="^70)
    println("✓ SCAN COMPLETE")
    println("="^70)
    println("Total results: $(length(results))")
    println("Output saved to: $output_dir/")
    
    # Final save
    if isempty(output_file)
        timestamp = Dates.format(Dates.now(), "yyyymmdd_HHMM")
        final_output = joinpath(output_dir, "right_boundary_L$(L)_lambda$(lambda_x)_$(timestamp)_final.json")
        open(final_output, "w") do f
            JSON.print(f, results, 4)
        end
        
        println("Final results: $final_output")
    end
    println()
    
    # Print summary statistics
    println("Summary of Rényi-2 Binder across scan:")
    println("  Min B = $(round(minimum(r["B"] for r in results), digits=4))")
    println("  Max B = $(round(maximum(r["B"] for r in results), digits=4))")
    M2_vals = [r["M2_bar"] for r in results]
    M4_vals = [r["M4_bar"] for r in results]
    println("  M2_bar range: [$(round(minimum(M2_vals), digits=4)), $(round(maximum(M2_vals), digits=4))]")
    println("  M4_bar range: [$(round(minimum(M4_vals), digits=4)), $(round(maximum(M4_vals), digits=4))]")
end

main()
