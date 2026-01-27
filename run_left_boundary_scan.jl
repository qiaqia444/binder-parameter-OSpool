#!/usr/bin/env julia

"""
Left Boundary Scan: Dephasing-Induced Phase Transition

Scan dephasing probability P from 0 to 1 at FIXED measurement strength λ.
This explores the left boundary of the phase diagram.

CRITICAL: Uses density matrix evolution (DiagonalStateMPS) instead of 
pure state trajectories, which is required for correct dephasing physics.
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

function parse_commandline()
    s = ArgParseSettings(description = "Left boundary scan: vary P at fixed λ")
    
    @add_arg_table! s begin
        "--L"
            help = "System size"
            arg_type = Int
            default = 12
        
        "--lambda_x"
            help = "X measurement strength (FIXED)"
            arg_type = Float64
            default = 0.3
        
        "--lambda_zz"
            help = "ZZ measurement strength (FIXED) - typically 0 for left boundary"
            arg_type = Float64
            default = 0.0
        
        "--P_min"
            help = "Minimum dephasing probability"
            arg_type = Float64
            default = 0.0
        
        "--P_max"
            help = "Maximum dephasing probability"
            arg_type = Float64
            default = 0.5
        
        "--P_steps"
            help = "Number of P values to scan"
            arg_type = Int
            default = 11
        
        "--P_zz_mode"
            help = "How to set P_zz: 'zero' (0), 'match' (=P_x), or 'fixed' (use --P_zz_value)"
            arg_type = String
            default = "zero"
        
        "--P_zz_value"
            help = "Fixed P_zz value if P_zz_mode='fixed'"
            arg_type = Float64
            default = 0.0
        
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
            default = "left_boundary_results"
        
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
    P_zz_mode = args["P_zz_mode"]
    P_zz_value = args["P_zz_value"]
    ntrials = args["ntrials"]
    maxdim = args["maxdim"]
    cutoff = args["cutoff"]
    seed = args["seed"]
    output_dir = args["output_dir"]
    output_file = args["output_file"]
    
    println("="^70)
    println("LEFT BOUNDARY SCAN: Dephasing-Induced Phase Transition")
    println("="^70)
    println("Physics: Varying dephasing P at fixed measurement λ")
    println("Method: Density matrix evolution (DiagonalStateMPS)")
    println()
    println("Parameters:")
    println("  System size L = $L")
    println("  X measurement λ_x = $lambda_x (FIXED)")
    println("  ZZ measurement λ_zz = $lambda_zz (FIXED)")
    println("  Dephasing scan: P_x ∈ [$P_min, $P_max] with $P_steps points")
    println("  P_zz mode: $P_zz_mode")
    if P_zz_mode == "fixed"
        println("    P_zz = $P_zz_value (fixed)")
    elseif P_zz_mode == "zero"
        println("    P_zz = 0 (no ZZ dephasing)")
    elseif P_zz_mode == "match"
        println("    P_zz = P_x (matched to X dephasing)")
    end
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
        # Determine P_zz based on mode
        if P_zz_mode == "zero"
            P_zz = 0.0
        elseif P_zz_mode == "match"
            P_zz = P
        elseif P_zz_mode == "fixed"
            P_zz = P_zz_value
        else
            error("Unknown P_zz_mode: $P_zz_mode. Use 'zero', 'match', or 'fixed'")
        end
        
        println("\n[$i/$P_steps] Running P_x = $(round(P, digits=3)), P_zz = $(round(P_zz, digits=3))")
        println("  (λ_x = $lambda_x, λ_zz = $lambda_zz)")
        println("-"^70)
        
        t_start = time()
        
        # Use density matrix evolution (CORRECT for dephasing)
        result = ea_binder_density_matrix(
            L;
            lambda_x = lambda_x,
            lambda_zz = lambda_zz,
            P_x = P,
            P_zz = P_zz,
            ntrials = ntrials,
            maxdim = maxdim,
            cutoff = cutoff,
            seed = seed + i  # Different seed for each P
        )
        
        t_elapsed = time() - t_start
        
        # Store results
        result_dict = Dict(
            "L" => L,
            "lambda_x" => lambda_x,
            "lambda_zz" => lambda_zz,
            "P_x" => P,
            "P_zz" => P_zz,
            "P_zz_mode" => P_zz_mode,
            "B" => result.B,
            "B_mean_of_trials" => result.B_mean_of_trials,
            "B_std_of_trials" => result.B_std_of_trials,
            "S2_bar" => result.S2_bar,
            "S4_bar" => result.S4_bar,
            "ntrials" => result.ntrials,
            "maxdim" => maxdim,
            "cutoff" => cutoff,
            "time_seconds" => t_elapsed
        )
        
        push!(results, result_dict)
        
        # Print results
        println("  B_EA = $(round(result.B, digits=4))")
        println("  B_mean = $(round(result.B_mean_of_trials, digits=4)) ± $(round(result.B_std_of_trials, digits=4))")
        println("  M₂² = $(round(result.S2_bar, digits=6))")
        println("  M₄² = $(round(result.S4_bar, digits=6))")
        println("  Time: $(round(t_elapsed, digits=1)) seconds")
        
        # Save intermediate results
        if !isempty(output_file)
            # Use specified filename for cluster jobs
            outpath = joinpath(output_dir, output_file)
        else
            # Use timestamped filename for local runs
            timestamp = Dates.format(now(), "yyyymmdd_HHMM")
            outpath = joinpath(output_dir, "left_boundary_L$(L)_lambda$(lambda_x)_$(timestamp).json")
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
        final_output = joinpath(output_dir, "left_boundary_L$(L)_lambda$(lambda_x)_$(timestamp)_final.json")
        open(final_output, "w") do f
            JSON.print(f, results, 4)
        end
        
        println("Final results: $final_output")
    end
    println()
    
    # Print summary statistics
    println("Summary of Binder parameter across scan:")
    println("  Min B = $(round(minimum(r["B"] for r in results), digits=4))")
    println("  Max B = $(round(maximum(r["B"] for r in results), digits=4))")
    println("  Range = $(round(maximum(r["B"] for r in results) - minimum(r["B"] for r in results), digits=4))")
end

main()
