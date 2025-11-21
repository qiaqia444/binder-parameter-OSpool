#!/usr/bin/env julia

"""
RBIM (Random Bond Ising Model) Simulation
Strong dephasing limit: λ_x = 0, only ZZ measurements survive

In this limit, the quantum circuit reduces to a classical RBIM due to gauge symmetry.
The correct order parameter is the EA (Edwards-Anderson) correlation:
    C_EA(i,j) = [⟨Z_i Z_j⟩²]_J

where the disorder average is over bond configurations.
"""

using ArgParse
using Random
using Statistics
using JSON
using ITensors
using ITensorMPS

# Import the BinderSim module for shared functions
include("src/BinderSim.jl")
using .BinderSim

include("src/BinderSimDephasing.jl")
using .BinderSimDephasing

function parse_commandline()
    s = ArgParseSettings(description = "RBIM simulation with strong dephasing (λ_x=0)")
    
    @add_arg_table! s begin
        "--L"
            help = "System size"
            arg_type = Int
            required = true
        "--P_x"
            help = "X dephasing probability"
            arg_type = Float64
            default = 0.5
        "--lambda_zz"
            help = "ZZ measurement strength"
            arg_type = Float64
            required = true
        "--ntrials"
            help = "Number of trials (disorder samples)"
            arg_type = Int
            default = 1000
        "--seed"
            help = "Random seed"
            arg_type = Int
            default = nothing
        "--sample"
            help = "Sample number (for parallel jobs)"
            arg_type = Int
            default = 1
        "--out_prefix"
            help = "Output file prefix"
            arg_type = String
            default = "rbim"
        "--maxdim"
            help = "Maximum MPS bond dimension"
            arg_type = Int
            default = 256
        "--cutoff"
            help = "SVD cutoff"
            arg_type = Float64
            default = 1e-12
    end
    
    return parse_args(s)
end

function main()
    args = parse_commandline()
    
    L = args["L"]
    P_x = args["P_x"]
    lambda_zz = args["lambda_zz"]
    ntrials = args["ntrials"]
    seed = args["seed"]
    sample = args["sample"]
    out_prefix = args["out_prefix"]
    maxdim = args["maxdim"]
    cutoff = args["cutoff"]
    
    # RBIM limit: λ_x = 0 (no X measurements)
    lambda_x = 0.0
    P_zz = P_x  # Same dephasing strength
    
    println("="^70)
    println("RBIM Simulation - Strong Dephasing Limit")
    println("="^70)
    println("System size:           L = $L")
    println("X dephasing:           P_x = $P_x")
    println("ZZ dephasing:          P_zz = $P_zz")
    println("X measurement:         λ_x = $lambda_x (RBIM: set to 0)")
    println("ZZ measurement:        λ_zz = $lambda_zz")
    println("Trials:                $ntrials")
    println("Sample:                $sample")
    println("Seed:                  $seed")
    println("Max bond dim:          $maxdim")
    println("Cutoff:                $cutoff")
    println("="^70)
    println()
    
    # Run simulation using the dephasing module
    # This will give us EA-type correlations automatically
    result = ea_binder_mc_dephasing(
        L;
        lambda_x=lambda_x,
        lambda_zz=lambda_zz,
        P_x=P_x,
        P_zz=P_zz,
        ntrials=ntrials,
        maxdim=maxdim,
        cutoff=cutoff,
        seed=seed
    )
    
    println("\nResults:")
    println("  EA Binder Parameter:     $(round(result.B, digits=6))")
    println("  Binder (mean of trials): $(round(result.B_mean_of_trials, digits=6))")
    println("  Binder (std of trials):  $(round(result.B_std_of_trials, digits=6))")
    println("  S2_bar:                  $(round(result.S2_bar, digits=6))")
    println("  S4_bar:                  $(round(result.S4_bar, digits=6))")
    println("  Trials completed:        $(result.ntrials)")
    
    # Save results
    output_data = Dict(
        "L" => L,
        "P_x" => P_x,
        "P_zz" => P_zz,
        "lambda_x" => lambda_x,
        "lambda_zz" => lambda_zz,
        "binder_parameter" => result.B,
        "binder_mean_of_trials" => result.B_mean_of_trials,
        "binder_std_of_trials" => result.B_std_of_trials,
        "S2_bar" => result.S2_bar,
        "S4_bar" => result.S4_bar,
        "ntrials" => result.ntrials,
        "seed" => seed,
        "sample" => sample,
        "maxdim" => maxdim,
        "cutoff" => cutoff,
        "measurement_type" => "rbim_strong_dephasing",
        "notes" => "RBIM limit: λ_x=0, EA correlations [⟨Z_i Z_j⟩²]_J"
    )
    
    # Create output directory if needed
    mkpath("output")
    
    # Save to JSON
    output_file = "output/$(out_prefix).json"
    open(output_file, "w") do f
        JSON.print(f, output_data, 2)
    end
    
    println("\n" * "="^70)
    println("Results saved to: $output_file")
    println("="^70)
    
    return 0
end

# Run main
exit(main())
