#!/usr/bin/env julia

"""
Dephasing channel simulation with P = 0.8
Run weak measurement + dephasing simulations and compute EA Binder parameter.
"""

using ArgParse
using Random
using JSON
using ITensors, ITensorMPS

include("src/BinderSim.jl")
using .BinderSim

include("src/BinderSimDephasing.jl")
using .BinderSimDephasing

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table! s begin
        "--L"
            help = "System size"
            arg_type = Int
            required = true
        "--P_x"
            help = "X dephasing probability"
            arg_type = Float64
            default = 0.8
        "--lambda_x"
            help = "X measurement strength"
            arg_type = Float64
            default = -1.0
        "--lambda_zz"
            help = "ZZ measurement strength"
            arg_type = Float64
            default = -1.0
        "--lambda"
            help = "Overall lambda (if >= 0, sets lambda_x=lambda, lambda_zz=1-lambda)"
            arg_type = Float64
            default = 0.5
        "--ntrials"
            help = "Number of Monte Carlo trials"
            arg_type = Int
            default = 2000
        "--seed"
            help = "Random seed"
            arg_type = Int
            default = 42
        "--sample"
            help = "Sample number for multiple runs"
            arg_type = Int
            default = 1
        "--out_prefix"
            help = "Output filename prefix"
            arg_type = String
            default = "dephasing_p08"
    end
    return parse_args(s)
end

function main()
    args = parse_commandline()
    
    L = args["L"]
    P_x = args["P_x"]
    lambda = args["lambda"]
    ntrials = args["ntrials"]
    seed_val = args["seed"]
    sample = args["sample"]
    out_prefix = args["out_prefix"]
    
    # Set lambda values
    if lambda >= 0
        lambda_x = lambda
        lambda_zz = 1.0 - lambda
    else
        lambda_x = args["lambda_x"]
        lambda_zz = args["lambda_zz"]
    end
    
    # P_zz hardcoded to 0.8 to match P_x
    P_zz = 0.8
    
    println("="^70)
    println("Dephasing Channel Simulation (P = 0.8)")
    println("="^70)
    println("System size L = $L")
    println("Lambda (X)    = $lambda_x")
    println("Lambda (ZZ)   = $lambda_zz")
    println("P_x           = $P_x")
    println("P_zz          = $P_zz")
    println("Trials        = $ntrials")
    println("Seed          = $seed_val")
    println("Sample        = $sample")
    println("="^70)
    
    # Run simulation
    result = ea_binder_mc_dephasing(
        L;
        lambda_x = lambda_x,
        lambda_zz = lambda_zz,
        P_x = P_x,
        P_zz = P_zz,
        ntrials = ntrials,
        seed = seed_val
    )
    
    # Prepare output
    output_data = Dict(
        "L" => L,
        "lambda_x" => lambda_x,
        "lambda_zz" => lambda_zz,
        "P_x" => P_x,
        "P_zz" => P_zz,
        "binder_parameter" => result.B,
        "binder_mean_of_trials" => result.B_mean_of_trials,
        "binder_std_of_trials" => result.B_std_of_trials,
        "S2_bar" => result.S2_bar,
        "S4_bar" => result.S4_bar,
        "ntrials" => result.ntrials,
        "seed" => seed_val,
        "sample" => sample,
        "measurement_type" => "dephasing_channel_p08"
    )
    
    # Create output directory if it doesn't exist
    mkpath("output")
    
    # Save results
    output_filename = "output/$(out_prefix).json"
    open(output_filename, "w") do f
        JSON.print(f, output_data, 2)
    end
    
    println("\nResults:")
    println("  EA Binder parameter: $(round(result.B, digits=6))")
    println("  Mean of trial Bs:    $(round(result.B_mean_of_trials, digits=6))")
    println("  Std of trial Bs:     $(round(result.B_std_of_trials, digits=6))")
    println("\nOutput saved to: $output_filename")
    println("="^70)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
