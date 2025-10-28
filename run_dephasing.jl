#!/usr/bin/env julia

using Pkg; Pkg.activate(".")
using ITensors, ITensorMPS
using JSON, Statistics, Random
using ArgParse

# Include our simulation modules
include("src/BinderSimDephasing.jl")
using .BinderSimDephasing

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table! s begin
        "--L"
            help = "System size"
            arg_type = Int
            default = 12
        "--P_x"
            help = "X dephasing probability"
            arg_type = Float64
            default = 0.0
        "--lambda_x"
            help = "Lambda_x parameter"
            arg_type = Float64
            default = 0.5
        "--lambda_zz"
            help = "Lambda_zz parameter"
            arg_type = Float64
            default = 0.5
        "--lambda"
            help = "Lambda value (sets lambda_x=lambda, lambda_zz=1-lambda)"
            arg_type = Float64
            default = -1.0
        "--maxdim"
            help = "Maximum bond dimension"
            arg_type = Int
            default = 256
        "--cutoff"
            help = "Cutoff for SVD truncation"
            arg_type = Float64
            default = 1e-12
        "--ntrials"
            help = "Number of Monte Carlo trials"
            arg_type = Int
            default = 1000
        "--chunk4"
            help = "Chunk size for 4-point correlators"
            arg_type = Int
            default = 50000
        "--seed"
            help = "Random seed"
            arg_type = Int
            default = 1234
        "--sample"
            help = "Sample index"
            arg_type = Int
            default = 1
        "--out_prefix"
            help = "Output file prefix"
            arg_type = String
            default = "binder"
        "--outdir"
            help = "Output directory"
            arg_type = String
            default = "output"
    end

    return parse_args(s)
end

function main()
    # Parse command line arguments
    args = parse_commandline()
    
    # Handle lambda parameter
    lambda_x = args["lambda_x"]
    lambda_zz = args["lambda_zz"]
    if args["lambda"] >= 0.0
        lambda_x = args["lambda"]
        lambda_zz = 1.0 - args["lambda"]
    end
    
    println("========================================")
    println("Dephasing Channel Binder Parameter Simulation")
    println("========================================")
    println("System size L = $(args["L"])")
    println("P_x = $(args["P_x"]), P_zz = 0.01")  # Fixed P_zz
    println("lambda_x = $(lambda_x), lambda_zz = $(lambda_zz)")
    println("Max bond dimension = $(args["maxdim"])")
    println("Number of trials = $(args["ntrials"])")
    println("Random seed = $(args["seed"])")
    println("Sample index = $(args["sample"])")
    println("========================================")
    
    # Set random seed
    Random.seed!(args["seed"])
    
    # Run the simulation with dephasing
    result = ea_binder_mc_dephasing(
        args["L"];
        lambda_x = lambda_x,
        lambda_zz = lambda_zz,
        P_x = args["P_x"],
        P_zz = 0.01,  # Fixed P_zz
        maxdim = args["maxdim"],
        cutoff = args["cutoff"],
        ntrials = args["ntrials"],
        chunk4 = args["chunk4"]
    )
    
    # Prepare output data
    output_data = Dict(
        "L" => args["L"],
        "lambda_x" => lambda_x,
        "lambda_zz" => lambda_zz,
        "P_x" => args["P_x"],
        "P_zz" => 0.01,  # Fixed P_zz
        "ntrials" => args["ntrials"],
        "seed" => args["seed"],
        "sample" => args["sample"],
        "maxdim" => args["maxdim"],
        "cutoff" => args["cutoff"],
        "binder_parameter" => result.B,
        "binder_mean_of_trials" => result.B_mean_of_trials,
        "binder_std_of_trials" => result.B_std_of_trials,
        "S2_bar" => result.S2_bar,
        "S4_bar" => result.S4_bar,
        "ntrials_completed" => result.ntrials,
        "measurement_type" => "dephasing_channel",
        "success" => true
    )
    
    # Create output directory if it doesn't exist
    mkpath(args["outdir"])
    
    # Write results to file with distinct naming
    output_filename = joinpath(args["outdir"], "dephasing_$(args["out_prefix"]).json")
    open(output_filename, "w") do f
        JSON.print(f, output_data, 2)
    end
    
    println("Success! Binder parameter = $(result.B)")
    println("Results saved to: $output_filename")
    
    return 0
end

# Run the main function
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
