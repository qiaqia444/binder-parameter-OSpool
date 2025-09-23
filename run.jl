#!/usr/bin/env julia

using Pkg; Pkg.activate(".")
using ITensors, ITensorMPS
using JSON, Statistics, Random
using ArgParse

# Include our simulation module
include("src/BinderSim.jl")
using .BinderSim

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table! s begin
        "--L"
            help = "System size"
            arg_type = Int
            default = 12
        "--lambda_x"
            help = "Lambda_x parameter"
            arg_type = Float64
            default = 0.5
        "--lambda_zz"
            help = "Lambda_zz parameter"
            arg_type = Float64
            default = 0.5
        "--lambda"
            help = "Lambda parameter (for consistency)"
            arg_type = Float64
            default = 0.5
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
    
    println("========================================")
    println("Standard BinderSim Quantum Trajectory Simulation")
    println("========================================")
    println("System size L = $(args["L"])")
    println("λₓ = $(args["lambda_x"]), λ_zz = $(args["lambda_zz"])")
    println("Max bond dimension = $(args["maxdim"])")
    println("Number of trials = $(args["ntrials"])")
    println("Random seed = $(args["seed"])")
    println("Sample index = $(args["sample"])")
    println("========================================")
    
    # Set random seed
    Random.seed!(args["seed"])
    
    # Run the simulation
    result = ea_binder_mc(
        args["L"];
        lambda_x = args["lambda_x"],
        lambda_zz = args["lambda_zz"],
        maxdim = args["maxdim"],
        cutoff = args["cutoff"],
        n_trials = args["ntrials"],
        chunk4 = args["chunk4"]
    )
    
    # Prepare output data
    output_data = Dict(
        "method" => "standard_binderSim",
        "L" => args["L"],
        "lambda_x" => args["lambda_x"],
        "lambda_zz" => args["lambda_zz"],
        "lambda" => args["lambda"],
        "maxdim" => args["maxdim"],
        "cutoff" => args["cutoff"],
        "ntrials" => args["ntrials"],
        "seed" => args["seed"],
        "sample" => args["sample"],
        "binder_parameter" => result["binder_parameter"],
        "chi_parameter" => result["chi_parameter"],
        "energy" => result["energy"],
        "correlation_length" => result["correlation_length"],
        "magnetization" => result["magnetization"],
        "susceptibility" => result["susceptibility"]
    )
    
    # Create output directory if it doesn't exist
    mkpath(args["outdir"])
    
    # Write results to file with distinct naming
    output_filename = joinpath(args["outdir"], "standard_$(args["out_prefix"]).json")
    open(output_filename, "w") do f
        JSON.print(f, output_data, 2)
    end
    
    println("\nSimulation completed successfully!")
    println("Results saved to: $output_filename")
    println("Binder parameter: $(round(result["binder_parameter"], digits=6))")
    println("Chi parameter: $(round(result["chi_parameter"], digits=6))")
    
    return 0
end

# Run the main function
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end