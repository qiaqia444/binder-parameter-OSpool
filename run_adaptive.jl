#!/usr/bin/env julia

using Pkg; Pkg.activate(".")
using ITensors, ITensorMPS, ITensorCorrelators
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
            help = "Sample number"
            arg_type = Int
            default = 1
        "--out_prefix"
            help = "Output file prefix"
            arg_type = String
            default = "result"
        "--outdir"
            help = "Output directory"
            arg_type = String
            default = "output"
    end

    return parse_args(s)
end

# Adaptive parameters for different system sizes
function get_adaptive_params(L)
    if L <= 12
        # Conservative parameters for small systems
        maxdim = 256
        cutoff = 1e-12
        chunk4 = 100000
    elseif L <= 16
        # Extremely conservative for L=16
        maxdim = 8
        cutoff = 1e-6
        chunk4 = 500
    else  # L >= 20
        # Minimal parameters for large systems
        maxdim = 4
        cutoff = 1e-4
        chunk4 = 100
    end
    
    return (maxdim=maxdim, cutoff=cutoff, chunk4=chunk4)
end

function main(args)
    # Parse arguments
    parsed_args = parse_commandline()
    
    L = parsed_args["L"]
    lambda_x = parsed_args["lambda_x"]
    lambda_zz = parsed_args["lambda_zz"]
    lambda = parsed_args["lambda"]
    maxdim = parsed_args["maxdim"]
    cutoff = parsed_args["cutoff"]
    ntrials = parsed_args["ntrials"]
    chunk4 = parsed_args["chunk4"]
    seed = parsed_args["seed"]
    sample = parsed_args["sample"]
    out_prefix = parsed_args["out_prefix"]
    outdir = parsed_args["outdir"]
    
    println("Starting Edwards-Anderson Binder parameter calculation...")
    println("L = $L, λₓ = $lambda_x, λ_zz = $lambda_zz")
    
    # Get adaptive parameters for this system size
    adaptive_params = get_adaptive_params(L)
    
    # Override with adaptive parameters if defaults were used
    if maxdim == 256 && L > 12
        maxdim = adaptive_params.maxdim
        println("Adaptive: Using maxdim = $maxdim for L = $L")
    end
    
    if cutoff == 1e-12 && L > 12
        cutoff = adaptive_params.cutoff
        println("Adaptive: Using cutoff = $cutoff for L = $L")
    end
    
    if chunk4 == 50000 && L > 12
        chunk4 = adaptive_params.chunk4
        println("Adaptive: Using chunk4 = $chunk4 for L = $L")
    end
    
    try
        # Run the simulation with adaptive parameters
        result = BinderSim.ea_binder_mc(L; 
            lambda_x=lambda_x, 
            lambda_zz=lambda_zz, 
            ntrials=ntrials, 
            maxdim=maxdim, 
            cutoff=cutoff, 
            chunk4=chunk4, 
            seed=seed
        )
        
        # Prepare output data
        output_data = Dict(
            "L" => L,
            "lambda_x" => lambda_x,
            "lambda_zz" => lambda_zz,
            "lambda" => lambda,
            "maxdim" => maxdim,
            "cutoff" => cutoff,
            "ntrials" => ntrials,
            "ntrials_completed" => result.ntrials,
            "chunk4" => chunk4,
            "seed" => seed,
            "sample" => sample,
            "out_prefix" => out_prefix,
            "binder" => result.B,
            "binder_mean_of_trials" => result.B_mean_of_trials,
            "binder_std_of_trials" => result.B_std_of_trials,
            "S2_bar" => result.S2_bar,
            "S4_bar" => result.S4_bar
        )
        
        # Create output directory if it doesn't exist
        mkpath(outdir)
        
        # Save results
        output_file = joinpath(outdir, "$(out_prefix).json")
        open(output_file, "w") do io
            JSON.print(io, output_data)
        end
        
        println("Success! Binder parameter = $(result.B)")
        println("Results saved to: $output_file")
        
    catch e
        println("ERROR: Simulation failed for L=$L")
        println("Error details: $e")
        
        # Create a failure record
        mkpath(outdir)
        failure_file = joinpath(outdir, "$(out_prefix)_FAILED.json")
        failure_data = Dict(
            "L" => L,
            "lambda_x" => lambda_x,
            "lambda_zz" => lambda_zz,
            "lambda" => lambda,
            "maxdim" => maxdim,
            "cutoff" => cutoff,
            "error" => string(e),
            "status" => "FAILED"
        )
        
        open(failure_file, "w") do io
            JSON.print(io, failure_data)
        end
        
        println("Failure record saved to: $failure_file")
        exit(1)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
end
