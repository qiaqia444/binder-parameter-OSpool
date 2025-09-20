#!/usr/bin/env julia

using ArgParse
using Pkg
Pkg.activate(".")

include("src/BinderSimForced.jl")
using .BinderSimForced
using JSON

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table! s begin
        "--L"
            help = "System size"
            arg_type = Int
            required = true
        "--lambda_x"
            help = "X measurement strength"
            arg_type = Float64
            required = true
        "--lambda_zz"
            help = "ZZ measurement strength"
            arg_type = Float64
            required = true
        "--lambda"
            help = "Overall lambda parameter"
            arg_type = Float64
            required = true
        "--ntrials"
            help = "Number of trials (should be 1 for forced measurements)"
            arg_type = Int
            default = 1
        "--seed"
            help = "Random seed (for reproducibility)"
            arg_type = Int
            default = 1234
        "--sample"
            help = "Sample number"
            arg_type = Int
            default = 1
        "--out_prefix"
            help = "Output file prefix"
            arg_type = String
            default = "forced_result"
        "--outdir"
            help = "Output directory"
            arg_type = String
            default = "output"
    end
    return parse_args(s)
end

function get_adaptive_params(L)
    if L <= 12
        # Conservative parameters for small systems
        maxdim = 256
        cutoff = 1e-12
    elseif L <= 16
        # Moderate parameters for medium systems
        maxdim = 128
        cutoff = 1e-10
    elseif L <= 20
        # More aggressive for larger systems
        maxdim = 64
        cutoff = 1e-8
    else  # L >= 24
        # Very aggressive for very large systems
        maxdim = 32
        cutoff = 1e-6
    end
    
    return (maxdim=maxdim, cutoff=cutoff)
end

function main()
    args = parse_commandline()
    
    println("Starting Edwards-Anderson Binder parameter calculation...")
    println("FORCED +1 MEASUREMENTS (Deterministic)")
    println("L = $(args["L"]), λₓ = $(args["lambda_x"]), λ_zz = $(args["lambda_zz"])")
    
    # Get adaptive parameters based on system size
    adaptive_params = get_adaptive_params(args["L"])
    println("Adaptive: Using maxdim = $(adaptive_params.maxdim) for L = $(args["L"])")
    println("Adaptive: Using cutoff = $(adaptive_params.cutoff) for L = $(args["L"])")
    
    # Create output directory if it doesn't exist
    if !isdir(args["outdir"])
        mkpath(args["outdir"])
    end
    
    try
        # Run forced +1 measurement simulation
        result = ea_binder_mc_forced(args["L"]; 
                                   lambda_x=args["lambda_x"], 
                                   lambda_zz=args["lambda_zz"],
                                   ntrials=args["ntrials"],
                                   maxdim=adaptive_params.maxdim,
                                   cutoff=adaptive_params.cutoff,
                                   seed=args["seed"])
        
        # Prepare output data
        output_data = Dict(
            "L" => args["L"],
            "lambda_x" => args["lambda_x"],
            "lambda_zz" => args["lambda_zz"],
            "lambda" => args["lambda"],
            "ntrials" => args["ntrials"],
            "seed" => args["seed"],
            "sample" => args["sample"],
            "maxdim" => adaptive_params.maxdim,
            "cutoff" => adaptive_params.cutoff,
            "binder_parameter" => result.B,
            "binder_mean_of_trials" => result.B_mean_of_trials,
            "binder_std_of_trials" => result.B_std_of_trials,
            "S2_bar" => result.S2_bar,
            "S4_bar" => result.S4_bar,
            "ntrials_completed" => result.ntrials,
            "measurement_type" => "forced_plus_one",
            "success" => true
        )
        
        # Save results
        output_file = joinpath(args["outdir"], "$(args["out_prefix"]).json")
        open(output_file, "w") do io
            JSON.print(io, output_data, 2)
        end
        
        println("Success! Binder parameter = $(result.B)")
        println("Results saved to: $output_file")
        
    catch e
        println("ERROR: Simulation failed for L=$(args["L"])")
        println("Error details: $e")
        
        # Save failure record
        failure_data = Dict(
            "L" => args["L"],
            "lambda_x" => args["lambda_x"],
            "lambda_zz" => args["lambda_zz"],
            "lambda" => args["lambda"],
            "seed" => args["seed"],
            "sample" => args["sample"],
            "maxdim" => adaptive_params.maxdim,
            "cutoff" => adaptive_params.cutoff,
            "measurement_type" => "forced_plus_one",
            "success" => false,
            "error" => string(e)
        )
        
        failure_file = joinpath(args["outdir"], "$(args["out_prefix"])_FAILED.json")
        open(failure_file, "w") do io
            JSON.print(io, failure_data, 2)
        end
        
        println("Failure record saved to: $failure_file")
        exit(1)
    end
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    main()
end