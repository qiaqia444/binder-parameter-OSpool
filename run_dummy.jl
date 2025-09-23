#!/usr/bin/env julia

using Pkg; Pkg.activate(".")
using ITensors, ITensorMPS
using JSON, Statistics, Random
using ArgParse

# Include our simulation module
include("src/BinderSimWithDummy.jl")
using .BinderSimWithDummy

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table! s begin
        "--L"
            help = "System size (physical sites, dummy will be added)"
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
        "--T_max"
            help = "Number of timesteps"
            arg_type = Int
            default = 20
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
            default = "binder_dummy"
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
    println("BinderSim with Dummy Site Simulation")
    println("========================================")
    println("Physical system size L = $(args["L"])")
    println("Total system size = $(args["L"] + 1) (includes dummy)")
    println("λₓ = $(args["lambda_x"]), λ_zz = $(args["lambda_zz"])")
    println("Max bond dimension = $(args["maxdim"])")
    println("Number of trials = $(args["ntrials"])")
    println("T_max = $(args["T_max"])")
    println("Random seed = $(args["seed"])")
    println("Sample index = $(args["sample"])")
    println("========================================")
    
    # Run multiple trials and collect results
    correlations = Float64[]
    
    for trial in 1:args["ntrials"]
        # Set unique seed for each trial
        trial_seed = args["seed"] + trial - 1
        
        # Run single simulation
        correlation = ea_binder_mc_with_dummy(
            args["L"];
            lambda_x = args["lambda_x"],
            lambda_zz = args["lambda_zz"],
            T_max = args["T_max"],
            maxdim = args["maxdim"],
            cutoff = args["cutoff"],
            seed = trial_seed
        )
        
        push!(correlations, correlation)
        
        if trial % 100 == 0
            println("Completed trial $trial/$(args["ntrials"])")
        end
    end
    
    # Calculate statistics
    mean_correlation = mean(correlations)
    std_correlation = std(correlations)
    
    # Calculate Binder-like parameter from correlation distribution
    corr2_mean = mean(correlations.^2)
    corr4_mean = mean(correlations.^4)
    binder_parameter = 1 - corr4_mean / (3 * corr2_mean^2)
    
    # Calculate chi-like parameter
    chi_parameter = args["L"] * var(correlations)
    
    # Prepare output data
    output_data = Dict(
        "method" => "dummy_site_binderSim",
        "L_physical" => args["L"],
        "L_total" => args["L"] + 1,
        "lambda_x" => args["lambda_x"],
        "lambda_zz" => args["lambda_zz"],
        "lambda" => args["lambda"],
        "maxdim" => args["maxdim"],
        "cutoff" => args["cutoff"],
        "ntrials" => args["ntrials"],
        "T_max" => args["T_max"],
        "seed" => args["seed"],
        "sample" => args["sample"],
        "mean_correlation" => mean_correlation,
        "std_correlation" => std_correlation,
        "binder_parameter" => binder_parameter,
        "chi_parameter" => chi_parameter,
        "all_correlations" => correlations[1:min(1000, length(correlations))]  # Save first 1000 for analysis
    )
    
    # Create output directory if it doesn't exist
    mkpath(args["outdir"])
    
    # Write results to file with distinct naming
    output_filename = joinpath(args["outdir"], "dummy_$(args["out_prefix"]).json")
    open(output_filename, "w") do f
        JSON.print(f, output_data, 2)
    end
    
    println("\nDummy site simulation completed successfully!")
    println("Results saved to: $output_filename")
    println("Mean correlation: $(round(mean_correlation, digits=6))")
    println("Binder parameter: $(round(binder_parameter, digits=6))")
    println("Chi parameter: $(round(chi_parameter, digits=6))")
    
    return 0
end

# Run the main function
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end