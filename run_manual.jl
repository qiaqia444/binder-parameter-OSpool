#!/usr/bin/env julia

# Manual correlators run script for L = [8, 10, 12] study WITHOUT ITensorCorrelators
# This uses manual tensor contractions to bypass ITensorCorrelators limitations

using Pkg
Pkg.activate(".")

include("src/BinderSimManual.jl")
using .BinderSimManual
using JSON
using Random

function parse_args()
    if length(ARGS) != 8
        println("Usage: julia run_manual.jl L lambda_x lambda_zz lambda ntrials seed sample out_prefix")
        println("Example: julia run_manual.jl 12 0.5 0.5 0.5 200 1001 1 manual_L12_lam0.5_s1")
        exit(1)
    end
    
    L = parse(Int, ARGS[1])
    lambda_x = parse(Float64, ARGS[2])
    lambda_zz = parse(Float64, ARGS[3])
    lambda = parse(Float64, ARGS[4])
    ntrials = parse(Int, ARGS[5])
    seed = parse(Int, ARGS[6])
    sample = parse(Int, ARGS[7])
    out_prefix = ARGS[8]
    
    return L, lambda_x, lambda_zz, lambda, ntrials, seed, sample, out_prefix
end

function get_manual_params(L::Int)
    # Conservative parameters for manual correlators
    if L <= 8
        return (maxdim=256, cutoff=1e-12, manual_maxdim=128, manual_cutoff=1e-10, chunk_size=50)
    elseif L <= 10
        return (maxdim=256, cutoff=1e-12, manual_maxdim=96, manual_cutoff=1e-10, chunk_size=40)
    elseif L <= 12
        return (maxdim=256, cutoff=1e-12, manual_maxdim=64, manual_cutoff=1e-10, chunk_size=30)
    else
        return (maxdim=256, cutoff=1e-12, manual_maxdim=48, manual_cutoff=1e-10, chunk_size=20)
    end
end

function main()
    L, lambda_x, lambda_zz, lambda, ntrials, seed, sample, out_prefix = parse_args()
    
    println("Starting Edwards-Anderson Binder parameter calculation with MANUAL correlators...")
    println("L = $L, λₓ = $lambda_x, λ_zz = $lambda_zz")
    
    # Get optimized parameters for this system size
    params = get_manual_params(L)
    println("Using maxdim = $(params.maxdim), manual_maxdim = $(params.manual_maxdim), chunk_size = $(params.chunk_size)")
    
    outdir = get(ENV, "OUTDIR", "output")
    mkpath(outdir)
    
    try
        # Run simulation with manual correlators
        result = ea_binder_mc_manual(L; 
                                   lambda_x=lambda_x, 
                                   lambda_zz=lambda_zz,
                                   ntrials=ntrials,
                                   maxdim=params.maxdim,
                                   cutoff=params.cutoff,
                                   manual_maxdim=params.manual_maxdim,
                                   manual_cutoff=params.manual_cutoff,
                                   chunk_size=params.chunk_size,
                                   seed=seed)
        
        # Prepare output data
        output_data = Dict(
            "L" => L,
            "lambda_x" => lambda_x,
            "lambda_zz" => lambda_zz,
            "lambda" => lambda,
            "ntrials" => ntrials,
            "ntrials_completed" => result.ntrials,
            "seed" => seed,
            "sample" => sample,
            "out_prefix" => out_prefix,
            "maxdim" => params.maxdim,
            "cutoff" => params.cutoff,
            "manual_maxdim" => params.manual_maxdim,
            "manual_cutoff" => params.manual_cutoff,
            "chunk_size" => params.chunk_size,
            "method" => "manual_correlators",
            "binder" => result.B,
            "binder_mean_of_trials" => result.B_mean_of_trials,
            "binder_std_of_trials" => result.B_std_of_trials,
            "S2_bar" => result.S2_bar,
            "S4_bar" => result.S4_bar,
            "status" => "SUCCESS"
        )
        
        # Save results
        output_file = joinpath(outdir, "$(out_prefix).json")
        open(output_file, "w") do f
            JSON.print(f, output_data, 2)
        end
        
        println("Success! Binder parameter = $(result.B) (manual correlators)")
        println("Results saved to: $output_file")
        
    catch e
        println("ERROR: Manual correlator simulation failed for L=$L")
        println("Error details: $e")
        
        # Save failure record
        failure_data = Dict(
            "L" => L,
            "lambda_x" => lambda_x,
            "lambda_zz" => lambda_zz,
            "lambda" => lambda,
            "maxdim" => params.maxdim,
            "manual_maxdim" => params.manual_maxdim,
            "method" => "manual_correlators",
            "status" => "FAILED",
            "error" => string(e)
        )
        
        output_file = joinpath(outdir, "$(out_prefix)_FAILED.json")
        open(output_file, "w") do f
            JSON.print(f, failure_data, 2)
        end
        
        println("Failure record saved to: $output_file")
        exit(1)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
