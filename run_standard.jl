#!/usr/bin/env julia

# Standard run script for L = [8, 10, 12] study using reliable ITensorCorrelators approach
# This uses the original BinderSim.jl module that works well for these system sizes

using Pkg
Pkg.activate(".")

include("src/BinderSim.jl")
using .BinderSim
using JSON
using Random

function parse_args()
    if length(ARGS) != 8
        println("Usage: julia run_standard.jl L lambda_x lambda_zz lambda ntrials seed sample out_prefix")
        println("Example: julia run_standard.jl 12 0.5 0.5 0.5 200 1001 1 std_L12_lam0.5_s1")
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

function get_standard_params(L)
    # Conservative parameters that work reliably for L ≤ 12
    if L <= 8
        maxdim = 256
        cutoff = 1e-12
        chunk4 = 100000
    elseif L <= 10
        maxdim = 256
        cutoff = 1e-12
        chunk4 = 80000
    else  # L = 12
        maxdim = 256
        cutoff = 1e-12
        chunk4 = 50000
    end
    
    return (maxdim=maxdim, cutoff=cutoff, chunk4=chunk4)
end

function main()
    L, lambda_x, lambda_zz, lambda, ntrials, seed, sample, out_prefix = parse_args()
    
    println("Starting Edwards-Anderson Binder parameter calculation...")
    println("L = $L, λₓ = $lambda_x, λ_zz = $lambda_zz")
    
    # Get optimized parameters for this system size
    params = get_standard_params(L)
    println("Using maxdim = $(params.maxdim), cutoff = $(params.cutoff), chunk4 = $(params.chunk4)")
    
    outdir = get(ENV, "OUTDIR", "output")
    mkpath(outdir)
    
    try
        # Run simulation with ITensorCorrelators (reliable for L ≤ 12)
        result = ea_binder_mc(L; 
                             lambda_x=lambda_x, 
                             lambda_zz=lambda_zz,
                             ntrials=ntrials,
                             maxdim=params.maxdim,
                             cutoff=params.cutoff,
                             chunk4=params.chunk4,
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
            "chunk4" => params.chunk4,
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
        
        println("Success! Binder parameter = $(result.B)")
        println("Results saved to: $output_file")
        
    catch e
        println("ERROR: Simulation failed for L=$L")
        println("Error details: $e")
        
        # Save failure record
        failure_data = Dict(
            "L" => L,
            "lambda_x" => lambda_x,
            "lambda_zz" => lambda_zz,
            "lambda" => lambda,
            "maxdim" => params.maxdim,
            "cutoff" => params.cutoff,
            "status" => "FAILED",
            "error" => string(e)
        )
        
        failure_file = joinpath(outdir, "$(out_prefix)_FAILED.json")
        open(failure_file, "w") do f
            JSON.print(f, failure_data, 2)
        end
        
        println("Failure record saved to: $failure_file")
        exit(1)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
