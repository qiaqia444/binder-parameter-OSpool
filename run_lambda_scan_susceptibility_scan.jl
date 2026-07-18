#!/usr/bin/env julia

"""
Lambda-Path Scan: Rényi-2 / Edwards-Anderson Overlap Susceptibilities
(ITensorCorrelators.jl)

Scan the interior phase-diagram path

    lambda_x  = delta * lambda
    lambda_zz = delta * (1 - lambda)
    q_x = q_zz = q  (FIXED)

from lambda = 0 (pure ZZ measurement, left edge) to lambda = 1 (pure X
measurement, right edge), at fixed dephasing q. This is the cluster entry
point sharing `lambda_scan_susceptibilities_itensorcorrelators.jl` with the
local convenience script `run_lambda_scan_susceptibility_simple.jl`.

Observables (both computed per Born-sampled trajectory, then trajectory-averaged):

    kappa_2  = [L + 2 sum_{i<j} <r_i r_j>] / L,   r_i = Z_i^bra Z_i^ket
               (Rényi-2 replica-overlap susceptibility, ITensorCorrelators.jl)

    kappa_EA = [L + 2 sum_{i<j} <Z_i Z_j>^2] / L
               (Edwards-Anderson susceptibility, ordinary trace expectations)
"""

using ArgParse
using JSON
using Dates
using Statistics

include("lambda_scan_susceptibilities_itensorcorrelators.jl")

function parse_commandline()
    s = ArgParseSettings(description = "Lambda-path scan: vary lambda at fixed delta, q; Rényi-2/EA overlap susceptibilities (ITensorCorrelators.jl)")

    @add_arg_table! s begin
        "--L"
            help = "System size"
            arg_type = Int
            default = 12

        "--lambda_min"
            help = "Minimum path parameter lambda (0 = pure ZZ, 1 = pure X)"
            arg_type = Float64
            default = 0.0

        "--lambda_max"
            help = "Maximum path parameter lambda"
            arg_type = Float64
            default = 1.0

        "--lambda_steps"
            help = "Number of lambda values to scan"
            arg_type = Int
            default = 11

        "--delta"
            help = "Path amplitude: lambda_x = delta*lambda, lambda_zz = delta*(1-lambda)"
            arg_type = Float64
            default = 0.7

        "--q"
            help = "Fixed dephasing strength (q_x = q_zz = q)"
            arg_type = Float64
            default = 0.1

        "--ntrials"
            help = "Number of Born-sampled Monte Carlo trajectories"
            arg_type = Int
            default = 20

        "--maxdim"
            help = "Maximum bond dimension for the dynamics"
            arg_type = Int
            default = 256

        "--cutoff"
            help = "SVD truncation cutoff for the dynamics"
            arg_type = Float64
            default = 1e-10

        "--T_max_factor"
            help = "Number of full layers is T_max_factor * L"
            arg_type = Int
            default = 2

        "--pair_batch_size"
            help = "Batch size for ITensorCorrelators.jl pair-correlator evaluation"
            arg_type = Int
            default = 4096

        "--skip_ea"
            help = "Skip the Edwards-Anderson susceptibility (kappa_EA) computation"
            action = :store_true

        "--seed"
            help = "Random seed"
            arg_type = Int
            default = 1234

        "--output_dir"
            help = "Output directory"
            arg_type = String
            default = "lambda_scan_susceptibility_results"

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
    lambda_min = args["lambda_min"]
    lambda_max = args["lambda_max"]
    lambda_steps = args["lambda_steps"]
    delta = args["delta"]
    q = args["q"]
    ntrials = args["ntrials"]
    maxdim = args["maxdim"]
    cutoff = args["cutoff"]
    T_max_factor = args["T_max_factor"]
    pair_batch_size = args["pair_batch_size"]
    compute_EA = !args["skip_ea"]
    seed = args["seed"]
    output_dir = args["output_dir"]
    output_file = args["output_file"]

    println("="^70)
    println("LAMBDA-PATH SCAN: Rényi-2 / Edwards-Anderson Overlap Susceptibilities")
    println("="^70)
    println("Physics: lambda_x = delta*lambda, lambda_zz = delta*(1-lambda), q_x=q_zz=q")
    println("Method: doubled-MPS Born sampling")
    println("Observables: trajectory-averaged kappa_2 (ITensorCorrelators.jl) and kappa_EA")
    println()
    println("Parameters:")
    println("  System size L = $L")
    println("  delta (path amplitude) = $delta")
    println("  q (FIXED dephasing) = $q")
    println("  lambda scan: [$lambda_min, $lambda_max] with $lambda_steps points")
    println("  Trajectories per point: $ntrials")
    println("  T_max = $T_max_factor * L = $(T_max_factor * L)")
    println("  Dynamics: maxdim=$maxdim, cutoff=$cutoff")
    println("  compute_EA = $compute_EA")
    println("="^70)
    println()

    # Create lambda values to scan
    lambda_values = lambda_steps <= 1 ? [lambda_min] : collect(range(lambda_min, lambda_max, length=lambda_steps))

    # Create output directory
    mkpath(output_dir)

    results = []

    for (i, lambda) in enumerate(lambda_values)
        lambda_x = delta * lambda
        lambda_zz = delta * (1.0 - lambda)

        println("\n[$i/$(length(lambda_values))] Running lambda = $(round(lambda, digits=4))")
        println("  (lambda_x = $(round(lambda_x, digits=4)), lambda_zz = $(round(lambda_zz, digits=4)), q = $q)")
        println("-"^70)

        t_start = time()

        result = ls_run_lambda_point(
            L, lambda;
            delta = delta,
            q = q,
            ntrials = ntrials,
            T_max = T_max_factor * L,
            maxdim = maxdim,
            cutoff = cutoff,
            pair_batch_size = pair_batch_size,
            seed = seed + i,  # Different seed for each lambda
            compute_EA = compute_EA,
        )

        t_elapsed = time() - t_start

        # Store results (trial-level arrays are intentionally not persisted)
        result_dict = Dict(
            "L" => L,
            "lambda" => lambda,
            "lambda_x" => result.lambda_x,
            "lambda_zz" => result.lambda_zz,
            "delta" => delta,
            "q" => q,
            "kappa2" => result.kappa2,
            "kappa2_se" => result.kappa2_se,
            "kappaEA" => result.kappaEA,
            "kappaEA_se" => result.kappaEA_se,
            "purity" => result.purity,
            "ntrials" => result.ntrials,
            "maxdim" => maxdim,
            "cutoff" => cutoff,
            "T_max" => result.T_max,
            "max_interphysical_linkdim" => result.max_interphysical_linkdim,
            "max_trace_error" => result.max_trace_error,
            "observable_method" => "itensorcorrelators",
            "time_seconds" => t_elapsed
        )

        push!(results, result_dict)

        # Print results
        println("  kappa2  = $(round(result.kappa2, digits=4)) ± $(round(result.kappa2_se, digits=4))")
        if compute_EA
            println("  kappaEA = $(round(result.kappaEA, digits=4)) ± $(round(result.kappaEA_se, digits=4))")
        end
        println("  Purity = $(round(result.purity, digits=4))")
        println("  Max trace error = $(result.max_trace_error)")
        println("  Time: $(round(t_elapsed, digits=1)) seconds")

        # Save intermediate results
        if !isempty(output_file)
            # Use specified filename for cluster jobs
            outpath = joinpath(output_dir, output_file)
        else
            # Use timestamped filename for local runs
            timestamp = Dates.format(now(), "yyyymmdd_HHMM")
            outpath = joinpath(output_dir, "lambda_scan_susceptibility_L$(L)_delta$(delta)_q$(q)_$(timestamp).json")
        end

        open(outpath, "w") do f
            JSON.print(f, results, 4)
        end
    end

    println("\n" * "="^70)
    println("✓ SCAN COMPLETE")
    println("="^70)
    println("Total results: $(length(results))")
    println("Results saved to: $output_dir/")
end

main()
