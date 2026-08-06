#!/usr/bin/env julia

"""
Right Boundary T_max Test Scan: Measurement-Induced Phase Transition with Rényi-2 Binder

Purpose:
Test whether the absence of a Binder crossing is a finite-time effect by scanning
T_max = T_max_factor * L for T_max_factor in {4, 8} at fixed right-boundary physics:
  lambda_x = 0.7, lambda_zz = 0.0, q_x = q_zz = q.

Output schema is intentionally aligned with `run_right_boundary_scan.jl` so existing
post-processing expectations remain compatible.
"""

using ArgParse
using JSON
using Dates
using Statistics

include("renyi2_right_boundary_core.jl")

function parse_commandline()
    s = ArgParseSettings(description = "Right boundary T_max test scan: vary q=P_x=P_zz with explicit T_max_factor")

    @add_arg_table! s begin
        "--L"
            help = "System size"
            arg_type = Int
            default = 16

        "--lambda_x"
            help = "X measurement strength (FIXED, must be 0.7 for the right edge)"
            arg_type = Float64
            default = 0.7

        "--lambda_zz"
            help = "ZZ measurement strength (FIXED, must be 0.0 for the right edge)"
            arg_type = Float64
            default = 0.0

        "--P_min"
            help = "Minimum dephasing strength q (q_x = q_zz = q)"
            arg_type = Float64
            default = 0.25

        "--P_max"
            help = "Maximum dephasing strength q (q_x = q_zz = q)"
            arg_type = Float64
            default = 0.45

        "--P_steps"
            help = "Number of q values to scan"
            arg_type = Int
            default = 5

        "--ntrials"
            help = "Number of Born-sampled Monte Carlo trajectories"
            arg_type = Int
            default = 100

        "--maxdim"
            help = "Maximum bond dimension for the dynamics"
            arg_type = Int
            default = 256

        "--cutoff"
            help = "SVD truncation cutoff for the dynamics"
            arg_type = Float64
            default = 1e-12

        "--T_max_factor"
            help = "Number of full layers is T_max_factor * L (typically 4 or 8 for this test set)"
            arg_type = Int
            default = 4

        "--obs_maxdim_factor"
            help = "Observable bond dimension is obs_maxdim_factor * maxdim"
            arg_type = Int
            default = 4

        "--obs_cutoff"
            help = "SVD truncation cutoff for the observable computation"
            arg_type = Float64
            default = 1e-14

        "--nboot"
            help = "Number of bootstrap resamples for the standard error of B"
            arg_type = Int
            default = 1000

        "--seed"
            help = "Random seed"
            arg_type = Int
            default = 42

        "--output_dir"
            help = "Output directory"
            arg_type = String
            default = "right_boundary_tmax_results"

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
    lambda_x = args["lambda_x"]
    lambda_zz = args["lambda_zz"]
    P_min = args["P_min"]
    P_max = args["P_max"]
    P_steps = args["P_steps"]
    ntrials = args["ntrials"]
    maxdim = args["maxdim"]
    cutoff = args["cutoff"]
    T_max_factor = args["T_max_factor"]
    obs_maxdim_factor = args["obs_maxdim_factor"]
    obs_cutoff = args["obs_cutoff"]
    nboot = args["nboot"]
    seed = args["seed"]
    output_dir = args["output_dir"]
    output_file = args["output_file"]

    @assert isapprox(lambda_x, 0.7; atol=1e-9) "Right-edge test requires lambda_x = 0.7."
    @assert isapprox(lambda_zz, 0.0; atol=1e-9) "Right-edge test requires lambda_zz = 0.0."
    @assert T_max_factor > 0 "T_max_factor must be positive."

    println("="^70)
    println("RIGHT BOUNDARY T_max TEST SCAN")
    println("="^70)
    println("Physics: lambda_x=0.7, lambda_zz=0.0, q_x=q_zz=q")
    println("Testing finite-time effects with T_max_factor=$T_max_factor")
    println()
    println("Parameters:")
    println("  System size L = $L")
    println("  q scan: q in [$P_min, $P_max] with $P_steps points")
    println("  Trajectories per point: $ntrials")
    println("  T_max = $T_max_factor * L = $(T_max_factor * L)")
    println("  Dynamics: maxdim=$maxdim, cutoff=$cutoff")
    println("  Observable: maxdim=$(obs_maxdim_factor * maxdim), cutoff=$obs_cutoff")
    println("="^70)
    println()

    q_values = P_steps <= 1 ? [P_min] : collect(range(P_min, P_max, length=P_steps))

    mkpath(output_dir)

    results = Dict{String,Any}[]

    for (i, q) in enumerate(q_values)
        println("\n[$i/$(length(q_values))] Running q = $(round(q, digits=4))")
        println("-"^70)

        t_start = time()

        result = rb_run_right_edge_point(
            L, Float64(q);
            lambda_x = lambda_x,
            lambda_zz = lambda_zz,
            ntrials = ntrials,
            T_max = T_max_factor * L,
            maxdim = maxdim,
            cutoff = cutoff,
            obs_maxdim = obs_maxdim_factor * maxdim,
            obs_cutoff = obs_cutoff,
            seed = seed + i,
            nboot = nboot,
        )

        t_elapsed = time() - t_start

        result_dict = Dict(
            "L" => L,
            "lambda_x" => lambda_x,
            "lambda_zz" => lambda_zz,
            "P_x" => q,
            "P_zz" => q,
            "T_max" => T_max_factor * L,
            "T_max_factor" => T_max_factor,
            "B" => result.B,
            "B_mean_of_trials" => result.B_mean_of_trials,
            "B_std_of_trials" => result.B_std_of_trials,
            "M2_bar" => result.M2_bar,
            "M4_bar" => result.M4_bar,
            "purity_bar" => result.purity_bar,
            "B_bootstrap_se" => result.B_bootstrap_se,
            "B_ci_low" => result.B_ci_low,
            "B_ci_high" => result.B_ci_high,
            "B2_ratio_of_mean_moments" => result.B2_ratio_of_mean_moments,
            "ntrials" => result.ntrials,
            "n_valid" => result.n_valid,
            "n_invalid" => result.n_invalid,
            "maxdim" => maxdim,
            "cutoff" => cutoff,
            "max_interphysical_linkdim" => result.max_interphysical_linkdim,
            "max_trace_error" => result.max_trace_error,
            "time_seconds" => t_elapsed
        )

        push!(results, result_dict)

        println("  B = $(round(result.B, digits=4)) +- $(round(result.B_bootstrap_se, digits=4))")
        println("  M2_bar = $(round(result.M2_bar, digits=6))")
        println("  M4_bar = $(round(result.M4_bar, digits=6))")
        println("  Purity = $(round(result.purity_bar, digits=4))")
        println("  Valid/Total = $(result.n_valid)/$(result.ntrials)")
        println("  Time: $(round(t_elapsed, digits=1)) seconds")

        if !isempty(output_file)
            outpath = joinpath(output_dir, output_file)
        else
            timestamp = Dates.format(now(), "yyyymmdd_HHMM")
            outpath = joinpath(output_dir, "right_boundary_tmax_L$(L)_Tf$(T_max_factor)_$(timestamp).json")
        end

        open(outpath, "w") do f
            JSON.print(f, results, 4)
        end
    end

    println("\n" * "="^70)
    println("SCAN COMPLETE")
    println("="^70)
    println("Total results: $(length(results))")
    println("Output saved to: $output_dir/")

    if isempty(output_file)
        timestamp = Dates.format(Dates.now(), "yyyymmdd_HHMM")
        final_output = joinpath(output_dir, "right_boundary_tmax_L$(L)_Tf$(T_max_factor)_$(timestamp)_final.json")
        open(final_output, "w") do f
            JSON.print(f, results, 4)
        end
        println("Final results: $final_output")
    end
end

main()
