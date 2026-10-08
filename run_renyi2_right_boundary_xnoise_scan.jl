#!/usr/bin/env julia

"""
Right Boundary Scan, Renyi-2 Binder measured AFTER THE LAST X-DEPHASING step

Same physics, parameters, CLI and output schema as `run_right_boundary_scan.jl`
(lambda_x = 0.7, lambda_zz = 0, q = P_x = P_zz scanned over [P_min, P_max], proposal-aligned
doubled-MPS Born sampling, replica-overlap Renyi-2 Binder, trajectory-averaged B), except that
the observable is evaluated on the state just BEFORE the last ZZ-dephasing layer, i.e. after the
X dephasing of step T_max - the point where the left-boundary EA code measures
(`ρ_after_X_noise`). The production `run_right_boundary_scan.jl` measures after the last ZZ
dephasing. See `renyi2_right_boundary_xnoise_core.jl`.

With the same --seed the trajectories are identical to those of `run_right_boundary_scan.jl`.

Smoke test:
  julia --project=. run_renyi2_right_boundary_xnoise_scan.jl --L 6 --P_min 0 --P_max 0.5 \\
      --P_steps 3 --ntrials 4 --maxdim 64 --nboot 100 --output_dir /tmp/renyi2_xnoise_smoke
"""

using ArgParse
using JSON
using Dates
using Statistics

include("renyi2_right_boundary_xnoise_core.jl")

function parse_commandline()
    s = ArgParseSettings(description = "Right boundary scan: Rényi-2 Binder after the last X dephasing, vs q=P_x=P_zz at fixed λ_x")

    @add_arg_table! s begin
        "--L"
            help = "System size"
            arg_type = Int
            default = 12

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
            default = 0.0

        "--P_max"
            help = "Maximum dephasing strength q (q_x = q_zz = q)"
            arg_type = Float64
            default = 0.5

        "--P_steps"
            help = "Number of q values to scan"
            arg_type = Int
            default = 11

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
            help = "Number of full layers is T_max_factor * L (observable after the X dephasing of the last layer)"
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
            default = "renyi2_right_boundary_xnoise_results"

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

    @assert isapprox(lambda_x, 0.7; atol=1e-9) "Proposal-aligned right edge requires lambda_x = 0.7 (delta=0.7, lambda=1)."
    @assert isapprox(lambda_zz, 0.0; atol=1e-9) "Proposal-aligned right edge requires lambda_zz = 0 (no weak ZZ measurement)."

    println("="^70)
    println("RIGHT BOUNDARY SCAN: Renyi-2 Binder after the LAST X-DEPHASING step")
    println("="^70)
    println("Physics: lambda=1, delta=0.7 => lambda_x=0.7, lambda_zz=0, q_x=q_zz=q")
    println("Method: doubled-MPS Born sampling (production dynamics, renyi2_right_boundary_core.jl)")
    println("Observable: trajectory-averaged replica-overlap Rényi-2 Binder,")
    println("            measured after the X dephasing of the last step (before the last ZZ dephasing)")
    println()
    println("Parameters:")
    println("  System size L = $L")
    println("  X measurement λ_x = $lambda_x (FIXED)")
    println("  ZZ measurement λ_zz = $lambda_zz (FIXED - no weak ZZ measurement)")
    println("  q scan: q_x = q_zz = q ∈ [$P_min, $P_max] with $P_steps points")
    println("  Trajectories per point: $ntrials")
    println("  T_max = $T_max_factor * L = $(T_max_factor * L)")
    println("  Dynamics: maxdim=$maxdim, cutoff=$cutoff")
    println("  Observable: maxdim=$(obs_maxdim_factor * maxdim), cutoff=$obs_cutoff")
    println("="^70)
    println()

    q_values = P_steps <= 1 ? [P_min] : collect(range(P_min, P_max, length=P_steps))

    mkpath(output_dir)

    results = []

    for (i, q) in enumerate(q_values)
        println("\n[$i/$(length(q_values))] Running q = $(round(q, digits=4)) (P_x = P_zz = q)")
        println("  (λ_x = $lambda_x, λ_zz = $lambda_zz)")
        println("-"^70)

        t_start = time()

        result = rb_xn_run_right_edge_point(
            L, Float64(q);
            lambda_x = lambda_x,
            lambda_zz = lambda_zz,
            ntrials = ntrials,
            T_max = T_max_factor * L,
            maxdim = maxdim,
            cutoff = cutoff,
            obs_maxdim = obs_maxdim_factor * maxdim,
            obs_cutoff = obs_cutoff,
            seed = seed + i,  # same per-q seed rule as run_right_boundary_scan.jl
            nboot = nboot,
        )

        t_elapsed = time() - t_start

        # Field names identical to run_right_boundary_scan.jl, plus the metadata below.
        result_dict = Dict(
            "L" => L,
            "lambda_x" => lambda_x,
            "lambda_zz" => lambda_zz,
            "P_x" => q,
            "P_zz" => q,
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
            "obs_maxdim" => obs_maxdim_factor * maxdim,
            "obs_cutoff" => obs_cutoff,
            "T_max" => T_max_factor * L,
            "max_interphysical_linkdim" => result.max_interphysical_linkdim,
            "max_trace_error" => result.max_trace_error,
            "time_seconds" => t_elapsed,
            "observable" => "renyi2_binder",
            "measurement_time" => "after the X dephasing of the last step t=T_max (before the last ZZ-dephasing layer)",
            "seed" => seed,
            "point_seed" => seed + i,
            "trial_seeds" => result.trial_seeds,
            "M2_per_trajectory" => result.M2_per_trajectory,
            "M4_per_trajectory" => result.M4_per_trajectory,
            "B2_per_trajectory" => result.B2_per_trajectory,
            "purity_per_trajectory" => result.purity_per_trajectory,
        )

        push!(results, result_dict)

        println("  B (trajectory-averaged) = $(round(result.B, digits=4)) ± $(round(result.B_bootstrap_se, digits=4)) (bootstrap SE)")
        println("  B (ratio of mean moments, diagnostic) = $(round(result.B2_ratio_of_mean_moments, digits=4))")
        println("  M₂_bar = $(round(result.M2_bar, digits=6))")
        println("  M₄_bar = $(round(result.M4_bar, digits=6))")
        println("  Purity = $(round(result.purity_bar, digits=4))")
        println("  Valid/Total = $(result.n_valid)/$(result.ntrials)")
        println("  Max trace error = $(result.max_trace_error)")
        println("  Time: $(round(t_elapsed, digits=1)) seconds")

        if !isempty(output_file)
            outpath = joinpath(output_dir, output_file)
        else
            timestamp = Dates.format(now(), "yyyymmdd_HHMM")
            outpath = joinpath(output_dir, "renyi2_xnoise_right_boundary_L$(L)_lambda$(lambda_x)_$(timestamp).json")
        end

        # Saved after every q so an interrupted scan keeps its finished points
        open(outpath, "w") do f
            JSON.print(f, results, 4)
        end
    end

    println("\n" * "="^70)
    println("✓ SCAN COMPLETE")
    println("="^70)
    println("Total results: $(length(results))")
    println("Output saved to: $output_dir/")
    println()

    println("Summary of Rényi-2 Binder (after last X dephasing) across scan:")
    println("  Min B = $(round(minimum(r["B"] for r in results), digits=4))")
    println("  Max B = $(round(maximum(r["B"] for r in results), digits=4))")
    M2_vals = [r["M2_bar"] for r in results]
    M4_vals = [r["M4_bar"] for r in results]
    println("  M2_bar range: [$(round(minimum(M2_vals), digits=4)), $(round(maximum(M2_vals), digits=4))]")
    println("  M4_bar range: [$(round(minimum(M4_vals), digits=4)), $(round(maximum(M4_vals), digits=4))]")
end

main()
