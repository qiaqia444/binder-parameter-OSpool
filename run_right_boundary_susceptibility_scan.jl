#!/usr/bin/env julia

"""
Right Boundary Scan (Rényi-2 Overlap Susceptibility, ITensorCorrelators.jl):
Measurement-Induced Phase Transition

Scan dephasing strength q = P_x = P_zz from 0 to 1/2 at FIXED measurement
strength λ_x = 0.7 (λ_zz = 0). This explores the right boundary of the
phase diagram.

Physics is the same doubled-MPS Born-sampling model validated in
`right_boundary_renyi2_proposal_aligned_v4.ipynb` and shared with
`run_right_boundary_itensorcorrelators_scan.jl`, but the observable here is
the Rényi-2 replica-overlap susceptibility

    chi2(m) = <Q^2>_m / L = L * M2(m),   Q = sum_i Z_i^bra Z_i^ket

computed per Born-sampled trajectory via ITensorCorrelators.jl n-point
correlators (see `renyi2_right_boundary_susceptibility_core.jl`), NOT the
Rényi-2 Binder ratio. The reported estimate is the trajectory average
chi2_bar = mean_m chi2(m); the connected susceptibility
[<Q^2> - <Q>^2]/L is kept as a diagnostic.
"""

using ArgParse
using JSON
using Dates
using Statistics

include("renyi2_right_boundary_susceptibility_core.jl")

function parse_commandline()
    s = ArgParseSettings(description = "Right boundary scan: vary q=P_x=P_zz at fixed λ_x, Rényi-2 overlap susceptibility (ITensorCorrelators.jl)")

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
            help = "Number of full layers is T_max_factor * L"
            arg_type = Int
            default = 4

        "--obs_maxdim_factor"
            help = "Unused by the ITensorCorrelators observable (kept for CLI compatibility)"
            arg_type = Int
            default = 4

        "--obs_cutoff"
            help = "Unused by the ITensorCorrelators observable (kept for CLI compatibility)"
            arg_type = Float64
            default = 1e-14

        "--nboot"
            help = "Number of bootstrap resamples for the standard error of chi2"
            arg_type = Int
            default = 1000

        "--seed"
            help = "Random seed"
            arg_type = Int
            default = 42

        "--output_dir"
            help = "Output directory"
            arg_type = String
            default = "right_boundary_susceptibility_results"

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
    println("RIGHT BOUNDARY SCAN: Rényi-2 Overlap Susceptibility (ITensorCorrelators.jl)")
    println("="^70)
    println("Physics: lambda=1, delta=0.7 => lambda_x=0.7, lambda_zz=0, q_x=q_zz=q")
    println("Method: doubled-MPS Born sampling (proposal-aligned)")
    println("Observable: trajectory-averaged Rényi-2 overlap susceptibility chi2 = <Q^2>/L")
    println()
    println("Parameters:")
    println("  System size L = $L")
    println("  X measurement λ_x = $lambda_x (FIXED)")
    println("  ZZ measurement λ_zz = $lambda_zz (FIXED - no weak ZZ measurement)")
    println("  q scan: q_x = q_zz = q ∈ [$P_min, $P_max] with $P_steps points")
    println("  Trajectories per point: $ntrials")
    println("  T_max = $T_max_factor * L = $(T_max_factor * L)")
    println("  Dynamics: maxdim=$maxdim, cutoff=$cutoff")
    println("="^70)
    println()

    # Create q values to scan (q_x = q_zz = q, matching the right-edge slice)
    q_values = P_steps <= 1 ? [P_min] : collect(range(P_min, P_max, length=P_steps))

    # Create output directory
    mkpath(output_dir)

    results = []

    for (i, q) in enumerate(q_values)
        println("\n[$i/$(length(q_values))] Running q = $(round(q, digits=4)) (P_x = P_zz = q)")
        println("  (λ_x = $lambda_x, λ_zz = $lambda_zz)")
        println("-"^70)

        t_start = time()

        result = rb_run_right_edge_susceptibility_point(
            L, Float64(q);
            lambda_x = lambda_x,
            lambda_zz = lambda_zz,
            ntrials = ntrials,
            T_max = T_max_factor * L,
            maxdim = maxdim,
            cutoff = cutoff,
            obs_maxdim = obs_maxdim_factor * maxdim,
            obs_cutoff = obs_cutoff,
            seed = seed + i,  # Different seed for each q
            nboot = nboot,
        )

        t_elapsed = time() - t_start

        # Store results
        result_dict = Dict(
            "L" => L,
            "lambda_x" => lambda_x,
            "lambda_zz" => lambda_zz,
            "P_x" => q,
            "P_zz" => q,
            "chi2" => result.chi,
            "chi2_mean_of_trials" => result.chi_mean_of_trials,
            "chi2_std_of_trials" => result.chi_std_of_trials,
            "chi2_connected_bar" => result.chi_connected_bar,
            "M2_bar" => result.M2_bar,
            "overlap_density_bar" => result.overlap_density_bar,
            "purity_bar" => result.purity_bar,
            "chi2_bootstrap_se" => result.chi_bootstrap_se,
            "chi2_ci_low" => result.chi_ci_low,
            "chi2_ci_high" => result.chi_ci_high,
            "ntrials" => result.ntrials,
            "n_valid" => result.n_valid,
            "n_invalid" => result.n_invalid,
            "maxdim" => maxdim,
            "cutoff" => cutoff,
            "T_max" => T_max_factor * L,
            "max_interphysical_linkdim" => result.max_interphysical_linkdim,
            "max_trace_error" => result.max_trace_error,
            "observable_method" => "itensorcorrelators",
            "time_seconds" => t_elapsed
        )

        push!(results, result_dict)

        # Print results
        println("  chi2 (proposal, trajectory-averaged) = $(round(result.chi, digits=4)) ± $(round(result.chi_bootstrap_se, digits=4)) (bootstrap SE)")
        println("  chi2_connected_bar = $(round(result.chi_connected_bar, digits=6))")
        println("  M2_bar = $(round(result.M2_bar, digits=6))")
        println("  Purity = $(round(result.purity_bar, digits=4))")
        println("  Valid/Total = $(result.n_valid)/$(result.ntrials)")
        println("  Max trace error = $(result.max_trace_error)")
        println("  Time: $(round(t_elapsed, digits=1)) seconds")

        # Save intermediate results
        if !isempty(output_file)
            # Use specified filename for cluster jobs
            outpath = joinpath(output_dir, output_file)
        else
            # Use timestamped filename for local runs
            timestamp = Dates.format(now(), "yyyymmdd_HHMM")
            outpath = joinpath(output_dir, "right_boundary_susceptibility_L$(L)_lambda$(lambda_x)_$(timestamp).json")
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
