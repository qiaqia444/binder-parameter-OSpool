#!/usr/bin/env julia

"""
Right Boundary Scan: dual Edwards-Anderson (EA) Binder parameter

Scan the dephasing strength q = P_x = P_zz over [P_min, P_max] (default [0, 0.5]) at
FIXED measurement strengths lambda_x = 0.7, lambda_zz = 0 (the right boundary).

Dynamics: the existing right-boundary doubled-MPS Born-sampling evolution
(`rb_evolve_right_edge_one_trial` in `renyi2_right_boundary_core.jl`), reused unchanged,
with the same per-point seed (`seed + i`) and per-trajectory seed derivation as
`run_right_boundary_scan.jl` - so the same seed analyses the same Born trajectories as
the Renyi-2 Binder scan.

Observable (see `dual_ea_right_boundary_core.jl`): open-chain dual endpoints a = 2..L
(N = L-1), X-strings W_ab, G2 = E(W_ab)^2, G4 = E(W_uv W_wx)^2 with ordinary normalised
trace expectations E(O) = Tr(rho O)/Tr(rho), measured on the final density matrix of
each trajectory (after the last ZZ-dephasing step of layer t = T_max).

Two Binder conventions are stored under separate names:
  B_mean_of_trials = mean_m [1 - M4_m/(3 M2_m^2)]
  B_ratio_of_means = 1 - mean(M4_m)/(3 mean(M2_m)^2)
There is intentionally NO bare "B", "S2_bar" or "S4_bar" key (those mean different
things in the left-boundary EA driver and the Renyi-2 right-boundary driver).

Smoke test:
  julia --project=. run_dual_ea_right_boundary_scan.jl --L 6 --P_min 0 --P_max 0.5 \\
      --P_steps 3 --ntrials 4 --maxdim 64 --nboot 100 --output_dir /tmp/dual_ea_smoke
Example scan:
  julia --project=. run_dual_ea_right_boundary_scan.jl --L 16 --P_min 0 --P_max 0.5 \\
      --P_steps 26 --ntrials 200 --maxdim 256 --cutoff 1e-12 --seed 42 \\
      --output_dir dual_ea_right_boundary_results
"""

using ArgParse
using JSON
using Dates
using Statistics

include("dual_ea_right_boundary_core.jl")

function parse_commandline()
    s = ArgParseSettings(description = "Right boundary scan: dual EA Binder parameter vs q = P_x = P_zz at fixed λ_x")

    @add_arg_table! s begin
        "--L"
            help = "System size"
            arg_type = Int
            default = 12

        "--lambda_x"
            help = "X measurement strength (FIXED; right boundary uses 0.7)"
            arg_type = Float64
            default = 0.7

        "--lambda_zz"
            help = "ZZ measurement strength (FIXED, must be 0.0 for the right boundary)"
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
            help = "Number of full layers is T_max_factor * L (observable measured after the last layer)"
            arg_type = Int
            default = 4

        "--nboot"
            help = "Bootstrap resamples for the standard error of B_ratio_of_means (0 disables)"
            arg_type = Int
            default = 1000

        "--seed"
            help = "Random seed"
            arg_type = Int
            default = 42

        "--verbose"
            help = "Print per-trajectory moments"
            action = :store_true

        "--output_dir"
            help = "Output directory"
            arg_type = String
            default = "dual_ea_right_boundary_results"

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
    nboot = args["nboot"]
    seed = args["seed"]
    output_dir = args["output_dir"]
    output_file = args["output_file"]

    @assert L >= 2 "Need L >= 2 (N = L-1 >= 1 dual endpoints)"
    @assert isapprox(lambda_zz, 0.0; atol=1e-9) "Right boundary requires lambda_zz = 0 (no weak ZZ measurement)."
    @assert 0.0 <= P_min <= P_max <= 0.5 "q must lie in [0, 0.5]"

    N = L - 1
    T_max = T_max_factor * L
    cost = dual_ea_cost_estimate(L, maxdim)

    println("="^70)
    println("RIGHT BOUNDARY SCAN: dual Edwards-Anderson Binder parameter")
    println("="^70)
    println("Dynamics: existing right-boundary doubled-MPS Born sampling (reused unchanged)")
    println("Observable: dual EA, open-chain endpoints a = 2..L (N = $N), trace expectations squared")
    println()
    println("Parameters:")
    println("  System size L = $L (N = $N dual endpoints)")
    println("  X measurement λ_x = $lambda_x, ZZ measurement λ_zz = $lambda_zz")
    println("  q scan: q_x = q_zz = q ∈ [$P_min, $P_max] with $P_steps points")
    println("  Trajectories per point: $ntrials")
    println("  T_max = $T_max_factor * L = $T_max (observable after the last ZZ-dephasing step)")
    println("  Dynamics: maxdim=$maxdim, cutoff=$cutoff (observable contraction is exact)")
    println("  Observable cost/trajectory (upper bound at chi=maxdim): " *
            "$(cost.n_sorted_quads) sorted 4-tuples ($(cost.n_ordered_quads) ordered), " *
            "~$(round(cost.approx_flops / 1e9, digits=2)) GFLOP")
    println("="^70)
    println()

    q_values = P_steps <= 1 ? [P_min] : collect(range(P_min, P_max, length=P_steps))

    mkpath(output_dir)
    if !isempty(output_file)
        outpath = joinpath(output_dir, output_file)
    else
        timestamp = Dates.format(now(), "yyyymmdd_HHMM")
        outpath = joinpath(output_dir, "dual_ea_right_boundary_L$(L)_lambda$(lambda_x)_$(timestamp).json")
    end

    results = []

    for (i, q) in enumerate(q_values)
        println("\n[$i/$(length(q_values))] Running q = $(round(q, digits=4)) (P_x = P_zz = q)")
        println("  (λ_x = $lambda_x, λ_zz = $lambda_zz)")
        println("-"^70)

        t_start = time()
        point_seed = seed + i   # same per-q seed rule as run_right_boundary_scan.jl

        r = dual_ea_run_right_edge_point(
            L, Float64(q);
            lambda_x = lambda_x,
            lambda_zz = lambda_zz,
            ntrials = ntrials,
            T_max = T_max,
            maxdim = maxdim,
            cutoff = cutoff,
            seed = point_seed,
            nboot = nboot,
            verbose = args["verbose"],
        )

        t_elapsed = time() - t_start

        result_dict = Dict(
            "observable" => "dual_EA_binder",
            "L" => L,
            "q" => q,
            "P_x" => q,
            "P_zz" => q,
            "lambda_x" => lambda_x,
            "lambda_zz" => lambda_zz,
            "boundary_conditions" => "open",
            "endpoint_convention" => "dual endpoints a=2..L, N=L-1; W_ab = prod_{k=min(a,b)}^{max(a,b)-1} X_k, W_aa = I",
            "N_endpoints" => N,
            "endpoint_min" => 2,
            "endpoint_max" => L,
            "expectation" => "E(O)=Tr(rho_m O)/Tr(rho_m); G2=E(W_ab)^2, G4=E(W_uv W_wx)^2 (no purity / HS normalisation)",
            "measurement_time" => "after the last ZZ-dephasing step of layer t=T_max",
            "T_max" => T_max,
            "T_max_factor" => T_max_factor,
            "maxdim" => maxdim,
            "cutoff" => cutoff,
            "seed" => seed,
            "point_seed" => point_seed,
            "trial_seeds" => r.trial_seeds,
            "ntrials_requested" => r.ntrials_requested,
            "n_completed" => r.n_completed,
            "n_finite_B" => r.n_finite_B,
            "B_mean_of_trials" => r.B_mean_of_trials,
            "B_std_of_trials" => r.B_std_of_trials,
            "B_se_of_trials" => r.B_se_of_trials,
            "B_ratio_of_means" => r.B_ratio_of_means,
            "B_ratio_of_means_bootstrap_se" => r.B_ratio_of_means_bootstrap_se,
            "nboot" => nboot,
            "M2_mean" => r.M2_mean,
            "M4_mean" => r.M4_mean,
            "M2_se" => r.M2_se,
            "M4_se" => r.M4_se,
            "M2_per_trajectory" => r.M2s,
            "M4_per_trajectory" => r.M4s,
            "B_per_trajectory" => r.Bs,
            "max_interphysical_linkdim" => r.max_interphysical_linkdim,
            "max_trace_error" => r.max_trace_error,
            "n_sorted_quads" => cost.n_sorted_quads,
            "n_ordered_quads" => cost.n_ordered_quads,
            "observable_seconds_mean" => r.obs_seconds_mean,
            "observable_seconds_max" => r.obs_seconds_max,
            "time_seconds" => t_elapsed,
        )

        push!(results, result_dict)

        println("  B_mean_of_trials = $(round(r.B_mean_of_trials, digits=4)) ± $(round(r.B_se_of_trials, digits=4)) (SE; std = $(round(r.B_std_of_trials, digits=4)))")
        println("  B_ratio_of_means = $(round(r.B_ratio_of_means, digits=4)) ± $(round(r.B_ratio_of_means_bootstrap_se, digits=4)) (bootstrap SE)")
        println("  M2_mean = $(round(r.M2_mean, digits=6)),  M4_mean = $(round(r.M4_mean, digits=6))")
        println("  Completed = $(r.n_completed)/$(r.ntrials_requested),  max χ = $(r.max_interphysical_linkdim),  max trace error = $(r.max_trace_error)")
        println("  Observable time/trajectory: mean $(round(r.obs_seconds_mean, digits=3)) s, max $(round(r.obs_seconds_max, digits=3)) s")
        println("  Time: $(round(t_elapsed, digits=1)) seconds")

        # Save after every q so an interrupted scan keeps its finished points
        open(outpath, "w") do f
            JSON.print(f, results, 4)
        end
    end

    println("\n" * "="^70)
    println("✓ SCAN COMPLETE")
    println("="^70)
    println("Total results: $(length(results))")
    println("Output saved to: $outpath")
    println()
    println("Summary of dual-EA Binder across scan:")
    for (name, key) in (("B_mean_of_trials", "B_mean_of_trials"), ("B_ratio_of_means", "B_ratio_of_means"))
        vals = [r[key] for r in results]
        println("  $name: min = $(round(minimum(vals), digits=4)), max = $(round(maximum(vals), digits=4))")
    end
end

main()
