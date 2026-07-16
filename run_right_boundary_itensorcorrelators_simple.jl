#!/usr/bin/env julia

"""
Right Boundary Scan (ITensorCorrelators.jl variant) - local convenience script

Scan q = P_x = P_zz at fixed λ_x = 0.7, λ_zz = 0.
Focus on measurement-induced entanglement transition using the proposal-aligned
Rényi-2 Binder, with the Rényi-2 moments computed via ITensorCorrelators.jl
n-point correlators instead of a global replica-overlap MPO.

Physics is the doubled-MPS Born-sampling model validated in
`right_boundary_renyi2_proposal_aligned_v4.ipynb`, shared with the cluster
entry point `run_right_boundary_itensorcorrelators_scan.jl` via
`renyi2_right_boundary_itensorcorrelators_core.jl`.
"""

using Dates
using JSON
using Statistics

include("renyi2_right_boundary_itensorcorrelators_core.jl")

# === EDIT THESE PARAMETERS ===
L = 8                    # System size
lambda_x = 0.7          # X measurement strength (FIXED)
lambda_zz = 0.0         # ZZ measurement strength (FIXED - no weak ZZ measurement)
P_steps = 11            # Number of q values
P_min = 0.0             # Minimum q
P_max = 0.5             # Maximum q

ntrials = 500           # Born-sampled Monte Carlo trajectories
maxdim = 256            # Max bond dimension (dynamics)
cutoff = 1e-12          # SVD truncation cutoff (dynamics)
T_max_factor = 4        # T_max = T_max_factor * L
obs_maxdim_factor = 4   # Unused by the ITensorCorrelators observable (kept for parity)
obs_cutoff = 1e-14      # Unused by the ITensorCorrelators observable (kept for parity)
nboot = 1000            # Bootstrap resamples for B standard error
seed = 42               # Random seed
output_dir = "right_boundary_itensorcorrelators_results"
# =============================

@assert isapprox(lambda_x, 0.7; atol=1e-9) "Proposal-aligned right edge requires lambda_x = 0.7."
@assert isapprox(lambda_zz, 0.0; atol=1e-9) "Proposal-aligned right edge requires lambda_zz = 0."

println("="^70)
println("RIGHT BOUNDARY SCAN (ITensorCorrelators.jl): Rényi-2 Binder Analysis")
println("="^70)
println("Physics: Measurement-induced transition (λ_x = 0.7, λ_zz = 0, q_x=q_zz=q)")
println()
println("Parameters:")
println("  L = $L")
println("  λ_x = $lambda_x (X measurement strength - FIXED)")
println("  λ_zz = $lambda_zz (no weak ZZ measurements)")
println("  q scan: [$P_min, $P_max] with $P_steps points (q_x = q_zz = q)")
println("  Trajectories per point: $ntrials")
println("  T_max = $T_max_factor * L = $(T_max_factor * L)")
println("  Focus: proposal-aligned, trajectory-averaged Rényi-2 Binder (ITensorCorrelators.jl)")
println("="^70)
println()

# Create output directory
mkpath(output_dir)

# Create q values
P_values = P_steps <= 1 ? [P_min] : range(P_min, P_max, length=P_steps)

results = []

for (i, q) in enumerate(P_values)
    println("\n[$i/$P_steps] q = $(round(q, digits=4)) (P_x = P_zz = q)")
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

    # Store results
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
        "max_interphysical_linkdim" => result.max_interphysical_linkdim,
        "max_trace_error" => result.max_trace_error,
        "observable_method" => "itensorcorrelators",
        "time_seconds" => t_elapsed
    )

    push!(results, result_dict)

    # Print results
    println("  B (proposal, trajectory-averaged) = $(round(result.B, digits=4)) ± $(round(result.B_bootstrap_se, digits=4)) (bootstrap SE)")
    println("  M2_bar = $(round(result.M2_bar, digits=4))")
    println("  M4_bar = $(round(result.M4_bar, digits=4))")
    println("  Purity = $(round(result.purity_bar, digits=4))")
    println("  Valid/Total: $(result.n_valid)/$(result.ntrials)")
    println("  Time: $(round(t_elapsed, digits=1))s")

    # Save intermediate results
    timestamp = Dates.format(Dates.now(), "yyyymmdd_HHMM")
    output_file = joinpath(output_dir, "right_boundary_itensorcorrelators_L$(L)_lambdax$(lambda_x)_$(timestamp).json")
    open(output_file, "w") do f
        JSON.print(f, results, 4)
    end
end

println("\n" * "="^70)
println("✓ SCAN COMPLETE")
println("="^70)
println("Results saved to: $output_dir/")

# Print summary
println("\nSummary:")
for (i, r) in enumerate(results)
    println("  P_x=$(round(r["P_x"], digits=2)): B=$(round(r["B"], digits=4)), M2=$(round(r["M2_bar"], digits=4)), M4=$(round(r["M4_bar"], digits=4))")
end
