#!/usr/bin/env julia

"""
Right Boundary Scan (Rényi-2 Overlap Susceptibility, ITensorCorrelators.jl)
- local convenience script

Scan q = P_x = P_zz at fixed λ_x = 0.7, λ_zz = 0.
Observable: trajectory-averaged Rényi-2 replica-overlap susceptibility
chi2 = <Q^2>/L, computed via ITensorCorrelators.jl n-point correlators
(see `renyi2_right_boundary_susceptibility_core.jl`).

Physics is the doubled-MPS Born-sampling model validated in
`right_boundary_renyi2_proposal_aligned_v4.ipynb`, shared with the cluster
entry point `run_right_boundary_susceptibility_scan.jl`.
"""

using Dates
using JSON
using Statistics

include("renyi2_right_boundary_susceptibility_core.jl")

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
nboot = 1000            # Bootstrap resamples for chi2 standard error
seed = 42               # Random seed
output_dir = "right_boundary_susceptibility_results"
# =============================

@assert isapprox(lambda_x, 0.7; atol=1e-9) "Proposal-aligned right edge requires lambda_x = 0.7."
@assert isapprox(lambda_zz, 0.0; atol=1e-9) "Proposal-aligned right edge requires lambda_zz = 0."

println("="^70)
println("RIGHT BOUNDARY SCAN: Rényi-2 Overlap Susceptibility Analysis")
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
println("  Focus: proposal-aligned, trajectory-averaged Rényi-2 overlap susceptibility (ITensorCorrelators.jl)")
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
        "max_interphysical_linkdim" => result.max_interphysical_linkdim,
        "max_trace_error" => result.max_trace_error,
        "observable_method" => "itensorcorrelators",
        "time_seconds" => t_elapsed
    )

    push!(results, result_dict)

    # Print results
    println("  chi2 (proposal, trajectory-averaged) = $(round(result.chi, digits=4)) ± $(round(result.chi_bootstrap_se, digits=4)) (bootstrap SE)")
    println("  chi2_connected_bar = $(round(result.chi_connected_bar, digits=4))")
    println("  M2_bar = $(round(result.M2_bar, digits=4))")
    println("  Purity = $(round(result.purity_bar, digits=4))")
    println("  Valid/Total: $(result.n_valid)/$(result.ntrials)")
    println("  Time: $(round(t_elapsed, digits=1))s")

    # Save intermediate results
    timestamp = Dates.format(Dates.now(), "yyyymmdd_HHMM")
    output_file = joinpath(output_dir, "right_boundary_susceptibility_L$(L)_lambdax$(lambda_x)_$(timestamp).json")
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
    println("  P_x=$(round(r["P_x"], digits=2)): chi2=$(round(r["chi2"], digits=4)), M2=$(round(r["M2_bar"], digits=4))")
end
